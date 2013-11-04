#! /usr/bin/env python

# Copyright (c) 2013 Genome Research Ltd. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Author: Iain Bancarz, ib5@sanger.ac.uk


"""
Wrapper for Plink diff functionality
Capture values from text output and reprocess into more useful formats

Three component classes:
- Parser for Plink output
- Data structure to hold raw results
- Interpreter to reprocess results and write to file

Key diff metrics:
- What is the intersection of samples/SNPs between datasets? Are there any 
  samples/SNPs in one dataset which are not present in the other?
- Are datasets equivalent? (Yes/no)
- Are datasets equivalent to within 'flip' (transpose of major/minor alleles)?
- Concordance on all shared calls (SNP/sample pairs)
- Concordance on calls which are both non-null
- Rate of calls with allele 'flip'
- Proportion of mismatches equivalent to within flip
- Breakdown of mismatches by SNP and sample

Note: Flips can be arbitrarily introduced by Plink operations such as recode, 
and are reported as mismatches by the Plink diff.
"""

import json, gzip, os, re, sys
from plink import PlinkValidator, PlinkToolsError

class PlinkDiffShared:
    """Define some shared constants for use as dictionary keys"""

    # Variables for both entire diff, and sample/SNP breakdowns
    COUNT_KEYS = ['MISMATCH', 'MISMATCH_NON_NULL', 'ALLELE_FLIP']
    MISMATCH_KEY = COUNT_KEYS[0]
    MISMATCH_NON_NULL_KEY = COUNT_KEYS[1]
    ALLELE_FLIP_KEY = COUNT_KEYS[2]

    # Variables for the diff as a whole
    TOTAL_KEYS = ['CALLS', 'CALLS_NON_NULL', 'SAMPLES', 'SNPS', 'CONCORDANT']
    CALLS_KEY = TOTAL_KEYS[0] # total shared (sample, snp) pairs
    CALLS_NON_NULL_KEY = TOTAL_KEYS[1] # as above, but both non-null
    SAMPLES_KEY = TOTAL_KEYS[2] # total shared samples
    SNPS_KEY = TOTAL_KEYS[3] # total shared snps
    CONCORDANT_KEY = TOTAL_KEYS[4] # concordance (as reported by Plink)

    # Variables for input datasets
    INPUT_KEYS = ['SNPS_1', 'SNPS_2', 'SAMPLES_1', 'SAMPLES_2', 
                  'NEW_SNPS', 'NEW_SAMPLES']
    SNPS_1_KEY = INPUT_KEYS[0] # size of SNP sets 1 & 2
    SNPS_2_KEY = INPUT_KEYS[1]
    SAMPLES_1_KEY = INPUT_KEYS[2] # size of sample sets 1 & 2
    SAMPLES_2_KEY = INPUT_KEYS[3]
    NEW_SNPS_KEY = INPUT_KEYS[4] # SNPs in 2 but not in 1
    NEW_SAMPLES_KEY = INPUT_KEYS[5] # Samples in 2 but not in 1


class PlinkDiffParser(PlinkDiffShared):
    """Run Plink's diff function and parse the output

    Requires horrible text munging
    Tested on Plink 1.0.7, may be broken by future Plink updates
    (Could replace by a class to parse .bed files and find diff results)

    Stores results in a PlinkDiffData object
    """

    def __init__(self):
        self.validator = PlinkValidator() # checks if Plink is available
        self.data = PlinkDiffData()
        self.verbose = False
        # compile regular expressions to parse Plink .log file
        logPatterns = [
            '\d+ markers to be included from', # first input dataset
            '\d+ individuals read from ',
            '\d+ markers to be merged from', # second input dataset
            '\d+ individuals merged from',
            'Of \d+ overlapping SNPs, \d+ were both genotyped',
            'and \d+ were concordant',
            # shared SNPs
            'Of these, \d+ are new, \d+ already exist in current data', 
            # shared samples
            'Of these, \d+ were new, \d+ were already in current data'
            ] 
        self.logExprs = []
        for p in logPatterns: self.logExprs.append(re.compile(p))

    def getData(self):
        return self.data

    def parseDiffFile(self, stem, removeInput=False):
        """Parse .diff output file from Plink, compress and find diff stats
        
        The .diff file has one line for each (snp, sample) pair that differs
        May be very large! Do not attempt to read entire file into memory
        Convert .diff file into compressed .csv
        """
        diffPath = stem+'.diff'
        outPath = stem+'.diff.csv.gz'
        csvHeader = 'SNP,FID,IID,NEW,OLD\n'
        if not (os.path.exists(diffPath)):
            raise PlinkToolsError("Plink .diff output file not found!")
        elif self.verbose:
            print "Plink .diff output file found."
        inFile = open(diffPath, 'r')
        outFile = gzip.open(outPath, 'w')
        outFile.write(csvHeader)
        firstLine = True
        while True:
            line = inFile.readline()
            if line=='':  # end of file
                break
            elif firstLine: # skip header line
                firstLine = False
                continue
            self.parseDiffLine(line, outFile)
        inFile.close()
        outFile.close()
        if self.verbose: print "Results read from Plink .diff output."
        if removeInput:
            os.remove(diffPath)
            if self.verbose: print "Removed original .diff file."

    def parseDiffLine(self, line, outFile):
        """Line of .diff file represents a mismatched (snp, sample) pair"""
        fields = re.split('\s+', line.strip())
        outFile.write(','.join(fields)+"\n")
        (snp, fid, iid, newCall, oldCall) = fields
        sample = (fid, iid) # sample ID = (family, individual)
        self.data.incrementGlobal(self.MISMATCH_KEY)
        self.data.incrementSample(sample, self.MISMATCH_KEY)
        self.data.incrementSnp(snp, self.MISMATCH_KEY)
        if newCall != "0/0" and oldCall != "0/0":
            # both calls are non-null
            self.data.incrementGlobal(self.MISMATCH_NON_NULL_KEY)
            self.data.incrementSample(sample, self.MISMATCH_NON_NULL_KEY)
            self.data.incrementSnp(snp, self.MISMATCH_NON_NULL_KEY)
        if newCall[::-1] == oldCall: 
            # calls are non-null, and differ only in allele reversal
            self.data.incrementGlobal(self.ALLELE_FLIP_KEY)
            self.data.incrementSample(sample, self.ALLELE_FLIP_KEY)
            self.data.incrementSnp(snp, self.ALLELE_FLIP_KEY)

    def parseLogFile(self, stem):
        """Munge the .log text from Plink to find totals for the dataset. 
        Unfortunately this is the only way to find these values from Plink 
        diff output."""
        logPath = stem+'.log'
        if not (os.path.exists(logPath)):
            raise PlinkToolsError("Plink .log output file not found!")
        elif self.verbose:
            print "Plink .log output file found."
        lines = open(logPath).readlines()
        # find: intersection, bothGenotyped, concordant, snps, samples
        logKeys = self.TOTAL_KEYS + self.INPUT_KEYS
        for key in logKeys:
            # set to None instead of 0 to detect parsing failure
            self.data.setGlobal(key, None)
        for i in range(len(lines)):
            line = lines[i]
            words = re.split('\s+', line.strip())
            if self.logExprs[0].match(line):
                self.data.setGlobal(self.SNPS_1_KEY, int(words[0]))
            elif self.logExprs[1].match(line):
                self.data.setGlobal(self.SAMPLES_1_KEY, int(words[0]))
            elif self.logExprs[2].match(line):
                self.data.setGlobal(self.SNPS_2_KEY, int(words[0]))
            elif self.logExprs[3].match(line):
                self.data.setGlobal(self.SAMPLES_2_KEY, int(words[0]))
            elif self.logExprs[4].match(line):
                self.data.setGlobal(self.CALLS_KEY, int(words[1]))
                self.data.setGlobal(self.CALLS_NON_NULL_KEY, int(words[4]))
            elif self.logExprs[5].match(line):
                self.data.setGlobal(self.CONCORDANT_KEY, int(words[1]))
            elif self.logExprs[6].match(line) and \
                    re.search('markers', lines[i-1]):
                self.data.setGlobal(self.NEW_SNPS_KEY, int(words[2]))
                self.data.setGlobal(self.SNPS_KEY, int(words[5]))
            elif self.logExprs[7].match(line) and \
                    re.search('individuals',lines[i-1]):
                self.data.setGlobal(self.NEW_SAMPLES_KEY, int(words[2]))
                self.data.setGlobal(self.SAMPLES_KEY, int(words[5]))
        parsedOK = True
        for key in logKeys:
            if self.data.getGlobal(key)==None:
                msg = "ERROR: Failed to parse "+key+" from Plink output\n";
                sys.stderr.write(msg)
                parsedOK = False
        if not parsedOK:
            raise PlinkToolsError("Failed to parse Plink diff .log file!")

    def run(self, stem1, stem2, outStem, cleanup=True, verbose=False):
        """Main method to run Plink diff and parse the output """
        self.verbose = verbose
        # if Plink diff output already present, do not run again
        if not (os.path.exists(outStem+".log") 
                and os.path.exists(outStem+".diff")):
            self.runBinaryDiff(stem1, stem2, outStem)
        self.parseDiffFile(outStem, cleanup)
        self.parseLogFile(outStem)
        return self.data

    def runBinaryDiff(self, stem1, stem2, outStem):
        """Run Plink executable with appropriate options to find diff 

        Inputs: Stems for two Plink binary input datasets, and for output
        Stem is Plink file path without .bed, .bim, .fam extension"""
        bmerge = "%s.bed %s.bim %s.fam" % (stem2, stem2, stem2)
        cmd = "plink --noweb --bfile "+stem1+" --bmerge "+bmerge+\
            " --make-bed --merge-mode 6 --out "+outStem
        if not self.verbose: cmd = cmd+" > /dev/null"
        status = os.system(cmd)
        if status!=0:
            msg = "Non-zero exit status from Plink executable!"
            raise PlinkToolsError(msg)

class PlinkDiffData(PlinkDiffShared):

    """Container for diff results:
    - Shared (snp, sample) pairs
    - Non-null shared pairs
    - Mismatched pairs
    - Mismatched non-null pairs
    - Mismatched pairs differing only by a major/minor allele flip

    Grand totals, and breakdowns by sample/SNP
    """

    def __init__(self):
        self.globalTotals = {}
        self.snpTotals = {}
        self.sampleTotals = {}
        self.inputs = {}

        self.congruentSNPs = True
        self.congruentSamples = True
        for key in self.COUNT_KEYS: self.globalTotals[key] = 0
        for key in self.TOTAL_KEYS: self.globalTotals[key] = 0
        for key in self.INPUT_KEYS: self.globalTotals[key] = 0
        self.globalKeys = self.COUNT_KEYS + self.TOTAL_KEYS + self.INPUT_KEYS

    def getGlobalDictionary(self):
        return self.globalTotals

    def getGlobal(self, key):
        return self.globalTotals[key]

    def setGlobal(self, key, value):
        if key not in self.globalKeys: 
            raise ValueError("Invalid key for global diff result: "+str(key))
        self.globalTotals[key] = value

    def incrementGlobal(self, key):
        self.globalTotals[key] += 1

    def addSample(self, sample):
        self.sampleTotals[sample] = {}
        for key in self.COUNT_KEYS:
            self.sampleTotals[sample][key] = 0

    def getSample(self, sample, key):
        return self.sampleTotals[sample][key]

    def setSample(self, sample, key, value):
        if key not in self.COUNT_KEYS:
            raise ValueError("Invalid key for sample diff: "+str(key))
        try:
            self.sampleTotals[sample][key] = value
        except KeyError:
            self.addSample(sample)
            self.sampleTotals[sample][key] = value

    def incrementSample(self, sample, key):
        try: 
            self.sampleTotals[sample][key] += 1
        except KeyError: 
            self.addSample(sample)
            self.sampleTotals[sample][key] += 1

    def addSnp(self, snp):
        self.snpTotals[snp] = {}
        for key in self.COUNT_KEYS:
            self.snpTotals[snp][key] = 0

    def getSnp(self, snp, key):
        return self.snpTotals[snp][key]

    def setSnp(self, snp, key, value):
        if key not in self.COUNT_KEYS:
            raise ValueError("Invalid key for SNP diff: "+str(key))
        try:
            self.snpTotals[snp][key] = value
        except KeyError:
            self.addSnp(snp)
            self.snpTotals[snp][key] = value

    def incrementSnp(self, snp, key):
        try: 
            self.snpTotals[snp][key] += 1
        except KeyError: 
            self.addSnp(snp)
            self.snpTotals[snp][key] += 1

    def getSnpList(self):
        return self.snpTotals.keys()

    def getSampleList(self):
        return self.sampleTotals.keys()

class PlinkDiffWriter(PlinkDiffShared):
    """Take diff data, calculate additional stats and write output"""

    def __init__(self, data, verbose=True):
        self.data = data
        self.verbose = verbose
        self.mismatchRate = 0.0
        self.mismatchNNRate = 0.0
        self.flipRate = 0.0
        self.flipMismatchRate = 0.0
        self.congruentSNPs = False
        self.congruentSamples = False
        self.findConcordance()
        self.findCongruence()
        self.equivalent = self.findEquivalence(allowFlip=False)
        self.equivalentWithinFlip = self.findEquivalence(allowFlip=True)

    def findConcordance(self):
        """ Find concordance on shared (sample, snp) pairs """
        calls = self.data.getGlobal(self.CALLS_KEY)
        callsNN = self.data.getGlobal(self.CALLS_NON_NULL_KEY)
        mismatch = self.data.getGlobal(self.MISMATCH_KEY)
        mismatchNN = self.data.getGlobal(self.MISMATCH_NON_NULL_KEY)
        flip = self.data.getGlobal(self.ALLELE_FLIP_KEY)
        # make sure we do not divide by zero
        if calls==0:
            sys.stderr.write("No shared calls in Plink data\n")
            self.mismatchRate = None
        else:
            self.mismatchRate = mismatch / float(calls)
        if callsNN==0:
            sys.stderr.write("No non-null shared calls in Plink data\n")
            self.mismatchNNRate = None
            self.flipRate = None
        else:
            self.mismatchNNRate = mismatchNN / float(callsNN)
            self.flipRate = flip / float(callsNN)
        if mismatch==0:
            sys.stderr.write("No mismatched calls in Plink data\n")
            self.flipMismatchRate = None
        else:
            self.flipMismatchRate = flip / float(mismatch)
        if self.verbose:
            print "Concordance: %.4f" % (1 - self.mismatchRate)
            print "Concordance on non-null calls: %.4f" % \
                (1-self.mismatchNNRate)
            print "Allele flip rate: %.4f" % self.flipRate
            print "Flips as proportion of mismatches: %.4f" % \
                self.flipMismatchRate

    def findCongruence(self):
        """Find congruence (or otherwise) of input SNP and samples sets"""
        if self.data.getGlobal(self.SNPS_1_KEY) == self.data.getGlobal(self.SNPS_2_KEY) and self.data.getGlobal(self.NEW_SNPS_KEY) == 0:
            self.congruentSNPs = True
        else:
            self.congruentSNPs = False
        if self.data.getGlobal(self.SAMPLES_1_KEY) == self.data.getGlobal(self.SAMPLES_2_KEY) and self.data.getGlobal(self.NEW_SAMPLES_KEY) == 0:
            self.congruentSamples = True
        else:
            self.congruentSamples = False
        if self.verbose:
            print "Congruent SNPs:", self.congruentSNPs
            print "Congruent samples:", self.congruentSamples

    def findEquivalence(self, allowFlip=False):
        """Find equivalence (or otherwise) of input data sets

        Requires no call mismatches, and congruent SNP/sample sets
        May be fully equivalent, or to within flip of major/minor alleles
        """
        if self.mismatchRate==0:
            callsEquivalent = True
        elif allowFlip and self.flipMismatchRate==1:
            # all mismatches are flips of major/minor allele
            callsEquivalent = True
        else:
            callsEquivalent = False
        if self.congruentSamples and self.congruentSnps and callsEquivalent:
            return True
        else:
            return False

    def writeSnpData(self, outPath):
        """Write a breakdown of stats for each SNP

        Plink diff doesn't break down how many non-null calls for each SNP
        So, only report mismatch rates with respect to all calls"""
        head = ["SNP", "MISMATCH", "MISMATCH_RATE", "FLIP", "FLIP_RATE"]
        out = open(outPath, 'w')
        out.write("\t".join(head)+"\n")
        snps = self.data.getSnpList()
        sampleTotal = float(self.data.getGlobal(self.SAMPLES_KEY))
        for snp in snps:
            mismatch = self.data.getSnp(snp, self.MISMATCH_KEY)
            mismatchRate = mismatch / sampleTotal
            flip = self.data.getSnp(snp, self.ALLELE_FLIP_KEY)
            flipRate = flip / sampleTotal
            fields = (snp, mismatch, mismatchRate, flip, flipRate)
            out.write("%s\t%d\t%.4f\t%d\t%.4f\n" % fields)
        out.close()

    def writeSampleData(self, outPath):
        """As with writeSnpData(), but with breakdown by sample"""
        head = ["FAMILY_ID", "INDIVIDUAL_ID", "MISMATCH", "MISMATCH_RATE", 
                "FLIP", "FLIP_RATE"]
        out = open(outPath, 'w')
        out.write("\t".join(head)+"\n")
        samples = self.data.getSampleList()
        snpTotal = float(self.data.getGlobal(self.SNPS_KEY))
        for sample in samples:
            (fid, iid) = sample # (family, individual) ID
            mismatch = self.data.getSample(sample, self.MISMATCH_KEY)
            mismatchRate = mismatch / snpTotal
            flip = self.data.getSample(sample, self.ALLELE_FLIP_KEY)
            flipRate = flip / snpTotal
            fields = (fid, iid, mismatch, mismatchRate, flip, flipRate)
            out.write("%s\t%s\t%d\t%.4f\t%d\t%.4f\n" % fields)
        out.close()

    def writeSummaryJson(self, outPath):
        """Write basic summary stats to file in .json format

        Write derived statistics, and values from raw data"""
        summary = {
            'MISMATCH': self.mismatchRate,
            'CONCORDANCE': 1 - self.mismatchRate,
            'MISMATCH_NON_NULL': self.mismatchNNRate,
            'CONCORDANCE_NON_NULL': 1 - self.mismatchNNRate,
            'FLIP_RATE': self.flipRate,
            'FLIP_MISMATCH_RATE': self.flipMismatchRate,
            'SAMPLES_CONGRUENT': self.congruentSamples,
            'SNPS_CONGRUENT': self.congruentSNPs,
            'EQUIVALENT': self.equivalent,
            'EQUIVALENT_WITHIN_FLIP': self.equivalentWithinFlip
            }
        dataGlobal = self.data.getGlobalDictionary()
        out = open(outPath, 'w')
        out.write(json.dumps([summary, dataGlobal]))
        out.close()

class PlinkDiffWrapper:
    """Class with 'main' method to run Plink diff and write results"""

    def __init__(self):
        self.snpSuffix = '_snps.txt'
        self.sampleSuffix = '_samples.txt'
        self.summaryJsonSuffix = '_summary.json'

    def run(self, stem1, stem2, outStem, cleanup, verbose):        
        data = PlinkDiffParser().run(stem1, stem2, outStem, cleanup, verbose)
        writer = PlinkDiffWriter(data, verbose)
        writer.writeSnpData(outStem+self.snpSuffix)
        writer.writeSampleData(outStem+self.sampleSuffix)
        writer.writeSummaryJson(outStem+self.summaryJsonSuffix)
