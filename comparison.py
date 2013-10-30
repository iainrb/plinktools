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
Required comparison functions:
- Are datasets equivalent?
- Are datasets equivalent to within major/minor allele flips?
- Are SNP sets congruent? If not, what is intersection?
- Are sample sets congruent? If not, what is intersection?
- Diff stats:
  - Shared SNPs and samples
  - Concordance on shared pairs, all calls
  - Concordance on shared pairs, calls where both are non-null
  - Breakdown by SNP and sample 


For diff stats:
A = Total shared pairs
B = Total shared pairs where both are non-null
P = Differing pairs, including null calls
Q = Differing pairs, excluding null calls
R = Differing pairs, excluding null calls and allele flips (??)

Concordance values:
1 - A/P
1 - B/Q
1 - B/R

P,Q,R from Plink diff output
A can be found from .bim, .fam files
B requires horrible munging of Plink log text, still better than re-implementing the Plink diff

Parse Plink diff for stats
Separate checks for equivalence of SNPs, samples, calls

Classes: PlinkDiffParser (run plink diff and parse output), PlinkComparator (main class with comparison functions), PlinkDiffResult (to read/write output)

Idea for compact Plink diff format: Header lists shared samples & SNPs
Let (i,j) be shared pair index
List differences in the form (i, j, old call, new call)
Or in JSON: sample -> snp -> (old, new)

Back-end: Horrible munging of Plink diff output text
Potentially to be replaced by parsing .bed files with Plinktools classes

Front-end: Data structure(s) of raw results, use to calculate summary stats

"""

import gzip, os, re

from plink import PlinkHandler, PlinkValidator, PlinkToolsError


class PlinkDiffShared:
    """Define some shared constants"""

    # Variables for entire dataset, and for sample/SNP breakdowns
    COUNT_KEYS = ['MISMATCH', 'MISMATCH_NON_NULL', 'ALLELE_FLIP']
    MISMATCH_KEY = COUNT_KEYS[0]
    MISMATCH_NON_NULL_KEY = COUNT_KEYS[1]
    ALLELE_FLIP_KEY = COUNT_KEYS[2]

    # Variables for the dataset as a whole
    TOTAL_KEYS = ['CALLS', 'CALLS_NON_NULL', 'SAMPLES', 'SNPS', 'CONCORDANT']
    CALLS_KEY = TOTAL_KEYS[0]
    CALLS_NON_NULL_KEY = TOTAL_KEYS[1]
    SAMPLES_KEY = TOTAL_KEYS[2]
    SNPS_KEY = TOTAL_KEYS[3]
    CONCORDANT_KEY = TOTAL_KEYS[4]


class PlinkDiffParser(PlinkDiffShared):

    def __init__(self):
        self.data = PlinkDiffData()

    def getData(self):
        return self.data

    def parseDiffFile(self, stem, removeInput=False, verbose=True):
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
        elif verbose:
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
        if verbose: print "Results read from Plink .diff output."
        if removeInput:
            os.remove(diffPath)
            if verbose: print "Removed original .diff file."

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

    def parseLogFile(self, stem, verbose=False):
        """Munge the .log text from Plink to find total number of shared 
        SNPs, and number which were both genotyped. Unfortunately this is 
        the only way to find the latter statistic from Plink diff output."""
        logPath = stem+'.log'
        if not (os.path.exists(logPath)):
            raise PlinkToolsError("Plink .log output file not found!")
        elif verbose:
            print "Plink .log output file found."
        lines = open(logPath).readlines()
        patterns = [
            'Of \d+ overlapping SNPs, \d+ were both genotyped',
            'and \d+ were concordant',
            # shared SNPs
            'Of these, \d+ are new, \d+ already exist in current data', 
            # shared samples
            'Of these, \d+ were new, \d+ were already in current data'
            ] 
        exprs = []
        for p in patterns: exprs.append(re.compile(p))
        # find: intersection, bothGenotyped, concordant, snps, samples
        for i in range(len(lines)):
            line = lines[i]
            words = re.split('\s+', line.strip())
            if exprs[0].match(line):
                self.data.setGlobal(self.CALLS_KEY, int(words[1]))
                self.data.setGlobal(self.CALLS_NON_NULL_KEY, int(words[4]))
            elif exprs[1].match(line):
                self.data.setGlobal(self.CONCORDANT_KEY, int(words[1]))
            elif exprs[2].match(line) and re.search('markers', lines[i-1]):
                self.data.setGlobal(self.SNPS_KEY, int(words[5]))
            elif exprs[3].match(line) and re.search('individuals',lines[i-1]):
                self.data.setGlobal(self.SAMPLES_KEY, int(words[5]))

    def run(self, stem1, stem2, outStem, verbose=False):
        """Main method to run Plink diff and parse the output """
        if not (os.path.exists(outStem+".log") 
                and os.path.exists(outStem+".diff")):
            self.runBinaryDiff(stem1, stem2, outStem, verbose)
        self.parseDiffFile(outStem, verbose=verbose)
        self.parseLogFile(outStem, verbose)
        print self.data.globalTotals
        return self.data

    def runBinaryDiff(self, stem1, stem2, outStem, verbose=False):
        """Run Plink executable with appropriate options to find diff 

        Inputs: Stems for two Plink binary input datasets, and for output
        Stem is Plink file path without .bed, .bim, .fam extension"""
        bmerge = "%s.bed %s.bim %s.fam" % (stem2, stem2, stem2)
        cmd = "plink --noweb --bfile "+stem1+" --bmerge "+bmerge+\
            " --make-bed --merge-mode 6 --out "+outStem
        if not verbose: cmd = cmd+" > /dev/null"
        status = os.system(cmd)
        if status!=0:
            msg = "Non-zero exit status from Plink executable!"
            raise PlinkToolsError(msg)

class PlinkDiffData(PlinkDiffShared):

    """Container for simple result counts:
    - Shared (snp, sample) pairs
    - Non-null shared pairs
    - Concordant pairs
    - Concordant non-null pairs
    - Concordant pairs within flip

    Grand totals, and breakdowns by sample/SNP
    """

    def __init__(self):
        self.globalTotals = {}
        self.snpTotals = {}
        self.sampleTotals = {}
        self.totalShared = 0
        self.totalSharedNonNull = 0       
        for key in self.COUNT_KEYS: self.globalTotals[key] = 0
        for key in self.TOTAL_KEYS: self.globalTotals[key] = 0

    def getGlobal(self, key):
        return self.globalTotals[key]

    def setGlobal(self, key, value):
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
        return self.snpTotals[key][snp]

    def setSnp(self, snp, key, value):
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

class PlinkDiffWriter(PlinkDiffShared):
    """Take diff data, calculate additional stats and write output"""

    def __init__(self, data):
        self.data = data
        self.findConcordance()

    def findConcordance(self):
        self.mismatch = self.data.getGlobal(self.MISMATCH_KEY) / float(self.data.getGlobal(self.CALLS_KEY))
        self.mismatchNonNull = self.data.getGlobal(self.MISMATCH_NON_NULL_KEY) / float(self.data.getGlobal(self.CALLS_NON_NULL_KEY))
        self.flipRate = self.data.getGlobal(self.ALLELE_FLIP_KEY) / float(self.data.getGlobal(self.CALLS_NON_NULL_KEY))
        print "Concordance:", 1 - self.mismatch
        print "Concordance on non-null calls:", 1 - self.mismatchNonNull
        print "Allele flip rate:", self.flipRate



def main():
    import sys
    stem1 = sys.argv[1]
    stem2 = sys.argv[2]
    out = sys.argv[3]
    data = PlinkDiffParser().run(stem1, stem2, out, True)    
    PlinkDiffWriter(data)

if __name__ == "__main__":
    main()
