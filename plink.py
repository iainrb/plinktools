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

# Originally based on PlinkHandler class from zCall

"""Process Plink format genotyping data
See http://pngu.mgh.harvard.edu/~purcell/plink/

A PlinkWriter subclass has been written for the zCall genotype caller
See https://github.com/wtsi-npg/zCall
"""


import json, os, re, struct, sys, time
from glob import glob
from checksum import ChecksumFinder

class PlinkHandler:

    """General purpose class to read/write plink data for zcall"""

    def __init__(self):
        self.validator = PlinkValidator()

    def blockToCalls(self, block, remainder):
        """Convert a binary 'block' of data to call codes

        Block represents all sample calls for a given SNP or vice versa
        Remainder is number of samples/snps left at end of block"""
        calls = []
        for i in range(len(block) - 1):
            parsed = self.parseBed(block[i])
            calls.extend(self.readGenotypes(parsed))
        lastCalls = self.readGenotypes(self.parseBed(block[-1]))
        # strip off null padding (if any)
        if remainder == 0:
            calls.extend(lastCalls)
        else:
            calls.extend(lastCalls[0:remainder])
        return calls

    def callsToBinary(self, calls):
        """Translate list of genotype calls to Plink binary

        4 genotype calls are packed into one byte of output
        Assume that null padding (if needed) has already been applied
        Returns a list of struct.pack strings corresponding to output bytes"""
        output = []
        i = 0
        while i < len(calls):
            byte = struct.pack('B', self.callsToByte(calls[i:i+4]))
            output.append(byte)
            i += 4
        return output

    def callsToByte(self, calls):
        """Convert list of 4 calls to an integer in Plink binary format 

        Create byte string of the form '01001101', convert to integer
        See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        """
        if len(calls) != 4:
            msg = "Must have exactly 4 calls for byte conversion!"
            raise PlinkToolsError(msg)
        byte = []
        for call in calls:
            if call==1: bcall = '00' # major homozygote, 'AA'
            elif call==2: bcall = '01' # heterozygote, 'AB'
            elif call==3: bcall = '11' # minor homozygote, 'BB'
            else: bcall = '10' # missing genotype, call=0
            byte.append(bcall)
        byteString = ''.join(byte)
        byteString = byteString[::-1] # reverse order of string characters
        return int(byteString, 2)    

    def findBlockSize(self, sampleTotal):
        """Find number of bytes for each SNP

        Plink encoding has up to 4 calls per byte; final byte may be padded with no-calls."""
        size = sampleTotal / 4 # integer division discards any fractional part
        if sampleTotal % 4 != 0: size += 1
        return size

    def parseBed(self, bed):
        """Parse a single byte from a .bed file"""
        parsed = bin(ord(bed))[2:]
        gap = 8 - len(parsed)
        if gap > 0: parsed = ''.join(['0']*gap)+parsed
        elif gap < 0: raise PlinkToolsError
        # parsed is now a string of the form '01101100'
        return parsed

    def readGenotypes(self, parsedBed):
        """Read a block of 4 genotypes from a parsed .bed string 

        Return in numeric format: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        parsedBed = parsedBed[::-1] # reverse character order
        i = 0
        gtypes = []
        while i < len(parsedBed):
            pair = parsedBed[i:i+2]
            if pair=='00': gtype = 1
            elif pair=='01': gtype = 2
            elif pair=='11': gtype = 3
            elif pair=='10': gtype = 0
            else: raise PlinkToolsError("Invalid genotype string")
            gtypes.append(gtype)
            i += 2
        return gtypes

    def validateHeader(self, inFile, snpMajor=True):
        """Read and validate Plink .bed header"""
        head = inFile.read(3)
        if ord(head[0])!=108 or ord(head[1])!=27:
            msg = "Header does not start with Plink 'magic number'!"
            raise PlinkToolsError(msg)
        if snpMajor:
            if ord(head[2])!=1:
                msg = "Plink .bed file not in SNP-major order."
                raise PlinkToolsError(msg)
        elif ord(head[2])!=0:
            msg = "Plink .bed file not in individual-major order."
            raise PlinkToolsError(msg)
        return head


class PlinkMerger(PlinkHandler):

    """Class to merge Plink datasets

    Merge for:
      - identical SNP sets, disjoint samples; or
      - disjoint SNP sets, identical samples
    Omits cross-checking in Plink, but faster, especially for multiple inputs
    """

    def convertChromosome(self, chrom):
        """Convert chromosome from string to numeric ID

        Assumes the chromosome in question is human!"""
        if re.search('\D+', chrom):
            if chrom=='X': 
                chrom = 23
            elif chrom=='Y': 
                chrom = 24
            elif chrom=='XY': 
                chrom = 25
            elif chrom=='MT': 
                chrom = 26
            else:
                raise ChromosomeNameError("Invalid chromosome name\n")
        else:
            chrom = int(chrom)
            if chrom < 0 or chrom > 26: 
                msg = "Numeric chromosome outside valid range\n"
                raise ChromosomeNameError(msg)
        return chrom

    def filesDisjoint(self, inPaths):
        """Check if files are disjoint, ie. no line occurs more than once

        Loads all inputs into memory -- may be problematic for large datasets"""
        lines = set()
        disjoint = True
        for inPath in inPaths:
            contents = open(inPath).readlines()
            for line in contents:
                if line in lines:
                    disjoint = False
                    break
                else:
                    lines.add(line)
        return disjoint

    def filesIdentical(self, path1, path2):
        """Compare md5 checksums of two files"""
        csf = ChecksumFinder()
        sum1 = csf.getMD5hex(path1)
        sum2 = csf.getMD5hex(path2)
        return sum1==sum2

    def findBedStems(self, prefix="."):
        """Find Plink binary stems (paths without .bed, .bim, .fam extension)"""
        bedList = glob(prefix+"*.bed")
        bedList.sort()
        plinkList = []
        for bed in bedList:
            (stem, ext) = os.path.splitext(bed)
            plinkList.append(stem)
            for suffix in (".bim", ".fam"):
                myFile = stem+suffix
                if not os.path.exists(myFile):
                    raise PlinkToolsError("Missing Plink data file "+myFile)
        return plinkList

    def merge(self, stems, outPrefix, bim=None, validate=True, verbose=False):
        snpTotal = len(open(stems[0]+".bim").readlines())
        if verbose:
            sys.stderr.write(time.asctime()+": Started.\n")
        congruentSNPs = self.pIdentical(stems, ".bim")
        if self.pIdentical(stems, ".bim") and self.pDisjoint(stems, ".fam"):
            if verbose:
                msg = "Found congruent SNPs, disjoint samples;"+\
                    " starting merge.\n"
                sys.stderr.write(msg)
            bedPaths = []
            samples = []
            for stem in stems:
                bedPaths.append(stem+".bed")
                samples.append(len(open(stem+".fam").readlines()))
            self.mergeBedCongruentSNPs(bedPaths, samples, outPrefix+".bed", 
                                       snpTotal, verbose)
            self.writeBimFamCongruentSNPs(stems, outPrefix)
        elif self.pDisjoint(stems, ".bim") and self.pIdentical(stems, ".fam"):
            if verbose:
                msg = "Found disjoint SNPs, congruent samples;"+\
                    " starting merge.\n"
                sys.stderr.write(msg)
            if bim!=None: stems = self.sortStemsInput(stems, bim, verbose)
            else: stems = self.sortStemsIlluminus(stems, verbose)
            self.mergeBedCongruentSamples(stems, outPrefix+".bed", verbose)
            self.writeBimFamCongruentSamples(stems, outPrefix)
        else:
            msg = "Merge conditions not satisfied; must have congruent SNPs"+\
                " and disjoint samples, or disjoint SNPs and congruent samples."
            raise PlinkToolsError(msg)
        if validate:
            ok = self.validator.run(outPrefix)
            if not ok:
                raise PlinkToolsError("Invalid output from Plink merge")
            elif verbose:
                msg = time.asctime()+": Plink output validated successfully.\n"
                sys.stderr.write(msg)
        if verbose: 
            sys.stderr.write(time.asctime()+": Finished.\n")

    def mergeBedCongruentSamples(self, stems, outPath, verbose=False):
        """Merge .bed files in SNP-major order, with identical samples

        Null-padding occurs at end of each SNP, depends on number of samples
        Therefore, inputs do not require decoding
        """
        maxSize = 50*(10**6) # read at most 50 MB
        outFile = open(outPath, 'w')
        outHead = False
        for stem in stems:
            bedFile = open(stem+".bed")
            head = self.validateHeader(bedFile)
            if not outHead: 
                outFile.write(head)
                outHead = True
            bedFile.seek(3) # beginning of data blocks
            while True:
                contents = bedFile.read(maxSize)
                if contents == '': break
                else: outFile.write(contents)
            bedFile.close()
        outFile.close()

    def mergeBedCongruentSNPs(self, inPaths, samples, outPath, snpTotal, 
                              verbose=True):
        """Merge .bed files in SNP-major order, with identical SNPs

        Arguments: paths to .bed files, and total numbers of samples per group
        Assumes that inputs contain disjoint samples, identical SNPs"""
        inputTotal = len(inPaths)
        if inputTotal!=len(samples):
            msg = "Mismatched lengths of input path and samples lists"
            raise PlinkToolsError(msg)
        inFiles = []
        for inPath in inPaths: inFiles.append((open(inPath, 'r')))
        outFile = open(outPath, 'w')
        # check Plink headers
        head = None
        for inFile in inFiles: head = self.validateHeader(inFile)
        outFile.write(head)
        # read byte blocks for each sample, recode (if needed), and output
        remainders = []
        sizes = []
        for total in samples: 
            remainders.append(total % 4)
            sizes.append(self.findBlockSize(total))
        for i in range(snpTotal): # for each snp
            recoded = False
            calls = []
            for j in range(inputTotal): # for each sample
                block = inFiles[j].read(sizes[j])
                if recoded==False and remainders[j]==0:
                    # special case -- can write input unchanged
                    outFile.write(block)
                else:
                    # must recode this and all subsequent blocks for the SNP
                    recoded = True
                    calls.extend(self.blockToCalls(block, remainders[j]))
            if recoded:
                if len(calls) % 4 != 0:
                    # null-pad output to an integer number of bytes
                    calls.extend([0]*(4 - (len(calls) % 4))) 
                outFile.write(''.join(self.callsToBinary(calls)))
            if verbose and ((i+1) % 10000 == 0 or i+1 == snpTotal):
                msg = "Merged samples for probe %s of %s\n" % (i+1, snpTotal)
                sys.stderr.write(msg)
        for i in range(inputTotal):
            if inFiles[i].read() != '':
                msg = "Bytes found after final expected sample in "+inPaths[i]
                raise PlinkToolsError(msg)
            else:
                inFiles[i].close()
        outFile.close()

    def pDisjoint(self, stems, suffix):
        """Check if plink files with given suffix are all disjoint"""
        plinkFiles = []
        for stem in stems: plinkFiles.append(stem+suffix)
        return self.filesDisjoint(plinkFiles)

    def pIdentical(self, stems, suffix):
        """Check if plink files with given suffix are all identical"""
        ident = True
        firstPlink = stems[0]+suffix
        for i in range(1, len(stems)):
            otherPlink = stems[i]+suffix
            if not self.filesIdentical(firstPlink, otherPlink):
                ident = False
                break
        return ident

    def sortStemsInput(self, stems, bimPath, verbose=False):
        """Sort a list of Plink file stems into same order as given .bim file

        Use to reproduce Plink's default sort order for testing"""
        snpOrder = {}
        lines = open(bimPath).readlines()
        for i in range(len(lines)):
            snp = re.split('\s+', lines[i])[1]
            snpOrder[snp] = i
        sortMap = {}
        for stem in stems:
            # order stems by position of initial SNP
            line = open(stem+".bim").readline()
            snp = re.split('\s+', line)[1]
            pos = snpOrder[snp]
            sortMap[pos] = stem
        sortKeys = sortMap.keys()
        sortKeys.sort()
        sortedStems = []
        for key in sortKeys:
            sortedStems.append(sortMap[key])
        if verbose:
            print "Inputs sorted in .bim order:"
            for stem in sortedStems: print stem
        return sortedStems

    def sortStemsIlluminus(self, stems, verbose=False):
        """Sort a list of Plink file stems into 'Illuminus' order

        Illuminus output .bim paths do not necessarily include SNP coordinates
        Instead, sort by filename
        Names of the form [prefix].part.[number].[chromosome].bim
        Chromosome is from (0-22, X, Y, XY, MT)

        If names are not in Illuminus form, return in lexical sort order
        """
        sortable = True
        sortMap = {}
        for stem in stems:
            terms = re.split('\.', stem)
            if len(terms)<3:
                if verbose:
                    sys.stderr.write("Insufficient terms for Illuminus sort\n")
                sortable = False; break
            else:
                part = terms[-2]
                chrom = terms[-1]
                if re.search('\D+', part): # contains non-numeric characters
                    if verbose:
                        msg = "Non-numeric characters in 'part' term: \""+\
                            part+"\"\n"
                        sys.stderr.write(msg)
                    sortable = False; break
                else:
                    part = int(part)
                try:
                    chrom = self.convertChromosome(chrom)
                except ChromosomeNameError:
                    if verbose:
                        sys.stderr.write("Cannot sort chromosomes\n")
                    sortable = False; break
                sortMap[(chrom, part)] = stem
        if sortable:
            keyList = sortMap.keys()
            keyList.sort()
            sortedStems = []
            for key in keyList: sortedStems.append(sortMap[key])
        else:
            if verbose:
                sys.stderr.write("Illuminus sort failed, using default.\n")
            sortedStems = stems
            sortedStems.sort()
        return sortedStems

    def writeBimFamCongruentSNPs(self, stems, outPrefix):
        """Write .bim, .fam files for congruent SNPs, disjoint samples"""
        status = os.system("cp %s %s" %  (stems[0]+".bim", outPrefix+".bim"))
        if status!=0: raise PlinkToolsError("Failed to write .bim file!")
        famPaths = []
        for stem in stems: famPaths.append(stem+".fam")
        status = os.system("cat "+" ".join(famPaths)+" > "+outPrefix+".fam")
        if status!=0: raise PlinkToolsError("Failed to write .fam file!")

    def writeBimFamCongruentSamples(self, sortedStems, outPrefix):
        """Write .bim, .fam files for congruent samples, disjoint SNPs

        Need to ensure stems are in same sort order as merged .bed files"""
        cmd = "cp %s %s" %  (sortedStems[0]+".fam", outPrefix+".fam")
        status = os.system(cmd)
        if status!=0: raise PlinkToolsError("Failed to write .fam file!")
        bimPaths = []
        for stem in sortedStems: bimPaths.append(stem+".bim")
        status = os.system("cat "+" ".join(bimPaths)+" > "+outPrefix+".bim")
        if status!=0: raise PlinkToolsError("Failed to write .bim file!")

class MafHetFinder(PlinkHandler):
    """Class to find autosome heterozygosity, split by MAF

    Assumes SNP-major input
    Use to apply het filter separately to high/low MAF, eg. for exome chips"""

    def mafSplitHetCounts(self, bedPath, samples, snpTotal, mafThreshold=0.01,
                          verbose=False):
        """Find het rates on SNPs above/below MAF threshold for each sample

        Populate an array: [low total, low het, high total, high het]"""
        counts = [None]*samples
        for i in range(samples): counts[i] = [0]*2
        highTotal = 0
        lowTotal = 0
        blockSize = self.findBlockSize(samples)
        remainder = samples % 4
        bed = open(bedPath)
        bed.seek(3) # skip headers
        for i in range(snpTotal):
            calls = self.blockToCalls(bed.read(blockSize), remainder)
            maf = self.findMAF(calls) # MAF for this SNP
            if maf <= mafThreshold: lowTotal += 1
            else: highTotal += 1
            for j in range(len(calls)):
                if calls[j]==2: # heterozygote
                    if maf <= mafThreshold: counts[j][0] += 1
                    else: counts[j][1] += 1
            if verbose and (i+1) % 10000 == 0:
                print "Found MAF for SNP", i+1, "of", snpTotal
        bed.close()
        return (lowTotal, highTotal, counts)

    def findMAF(self, calls):
        """Find minor allele frequency from a given set of calls"""
        totalCalls = 0 # total non-null calls for this SNP
        minorCalls = 0 # total 'B' allele calls
        for call in calls:
            if call==0: continue
            totalCalls += 2 # count both alleles towards total
            if call==2: minorCalls+=1 # minor het
            elif call==3: minorCalls+=2 # minor hom
        try: maf = minorCalls/float(totalCalls)
        except ZeroDivisionError: maf = 0
        if maf > 0.5: maf = 1 - maf # by definition, MAF is less frequent
        return maf

    def readSampleNames(self, famPath):
        """Read sample individual names from Plink .fam file

        First field in .fam line is family name, second is individual name"""
        famLines = open(famPath).readlines()
        names = []
        for line in famLines:
            names.append(re.split('\s+', line.strip()).pop(1))
        return names

    def runJson(self, outPath, bedPath, famPath, 
                snpTotal, mafThreshold=0.01, verbose=False, digits=6):
        """Run MAF/het calculation and write JSON output """
        sampleNames = self.readSampleNames(famPath)
        countInfo = self.mafSplitHetCounts(bedPath, len(sampleNames), 
                                           snpTotal, mafThreshold, verbose)
        self.writeJson(outPath, countInfo, sampleNames, snpTotal, mafThreshold, 
                       verbose)

    def runText(self, outPath, bedPath, famPath, 
                snpTotal, mafThreshold=0.01, verbose=False, digits=6):
        """Run MAF/het calculation and write plain text output """
        samples = len(self.readSampleNames(famPath))
        countInfo = self.mafSplitHetCounts(bedPath, samples, 
                                           snpTotal, mafThreshold, verbose)
        self.writeText(outPath, countInfo, samples, snpTotal, mafThreshold, 
                       verbose)

    def writeJson(self, outPath, countInfo, sampleNames, 
                  snpTotal, mafThreshold=0.01, verbose=False, digits=6):
        """Find heterozygosity for high/low MAF and write to .json

        .json format is similar to qc_results.json in genotyping PL"""
        samples = len(sampleNames)
        (lowTotal, highTotal, counts) = countInfo
        lowTotal = float(lowTotal)
        highTotal = float(highTotal)
        output = {}
        for i in range(samples):
            try: lmh = counts[i][0]/lowTotal # low MAF het
            except ZeroDivisionError: lmh = 0
            try: hmh = counts[i][1]/highTotal # high MAF het
            except ZeroDivisionError: hmh = 0
            result = { 'low_maf_het':[1, round(lmh, digits)], 
                       'high_maf_het':[1, round(hmh, digits)]}
            output[sampleNames[i]] = result
        out = open(outPath, 'w')
        json.dump(output, out)
        out.close()
        
    def writeText(self, outPath, countInfo, samples, snpTotal, 
                  mafThreshold=0.01, verbose=False, digits=6):
        """Find heterozygosity for high/low MAF and write to text"""
        out = open(outPath, 'w')
        out.write("#MAF_THRESHOLD:%s\n" % (mafThreshold,))
        (lowTotal, highTotal, counts) = countInfo
        out.write("#LOW_MAF_TOTAL:%s\n" % (lowTotal,))
        out.write("#HIGH_MAF_TOTAL:%s\n" % (highTotal,))
        out.write("#LOW_MAF_HET\tHIGH_MAF_HET\tALL_HET\n")
        lowTotal = float(lowTotal)
        highTotal = float(highTotal)
        snpTotal = float(snpTotal)
        for i in range(samples):
            try: lmh = counts[i][0]/lowTotal # low MAF het
            except ZeroDivisionError: lmh = 0
            try: hmh = counts[i][1]/highTotal # high MAF het
            except ZeroDivisionError: hmh = 0
            allHet = (counts[i][0] + counts[i][1])/snpTotal
            out.write("%s\t%s\t%s\n" % (round(lmh, digits), 
                                        round(hmh, digits),
                                        round(allHet, digits)))
        out.close()



class PlinkToolsError(Exception):
    """Exception class for generic Plinktools errors"""
    pass

class ChromosomeNameError(PlinkToolsError):
    """Exception class to identify parse errors for chromosome names"""
    pass

class PlinkValidator:
    """Simple class to run Plink executable and validate data"""
    def __init__(self):
        if not self.plinkAvailable():
            raise PlinkToolsError("Plink executable not found")

    def checkStem(self, stem):
        """Check if given file stem is a valid binary or non-binary datset

        If both types are present, binary takes precedence
        Only checks for existence of files, not readability or correctness"""
        binarySuffixes = ['.bed', '.bim', '.fam']
        nonBinarySuffixes = ['.ped', '.map']
        valid = True
        binary = True
        for suffix in binarySuffixes: 
            if not os.path.exists(stem+suffix):
                binary = False
                break
        if not binary:
            for suffix in nonBinarySuffixes:
                if not os.path.exists(stem+suffix):
                    valid = False # neither binary nor non-binary exists
                    break
        return (valid, binary)

    def plinkAvailable(self):
        """Check that Plink executable is available on PATH"""
        status = os.system('which plink > /dev/null')
        if (status!=0): return False
        else: return True

    def run(self, stem, binary=True, cleanup=True):
        """Check if Plink runs with non-zero exit status"""
        if binary: cmd = 'plink --bfile '+stem+' > /dev/null'
        else: cmd = 'plink '+stem+' > /dev/null'
        status = os.system(cmd)
        if cleanup:
            plinkFiles = ['.pversion', 'plink.hh', 'plink.log', 'plink.nof', 
                          'plink.nosex', 'plink.nof']
            for pFile in plinkFiles: 
                if os.path.exists(pFile): os.remove(pFile)
        if (status!=0): return False
        else: return True
