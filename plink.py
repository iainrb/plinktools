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


import os, re, struct, sys, time
from glob import glob
from checksum import ChecksumFinder
import json # used for temporary troubleshooting in equivalence test

class PlinkHandler:

    """General purpose class to read/write plink data for zcall"""

    def __init__(self):
        pass

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
            raise ValueError("Must have exactly 4 calls for byte conversion!")
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
        elif gap < 0: raise ValueError
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
            else: raise ValueError("Invalid genotype string")
            gtypes.append(gtype)
            i += 2
        return gtypes

    def validateHeader(self, inFile, snpMajor=True):
        """Read and validate Plink .bed header"""
        head = inFile.read(3)
        if ord(head[0])!=108 or ord(head[1])!=27:
            msg = "Header does not start with Plink 'magic number'!"
            raise ValueError(msg)
        if snpMajor:
            if ord(head[2])!=1:
                raise ValueError("Plink .bed file not in SNP-major order.")
        elif ord(head[2])!=0:
            raise ValueError("Plink .bed file not in individual-major order.")
        return head

class PlinkEquivalenceTester(PlinkHandler):
    """Test equivalence of two Plink datasets

    Plink executable can do diff on SNP calls, but this doesn't check 
    congruence of SNP/sample sets, or deal with the flip of major/minor
    alleles introduced by certain plink actions (eg. merge)

    This class can compare binary datasets or non-binary .ped files"""

    def bedEquivalent(self, bedPath1, bedPath2, flip, samples, 
                      verbose=True):
        """Check equivalence of calls in .bed files

        ASSUMPTIONS:
        - SNP sets are identical to within allele flips
        - Sample sets are identical
        flip argument is array of 0 for no flip, 1 for flip on each SNP value
        snpTotal, samples arguments are total number of snps, samples"""
        blockSize = self.findBlockSize(samples)
        bed1 = open(bedPath1)
        bed2 = open(bedPath2)
        bed1.seek(3) # skip headers
        bed2.seek(3)
        equiv = True
        remainder = samples % 4
        for i in range(len(flip)):
            calls1 = self.blockToCalls(bed1.read(blockSize), remainder)
            calls2 = self.blockToCalls(bed2.read(blockSize), remainder)
            for j in range(len(calls1)):
                flipped = flip[j]!=0
                callsOK = self.callsEquivalent(calls1[j], calls2[j], flipped)
                if not callsOK:
                    equiv = False
                    if verbose:
                        msg = "Non-equivalent calls at SNP %s, sample %s\n" % \
                            (i, j)
                        sys.stderr.write(msg)
                        break
        bed1.close()
        bed2.close()
        return equiv

    def bimEquivalent(self, bimPath1, bimPath2, verbose=True):
        """Check SNP equivalence to within flip of major/minor alleles

        Example of flip: AG -> GA"""
        lines1 = open(bimPath1).readlines()
        lines2 = open(bimPath2).readlines()
        snpTotal = len(lines1)
        equiv = True
        flip = [0]*snpTotal
        if snpTotal!=len(lines2):
            equiv = False
            flip = None
            if verbose:
                sys.stderr.write("SNP totals in .bim files do not match!\n")
        else:
            for i in range(snpTotal):
                (name1, a1a, a1b) = self.parseBimLine(lines1[i])
                (name2, a2a, a2b) = self.parseBimLine(lines2[i])
                if name1!=name2:
                    equiv = False
                    if verbose:
                        msg = "Mismatched SNP names at index "+str(i)+"\n"
                        sys.stderr.write(msg)
                elif a1a==a2a and a1b==a2b:
                    pass # flip[i] kept as 0
                elif a1a==a2b and a1b==a2a:
                    flip[i] = 1
                else:
                    equiv = False
                    if verbose:
                        msg = "Allele values for SNP index "+str(i)+\
                            "incompatible with flip!\n"
                        sys.stderr.write(msg)
        return (equiv, flip)

    def callsEquivalent(self, c1, c2, flip):
        """Compare call codes for equivalence

        N = 0, AA = 1, AB = 2, BB = 3
        flip = transposition of major (A) and minor (B) alleles
        """
        equiv = False
        if c1==0:
            if c2 == 0: equiv = True
        elif c1==1:
            if flip and c2 == 3: equiv = True
            elif not flip and c2 == 1: equiv = True
        elif c1==2:
            if c2==2: equiv = True
        elif c1==3:
            if flip and c2 == 1: equiv = True
            elif not flip and c2 == 3: equiv = True
        return equiv

    def compareBinary(self, stem1, stem2, verbose=False):
        """Compare two Plink binary datasets with given prefixes"""
        (famOK, samples) = self.famEquivalent(stem1+".fam", stem2+".fam", 
                                              verbose)
        if not famOK: 
            if verbose: sys.stderr.write("Mismatched .fam files!\n")
            return False
        (bimOK, flip) = self.bimEquivalent(stem1+".bim", stem2+".bim", verbose)
        if not bimOK:
            if verbose: sys.stderr.write("Mismatched .bim files!\n")
            return False
        bedOK = self.bedEquivalent(stem1+".bed", stem2+".bed", 
                                   flip, samples, verbose)
        if not bedOK:
            if verbose: sys.stderr.write("Mismatched .bed files!\n")
            return False
        else:
            return True

    def comparePed(self, pedPath1, pedPath2, flip=True, verbose=False):
        """Compare two .ped files

        If flip==True, assume major/minor alleles in one input are swapped
        (Unlike compareBed, flip status is an 'all-or-nothing' option)
        Note that .ped files are sample-major by default
        """
        in1 = open(pedPath1)
        in2 = open(pedPath2)
        equiv = True
        i = 0
        while True:
            calls1 = self.parsePedLine(in1.readline())
            calls2 = self.parsePedLine(in2.readline())
            if calls1==None:
                if calls2!=None: raise ValueError("Inputs of unequal length!")
                else: break
            elif len(calls1)!=len(calls2):
                raise ValueError("SNP sets of unequal length!")
            for j in range(len(calls1)):
                c1 = calls1[j]
                c2 = calls2[j]
                if flip: c2 = c2[1]+c2[0]
                if c1!=c2:
                    equiv = False
                    if verbose:
                        msg = "Non-equivalent genotypes at sample index "+\
                            str(i)+", SNP index "+str(j)+": "+c1+", "+c2+"\n"
                        sys.stderr.write(msg)
                    break
            if not equiv: break
            i += 1
            if verbose and i % 20 == 0:
                print "Read .ped line "+str(i)
                sys.stdout.flush()
        in1.close()
        in2.close()
        return equiv

    def famEquivalent(self, famPath1, famPath2, verbose=True):
        """Check sample info in .fam files"""
        genders1 = self.readFamGenders(famPath1)
        genders2 = self.readFamGenders(famPath2)
        [keys1, keys2] = [genders1.keys(), genders2.keys()]
        equiv = True
        samples = len(keys1)
        if samples != len(keys2):
            equiv = False
            if verbose: 
                sys.stderr.write(".fam files of unequal length!\n")
        elif set(keys1) != set(keys2):
            equiv = False
            if verbose: 
                sys.stderr.write("Mismatched sample IDs in .fam files!\n")
        else:
            for key in keys1:
                if genders1[key] != genders2[key]:
                    equiv = False
                    if verbose:
                        msg = "Mismatched sample genders in .fam files!\n"
                        sys.stderr.write(msg)
        return (equiv, samples)

    def readFamGenders(self, famPath):
        """Read sample identifiers and genders from .fam file"""
        lines = open(famPath).readlines()
        genders = {}
        for line in lines:
            (name, family, gender) = self.parseFamLine(line)
            genders[(name, family)] = gender
        return genders

    def parseBimLine(self, bimLine):
        """Get SNP name and allele values from .bim line"""
        words = re.split("\s+", bimLine.strip())
        name = words[1]
        allele1 = words[4]
        allele2 = words[5]
        return (name, allele1, allele2)

    def parseFamLine(self, famLine):
        """Get name, family, gender from .fam line"""
        words = re.split('\s+', famLine)
        name = words[0]
        family = words[1]
        gender = words[4]
        return (name, family, gender)

    def parsePedLine(self, line):
        """Parse single line from a .ped file into allele pairs

        Eg. [AG, GC, AA, TT, ...]"""
        if line=='': return None
        words = re.split('\s+', line.strip())
        i = 6 # omit header information
        calls = []
        while i < len(words):
            calls.append(words[i]+words[i+1])
            i += 2
        return calls
            

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
                    raise ValueError("Missing Plink data file "+myFile)
        return plinkList

    def merge(self, stems, outPrefix, bim=None, verbose=False):
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
            raise ValueError(msg)
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
            raise ValueError(msg)
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
                raise ValueError(msg)
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
        cmd = "cp %s %s" %  (stems[0]+".bim", outPrefix+".bim")
        os.system(cmd)
        famPaths = []
        for stem in stems: famPaths.append(stem+".fam")
        cmd = "cat "+" ".join(famPaths)+" > "+outPrefix+".fam"
        os.system(cmd)

    def writeBimFamCongruentSamples(self, sortedStems, outPrefix):
        """Write .bim, .fam files for congruent samples, disjoint SNPs

        Need to ensure stems are in same sort order as merged .bed files"""
        cmd = "cp %s %s" %  (sortedStems[0]+".fam", outPrefix+".fam")
        os.system(cmd)
        bimPaths = []
        for stem in sortedStems: bimPaths.append(stem+".bim")
        cmd = "cat "+" ".join(bimPaths)+" > "+outPrefix+".bim"
        os.system(cmd)

class ChromosomeNameError(Exception):
    """Exception class to identify parse errors for chromosome names"""

    def __init__(self, message, Errors):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)
