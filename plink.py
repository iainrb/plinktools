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

class PlinkHandler:

    """Read/write plink data for zcall"""

    def __init__(self, snpTotal):
        """Assume snp-major input; different samples, same SNP set"""
        self.snpTotal = snpTotal

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

    def __init__(self):
        """Overrides parent class -- snp total is not an instance variable"""
        pass

    def bedEquivalent(self, bedPath1, bedPath2, flip, snpTotal, samples, 
                      verbose=True):
        """Check equivalence of calls in .bed files

        ASSUMPTIONS:
        - SNP sets are identical to within allele flips
        - Sample sets are identical
        flip argument is array of 0 for no flip, 1 for flip
        snpTotal, samples arguments are total number of snps, samples"""
        blockSize = self.findBlockSize(samples)
        bed1 = open(bedPath1)
        bed2 = open(bedPath2)
        bed1.seek(3) # skip headers
        bed2.seek(3)
        equiv = True
        remainder = samples % 4
        for i in range(snpTotal):
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
                    pass # flip[i] = 0
                elif a1a==a2b and a1b==a2a:
                    flip[i] = 1
                else:
                    equiv = False
                    if verbose:
                        msg = "Allele values for SNP index "+str(i)+\
                            "incompatible with flip!\n"
                        sys.stderr.write(msg)
        return (equiv, flip, snpTotal)

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

    def compare(self, stem1, stem2, verbose=False):
        """Compare two Plink binary datasets with given prefixes"""
        (famOK, samples) = self.famEquivalent(stem1+".fam", stem2+".fam", 
                                              verbose)
        if not famOK: 
            if verbose: sys.stderr.write("Mismatched .fam files!")
            return False
        (bimOK, flip, snpTotal) = self.bimEquivalent(stem1+".bim", stem2+".bim",
                                                     verbose)
        if not bimOK:
            if verbose: sys.stderr.write("Mismatched .bim files!")
            return False
        bedOK = self.bedEquivalent(stem1+".bed", stem2+".bed", 
                                   flip, snpTotal, samples,
                                   verbose)
        if not bedOK:
            if verbose: sys.stderr.write("Mismatched .bed files!")
            return False
        else:
            return True

    def comparePed(self, pedPath1, pedPath2, flip=True, verbose=True):
        """Compare two .ped files

        If flip==True, assume major/minor alleles in one input are swapped
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
                if calls2!=None: raise ValueError
                else: break
            elif len(calls1)!=len(calls2):
                raise ValueError
            for j in range(len(calls1)):
                c1 = calls1[j]
                c2 = calls2[j]
                if flip: c2 = c2[1]+c2[0]
                if c1!=c2:
                    equiv = False
                    break
            i += 1
            if verbose and i % 10 == 0:
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
            words = re.split('\s+', line)
            name = words[0]
            family = words[1]
            gender = words[4]
            genders[(name, family)] = gender
        return genders

    def parseBimLine(self, bimLine):
        """Get SNP name and allele values from .bim line"""
        words = re.split("\s+", bimLine.strip())
        name = words[1]
        allele1 = words[4]
        allele2 = words[5]
        return (name, allele1, allele2)

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

    Merge for identical SNP sets, disjoint samples
    Omits cross-checking in Plink, but faster, especially for multiple inputs
    """

    def __init__(self):
        self.timeFormat = "%Y-%m-%d_%H:%M:%S"

    def filesIdentical(self, path1, path2):
        cmd = "diff -q %s %s &> /dev/null" % (path1, path2)
        status = os.system(cmd)
        if status==0: return True
        else: return False

    def mergeBedSnpMajor(self, inPaths, samples, outPath, snpTotal, 
                         log=None):
        """Merge two or more Plink .bed files in SNP-major order

        Arguments: paths to .bed files, and total numbers of samples per group
        Inputs contain disjoint samples, identical SNPs"""
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
            if log!=None and ((i+1) % 10000 == 0 or i+1 == snpTotal):
                msg = "Merged samples for probe %s of %s" % (i+1, snpTotal)
                t = time.strftime(self.timeFormat, time.localtime())
                log.write(t+" "+msg+"\n")
                log.flush()
        for i in range(inputTotal):
            if inFiles[i].read() != '':
                msg = "Bytes found after final expected sample in "+inPaths[i]
                raise ValueError(msg)
            else:
                inFiles[i].close()
        outFile.close()

    def merge(self, stems, outPrefix, logPath=None):
        bimPath = stems[0]+".bim"
        snpTotal = len(open(bimPath).readlines())
        if logPath!=None: 
            log = open(logPath, 'w')
            t = time.strftime(self.timeFormat, time.localtime())
            log.write(t+" Started merge.\n")
        else: 
            log = None
        for i in range(1, len(stems)):
            otherBim = stems[i]+".bim"
            if not self.filesIdentical(bimPath, otherBim):
                raise ValueError("Non-identical .bim files!")
        if log!=None:
            t = time.strftime(self.timeFormat, time.localtime())
            log.write(t+" Checked .bim files.\n")
        bedPaths = []
        samples = []
        for stem in stems:
            bedPaths.append(stem+".bed")
            samples.append(len(open(stem+".fam").readlines()))
        self.mergeBedSnpMajor(bedPaths, samples, outPrefix+".bed", 
                              snpTotal, log)
        # Write appropriate .bim, .fam files
        cmd = "cp %s %s" %  (stems[0]+".bim", outPrefix+".bim")
        os.system(cmd)
        famPaths = []
        for stem in stems: famPaths.append(stem+".fam")
        cmd = "cat "+" ".join(famPaths)+" > "+outPrefix+".fam"
        os.system(cmd)
        if log!=None: 
            t = time.strftime(self.timeFormat, time.localtime())
            log.write(t+" Finished.\n")
            log.close()
