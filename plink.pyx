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


"""Process Plink format genotyping data

See http://pngu.mgh.harvard.edu/~purcell/plink/
"""

import json, struct

class PlinkHandler:

    def __init__(self, snpTotal):
        """Assume snp-major input; different samples, same SNP set"""
        self.snpTotal = snpTotal

        #self.sampleTotal = sampleTotal
        #self.remainder = sampleTotal % 4
        #self.blockSize = self.findBlockSize(self.sampleTotal)
        #print "REMAINDER", self.remainder
        #print "BLOCK_SIZE", self.blockSize

    def callsToBinary(self, calls):
        """Translate list of genotype calls to Plink binary

        4 genotype calls are packed into one byte of output
        Assume that null padding (if needed) has already been applied
        Returns a list of struct.pack strings corresponding to output bytes"""
        output = []
        i = 0
        while i < len(calls):
            byte = struct.pack('B', self.callsToByte(calls[i:i+4]))
            #print "\t", calls[i:i+4], ord(byte)
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

    def mergePair(self, inPaths, samples, outPath):
        """Merge two Plink .bed files in SNP-major order"""
        inFiles = (open(inPaths[0], 'r'), open(inPaths[1], 'r'))
        outFile = open(outPath, 'w')
        # check Plink headers
        head = None
        for inFile in inFiles:
            head = inFile.read(3)
            if ord(head[0])!=108 or ord(head[1])!=27:
                msg = "Header does not start with Plink 'magic number'!"
                raise ValueError(msg)
            elif ord(head[2])!=1:
                raise ValueError("Plink .bed files must be in SNP-major order.")
        outFile.write(head)
        # read byte blocks for each sample
        remainders = (samples[0] % 4, samples[1] % 4)
        blockSize = (self.findBlockSize(samples[0]), 
                     self.findBlockSize(samples[1]))
        for i in range(self.snpTotal):
            blocks = (inFiles[0].read(blockSize[0]),
                      inFiles[1].read(blockSize[1]))
            if remainders[0]==0 and remainders[1]==0:
                # write both blocks without change
                for j in (0,1): outFile.write(blocks[j])
            else:
                calls = []
                for j in (0,1):
                    for k in range(blockSize[j] - 1):
                        parsed = self.parseBed(blocks[j][k])
                        calls.extend(self.readGenotypes(parsed))
                    lastCalls = self.readGenotypes(self.parseBed(blocks[j][-1]))
                    if remainders[j]!=0: 
                        calls.extend(lastCalls[0:remainders[j]])
                    else:
                        calls.extend(lastCalls)
                calls.extend([0]*(len(calls) % 4)) # null padding, if needed
                #print i, calls, len(calls), self.callsToBinary(calls)
                outFile.write(''.join(self.callsToBinary(calls)))
        for i in (0,1):
            result = inFiles[i].read()
            if result!='':
                msg = "Error: "+len(result)+" bytes found after final "+\
                    "expected sample in "+inPaths[i]
                raise ValueError(msg)
        for myFile in (inFiles[0], inFiles[1], outFile): myFile.close()

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
