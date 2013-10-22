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

import argparse, json
from plink import PlinkEquivalenceTester

def main():
    """Method to run as script from command line"""
    description = "Compare two Plink datasets, in binary or plain-text format."
    parser = argparse.ArgumentParser(description=description)   
    parser.add_argument('--stem', required=True, action='append', 
                        metavar="PATH", help="Prefix for Plink datasets. Must be specified exactly twice. Input datasets must be either both binary, or both non-binary.")
    parser.add_argument('--flip', required=False, action='store_true',
                        help="For non-binary datasets only; indicates major/minor alleles in one input are reversed with respect to the other. Optional.")
    parser.add_argument('--out', required=False, metavar="PATH", 
                        help="Path for .json result output. Optional.")
    parser.add_argument('--silent', required=False, action='store_true',
                        help="Do not print results to standard output. Overrides --verbose.")
    parser.add_argument('--verbose', required=False, action='store_true',
                        help="Run comparison methods in verbose mode.")
    args = vars(parser.parse_args())
    equiv = False
    pet = PlinkEquivalenceTester()
    stems = args['stem']
    if len(stems)!=2:
        raise ValueError("Must specify exactly two Plink input stems")
    binaryType = []
    verbose = False
    if args['silent']==False and args['verbose']==True: verbose = True
    for stem in stems:
        (valid, binary) = pet.checkStem(stem)
        if not valid: raise ValueError("Invalid Plink stem: "+stem)
        binaryType.append(binary)
    if binaryType[0]!=binaryType[1]:
        raise ValueError("Mismatched input formats (binary and non-binary)")
    elif binaryType[0]==True:
        if verbose: print "Found binary datasets; starting comparison."
        equiv = pet.compareBinary(stems[0], stems[1], verbose)
    else:
        if verbose: print "Found non-binary datasets; starting comparison."
        equiv = pet.comparePed(stems[0]+".ped", stems[1]+".ped", args['flip'],
                               verbose)
    if not args['silent']:
        if equiv: print "OK: Plink datasets match."
        else: print "NOT_OK: Plink datasets do not match."
    if args['out']!=None:
        out = open(args['out'], 'w')
        result = [stems[0], stems[1], equiv]
        json.dump(result, out)
        out.close()

if __name__ == "__main__":
    main()
