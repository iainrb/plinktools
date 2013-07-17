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

"""Find heterozygosity separately for SNPs with high/low MAF"""

import argparse, os
from plink import MafHetFinder

def main():
    """Method to run as script from command line"""
    description = "Find heterozygosoity separately for SNPs with high/low minor allele frequency (MAF)."
    outDefaultJson = "het_by_maf.json"
    outDefaultText = "het_by_maf.txt"
    mafThreshold = 0.01
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--in', required=True, metavar="PATH", 
                        help="Plink binary stem (path without .bed, .bim, .fam extension) for input. Required.")
    parser.add_argument('--out', metavar="PATH", default=outDefaultJson,
                        help="Path for output. Optional, defaults to "+outDefaultJson+" in current working directory (unless --text is specified).")
    parser.add_argument('--threshold', metavar="NUM", default=mafThreshold, 
                        type=float, help="Threshold between 0 and 1 for the boundary between high and low minor allele frequency groups. Optional, defaults to "+str(mafThreshold)+".")
    parser.add_argument('--text', action="store_true", 
                        help="Write output in plain text format instead of default .json; changes default output to "+outDefaultText+".")
    parser.add_argument('--verbose', action="store_true", 
                        help="Print additional information to standard output/error.")
    args = vars(parser.parse_args())
    # validate arguments
    if args.has_key('threshold'):
        t = args['threshold']
        if t<0 or t>1:
            raise ValueError("Threshold must be between 0 and 1")
        else:
            mafThreshold = t
    for suffix in ('.bed', '.bim', '.fam'):
        inPath = args['in']+suffix
        if not os.access(inPath, os.R_OK):
            raise OSError("Cannot read input "+inPath)
    outDir = os.path.dirname(os.path.abspath(args['out']))
    if not os.access(outDir, os.W_OK):
         raise OSError("Cannot write to output directory "+outDir)
    # run MAF and heterozygosity computation
    stem = os.path.abspath(args['in'])
    snpTotal = len(open(stem+".bim").readlines())
    mhf = MafHetFinder()
    if args['text']:
        if args['out'] == outDefaultJson: args['out'] = outDefaultText 
        mhf.runText(args['out'], stem+".bed", stem+".fam", snpTotal,
                    mafThreshold, args['verbose'])
    else:
        mhf.runJson(args['out'], stem+".bed", stem+".fam", snpTotal,
                    mafThreshold, args['verbose'])

if __name__ == "__main__":
    main()
