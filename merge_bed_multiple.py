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

"""Script to merge multiple Plink .bed files"""

import argparse, json, sys
from plink import PlinkMerger

def main():
    """Method to run as script from command line"""
    description = "Merge multiple Plink binary datasets"
    epilog = "Default behaviour is to merge all plink binary datasets in the current working directory. Can change directory using --dir or specify a list of datasets using --stems. Specifying both --dir and --stems is not permitted."
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('--stems', required=False, metavar="PATH", 
                        help="Path to .json file listing Plink stems (paths without .bed, .bim, .fam extension). Optional, see below.")
    parser.add_argument('--dir', required=False, metavar="PATH", default=None,
                        help="Path to directory containing Plink binary datasets to merge. Optional, see below.")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Plink stem (path without .bed, .bim, .fam extension) for merged output. Required.")
    parser.add_argument('--bim', metavar="PATH", 
                        help="Path to .bim file for SNP sort order, used in case with congruent SNPs and disjoint samples. Optional, uses default sort if absent.")
    args = vars(parser.parse_args())
    pm = PlinkMerger()
    if args['stems']!=None:
        if args['dir']!=None:
            raise ValueError("Can specify at most one of --stems, --dir")
        else:
            stemList = json.loads(open(args['stems']).read())
    else:
        if args['dir']==None: inputDir = "."
        else: inputDir = args['dir']
        stemList = pm.findBedStems(inputDir)
    pm.merge(stemList, args['out'], args['bim'])

if __name__ == "__main__":
    main()
