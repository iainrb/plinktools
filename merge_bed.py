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

"""Script to merge two or more Plink binary datasets"""

import argparse, json, sys
from plink import PlinkMerger

def main():
    """Method to run as script from command line"""
    description = "Merge two or more Plink binary datasets. Inputs must have either disjoint samples and congruent SNPs, or congruent SNPs and disjoint samples."
    epilog = "Default behaviour is to merge all plink binary datasets in the current working directory. Can change directory and file prefix using --prefix or specify a list of datasets using --stems. Specifying both --prefix and --stems is not permitted."
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('--list', required=False, metavar="PATH", 
                        help="Path to .json file listing Plink stems (paths without .bed, .bim, .fam extension). Optional, see below.")
    parser.add_argument('--prefix', required=False, metavar="PATH", 
                        default=None,
                        help="Prefix for Plink binary datasets. Merge will glob all datasets of the form [prefix]*.bed.  Optional, see below.")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Plink stem (path without .bed, .bim, .fam extension) for merged output. Required.")
    parser.add_argument('--bim', metavar="PATH", 
                        help="Path to .bim file for SNP sort order, used in case with congruent SNPs and disjoint samples. Optional, uses default sort if absent.")
    parser.add_argument('--no-validate',  action="store_true", 
                        help="Do not attempt to run PLINK to validate final output. Allows merging without use of the Plink executable.")
    parser.add_argument('--verbose', action="store_true", 
                        help="Print additional information to standard output/error.")
    args = vars(parser.parse_args())
    pm = PlinkMerger()
    if args['list']!=None:
        if args['prefix']!=None:
            raise ValueError("Can specify at most one of --list, --prefix")
        else:
            stemList = json.loads(open(args['list']).read())
    else:
        if args['prefix']==None: inputPrefix = "."
        else: inputPrefix = args['prefix']
        stemList = pm.findBedStems(inputPrefix)
    if args['no-validate']==True: validate = False
    else: validate = True
    pm.merge(stemList, args['out'], args['bim'], validate, args['verbose'])

if __name__ == "__main__":
    main()
