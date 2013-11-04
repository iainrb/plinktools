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

"""Wrapper script for Plink's diff capability

More straightforward interface and additional statistics"""

import argparse
from plink import PlinkValidator, PlinkToolsError
from comparison import PlinkDiffWrapper

def main():
     """Method to run as script from command line"""
     description = "Compute differences between two Plink binary datasets."
     parser = argparse.ArgumentParser(description=description)   
     parser.add_argument('--in', required=True, action='append', 
                         metavar="PATH", help="Prefix for Plink datasets. Must be specified exactly twice. Input datasets must both be in binary format.")
     parser.add_argument('--out', required=False, metavar="PATH", 
                         help="Prefix for Plink output. Optional.")
     parser.add_argument('--uncompressed', required=False, 
                         action='store_true',
                          help="Retain the raw, uncompressed .diff output from the Plink executable. (The uncompressed file may be very large.)")
     parser.add_argument('--verbose', required=False, action='store_true',
                         help="Print additional information to stdout.")
     args = vars(parser.parse_args())
     stems = args['in']
     validator = PlinkValidator()
     if len(stems)!=2:
         raise ValueError("Must specify exactly two Plink input stems")
     for stem in stems:
         (valid, binary) = validator.checkStem(stem)
         if not valid: 
              raise PlinkToolsError("Invalid Plink stem: "+stem)
         elif not binary: 
              raise PlinkToolsError("Non-binary Plink stem: "+stem)
     if args['out']==None: outStem = 'plinktools'
     else: outStem = args['out']
     cleanup = not args['uncompressed'] # clean up the raw output?
     verbose = args['verbose']
     PlinkDiffWrapper().run(stems[0], stems[1], outStem, cleanup, verbose)


if __name__ == "__main__":
    main()
