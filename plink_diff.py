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

import argparse, sys
from plink import PlinkValidator, PlinkToolsError
from comparison import PlinkDiffWrapper

# want options:
# retain NO detailed output (--brief)
# retain gzipped output (--gzip)
# otherwise, retain uncompressed output

def main():
     """Method to run as script from command line"""
     description = "Wrapper for the diff function of the Plink executable. Compute differences between two Plink binary datasets."
     epilog = "The 'prefix' for input is the path to a Plink binary dataset without filename extensions (.bed, .bim, .fam)."
     parser = argparse.ArgumentParser(description=description, epilog=epilog)
     parser.add_argument('--in1', required=True, 
                         metavar="PATH", help="First input stem. Prefix for a Plink binary dataset.")
     parser.add_argument('--in2', required=True, 
                         metavar="PATH", help="Second input stem. Similar to --in1 above.")
     parser.add_argument('--out', required=False, metavar="PATH", 
                         help="Prefix for output files. Optional.")
     parser.add_argument('--brief', required=False, 
                         action='store_true',
                         help="Report only whether the two datasets differ; delete all output from Plink and do not write any additional files.")
     parser.add_argument('--gzip', required=False, 
                         action='store_true',
                         help="Create a gzip compressed copy of the .diff output from the Plink executable. Only relevant if --brief is not enabled. If neither --brief nor --gzip is enabled, the uncompressed Plink .diff output (which may be quite large) will be retained.")
     parser.add_argument('--verbose', required=False, action='store_true',
                         help="Print additional information to stdout.")
     # Have --in1 and --in2 because Percolate workflows require all command line option names to be distinct; so we can't specify the --in option twice.
     args = vars(parser.parse_args())
     stems = (args['in1'], args['in2'])
     validator = PlinkValidator()
     for stem in stems:
          if stem==None: valid = False
          else: (valid, binary) = validator.checkStem(stem)
          if not valid: 
               raise PlinkToolsError("Invalid Plink stem: "+str(stem))
          elif not binary: 
               raise PlinkToolsError("Non-binary Plink stem: "+str(stem))
     if args['out']==None: outStem = 'plinktools'
     else: outStem = args['out']
     brief = args['brief']
     gzip = args['gzip']
     verbose = args['verbose']
     equivalent = PlinkDiffWrapper().run(stems[0], stems[1], outStem, 
                                         brief, gzip, verbose)
     if not equivalent:
          sys.stderr.write("Plink binary datasets are not equivalent\n")
          exit(1)
     else:
          exit(0)


if __name__ == "__main__":
    main()
