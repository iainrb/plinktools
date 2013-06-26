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

"""Simple front-end script to compare plink (non-binary) .ped files

Assumes that all alleles in one file are flipped with respect to the other"""

import sys
from plink import PlinkEquivalenceTester

if len(sys.argv)!=3:
    print "Usage: "+sys.argv[0]+" [.ped path 1] [.ped path 2]"
    sys.exit(0)

pet = PlinkEquivalenceTester()
ok = pet.comparePed(sys.argv[1], sys.argv[2], flip=True, verbose=True)
if ok: print "OK: .ped files match"
else: print "NOT_OK: .ped files do not match"
