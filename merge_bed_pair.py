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

"""Script to merge two Plink .bed files"""


import os, sys
from plink import PlinkMerger

stem1 = sys.argv[1]
stem2 = sys.argv[2]
out = sys.argv[3]

snps1 = len(open(stem1+".bim").readlines())
snps2 = len(open(stem2+".bim").readlines())
if snps1 != snps2:
    raise ValueError("Mismatched SNP set lengths!")
samples1 = len(open(stem1+".fam").readlines())
samples2 = len(open(stem2+".fam").readlines())

pm = PlinkMerger()
pm.merge((stem1, stem2), out)

