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
from plink import PlinkEquivalenceTester

stem1 = sys.argv[1]
stem2 = sys.argv[2]

snps1 = len(open(stem1+".bim").readlines())
snps2 = len(open(stem2+".bim").readlines())
if snps1 != snps2:
    raise ValueError("Mismatched SNP set lengths!")

pet = PlinkEquivalenceTester(snps1)

if pet.famEquivalent(stem1+".fam", stem2+".fam"): print ".fam files match."
else: raise ValueError
samples = len(open(stem1+".fam").readlines())

(equiv, flip) = pet.bimEquivalent(stem1+".bim", stem2+".bim")
if equiv: print ".bim files match."
else: raise ValueError
print "%s of %s allele pairs flipped." % (sum(flip), len(flip))

ok = pet.bedEquivalent(stem1+".bed", stem2+".bed", flip, snps1, samples)

if ok: print "OK: .bed files match."
else: print "NOT_OK: .bed files do not match!"
