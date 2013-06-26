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

import os, sys, unittest
from hashlib import md5
from tempfile import mkdtemp
from plink import PlinkEquivalenceTester, PlinkMerger

class TestPlink(unittest.TestCase):

    def setUp(self):
        self.dataDir = '/nfs/gapi/data/genotype/plinktools_test'
        self.outDir = os.path.abspath(mkdtemp(dir='.'))

    def getMD5hex(self, inPath):
        """Get MD5 checksum for contents of given file, in hex format"""
        m = md5()
        m.update(open(inPath).read())
        checksum = m.hexdigest()
        return checksum

    def test_equivalence(self):
        """Test equivalence check on pairs of datasets"""
        stem1 = os.path.join(self.dataDir, 'samples_part_100_small')
        stem2 = os.path.join(self.dataDir, 'samples_part_100_small')
        pet = PlinkEquivalenceTester()
        match = pet.compare(stem1, stem2)
        self.assertTrue(match)
        stem3 = os.path.join(self.dataDir, 'samples_part_101_small')
        match = pet.compare(stem1, stem3)
        self.assertFalse(match)

    def test_merge(self):
        """Merge two .bed datasets with 100 and 81 samples respectively

        Size of second dataset is not divisible by 4 -- tests recoding"""
        stem0 = os.path.join(self.dataDir, 'samples_part_000.smajor')
        stem1 = os.path.join(self.dataDir, 'samples_part_028.smajor')
        prefix = 'merge_test'
        out = os.path.join(self.outDir, 'merge_test')
        PlinkMerger().merge((stem0, stem1), out)
        self.assertTrue(os.path.exists(out+".bed"))
        md5 = self.getMD5hex(out+".bed")
        self.assertEqual(md5, 'db7dc2a92d339817d6d18974b0add3c7')
        # run plink on output
        startDir = os.getcwd()
        os.chdir(self.dataDir)
        os.chdir(self.outDir)
        self.assertEqual(os.system('plink --bfile '+prefix+' > /dev/null'), 0)
        os.chdir(startDir)

if __name__ == "__main__":
    unittest.main(verbosity=2)
