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

import os, pyximport, sys, unittest
pyximport.install()
from hashlib import md5
from tempfile import mkdtemp
from plink import PlinkHandler

class TestPlink(unittest.TestCase):

    def setUp(self):
        self.dataDir = 'data'
        self.outDir = os.path.abspath(mkdtemp(dir=self.dataDir))

    def getMD5hex(self, inPath):
        """Get MD5 checksum for contents of given file, in hex format"""
        m = md5()
        m.update(open(inPath).read())
        checksum = m.hexdigest()
        return checksum

    def test_merge(self):
        bed0 = os.path.join(self.dataDir, 'merge0.bed')
        bed1 = os.path.join(self.dataDir, 'merge1.bed')
        bed2 = os.path.join(self.dataDir, 'merge2.bed')
        prefix = 'merge_test'
        out = os.path.join(self.outDir, 'merge_test.bed')
        snps = 538448
        PlinkHandler(snps).mergeBed((bed0, bed1, bed2), (11, 11, 11), out)
        self.assertTrue(os.path.exists(out))
        md5 = self.getMD5hex(out)
        self.assertEqual(md5, '9bab68c95b13de8200b350e1db00fa04')
        # run plink on output
        startDir = os.getcwd()
        os.chdir(self.dataDir)
        mergedFam = os.path.join(self.outDir, prefix+'.fam')
        os.system('cat merge0.fam merge1.fam merge2.fam > '+mergedFam)
        os.system('cp merge0.bim '+os.path.join(self.outDir, prefix+'.bim'))
        os.chdir(self.outDir)
        self.assertEqual(os.system('plink --bfile '+prefix+' > /dev/null'), 0)
        os.chdir(startDir)

if __name__ == "__main__":
    unittest.main(verbosity=2)
