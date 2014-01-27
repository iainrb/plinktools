Plinktools: Python code to process Plink genotyping files
=========================================================

Author: Iain Bancarz, ib5@sanger.ac.uk

Introduction
------------

Plinktools is a Python package for processing Plink format genotyping data. 
It was developed for the following applications not supported by the standard 
Plink executable. Plinktools functions and respective executable scripts are:

1. Fast merge of binary datasets: `merge_bed.py`
2. Heterozygosity calculation by high/low MAF: `het_by_maf.py`
3. Extended diff functionality: `plink_diff.py`

Run any script with --help for more information.

Prerequisites are Python >= 2.7 and Plink 1.0.7. The `plink` executable should 
be visible in the PATH variable. The `plink_diff.py` script parses the log 
text from the Plink executable, and is not guaranteed to work with any 
version of Plink other than 1.0.7.

Fast binary merge
-----------------

The standard Plink merge function performs extended cross-checking and 
recoding of data. This is rather slow, and does not scale well with an 
increasing number of inputs.

To address this issue, Plinktools implements a fast binary merge for two 
important special cases: Congruent SNP sets with disjoint samples, and 
congruent samples with disjoint SNPs. In both cases, input and output is in 
the default SNP-major format. In the case of congruent samples, the fast 
merge strips off headers and concatenates the input .bed files with no need 
for recoding. For congruent SNPs, a fast merge can be done if the number of 
samples in each input is divisible by 4; otherwise recoding is necessary and 
the merge will be slowed, although still somewhat faster than Plink.

Heterozygosity calculation by high/low MAF
------------------------------------------

As part of the Wellcome Trust Oxford SOP for exome chip QC, heterozygosity 
for each sample is calculated separately for SNPs with minor allele frequency 
above and below 1%. This test has been implemented in Plinktools.

Extended diff
-------------

A wrapper for Plink's capability to compare two datasets. It constructs and 
executes the appropriate command, computes some additional statistics, and 
records the results in machine-readable JSON format. It also supports a 
simple equivalence test on two datasets, via the --brief option.

References
----------

The Plink data format was created by Shaun Purcell et al. See: 
http://pngu.mgh.harvard.edu/~purcell/plink

Plinktools was created to support the WTSI genotyping pipeline: 
https://github.com/wtsi-npg/genotyping

Plinktools is a prerequisite for the WTSI extension of the zCall genotype 
caller: https://github.com/wtsi-npg/zCall
