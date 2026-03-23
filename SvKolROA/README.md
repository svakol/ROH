# Readme Runs of Autozygosity
This is a program that implements the runs of autozygosity finding method found in the paper Genome-wide autozygosity mapping in human populations (Wang et al. 2009). 

Version: 1.0.0\
Author: Svante Koldenius\
Date: 22-03-2026

## Summary of algorthim 

1) Calculate a LOD-score for a window size as per Paper 2. This LOD-score takes genotype and allele frequency into account, and predicts how likely a window is to be autozygous.

2) If these LOD-scores have a bimodal distribution, the second peak is thought to be the autozygous windows. The threshold for autozygosity is the minima between these two peaks. The program tries to find these peaks automatically; if this is not possible a config file with peaks for each population has to be submitted. If there is no bimodal distribution for a populations, this populations should be skipped. 

3) Concatenate the autozygous windows and cluster them into three size categories. If specified by the -r flag, plot the results. 


## Usage: 
SvKolROA/execute.sh -i [input filepath] -o [output filename template]\
\
Optional Parameters: -p [population config filepath] -t [threshold config filepath] -w [window size] -r [plotting parameters]\
\
Plotting options:\
    p -- plot by population\
    c -- plot by chromosome\
    i -- plot by individual\

### If threshold config file is needed: 
If the threshold config file is needed, simply make one and rerun the program exactly as before but with the additional -t flag. If same -o name is used the program will check which files exist, and not make them if they already exist.

## Input Files
Required:\
A vcf file with snp-data. 
\
\
Potential:\
If the peak finding algorithm cannot identify two clear peaks, or if you disagree with the peaks created by the program, a config file with Population-Peak information. 
\
\
Optional:\
A config file with Sample-Population information, if the vcf-file contains multiple populations. 

## Output Files

1. LOD-scores for each window across all individuals in the input vcf file.
2. Stretches of the genome that the program has evaluated to be autozygous, categorized by size.
3) Density plots for LOD-scores of unconcatenated window segments for each population.
4) If specified: plots of autozygous windows, clustered by size category and population/chromosome/individual.

## Contigencies
For this program to run, certain R packages and Python modules need to be available in the environment.\
\
Python (3.13.5):
numpy (2.3.4)\
matplotlib \
scipy (1.17.1)\
\
R (4.5.2):\
mclust (6.1.2)\
vioplot (0.5.1)\
tidyverse (2.0.0)\

## Population config file 
The population config file should be a tab delimeted file with two columns, the first column should be individual ids identical to those in the vcf file, and the second column their corresponding populations. 


## Theshold config file
The threshold config file should be a tab delimeted file with two to three columns. The first column should be the population id, same as in the population config file. The second column should be the x-coordiate for the first peak in the bimodal distribution, the third column should be the second peak in the bimodal distribution. If the distribution for a certain population is not bimodal, or there is another reason for it to be skipped, the second column should contain the string "skip" (without quotes) and the third column should be left blank. If your input file has only a single population, the file should have a single entry where the first column should have the string "running_one_population_only" (without quotes). 

## Cited Papers

Cited papers:\
&nbsp;&nbsp;&nbsp;&nbsp;Paper 1: Wang, S., Haynes, C., Barany, F., & Ott, J. (2009).\
&nbsp;&nbsp;&nbsp;&nbsp;Genome-wide autozygosity mapping in human populations. Genet. Epidemiol. 33, 172-180\
&nbsp;&nbsp;&nbsp;&nbsp;Paper 2: Pemberton, T. J., Absher, D., Feldman. M. W., Myers, R. M., Rosenberg, N. A., Jun Z. L. (2012).\
&nbsp;&nbsp;&nbsp;&nbsp;Genomic patterns of homozygosity in worldwide human populations. Am. J. Hum. Genet. 91, 275–292\