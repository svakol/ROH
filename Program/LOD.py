#!/usr/bin/env python3
"""
LOD.py

Description: This script reads a vcf-file, calculates a lod score (probability of
    autozygosity) for each segment along the genome. It is part of the program SvKolROH.
    

Procedure: Read the config files. For each population, run through the vcf-file and calculate
    (using resampling as in paper 1) allele frequencies for the population, then LOD-scores for each individual,
    then write the result to the output file. 


User defined functions:
    calculate_lod(): Calculates a LOD-score as per paper 2 for a genotype and allele frequency. 

Cited papers: 
    Paper 1: Wang, S., Haynes, C., Barany, F., & Ott, J. (2009). 
    Genome-wide autozygosity mapping in human populations. Genet. Epidemiol. 33, 172-180
    Paper 2: Pemberton, T. J., Absher, D., Feldman. M. W., Myers, R. M., Rosenberg, N. A., Jun Z. L. (2012). 
    Genomic patterns of homozygosity in worldwide human populations. Am. J. Hum. Genet. 91, 275–292

Version: 1.0.0
Author: Svante Koldenius
Date: 22-03-2026
"""

# Import necessary modules
import sys
import random
import re
import math

# Take input variables. pop_config stores which population each individual belongs to. 
input_file = sys.argv[1]
output_file = sys.argv[2]
window_size = int(sys.argv[3])
if len(sys.argv) > 4:
    pop_config = sys.argv[4]

    #current_seed = sys.argv[4]

# Set seed
#random.seed(current_seed)

# Define Error Rate
e = 0.001
'''
Define a function that calculates the lod score, taking genotype and allele frequencies as input. 
pK1 is the probability that the genotype is autozygous, and pK0 is the probability that it is not.  
'''
def calculate_lod(genotype,allele_freq):    
    
    # If any genotype is missing, p of both autozygozity and not autozygosity is set to zero. 
    if genotype[0] == "." or genotype[1] == ".":
        pK1 = 1
        pK0 = 1
        
    # If heterozygous, calculate probabilities according to paper 2.
    elif genotype[0] == genotype[1]:
        pK1 = (1-e)*allele_freq[genotype[0]] + e*allele_freq[genotype[0]]*allele_freq[genotype[0]]
        pK0 = allele_freq[genotype[0]]*allele_freq[genotype[0]]
    
    # If not (homozygous), calculate probablities according to paper 2. 
    else:
        pK1 = 2*e*allele_freq[genotype[0]]*allele_freq[genotype[1]]
        pK0 = 2*allele_freq[genotype[0]]*allele_freq[genotype[1]]
    
    # Return a lod score based on these probabilites, calculated according to paper 2. 
    return(math.log((pK1/pK0),10))


'''
Define empty variables.
'''
lod_list = []
lod_scores = []
pop_dict = dict()
individuals_ids = []
individuals_dict = dict()

'''
Processing of the population config file. 
Create a dictionary with the populations as keys, and a list of individuals as the answers.
If there is no population config file, create one with one population for all individuals, 
as to not have to duplicate code later. 
'''
if pop_config:
    with open(pop_config) as f:
        for line in f:
            if line.split()[1] in pop_dict:
                pop_dict[line.split()[1]].append(line.split()[0])
            else:
                pop_dict[line.split()[1]] = [line.split()[0]]
else:
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#Chrom"):
                ind_names = line.split()[9:]
                for i in ind_names:
                    pop_dict["dummy_population"] = i
                break

# Open a vcf file as input, and an output file as output
with open(input_file, "r") as f, open(output_file, "w") as o:
    
    '''
    Run the analysis for each population in the config dictionary, if there is
    no population config file, there should only be one pop. Note: if pop_dict: remains
    as it is a pain in the ass to fix the indentation. Fix at a later date. 
    '''
    if pop_dict:
        for pop in pop_dict:
            f.seek(0)
            for line in f:
               
                # Process the header: 
                if line.startswith("#CHROM"):
        
                    # Create a list of each individuals, which will be used to only calculate for those individuals in the correct population. 
                    ind_names = line.split()[9:]
                    '''
                    Create a list of individual ids for the correct population, 
                    and a dictionary dictionary, individual_dict which will store the LOD score, 
                    position and chromosome for each individual -- snp combination. 
                    '''
                    individuals_ids = []
                    for i in line.split()[9:]:
                        if i in pop_dict[pop]:
                            individuals_ids.append(i)    
                    
                    for i in individuals_ids:
                        # List of snp specific lods, list of positions, chromosome. 
                        individuals_dict[i] = [[],[],[]]
                    
                    # Write a header to the output file. 
                    o.write("#LODscore\tstart\tend\tchromosome\tindividual\n")
                
                    
                # For each line in the file that is not a header: 
                if not line.startswith("#"):
                    
                    
                    
                    '''
                    Save the chromosome, the position in the dna, and the format of the snps, for each line.
                    '''
                    chrom = line.split()[0]
                    pos = line.split()[1]
                    format_ = line.split()[8].split(":")
        
                    ###--- Part 1: Estimate allele frequencies using resampling. -----------
                    # samples_unproccessed: all samples/individuals from the vcf for that snp. Samples: samples from correct population. 
                    samples_unprocessed = line.split()[9:]
                    samples = []
                    
                    # Keep samples from correct population. Loop over the length of the unproccessed samples. 
                    for i in range(len(samples_unprocessed)):
                        
                        # If the individual is in the correct population, add it to the list of samples. 
                        if ind_names[i] in pop_dict[pop]:
                            samples.append(samples_unprocessed[i])
                    
                    # Reset variables. Haplotypes are empty as to skip an snp if all individuals are missing for that snp. allele_freq is used to store allele frequencies.
                    haplotypes = []
                    allele_freq = dict()
                    
                    # For each individual in the list of samples, add its genotype to the list of haplotypes. 
                    for ind in samples:
                        
                        # Get the genotype (same position as GT in format_ variable, as per vcf-format)
                        genotype = ind.split(":")[format_.index("GT")]
                       
                        # Split it according to delimeter, both phased and unphased data works. 
                        genotype = re.split(r'[/|]', genotype)
                        
                        # Remove all missing values and add the resulting list to the list of haplotypes. 
                        genotype = [x for x in genotype if x != "."]
                        haplotypes.extend(genotype)
                    
                    # Only continue if there is at least one haplotype. 
                    if haplotypes:
                                                
                        # Resample 40 times as per paper 1. 
                        resampled = random.choices(haplotypes, k = 40)
                        
                        # Calculate allele frequencies, and make a dictionary with the allele frequencies as results and alleles as keys. 
                        for n in resampled:
                            if n not in allele_freq: 
                                allele_freq[n] = 1
                            else:
                                allele_freq[n] += 1
                        
                        for i in allele_freq:
                            allele_freq[i] = allele_freq[i]/40
            
                        ### ------------- Part 2: Calculate LOD scores -----------------
            
                        '''
                        For each individual, calculate a lod score and add it to the 
                        dictionary individuals_dict, which stores information about each
                        snp. 
                        '''
                        
                        for ind in individuals_ids:

                            # Get the sample from list of correct samples. 
                            sample = samples[individuals_ids.index(ind)]  
                            
                            # Get the genotype (same position as GT in format_ variable, as per vcf-format)
                            genotype = sample.split(":")[format_.index("GT")]
                            
                            # Split it according to delimeter, both phased and unphased data works. 
                            genotype = re.split(r'[/|]', genotype)
                            
                            # If an allele is not in the allele frequency index due to resampling, consider it to be missing. 
                            for i in genotype:
                                if i not in allele_freq:
                                    genotype[genotype.index(i)] = "."
                           
                            # Calculate a lod score using the genotype for this individual, and the allele frequencies for this snp
                            lod = calculate_lod(genotype, allele_freq)
                            
                            # Add the results, plus the position and the chromosome, to the dictionary connecting snps to individuals
                            individuals_dict[ind][0].append(lod)
                            individuals_dict[ind][1].append(pos)
                            individuals_dict[ind][2].append(chrom)
                            
                            # If the length of the list is longer than the window size, remove the first element of each list as to loop over the genome. 
                            if len(individuals_dict[ind][0]) > window_size:
                                del individuals_dict[ind][0][0]
                                del individuals_dict[ind][1][0]
                                del individuals_dict[ind][2][0]
                        
                        
                        ### --------- Part 3: write to output -----------------------
                        
                        # For each individual, write to output. 
                        for ind in individuals_ids:
                            
                            # Get the most recent chromosome. 
                            current_chromosome = individuals_dict[ind][2][-1]
                            
                            '''
                            Get the position of that chromosome in the list of chromosomes,
                            if we are at boundaries between chromosomes this will 
                            make sure that we do not calculade lod over chromosomal boundaries.
                            '''
                            
                            chrom_start = individuals_dict[ind][2].index(current_chromosome)
                            
                            # Sum the LOD scores over the window size, as per paper 1. 
                            out_lod = sum(individuals_dict[ind][0][chrom_start:])
                            
                            # Save starting and endin position of this window, along with current chromosome, individual, and 
                            start_pos = individuals_dict[ind][1][chrom_start]
                            end_pos = individuals_dict[ind][1][-1]
                            o.write(f'{out_lod}\t{start_pos}\t{end_pos}\t{current_chromosome}\t{ind}\n')

