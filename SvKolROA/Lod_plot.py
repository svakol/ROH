#!/usr/bin/env python3
"""
Lod_plot.py 

Description: This script takes a tab delimeted file containing LOD scores (likelihood
    of autozygosity), tries to find the two peaks of a bimodal distribution, and categorizes
    a segment as either autozygous or not based on which peak it belongs to. Then it concatenates
    neighbouring autozygous segments.

Procedure: 

Cited papers: 
    Paper 1: Wang, S., Haynes, C., Barany, F., & Ott, J. (2009). 
    Genome-wide autozygosity mapping in human populations. Genet. Epidemiol. 33, 172-180
    Paper 2: Pemberton, T. J., Absher, D., Feldman. M. W., Myers, R. M., Rosenberg, N. A., Jun Z. L. (2012). 
    Genomic patterns of homozygosity in worldwide human populations. Am. J. Hum. Genet. 91, 275–292

Version: 1.0.0
Author: Svante Koldenius
Date: 22-03-2026
"""

# Import necessary modules. Scipy, matplotlib, and numpy are not standard modules and needs to be installed if not on your computer. 
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks


# Take input variables. 
input_file = sys.argv[1]
output_file = sys.argv[2]
pop_config = ""
pop = sys.argv[3] 
pop_config = sys.argv[4]
if len(sys.argv) > 5:
    low_peak = float(sys.argv[5])
    high_peak = float(sys.argv[6])

## ----------- part 1: Find peaks ---------

'''
Processing of the population config file. 
Create a dictionary with the populations as keys, and a list of individuals as the answers.
If there is no population config file, create one with one population for all individuals, 
as to not have to duplicate code later. 
'''

pop_dict = dict()
if pop_config:
    with open(pop_config) as f:
        for line in f:
            if line.split()[1] in pop_dict:
                pop_dict[line.split()[1]].append(line.split()[0])
            else:
                pop_dict[line.split()[1]] = [line.split()[0]]

# Read the lod scores and add to a list. Only add them if they come from the correct population. If there is no pop_dict, add all of them. 
lod_scores = []
with open(input_file) as f:
    next(f)
    for line in f:
        if not pop_dict:
            lod_scores.append(float(line.split()[0]))
        elif line.split()[-1] in pop_dict[pop]:
            lod_scores.append(float(line.split()[0]))
            
# Do gaussian kernel density estimation as per paper 1. 
lod_array = np.array(lod_scores, dtype=np.float64)
kde = stats.gaussian_kde(lod_array)

# Find the peaks of this distribution
array = np.linspace(lod_array.min(), lod_array.max(), 50)
evaluated = kde.evaluate(array)
peaks, _ = find_peaks(evaluated, prominence = 0.0001)

# If there are two peaks, take the minima between these two peaks
if len(peaks) == 2:
    threshold = minimize_scalar(kde, bounds=(array[peaks[0]], array[peaks[1]]), method='bounded')

# Else, use the peaks given by user
elif len(sys.argv) > 5:
    threshold = minimize_scalar(kde, bounds=(low_peak, high_peak), method='bounded')

#Else, plot the peaks and exit with feedback to user. 
else:
    print("Wrong amount of peaks, look at plot and decide between which areas to check for values.")
    
    plt.plot(array, evaluated, label="KDE")
    plt.xlabel("LOD-Scores")
    plt.ylabel("Density")
    plt.xticks(np.arange(round(min(array))-1, round(max(array))+1, 5))
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tick_params(axis='both', which='minor', labelsize=5)
    plt.xticks(rotation=90)
    plt.scatter(array[peaks], evaluated[peaks], color="red", label="Peaks")
    plt.savefig(f"{pop}.density.pdf", format="pdf")  
    sys.exit()

# Plot peaks and threshold. 
plt.plot(array, evaluated, label="KDE")
plt.xlabel("LOD-Scores")
plt.ylabel("Density")
plt.xticks(np.arange(round(min(array))-1, round(max(array))+1, 5))
plt.tick_params(axis='both', which='major', labelsize=6)
plt.tick_params(axis='both', which='minor', labelsize=5)
plt.xticks(rotation=90)
plt.scatter(array[peaks], evaluated[peaks], color="red", label="Peaks")
plt.scatter(threshold.x, threshold.fun, color="blue", label="Threshold")
plt.savefig(f"{pop}.density.pdf", format="pdf")  


## ------------------  Part 2: concatenate segments ---------

# Create an empty list to store significan LOD-segments
sig_lod = []

# Add all information about the significant window segments to the list, each window being a sublist. 
with open(input_file) as f:
    for line in f:
        if line[0] != "#":
            if (float(line.split()[0]) >= threshold.fun) and (line.split()[4] in pop_dict[pop]):
                # Start, Stop, Chromosome, Individual
                sig_lod.append([int(line.split()[1]), int(line.split()[2]), line.split()[3], line.split()[4]])

# Define functions that returns starting position, chromosome, and individual from the above list.
def start_pos(ind): return ind[0]
def chrom(ind): return ind[2]
def individual(ind): return ind[3]

# Use these functions with lambda to sort the list according to chromosome, individual, and starting position.
sorted_sid_lod = sorted(sig_lod, key = lambda x: (chrom(x), individual(x), start_pos(x)))

'''
Loop backwards over the segments, change the end of a segment segment if the chromosomes are the same, 
the individuals are the same, and the end positon and start positions overlap. Then delete the compared segment.
Looping backwards ensures that we don't delete segements we haven't compared.
'''
for i in range(len(sorted_sid_lod)-1,0, -1):
    if sorted_sid_lod[i][2] == sorted_sid_lod[i-1][2] and sorted_sid_lod[i][3] == sorted_sid_lod[i-1][3] and sorted_sid_lod[i][0] <= sorted_sid_lod[i-1][1]:
        sorted_sid_lod[i-1][1] = sorted_sid_lod[i][1]
        del sorted_sid_lod[i]

# Write this to output.
with open(output_file, "w") as o:
    o.write("#Start\tStop\tChrom\tInd\tLength\n")
    for i in sorted_sid_lod:
        length = i[1]-i[0]+1
        o.write(f"{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\t{length}\n")
