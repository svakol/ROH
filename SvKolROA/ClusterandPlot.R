#!/usr/bin/env Rscript
# ClusterandPlot.R

#Description: This script reads a tsv-separated file 

#Version: 1.0.0
#Author: Svante Koldenius
#Date: 22-03-2026



# Load packages
library(mclust)
library(vioplot)
library(tidyverse)

# Set environment
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

# Get inputs
input_file = args[1]
pop_config = args[2]
output_plot_template_name = args[3]
plotting_params = strsplit(args[4], "")[[1]]

# Read data
pop_config = read_tsv(pop_config, col_names = F)
data = read_tsv(input_file)

## ----------------- Step 1: Process Data --------------------

# Create an empty dataframe with the correct headers to store the data. 
processed_data = data.frame(start = double(), stop = double(), chrom = character(), ind = character(), length = double(), Model_categories = double(), Size_classifications = character(), pop = character())

# Rename the header to make it easier to handle.
data = data |> rename('Start' = `#Start`)

# Remove populations from the config file that aren't in the input
pop_config = pop_config[pop_config$X1 %in% unique(data$Ind),]

# For each population.
for (i in unique(pop_config$X2)){
  
  # Create a dataframe where only this population is kept.
  df <- data.frame(
    start = data$Start[data$Ind %in% pop_config$X1[pop_config$X2 == i]],
    stop = data$Stop[data$Ind %in% pop_config$X1[pop_config$X2 == i]],
    chrom = data$Chrom[data$Ind %in% pop_config$X1[pop_config$X2 == i]],
    ind = data$Ind[data$Ind %in% pop_config$X1[pop_config$X2 == i]],
    length = data$Length[data$Ind %in% pop_config$X1[pop_config$X2 == i]],
    pop = i
  )
  
  # Cluster using Mclust, and add these categories to the dataframe.
  model = Mclust(df$length, G = 3)
  df["Model_categories"] = model$classification

  # Create new size boundaries.
  B1 = (max(df$length[df$Model_categories == 1])+min(df$length[df$Model_categories == 2])/2)
  B2 = (max(df$length[df$Model_categories == 2])+min(df$length[df$Model_categories == 3])/2)
  
  # Add these size boundaries to the dataframe.
  df["Size_classifications"] = "B"
  df$Size_classifications[df$length <= B1] = "A"
  df$Size_classifications[df$length > B2] = "C"

  # Store this dataframe in the dataframe with the processed data. 
  processed_data = rbind(processed_data, df)
  
}

# Write it to the input file, as to save the clustering even if no plotting is done. 
write.table(processed_data, input_file, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)

### ------------------- Step 2: Plot -------------------------

## Plot by Population
if ("p" %in% c("p")){
# Create colors:
cols <- colorRampPalette(c("lightblue", "darkblue"))(length(unique(processed_data$pop)))

pop_plot = paste(output_plot_template_name,".pop",".png", sep = "")
png(filename=pop_plot)

# plotting parameters 
par(mfrow = c(2,2))
par(mgp=c(4,1,0))
par(las=1)
par(mar = c(5,7,4,2) + 0.1)
fig.dim = c(12, 6)

vioplot(
  processed_data$length[processed_data$Size_classifications == "A"] ~ processed_data$pop[processed_data$Size_classifications == "A"],
  col = cols,
  main = "Small Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Population"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "B"] ~ processed_data$pop[processed_data$Size_classifications == "B"],
  col = cols,
  main = "Medium Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Population"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "C"] ~ processed_data$pop[processed_data$Size_classifications == "C"],
  col = cols,
  main = "Large Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Population"
)
vioplot(
  processed_data$length ~ processed_data$pop,
  col = cols,
  main = "All Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Population"
)
dev.off()
}


## Plot by Chromosome

if ("c" %in% plotting_params){
# Create colors:
cols <- colorRampPalette(c("pink", "darkred"))(length(unique(processed_data$chrom)))
chrom_plot = paste(output_plot_template_name,".chrom",".png", sep = "")

png(filename=chrom_plot)

# plotting parameters 
par(mfrow = c(2,2))
par(mgp=c(4,1,0))
par(las=1)
par(mar = c(5,7,4,2) + 0.1)
fig.dim = c(12, 6)

vioplot(
  processed_data$length[processed_data$Size_classifications == "A"] ~ processed_data$chrom[processed_data$Size_classifications == "A"],
  col = cols,
  main = "Small Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Chromosome"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "B"] ~ processed_data$chrom[processed_data$Size_classifications == "B"],
  col = cols,
  main = "Medium Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Chromosome"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "C"] ~ processed_data$chrom[processed_data$Size_classifications == "C"],
  col = cols,
  main = "Large Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Chromosome"
)
vioplot(
  processed_data$length ~ processed_data$chrom,
  col = cols,
  main = "All Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Chromosome"
)
dev.off
}

## Plot by individual

if ("i" %in% plotting_params){
# Create colors:
cols <- colorRampPalette(c("lightyellow", "darkred"))(length(unique(processed_data$ind)))
ind_plot = paste(output_plot_template_name,".ind",".png", sep = "")

png(filename=ind_plot)

# plotting parameters 
par(mfrow = c(2,2))
par(mgp=c(4,1,0))
par(las=1)
par(mar = c(5,7,4,2) + 0.1)
fig.dim = c(12, 6)

vioplot(
  processed_data$length[processed_data$Size_classifications == "A"] ~ processed_data$ind[processed_data$Size_classifications == "A"],
  col = cols,
  main = "Small Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Individual"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "B"] ~ processed_data$ind[processed_data$Size_classifications == "B"],
  col = cols,
  main = "Medium Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Individual"
)
vioplot(
  processed_data$length[processed_data$Size_classifications == "C"] ~ processed_data$ind[processed_data$Size_classifications == "C"],
  col = cols,
  main = "Large Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Individual"
)
vioplot(
  processed_data$length ~ processed_data$ind,
  col = cols,
  main = "All Segments",
  ylab = "Size of ROH (bp)",
  xlab = "Individual"
)

dev.off
}














