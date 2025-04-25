#!/usr/bin/Rscript
require(pheatmap)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the data from the csv file
data <- read.csv(input_file)

# Process the data as before
rownames(data) <- data[,1]
data <- data[,3:ncol(data)-1]
data[is.na(data)] <- 0
data_matrix <- as.matrix(data)
data_matrix[upper.tri(data_matrix)] <- t(data_matrix)[upper.tri(data_matrix)]

# Create the heatmap and print to pdf
pdf(output_file)
heatmap_plot <- pheatmap(as.data.frame(data_matrix))
dev.off()
