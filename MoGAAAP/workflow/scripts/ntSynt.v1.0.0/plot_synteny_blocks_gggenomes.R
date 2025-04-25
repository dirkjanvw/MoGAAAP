#!/usr/bin/env Rscript
library(argparse)
library(dplyr)
library(gggenomes)
library(gtools)
library(scales)

# Example script for generating ntSynt synteny ribbon plots using gggenomes

# Parse the input arguments
parser <- ArgumentParser(description = "Plot the ntSynt synteny blocks using gggenomes")
parser$add_argument("-s", "--sequences", help = "Input sequence lengths TSV", required = TRUE)
parser$add_argument("-l", "--links", help = "Synteny block links", required = TRUE)
parser$add_argument("--minlen", help = "Minimum length of a sequence in bases to be included (default 1 Mbp)",
                    default = 1e6, required = FALSE, type = "double")
parser$add_argument("--scale", help = "Length of scale bar in bases (default 1 Gbp)", default = 1e9,
                    required = FALSE, type = "double")
parser$add_argument("-p", "--prefix",
                    help = "Output prefix for PNG image (default synteny_gggenomes_plot)", required = FALSE,
                    default = "synteny_gggenomes_plot")

args <- parser$parse_args()


# Read in and prepare sequences
sequences <- read.csv(args$sequences, sep = "\t", header = TRUE)
sequences <- sequences[which(sequences$length > args$minlen),]

# https://stackoverflow.com/questions/32378108/using-gtoolsmixedsort-or-alternatives-with-dplyrarrange
mixedrank <- function(x) order(gtools::mixedorder(x))
sequences <- sequences %>%
  arrange(mixedrank(bin_id))


# Read in and prepare synteny links
links_ntsynt <- read.csv(args$links,
                         sep = "\t", header = TRUE)
links_ntsynt$seq_id <- factor(links_ntsynt$seq_id,
                              levels = mixedsort(unique(links_ntsynt$seq_id)))
links_ntsynt <- links_ntsynt[mixedorder(links_ntsynt$seq_id), ]
links_ntsynt$seq_id2 <- as.character(links_ntsynt$seq_id2)
links_ntsynt$colour_block <- as.factor(links_ntsynt$colour_block)

# Make the ribbon plot - these layers can be fully customized as needed!
make_plot <- function(links, sequences, add_scale_bar = FALSE) {
  p <- gggenomes(seqs = sequences, links = links)
  plot <- p + theme_gggenomes_clean(base_size = 15) +
    geom_link(aes(alpha = 0.5), offset = 0) +
    geom_seq(size = 2, colour = "grey") + # draw contig/chromosome lines
    geom_bin_label(aes(label = bin_id), size = 6, hjust = 0.9) + # label each bin
    geom_seq_label(aes(label = seq_id), vjust = 1.1, size = 4) + # Can add seq labels if desired
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # Remove x axis
    theme(legend.position="none")
  return(plot)
}

synteny_plot <- make_plot(links_ntsynt, sequences, add_scale_bar = FALSE)

# Save the ribbon plot
ggsave(paste(args$prefix, ".png", sep = ""), synteny_plot,
       units = "cm", width = 50, height = 20, bg = "white")

cat(paste("Plot saved:", paste(args$prefix, ".png", sep = ""), "\n", sep = " "))
