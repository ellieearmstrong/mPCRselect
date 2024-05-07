#!/usr/bin/env Rscript

# Simple plot of PLINK2 PC1 versus 2

# Set up arguments
pca_args = commandArgs(trailingOnly=TRUE) # Get input arguments
pca <- read.table(pca_args) # PLINK2 eigenvec file
eigenval <- read.table(gsub(".eigenvec",".eigenval", pca_args)) #PLINK2 eigenval file
pc1_percent <- round(eigenval[1,]/sum(eigenval)*100, digits = 2)
pc2_percent <- round(eigenval[2,]/sum(eigenval)*100, digits = 2)

pdf(gsub(".eigenvec", ".pca.pdf", pca_args), width = 10, height = 10)
	pc1lab = paste("PC1:", pc1_percent, "%")
	pc2lab = paste("PC2:", pc2_percent, "%")
	plot(pca$V2, pca$V3, xlab = pc1lab, ylab = pc2lab)
dev.off()	

