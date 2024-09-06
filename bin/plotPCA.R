#!/usr/bin/env Rscript

# plotPCA.R from mPCRselect version 0.3.1
# Michael G. Campana, 2024
# Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to mPCRselect;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# We politely request that this work be cited as:
# Armstrong EE, Li C, Campana MG, Ferrari T, Kelley JL, Petrov DA, Solari KA, Mooney JA.
# In prep. Recommendations for population and individual diagnostic SNP selection in non-
# model species.

# Simple plot of PLINK2 PC1 versus 2

# Set up arguments
pca_args = commandArgs(trailingOnly=TRUE) # Get input arguments
pca <- read.table(pca_args[1]) # PLINK2 eigenvec file
pops <- read.csv(pca_args[2]) # Population assignment file
pop_list <- c()
for (i in 1:length(pca$V1)) {if (pca$V1[i] %in% pops$Sample) {pop_list <- c(pop_list, pops$Population[i])} else {pop_list <- c(pop_list, 'NoPopAssignment')}} # Handling for samples in VCF without pop assignments
eigenval <- read.table(gsub(".eigenvec",".eigenval", pca_args[1])) #PLINK2 eigenval file
total_variance <- sum(read.table(gsub(".eigenvec",".rel.diag", pca_args[1])))
pc1_percent <- round((eigenval[1,]/total_variance)*100, digits = 2)
pc2_percent <- round((eigenval[2,]/total_variance)*100, digits = 2)

png(file=gsub(".eigenvec", ".pca.png", pca_args[1]), width = 1000, height = 1000)
	pc1lab = paste("PC1: ", pc1_percent, "%", sep = "")
	pc2lab = paste("PC2: ", pc2_percent, "%", sep = "")
	poplevels <- length(levels(as.factor(pop_list)))
	plot(pca$V2, pca$V3, xlab = pc1lab, ylab = pc2lab, col=as.factor(pop_list), pch=c(1:poplevels)[as.factor(pop_list)])
	legend("bottom", legend=c(levels(as.factor(pop_list))), col=c(1:poplevels), pch=c(1:poplevels))
dev.off()	

