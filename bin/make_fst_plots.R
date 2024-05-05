#!/usr/bin/env Rscript

# Input is PLINK .raw format encoded using --recode A

library(tidyverse)
library(caret)
library(ggplot2)

fst_args = commandArgs(trailingOnly=TRUE) # Get input arguments
pop_file1 <- fst_args[1] # Raw file of genotypes for population 1
pop_file2 <- fst_args[2] # Raw file of genotypes for population 2
outstem <- fst_args[3] # Output file stem
dcontrol <- as.numeric(fst_args[4]) # Number of control replicates
resamples <- as.numeric(fst_args[5]) # Number of individuals to resample data to
requested_snps <- as.numeric(fst_args[6]) # Maximum number of SNPs to evaluate
snp_interval <- as.numeric(fst_args[7]) # Interval at which to plot SNP accuracy

num_individuals1 <- rep(totalind1, each = 2 * (requested_snps/snp_interval))
num_individuals2 <- rep(totalind2, each = 2 * (requested_snps/snp_interval))
nmarks <- c()
for(i in 1:(requested_snps/snp_interval)){
	nmarks  <- c(nmarks, i * snp_interval)
}
num_markers <- rep(nmarks, times = 2)
type <- rep(c("Highest_Fst", "Random"), each = (requested_snps/snp_interval))
accuracy <- rep(NA, 2 * (requested_snps/snp_interval))  
result <- data.frame(`Num of Individuals Pop1` = num_individuals1,`Num of Individuals Pop2` = num_individuals2, `Num of Markers` = num_markers,`Marker Selection` = type, `Accuracy` = accuracy)

dfPop1 = read_delim(pop_file1, delim = "\t") 
dfPop2 = read_delim(pop_file2, delim = "\t")
totalind1 <- nrow(dfPop1) - 1
totalind2 <- nrow(dfPop2) - 1
dfPop1 <- dfPop1[1:totalind1, 7:ncol(dfPop1)]
dfPop2 <- dfPop2[1:totalind2, 7:ncol(dfPop2)]
dfPop1 <- dfPop1[, colSums(is.na(dfPop1)) == 0]
dfPop2 <- dfPop2[, colSums(is.na(dfPop2)) == 0]

# Resample populations to account for varying sample sizes
dfPop1 <- sample_n(dfPop1, resamples, replace = TRUE)
dfPop2 <- sample_n(dfPop2, resamples, replace = TRUE)

common_columns_ab = intersect(names(dfPop1), names(dfPop2))
dfPop1 <- dfPop1[, common_columns_ab]
dfPop2 <- dfPop2[, common_columns_ab]
totalsnps <- ncol(dfPop1) # Adjust maximum number of SNPs to those that can be used
if (requested_snps > totalsnps) {requested_snps <- totalsnps} # Fix bad values
if (snp_interval > requested_snps) {snp_interval <- requested_snps} # Fix bad values

markers <- colnames(dfPop1)
frequencies <- colSums(dfPop1) / (2 * resamples)
frequencies <- as.numeric(frequencies)
freq_Pop1 <- data.frame(markers, frequencies)
frequencies2 <- colSums(dfPop2) / (2 *  resamples)
frequencies2 <- as.numeric(frequencies2)
freq_Pop2 <- data.frame(markers, frequencies2)

Fst_list <- vector("numeric", length = totalsnps) 
for(i in 1:totalsnps){
  p1 = freq_Pop1[i,2]
  q1 = 1-p1
  p2 = freq_Pop2[i,2]
  q2 = 1-p2
  pbar = (p1+p2)/2
  qbar = (q1+q2)/2
  varp = ((pbar-p1)^2 + (pbar-p2)^2)/2
  Fst <- varp/(pbar * qbar)
  Fst_list[i] <- Fst
}

nan_positions <- which(is.nan(Fst_list))
Fst_list[nan_positions] <- 0
top_indices<- order(unlist(Fst_list), decreasing = TRUE)[1:requested_snps]
freq_Pop1$frequencies[freq_Pop1$frequencies == 0] <- 0.00001
freq_Pop2$frequencies2[freq_Pop2$frequencies2 == 0] <- 0.00001
freq_Pop1$frequencies[freq_Pop1$frequencies == 1] <- 0.99999
freq_Pop2$frequencies2[freq_Pop2$frequencies2 == 1] <- 0.99999

# Output list of highest Fst sites
write.table(colnames(dfPop1)[top_indices], file = paste(outstem,".highFst.csv", sep = ''), row.names = FALSE, col.names = FALSE)

#compute delta

deltaPop1 = data.frame(resamples, totalsnps)
for(i in 1:resamples){
	for(j in 1:totalsnps){
		if(dfPop1[i,j] == 1){
			deltaPop1[i,j] = 0
		}else{
			deltaPop1[i,j] = 1
    	}
 	 }
}
deltaPop2 = data.frame(resamples, totalsnps)
for(i in 1:resamples){
	for(j in 1:totalsnps){
    	if(dfPop2[i,j] == 1){
      		deltaPop2[i,j] = 0
    	}else{
      		deltaPop2[i,j] = 1
    	}
    }
}
for(x in 1:(requested_snps/snp_interval)) {
	m <- x*snp_interval
    print(paste("Highest Fst: Num of Markers:", m))
    ppl1value = numeric(resamples)
	for(i in 1:resamples){
    	ppl1value[i] = 1
	}
  	for(i in 1:resamples){
    	for(j in 1:m){
      		if(dfPop1[i,top_indices[j]] == 2){
        		ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,top_indices[j]])*freq_Pop1[top_indices[j], 2]*freq_Pop1[top_indices[j], 2]
      		}else if(dfPop1[i,top_indices[j]]== 0){
        		ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,top_indices[j]])*(1-freq_Pop1[top_indices[j], 2])*(1-freq_Pop1[top_indices[j], 2])
      		}else {
        		ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,top_indices[j]])*freq_Pop1[top_indices[j], 2]*(1-freq_Pop1[top_indices[j], 2])
      		}
    	}
	}
    ppl2value = numeric(resamples)
	for(i in 1:resamples){
		ppl2value[i] = 1
  	}
  	for(i in 1:resamples){
    	for(j in 1:m){
      		if(dfPop1[i,top_indices[j]] == 2){
        		ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,top_indices[j]])*freq_Pop2[top_indices[j], 2]*freq_Pop2[top_indices[j], 2]
      		}else if(dfPop1[i,top_indices[j]] == 0){
    			ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,top_indices[j]])*(1-freq_Pop2[top_indices[j], 2])*(1-freq_Pop2[top_indices[j], 2])
      		}else {
				ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,top_indices[j]])*freq_Pop2[top_indices[j], 2]*(1-freq_Pop2[top_indices[j], 2])
      		}
    	}
  	}
  	matrix1 <- matrix(1, nrow = 1, ncol = resamples)
  	for(i in 1:resamples){
    	if(ppl1value[i] >= ppl2value[i]){
      		matrix1[i] = 1
    	}else{
      		matrix1[i] = 2
  		}
  	}
  	ppl1value = numeric(resamples)
  	for(i in 1:resamples){
      	ppl1value[i] = 1
  	}
  	for(i in 1:resamples){
   	 	for(j in 1:m){
      		if(dfPop2[i,top_indices[j]] == 2){
        		ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,top_indices[j]])*freq_Pop1[top_indices[j], 2]*freq_Pop1[top_indices[j], 2]
      		}else if(dfPop2[i,top_indices[j]]== 0){
        		ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,top_indices[j]])*(1-freq_Pop1[top_indices[j], 2])*(1-freq_Pop1[top_indices[j], 2])
      		}else {
        		ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,top_indices[j]])*freq_Pop1[top_indices[j], 2]*(1-freq_Pop1[top_indices[j], 2])
      		}
    	}
  	}
  	ppl2value = numeric(resamples)
  	for(i in 1:resamples){
      	ppl2value[i] = 1
  	}
  	for(i in 1:resamples){
    	for(j in 1:m){
      		if(dfPop2[i,top_indices[j]] == 2){
        		ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,top_indices[j]])*freq_Pop2[top_indices[j], 2]*freq_Pop2[top_indices[j], 2]
     		}else if(dfPop2[i,top_indices[j]] == 0){
        		ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,top_indices[j]])*(1-freq_Pop2[top_indices[j], 2])*(1-freq_Pop2[top_indices[j], 2])
      		}else {
        		ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,top_indices[j]])*freq_Pop2[top_indices[j], 2]*(1-freq_Pop2[top_indices[j], 2])
      		}
    	}
	}
  	matrix2 <- matrix(0, nrow = 1, ncol = resamples)
  	for(i in 1:resamples){
    	if(ppl1value[i]>=ppl2value[i]){
      		matrix2[i] = 1
    	}else{
      		matrix2[i] = 2
  		}
  	}
  	if (all(matrix1 == 2)) {
  		accuracy <- 0.5
	}else if (all(matrix2 == 1)) {
  		accuracy <- 0.5
	}
	else{
  		actural_label <- matrix(c(rep(1, resamples), rep(2, resamples)), nrow = 2, ncol = resamples, byrow = TRUE)
  		combined_matrix <- rbind(matrix1, matrix2)
 		predicted_classes <- combined_matrix
  		confusion_matrix<- confusionMatrix(table(predicted_classes, actural_label))
  		accuracy <- confusion_matrix$overall["Accuracy"]
	}
	l <- 0
  	l <- x
  	result[l,5] <- accuracy
}
for(y in 1:(requested_snps/snp_interval)){
	m <- y*snp_interval
   	print(paste("Random: Num of Markers:", m))
    sum = 0
    for(d in 1:dcontrol){
    	set.seed(d)
    	random_indices <- sample(1:totalsnps, requested_snps+40, replace = FALSE)
    	ppl1value = numeric(resamples)
  		for(i in 1:resamples){
    		ppl1value[i] = 1
		}
  		for(i in 1:resamples){
    		for(j in 1:m){
      			if(dfPop1[i,random_indices[j]] == 2){
        			ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,random_indices[j]])*freq_Pop1[random_indices[j], 2]*freq_Pop1[random_indices[j], 2]
      			}else if(dfPop1[i,random_indices[j]] == 0){
        			ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,random_indices[j]])*(1-freq_Pop1[random_indices[j], 2])*(1-freq_Pop1[random_indices[j], 2])
      			}else{
        			ppl1value[i] = ppl1value[i]*(2-deltaPop1[i,random_indices[j]])*freq_Pop1[random_indices[j], 2]*(1-freq_Pop1[random_indices[j], 2])
      			}
    		}
  		}
     	ppl2value = numeric(resamples)
  		for(i in 1:resamples){
    		ppl2value[i] = 1
  		}
  		for(i in 1:resamples){
    		for(j in 1:m){
      			if(dfPop1[i,random_indices[j]] == 2){
        			ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,random_indices[j]])*freq_Pop2[random_indices[j], 2]*freq_Pop2[random_indices[j], 2]
      			}else if(dfPop1[i,random_indices[j]] == 0){
        			ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,random_indices[j]])*(1-freq_Pop2[random_indices[j], 2])*(1-freq_Pop2[random_indices[j], 2])
      			}else{
        			ppl2value[i] = ppl2value[i]*(2-deltaPop1[i,random_indices[j]])*freq_Pop2[random_indices[j], 2]*(1-freq_Pop2[random_indices[j], 2])
      			}
    		}
  		}

  		matrix1 <- matrix(1, nrow = 1, ncol = resamples)
  		for(i in 1:resamples){
    		if(ppl1value[i] > ppl2value[i]){
      			matrix1[i] = 1
    		}else{
      			matrix1[i] = 2
  			}
  		}
  		ppl1value = numeric(resamples)
  		for(i in 1:resamples){
      		ppl1value[i] = 1
  		}
  		for(i in 1:resamples){
    		for(j in 1:m){
      			if(dfPop2[i,random_indices[j]] == 2){
        			ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,random_indices[j]])*freq_Pop1[random_indices[j], 2]*freq_Pop1[random_indices[j], 2]
      			}else if(dfPop2[i,random_indices[j]]== 0){
       				ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,random_indices[j]])*(1-freq_Pop1[random_indices[j], 2])*(1-freq_Pop1[random_indices[j], 2])
      			}else{
        			ppl1value[i] = ppl1value[i]*(2-deltaPop2[i,random_indices[j]])*freq_Pop1[random_indices[j], 2]*(1-freq_Pop1[random_indices[j], 2])
      			}
    		}
  		}
  		ppl2value = numeric(resamples)
  		for(i in 1:resamples){
      		ppl2value[i] = 1
  		}
  		for(i in 1:resamples){
    		for(j in 1:m){
      			if(dfPop2[i,random_indices[j]] == 2){
        			ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,random_indices[j]])*freq_Pop2[random_indices[j], 2]*freq_Pop2[random_indices[j], 2]
      			}else if(dfPop2[i,random_indices[j]]== 0){
        			ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,random_indices[j]])*(1-freq_Pop2[random_indices[j], 2])*(1-freq_Pop2[random_indices[j], 2])
      			}else{
        			ppl2value[i] = ppl2value[i]*(2-deltaPop2[i,random_indices[j]])*freq_Pop2[random_indices[j], 2]*(1-freq_Pop2[random_indices[j], 2])
      			}
    		}
  		}
  		matrix2 <- matrix(0, nrow = 1, ncol = resamples)
  		for(i in 1:resamples){
    		if(ppl1value[i]>ppl2value[i]){
      			matrix2[i] = 1
    		}else{
      			matrix2[i] = 2
  			}
  		}
  
  		actural_label <- matrix(c(rep(1, resamples), rep(2, resamples)), nrow = 2, ncol = resamples, byrow = TRUE)
  		combined_matrix <- rbind(matrix1, matrix2)
  		predicted_classes <- combined_matrix
  		confusion_matrix<- confusionMatrix(table(predicted_classes, actural_label))
  		accuracy <- confusion_matrix$overall["Accuracy"]
  		sum = sum + accuracy
   	}
	sum = sum/dcontrol
   	l <- 0
  	l <- (requested_snps/snp_interval) + y
  	result[l,5] <- sum
}

write.csv(result, file = paste(outstem,".csv", sep = ''), row.names = FALSE) # At this point Julie reran this analysis multiple times to get a distribution and generate SDs

pdf(NULL)
ggplot(result, aes(x = `Num.of.Markers`, y = Accuracy, shape = `Marker.Selection`, color = `Marker.Selection`)) +
	geom_point(size = 3) +
	labs(x = "Number of Markers", y = "Accuracy", shape = "Marker Selection") +
	scale_shape_manual(name = "Marker Selection", values = c(21, 16), labels = c("Random", "Highest Fst")) + 
	scale_colour_manual(name = "Marker Selection", values = c("blue","red"), labels = c("Random", "Highest Fst"))
	theme_bw() +
	theme(panel.background = element_rect(fill = "white"),text = element_text(size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 10))
dev.off()	

ggsave(paste(outstem,".png", sep = ''), width = 10, height = 8, dpi = 300)

dev.cur()