# CHANGELOG  

### Versions 0.3.0  
Added NULL option bypass to all filter processes  
Removed PCA steps and associated scripts  
Added PLINK LD thinning  
Subspecies now called populations  
Removed '.recode' from VCF names  
splitPopulations converts to PLINK2 raw format  
get_best_snps uses PLINK2 SNP ID  
finalSNPs now fstFinalSNPs and piFinalSNPs as Fst and Pi SNPs are chosen separately  
concatFinalSNPs now merges Fst and Pi SNPs  

### Version 0.2.0  
Added dummy process to generate mPCR primers  
Added process to generate BaitsTools baits  
Removed non-functional smartPCA process  
Updated to nextflow DSL2  
More efficient VCFtools output compression  
Removed params.mPCRbindir  
$launchDir instead of $baseDir  
Rewrote plinkPCA process to use plink2  

### Version 0.1.0  
First functional version of the mPCR SNP selection pipeline  
