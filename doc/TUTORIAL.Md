# mPCRselect Pipeline Tutorial  

__Ellie E. Armstrong, Chenyang Li, Katherine A. Solari, Jazlyn A. Mooney, Michael G. Campana, 2022-2024__  
Stanford University  
University of Southern California  
Smithsonian Institution  

1. Prepare required input files  
mPCRselect's primary input is a multi-sample variant call format (VCF) file of genotypes, such as that produced by BCFtools mpileup/call [1] or the Genome Analysis Toolkit [2] GenotypeGVCFs command. Please ensure that chromosome names are compatible with [PLINK 2](https://www.cog-genomics.org/plink/2.0/) [3] before running mPCRselect.  

mPCRselect also requires a two-column headered comma-separated values (CSV) format file assigning individuals in the VCF to populations, e.g.:

Sample,Population  
\<samp1\>,\<pop1\>  
\<samp2\>,\<pop1\>  
\<samp3\>,\<pop2\>  
\<samp4\>,\<pop2\>  
\<samp5\>,\<pop3\>  
\<samp6\>,\<pop3\>  

*NB: Samples in the input VCF that do not appear in the population file will be marked as 'NoPopAssignment' in the final PCA output. They will not be used for population-specific SNP selection.*  

2. Prepare optional input files  
mPCRselect can perform optional filters using the following files:  
2a. Sample list  
mPCRselect can automatically remove unwanted samples from the original provided VCF using a sample list. This is a simple headerless text file listing the samples to remove, e.g.:  

\<samp1\>  
\<samp3\>  

2b. mPCRselect can remove low-mappability and other unwanted regions using a browser extensible data (BED) format file specifying the genomic coordinates of regions to remove.  

2c. mPCRselect can automatically filter variants to retain only polymorphisms on specified chromosomes. Chromosomes are specified in a headerless text file listing chromosomes to retain, e.g:  

\<chr1\>  
\<chr2\>  

2d. To produce multiplex primers or hybridization capture baits from the selected sites, mPCRselect requires the reference genome sequence in FASTA format.  

3. Prepare the configuration file. The simplest way is to copy and modify the `nextflow.config` file included in the mPCRselect repository:  

3a. Copy the configuration file:  
`git clone https://github.com/ellieearmstrong/mPCRselect`  
`cp mPCRselect/nextflow.config <your_name>.config`  

3b. In a text editor, modify your configuration file for your analysis. Most run parameter default values will suffice for most analyses, but the file path names for the input files MUST be modified and specified. For file path shortcuts (such as the $launchDir variable), please consult the Nextflow documentation. Parameters are described below:

*Input parameters that MUST be modified in the configuration file*  
`vcf`: Path to the gzipped input VCF.  
`populations`: Path to the population-assignment CSV.  
`samples`: Path to the list of samples to remove. Set to "NULL" to bypass filter.  
`map_bed`: Path to the BED file of regions to remove. Set to "NULL" to bypass filter.
`chr_file`: Path to the list of chromosomes to retain. Set to "NULL" to bypass filter.
`refseq`: Path to the reference sequence (FASTA format). Can be ignored if not generating primers or capture baits.  
`haploidN`: Number of retained haploid chromosomes in the genome assembly as specified in the PLINK 2 documentation. Default is `18`.  

*Other SNP-selection pipeline parameters*  
`outdir`: Path to the output directory. Default is `mPCR_results`.  
`minGQ`: Minimum genotype quality (GQ). Default is `20`. Set to "NULL" to bypass filter.  
`indvMissing`: Maximum individual missingness to retain individual. Default is `0.50`. Set to "NULL" to bypass filter.  
`cull`: Distance in bp to cull SNPs (remove both SNPs). Set to "NULL" to bypass filter.  
`siteMissing`: Maximum site missingness (as specified in VCFtools [4]). Default is `0.8`. Set to 0 to ignore site-missingness.  
thin = 100000 // Distance to randomly thin for LD avoidance. Set to "NULL" to bypass filter.
`plinkLD_indep_pairwise`: Parameter string for PLINK 2 linkage disequilbrium (LD)  `-indep-pairwise` parameter. Default is `'500kb 0.2'`. Set to "NULL" to ignore PLINK LD filter.  
`minMAF`: Minimum minor allele frequency (MAF) for PLINK 2 LD filter. Default is `0.05`.  
`maxPiSNPs`: Maximum number of highest pi SNPs to retain. Default is `1000`.  
`minPi`: Minimum site-pi (as calculated by VCFtools) for retention. Default is `0.1`.  
`maxPiPC`: Maximum proportion of pi SNPs to retain. Default is `0.1`.  
`maxFstSNPs`: Maximum number of Fst SNPs to retain for VCFtools analysis. Default is `1000`.  
`minFst`: Minimum site Fst for retention using VCFtools. Default is `0.1`.  
`maxFstPC`: Maximum proportion of Fst SNPs to retain using VCFtools.  Default is `0.1`.  
`Fst_plot_repet`: Number of analysis repetitions per population combination for Fst greedy algorithm and plots. Default is `20`.  
`Dcontrol`: Number of control (randomly selected SNPs) replicates for Fst greedy algorithm. Default is `20`.  
`FstIndResamples`: Number of individuals to resample per population for Fst greedy algorithm. Default is `20`.  
`MaxFstEval`: Maximum number of SNPs to evaluate for Fst greedy algorithm. Default is `1000`.  
`FstSNPInterval`: Interval at which to plot evaluated SNPs for Fst greedy algorithm. Default is `50`.  
	
*NGS-PrimerPlex parameters*  
`makePrimers`: Switch controlling whether to make mPCR primers (O = no, 1 = yes). Default is `0`.  
`primerSeq1`: NGS-PrimerPlex forward adapter construct sequence for NGS. Default is `'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'`.  
`primerSeq2`: NGS-PrimerPlex reverse adapter construct sequence for NGS. Default is `'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'`.  
`NPP_params`: Parameter string for NGS-PrimerPlex. See [NGS-PrimerPlex documentation](https://github.com/aakechin/NGS-PrimerPlex) for details. Default is `'-minampllen 26 -optampllen 50 -maxampllen 120 -maxprimerlen 26 -minprimerlen 18 -blast -skip -primernum1 25 -tries 500'`.  
	
*BaitsTools parameters*  
`makeBaits`: Switch controlling whether to make capture baits (O = no, 1 = yes). Default is `0`.  
`baitsparams`: Parameter string for BaitsTools vcf2baits. See [BaitsTools documentation](https://github.com/campanam/BaitsTools) fro details. Default is `'-w -c -N -n -x -J -y -B'`.  
*NB: PLINK 2 VCF output does not retain chromosome information in the header, so some BaitsTools options (such as `-l`) will not work.*  

3c. Optionally, add a custom profile for your computer architecture to the profiles section of the config file. Consult the [Nextflow documentation](https://www.nextflow.io/) and your computing staff for optimizing these parameters.  

4. Execute the pipeline:  
`nextflow run ellieearmstrong/mPCRselect -r main -c <your_name>.config`.  

Add the `-bg` command-line option to run the pipeline in the background and/or the `-resume` option to restart an incomplete prior run.  

5. Output files  
Within the specified output directory, mPCRselect will generate subdirectories containing the output from each process. Most processes will produce a (b)gzipped VCF file and a log of the process results. Subdirectories include:  
`01_RemoveSamples`: VCF after removal of unwanted samples. Directory will be omitted if process bypassed.  
`02_GQSNPs`: VCF after filtering for GQ.  
`03_SingletonSNPs`: VCF after singleton and individual-unique doubleton removal. Also will include a '.singleton' file listing singleton sites removed.  
`04_MissingIndvSNPs`: VCF after filtering individuals for missingness.  
`05_CullSNPs`: VCF after culling SNPs for having a nearby SNP.  
`06_MapSNPs`: VCF after removing low-mappability and unwanted regions. Directory will be omitted if process bypassed.  
`07_SiteFilterSNPs`: VCF after removing indels, non-biallelic SNPs, sites with too much missing data, and sites that do not pass quality filters.  
`08_ChrSNPs`: VCF after removing unwanted chromosomes and contigs. Directory will be omitted if process bypassed.  
`09_ThinSNPs`: VCF after thinning SNPs by physical distance.  
`10_PlinkLD`: VCF after filtering by linkage-disequilibrium in PLINK 2. Directory will also include PLINK 2 files ('.bed', '.bim', '.fam').  
`11_PopulationsSNPs`: VCFs after splitting by population assignment. Directory will also contain these files converted to RAW format.  
`12_PopulationsPi`: VCFs of highest pi sites for each population.  
`13_FstSNPs`: VCF of highest Fst sites (as calculated by VCFtools).  
`14_FstPlots`: Output from the greedy algorithm. For each replicate run, there will be a CSV ('repN.csv') and plot (PNG format) comparing the accuracy of highest Fst SNPs to randomly selected SNPs. There is also a list of the highest Fst sites for that replicate ('repN.highFst.csv').  
`15_FstFinalSNPs`: VCF of the final selected Fst SNPs.  
`16_PiFinalSNPs`: VCF of the final selected pi SNPs. There will also be a text file ('.RMP.txt') giving the highest and lowest calculated RMP values and the pi site allele frequencies as calculated by VCFtools.  
`17_Fst_Pi_SNPs`: VCF of the concatenated Fst and pi SNPs.  
`18_mPCRPrimers`: Output of NGS-PrimerPlex. Directory will be omitted if primers are not made.  
`19_Baits`: Output of BaitsTools. Directory will be omitted if baits are not made.  
`20_plinkPCA`: PLINK 2 PCA output and plots for the primary stages of the pipeline.  

6. Optimization of the SNP panel  
It is difficult to predict the number of SNPs that will successfully yield informative mPCR primers or baits in advance since each genome has its own specific characteristics (e.g. GC content, repetitive elements, etc.). End-users may need to optimize their mPCRselect configurations analytically to recover the optimal panel. Since it may not be possible to generate usable primers or baits for all SNPs in the mPCRselect panel, we recommend that end-users generate SNP panels that are ~10-fold larger than the expected number of required SNPs. We found that requesting ~1000 SNPs worked well for panels where only ~100 SNPs were needed. We note that there are diminishing returns with increased panel size in terms of population assignment and random match probability. The graphic output from the `14_FstPlots` output can help end-users determine the number of SNPs they need to retain from the designed panel. The 'highFst.csv' lists in this directory are sorted from highest to lowest Fst, so we recommend selecting SNPs that are higher in the list first. The random match probability of a subselected panel can be recalculated by the end-user using the `RMP_calc.R` script included in the `bin` directory of this repository.  

6a. Make sure you have `tidyverse` [5] installed in your R environment.  

6b. Subset the '.RMP.txt' file from `16_PiFinalSNPs` directory to only include the frequency lines of the selected SNPs. Remove all headers and RMP calculation lines.  

6c. Run the `RMP_calc.R` script on the subselected file:  
`Rscript RMP_calc.R <subselected file>`  
## References  
1. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
2.  McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
3. Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M., Lee, J.J. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, __4__, s13742-015-0047-8. DOI: [10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8).  
4. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
6. Wickham, H., Averick, M., Bryan, J., Chang, W., D'Agostino McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T.L., Miller, E., Bache, S.M., Müller, K., Ooms, J., Robinson, D., Seidel, D.P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. Yutani, H. (2019). Welcome to the Tidyverse. *J Open Source Softw*, __4__, 1686. DOI: [10.21105/joss.01686](https://joss.theoj.org/papers/10.21105/joss.01686).  
