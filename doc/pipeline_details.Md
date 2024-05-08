# mPCRselect Pipeline Details  

__Ellie E. Armstrong, Chenyang "Julie" Li, Katherine A. Solari, Jazlyn A. Mooney, Michael G. Campana, 2022-2024__  
Stanford University  
University of Southern California  
Smithsonian Institution  

The following describes the mPCRselect pipeline processes in the order in which they are executed.  

1. `removeSamples`: Using VCFtools [1], the `removeSamples` process removes unwanted samples from the original input VCF file. This process can be bypassed by setting the params.samples variable to 'NULL'.  
2. `filterGQ`: Using VCFtools, the `filterGQ` process removes sites with a genotype quality (GQ) below the specified threshold. This filter can be ignored by setting the params.minGQ variable to 'NULL'.  
3. `removeSingletons`: Using VCFtools, the `removeSingletons` process removes singleton and individual-unique doubleton sites.  
4. `removeMissingIndv`: Using VCFtools and the custom Ruby script `indvmiss.rb`, the `removeMissingIndv` process removes samples below the VCFtools completion threshold (i.e., between 0 for allowing 100% missing data and 1 for requiring 100% complete data). This filter can be ignored by setting the params.indvMissing variable to 'NULL'.  
`indvmiss.rb` is executed: `indvmiss.rb <VCFtools idepth file> <missingness>`  
5. `cullSNPs`: Using a custom Python script (`Culling.py`), this process removes all sites for which there is another site within the specified window (in bp) up or downstream of it. This process reduces the chance of misprimingmwhen developing multiplex PCR primers. This filter can be ignored by setting the params.cull variable to 'NULL'.  
`Culling.py` is executed: `Culling.py <input uncompressed VCF> <cull distance in bp>`  
6. `filterMappability`: The `filterMappability` process uses BEDtools [2] intersect to remove low-mappability regions specified in the input BED file. This process can be bypassed by setting the params.map_bed variable to 'NULL'.  
7. `filterSites`: The `filterSites` process using VCFtools to remove non-biallelic sites, indels, filtered genotypes, and those sites with too much missing data.  
8. `filterChr`: Using VCFtools, the `filterChr` process retains only sites assigned to the chromosomes listed in the input chromosome list. This process can be bypassed by setting the params.chr_file variable to 'NULL'.  
9. `thinSNPs`: The `thinSNPs` process uses VCFtools to thin SNPs by physical distance in bp. This filter can be ignored by setting the params.thin variable to 'NULL'.  
10. `plinkLD`: The `plinkLD` process uses PLINK 2 [3] to thin SNPs by linkage-disequilibrium. This filter can be ignored by setting the params.plinkLD_indep_pairwise variable to 'NULL'.  
11. `splitPopulations`: The `splitPopulations` process uses VCFtools to subset the filtered VCF by population to produce individual VCFs for each population in the dataset. These data are converted to RAW format using PLINK 2.  
12. `optimizePi`: The `optimizePi` process calculates per-site pi values for each population-specific VCF using VCFtools `--site-pi`. The process then uses the custom Ruby `get_pi_sites.rb` script to identify the sites with the highest pi values. The process then uses VCFtools to subselect the highest pi sites.  
`get_pi_sites` is executed: `get_pi_sites.rb <input VCF> <minimum-Pi> <hard-maximum-number-of-sites> <soft-percent-sites-to-retain>`  
13. `fstSNPs`: The `fstSNPs` process chooses sites that have highest Fst values between populations using VCFtools. It uses the custom Ruby script `split_pops.rb` to generate population lists for Fst comparisons and the `get_pi_sites.rb` to identify the sites with the highest Fst values. The process then uses VCFtools to subselect the highest Fst sites.  
`split_pops.rb` is executed: `split_pops.rb <population assignment CSV>`  
14. `makeFstPlots`: Using the custom R script `make_fst_plots.R`, the `makeFstPlots` process performs the novel greedy algorithm to identify the highest Fst sites and plots accuracy of the Fst-chosen SNPs against randomly chosen ones.  
`make_fst_plots.R` is executed: `make_fst_plots.R <population1 RAW> <population2 RAW> <output file stem> <Number of random individual resamples> <maximum number of SNPs> <plotting interval of SNPs>`  
15. `fstFinalSNPs`: The `fstFinalSNPs` process uses the custom Ruby script `get_best_snps.rb` to find the sites that appear the most frequently across the maximum Fst runs (both greedy-algorithm and using VCFtools). The process then uses VCFtools to select the final Fst site VCF.  
`get_best_snps.rb` is executed: `get_best_snps.rb <maximum Fst SNPs> <input VCF1> ... <input VCFn>`  
16. `piFinalSNPs`: The `piFinalSNPs` process uses the `get_best_snps.rb` script to identify the sites that appear most frequently in the selected highest-pi datasets. The process then uses VCFtools to select the final pi site VCF.  
17. `concatFinalSNPs`: The `concatFinalSNPs` process uses BCFtools [4] concat to merge the final Fst and pi selected SNPs.  
18. `indexRef`: The `indexRef` process uses BWA [5,6] index to index the reference genome for downstream NGS-PrimerPlex [7] primer design.  
19. `makePrimers`: The `makePrimers` process uses the custom Ruby script `vcf_2_ngsprimerplex.rb` to convert the final selected site VCF to NGS-PrimerPlex input format. It then uses the `NGS_primerplex.py` script from NGS-PrimerPlex to generate multiple PCR primers.  
`vcf_2_ngsprimerplex.rb` is executed: `vcf_2_ngsprimerplex.rb <input VCF>`
20. `makeBaits`: The `makeBaits` process uses BaitsTools [8] vcf2baits to generate bait sequences from the final selected site VCF.  
21. `plinkPCA`: The `plinkPCA` process performs principal component analysis (PCA) on the original VCF, the filtered VCF, the final Fst VCF, the final pi VCF, and final concatenated Fst and pi VCF to ensure that the SNP panel captured the observed variability. First, the process uses VCFtools to retain only biallelic sites. Next, PCA is performed using PLINK 2. Finally, PCAs are plotted using the custom R script `plotPCA.R`.  
`plotPCA.R` is executed: `plotPCA.R <PLINK 2 eigenvec file>`  

## References  
1. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
2. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, DOI: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
3. Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M., Lee, J.J. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, __4__, s13742-015-0047-8. DOI: [10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8).  
4. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
5. Li, H., Durbin, R. (2009) Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, __25__, 1754-1760. DOI: [10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324).  
6. Li, H., Durbin, R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. *Bioinformatics*, __26__, 589-595. DOI: [10.1093/bioinformatics/btp698](https://doi.org/10.1093/bioinformatics/btp698).  
7. Kechin, A., Borobova, V., Boyarskikh, U., Khrapov, E., Subbotin, S., Filipenko, M. (2020) NGS-PrimerPlex: high-throughput primer design for multiplex polymerase chain reactions. *PLoS Comput Biol*, __16__, e1008468. DOI: [10.1371/journal.pcbi.1008468](https://doi.org/10.1371/journal.pcbi.1008468).  
8. Campana, M.G. (2018) BaitsTools: software for hybridization capture bait design. *Mol Ecol Resour*, __18__, 356-361. DOI: [10.1111/1755-0998.12721](https://doi.org/10.1111/1755-0998.12721).  




