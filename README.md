# mPCRselect  
Pipeline for selecting ancestry and individual identification informative SNPs from variant call files (VCF).  

## Dependencies  
mPCRselect depends on the following software:  

Programming Languages:  
* Nextflow v. >= 23.10.0 [1]  
* Ruby v. >= 3.2.2 [2]  
* R v. >= 4.2.3 [3]  
* Python v. >= 3.0 [4]  

POSIX Utilities:
* gzip  
* awk  

Required Packages:  
* VCFtools v. 0.1.16 [5]  
* BEDtools v. >= 2.31.0 [6]  
* PLINK2 v. >= 2.00a5.10 [7]  
* tidyverse v. 1.3.1 [8]  
* caret v. 6.0.94 [9]  
* ggplot2 v. 3.4.0 [10]  
* bgzip from HTSlib v. >= 1.18 [11]  
* BCFtools v. >= 1.18 [11]  

Optional Packages:  
For multiplex primer design:  
* NGS_primerplex.py from NGS-PrimerPlex [8]  
* BWA v. 0.7.17 [9] and NGS-PrimerPlex dependencies  

For hybridization capture bait design:  
* BaitsTools v. 1.8.1 [10]  


## Installation:  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. Ruby: A Programmer's Best Friend (2024) (http://www.ruby-lang.org). Accessed 7 May 2024.  
3. R Core Team (2020) *R: A language and environment for statistical computing.* R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/.  
4. python (2024) https://www.python.org/. Accessed 7 May 2024.  
5. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
6. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, DOI: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
7. 
8. Wickham, H., Averick, M., Bryan, J., Chang, W., D'Agostino McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T.L., Miller, E., Bache, S.M., Müller, K., Ooms, J., Robinson, D., Seidel, D.P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. Yutani, H. (2019). Welcome to the Tidyverse. *J Open Source Softw*, __4__, 1686. DOI: [10.21105/joss.01686](https://joss.theoj.org/papers/10.21105/joss.01686).  
9. Kuhn, M. (2008) Building Predictive Models in R Using the caret package. *J Stat Soft*, _28_, 1–26. DOI: [10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05).  
10. Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag, New York, USA.  
11. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
