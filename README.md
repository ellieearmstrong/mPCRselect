# mPCRselect  
Pipeline for selecting ancestry and individual identification informative SNPs from variant call files (VCF).  

__Ellie E. Armstrong, Chenyang Li, Katherine A. Solari, Jazlyn A. Mooney, Michael G. Campana, 2022-2024__  
Stanford University  
University of Southern California  
Smithsonian Institution  

## Introduction  
mPCRselect is a Nextflow [1] DSL2 pipeline for selecting and optimizing single-nucleotide polymorphisms (SNPs) for panels that reflect known population structure and identify individuals through minimized random match probability (RMP). Optionally, the pipeline can produce either massively multiplex polymerase chain reaction (PCR) primers or hybridization capture baits. Detailed descriptions of the pipeline processes and scripts are available in the [pipeline documentation](doc/pipeline_details.Md). A diagram of the complete pipeline is available [here](doc/mPCRselect.mmd). Please see the [tutorial](doc/TUTORIAL.Md) for instructions on how to configure and run the pipeline.  

mPCRselect is currently primarily designed for autosomal SNPs as it is not yet compatible with non-diploid chromosomes (e.g. sex chromosomes). Non-diploid and sex chromosomes should be removed before running the pipeline or via the `chr_file` parameter in the configuration file.  

## License  
The script `make_fst_plots.R` is derived from from code written by Chenyang Li (2024) available at: https://github.com/ChenyangLi6/SNP-panel. Original and unmodified code are subject to copyright under the University of Southern California Research License 2.0 (USC-RL v2.0), whose terms are available in the [USC-RLv2.0.txt](licenses/USC-RLv2.0.txt) file.  

![image](https://user-images.githubusercontent.com/19614608/118704084-bf02f280-b7e4-11eb-8d59-0ce648313d9e.png)  
With the exception of the copyrighted original and unmodified code in `make_fst_plots.R`, to the extent possible under law, the Smithsonian Institution and Stanford University have waived all copyright and related or neighboring rights to mPCRselect; this work is published from the United States. Please see the terms in the [CC0.txt](licenses/CC0.txt) file.    

## Installation:  
All [dependencies](#dependencies) can be installed manually following the instructions included in their external documentation. While Nextflow and the optional dependencies must be installed manually, the remaining dependencies can be installed using Conda/Mamba through the 'conda' profile included in the default mPCRselect configuration.  

Primary pipeline installation:  
1. Install Nextflow: `curl -s https://get.nextflow.io | bash`  
2. Install Conda/Mamba (Recommended): See installation instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [here](https://mamba.readthedocs.io/en/latest/installation.html#installation)  
3. Install the pipeline: `nextflow pull ellieearmstrong/mPCRselect -r <version>`, where `version` is the needed release. Setting the version as `main` will get the latest primary release.  

Install `NGS_primerplex.py` script:  
1. Clone the NGS-PrimerPlex repository: `git clone https://github.com/aakechin/NGS-PrimerPlex`  
2. Install the Python dependencies: `bash install_for_linux.sh`  
3. Make the script executable: `chmod +x NGS_primerplex.py`  
4. Modify the shebang line of the `NGS_primerplex.py` script: Change the first line from `#!/usr/bin/python3` to `#!/usr/bin/env python3` 
5. Move the `NGS_primerplex.py` script to a directory in your PATH variable.  

Install BaitsTools:  
1. Install Ruby and RubyGems: See installation instructions [here](http://www.ruby-lang.org).  
2. Build and install the BaitsTools gem following the instructions [here](https://github.com/campanam/BaitsTools).  

## Dependencies  
mPCRselect depends on the following software to perform the SNP selection and optimization pipeline:  

Programming Languages:  
* [Nextflow](https://www.nextflow.io/) v. >= 23.10.0 [1]  
* [Ruby](http://www.ruby-lang.org) v. >= 3.2.2 [2]  
* [R](https://www.r-project.org/) v. >= 4.2.3 [3]  
* [Python](https://www.python.org/) v. >= 3.0 [4]  

Standard Unix/Linux/POSIX Utilities:
* gzip  
* awk  
* Bash

Required Bioinformatics Packages:  
* [VCFtools](https://vcftools.github.io/index.html) v. 0.1.16 [5]  
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) v. >= 2.31.0 [6]  
* [PLINK 2](https://www.cog-genomics.org/plink/2.0/) v. >= 2.00a5.10 [7]  
* [bgzip](http://www.htslib.org/) from HTSlib v. >= 1.18 [8]  
* [BCFtools](http://www.htslib.org/) v. >= 1.18 [9]  

Required R Packages:
* [tidyverse](https://www.tidyverse.org/) v. 1.3.1 [9]  
* [caret](https://topepo.github.io/caret/) v. 6.0.94 [10]  
* [ggplot2](https://ggplot2.tidyverse.org/) v. 3.4.0 [11]  

## Optional Dependencies:  
The following packages are optional, but are required for multiplex primer and/or hybridization capture bait design.  

For multiplex primer design:  
* `NGS_primerplex.py` from [NGS-PrimerPlex](https://github.com/aakechin/NGS-PrimerPlex) [12]  
* [BWA](http://bio-bwa.sourceforge.net/) v. 0.7.17 [13,14]
* NGS-PrimerPlex python dependencies installed using the NGS-PrimerPlex `install_for_linux.sh` script  

For hybridization capture bait design:  
* [BaitsTools](https://github.com/campanam/BaitsTools) v. 1.8.1 [15]  

## Pipeline Configuration  
A standard local configuration profile that installs the required dependencies through Conda/Mamba is included in the `nextflow.config` file under the `conda` profile. The maximum number of processors used can be specified by the modifying the maxForks parameter (default = 32 in the `conda` profile).  

A basic explanation of configuring the software parameters is available in the [tutorial](doc/TUTORIAL.Md).  

Given the wide variety of computing architectures, we cannot provide detailed configuration settings for all software processes. Please consult your computing staff and the Nextflow documentation to generate custom profiles for your system.  

## Running the Pipeline  
Enter `nextflow run ellieearmstrong/mPCRselect -r <version> -c <config_file> -profile conda` to run the pipeline, where `version` is the installed mPCRselect release. Further details on running Nextflow pipelines are available in the official Nextflow documentation.  

## Running the Example Dataset  
A small example dataset is available in [Dryad](https://dx.doi.org/10.5061/dryad.0k6djhb96). This dataset can be run using the defaults in the `nextflow.config` file in the mPCRselect repository.  
1. Download and unzip the data.  
2. Rename the unzipped data directory to `RawData`.  
3. Run the pipeline `nextflow run ellieearmstrong/mPCRselect -r main -profile conda`.  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. Ruby: A Programmer's Best Friend (2024) http://www.ruby-lang.org. Accessed 7 May 2024.  
3. R Core Team (2020) *R: A language and environment for statistical computing.* R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/.  
4. python (2024) https://www.python.org/. Accessed 7 May 2024.  
5. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
6. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, DOI: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
7. Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M., Lee, J.J. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, __4__, s13742-015-0047-8. DOI: [10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8).  
8. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
9. Wickham, H., Averick, M., Bryan, J., Chang, W., D'Agostino McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T.L., Miller, E., Bache, S.M., Müller, K., Ooms, J., Robinson, D., Seidel, D.P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. Yutani, H. (2019). Welcome to the Tidyverse. *J Open Source Softw*, __4__, 1686. DOI: [10.21105/joss.01686](https://joss.theoj.org/papers/10.21105/joss.01686).  
10. Kuhn, M. (2008) Building Predictive Models in R Using the caret package. *J Stat Soft*, __28__, 1–26. DOI: [10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05).  
11. Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag, New York, USA.  
12. Kechin, A., Borobova, V., Boyarskikh, U., Khrapov, E., Subbotin, S., Filipenko, M. (2020) NGS-PrimerPlex: high-throughput primer design for multiplex polymerase chain reactions. *PLoS Comput Biol*, __16__, e1008468. DOI: [10.1371/journal.pcbi.1008468](https://doi.org/10.1371/journal.pcbi.1008468).  
13. Li, H., Durbin, R. (2009) Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, __25__, 1754-1760. DOI: [10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324).  
14. Li, H., Durbin, R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. *Bioinformatics*, __26__, 589-595. DOI: [10.1093/bioinformatics/btp698](https://doi.org/10.1093/bioinformatics/btp698).  
15. Campana, M.G. (2018) BaitsTools: software for hybridization capture bait design. *Mol Ecol Resour*, __18__, 356-361. DOI: [10.1111/1755-0998.12721](https://doi.org/10.1111/1755-0998.12721).  
