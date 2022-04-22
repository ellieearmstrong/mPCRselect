File description:
```
1. Unimputed vcf file with all snp called tigers: genotypes.highcov.vcf.gz
2. Imputed vcf files from gencove: bengal-data-imputed.vcf.gz, southchina_for_merge.vcf.gz
```

First download all files from gencove 
```
gencove download --samples ids.txt --
```

Merge vcfs 
```
bcftools merge -O z -o bengal-data-imputed.vcf.gz *.vcf.gz
```

Change scaffold names from imputed vcfs so they match original files
```
zcat bengal-data-imputed.vcf.gz | sed 's/-HRSCAF/;HRSCAF/g' > bengal-data-imputed.eh.vcf
bgzip bengal-data-imputed.eh.vcf
tabix -p vcf bengal-data-imputed.eh.vcf.gz

zcat southchina_for_merge.vcf.gz | sed 's/-HRSCAF/;HRSCAF/g' > southchina_for_merge.eh.vcf 
bgzip southchina_for_merge.eh.vcf 
tabix -p vcf southchina_for_merge.eh.vcf.gz
```

Merge with other files 
```
bcftools merge -O z -o all-merge-mpcr.vcf.gz bengal-data-imputed.eh.vcf.gz genotypes.highcov.vcf.gz southchina_for_merge.vcf.gz
```

Remove duplicates
```
bcftools view -S \^captives_duplicates.txt all-merge-mpcr.vcf.gz -O z -o all-merge-dups-generics-removed-mpcr.vcf.gz
```
Filter out SNPs with very low genotype quality (removing SNPs that we think are probably not real)
```
vcftools --gzvcf all-merge-dups-generics-removed-mpcr.vcf.gz --out all-merge-dups-generics-removed-GQ20-mpcr --minGQ 20 --recode
```
Remove singletons, including doubletons unique to one individual 
```
vcftools --vcf all-merge-dups-generics-removed-GQ20-mpcr.recode.vcf --singletons --out all-merge-dups-generics-removed-GQ20-mpcr_singletons
cut -f 1,2 all-merge-dups-generics-removed-GQ20-mpcr_singletons.singletons > singletons2remove.txt
vcftools --vcf all-merge-dups-generics-removed-GQ20-mpcr.recode.vcf --exclude-positions singletons2remove.txt --recode --out all-merge-dups-generics-removed-GQ20-NoSing-mpcr
```

