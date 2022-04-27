Load a conda environment or load these modules however you please
```
module load biology bcftools
module load biology bedops
module load biology vcftools
```

First pull coordinates of indels
```
vcftools --gzvcf genotypes.highcov.vcf.gz --keep-only-indels --recode --recode-INFO-all --out genotypes.highcov.indels
```

Convert coordinates to bed format
```
vcf2bed < genotypes.highcov.indels.recode.vcf > genotypes.highcov.indels.bed
```

Restrict only to columns we need and add headers
```
cut -f1,2,3 genotypes.highcov.indels.bed > genotypes.highcov.indels.red.bed
```

Subtract/add proper amounts to each column, here we start w/25 bp around each indel
Bed files are 0 based which is a real bitch, so for coordinate 1, [0,1] would be the start/stop. 
pretty sure this is right though

```
awk '{$2=$2-25; $3=25+$3} {print}' genotypes.highcov.indels.red.bed > genotypes.highcov.indels.red.25bp.bed
```

Okay now we need to make sure of two things 1, that we are not running off our coordinate system (the genome) and 2, that we are actually in bed format
First we can just covert any lines that are negative to 0 and also remove any places where the chromEnd coordinate is longer than our scaffold
```
awk '$2 < 0 {$2=0} {print}' genotypes.highcov.indels.red.25bp.bed > genotypes.highcov.indels.red.25bp.start.bed
```

Before I do this second part, I am going to reduce to only the autosomes because who cares about the rest of the genome anyways

Note to ellie: move this wayyyyy up, no need to process all the junk regions

REMIND ME TO UPDATE THIS TO HAVING X + Y (although this prob wont work for y...)

```
sed -e 's/Scaffold_120;HRSCAF_207/chrA1/g' | sed -e 's/Scaffold_13;HRSCAF_59/chrA2/g' | sed -e 's/Scaffold_11;HRSCAF_30/chrA3/g' | sed -e 's/Scaffold_1;HRSCAF_1/chrB1/g' | sed -e 's/Scaffold_7;HRSCAF_14/chrB2/g' | sed -e 's/Scaffold_3;HRSCAF_9/chrB3/g' | sed -e 's/Scaffold_430;HRSCAF_546/chrB4/g' | sed -e 's/Scaffold_22;HRSCAF_95/chrC1/g' | sed -e 's/Scaffold_6;HRSCAF_13/chrC2/g' | sed -e 's/Scaffold_9;HRSCAF_19/chrD1/g' | sed -e 's/Scaffold_8;HRSCAF_18/chrD2/g' | sed -e 's/Scaffold_14;HRSCAF_68/chrD3/g' | sed -e 's/Scaffold_10;HRSCAF_21/chrD4/g' | sed -e 's/Scaffold_119;HRSCAF_206/chrE1/g' | sed -e 's/Scaffold_90;HRSCAF_172/chrE2/g' | sed -e 's/Scaffold_18;HRSCAF_89/chrE3/g' | sed -e 's/Scaffold_15;HRSCAF_71/chrF1/g' | sed -e 's/Scaffold_2;HRSCAF_5/chrF2/g' genotypes.highcov.indels.red.25bp.start.bed > genotypes.highcov.indels.red.25bp.start.pcc.bed
```

Okay now remove shit based on scaffold length

Now add in bed file headers or everything will yell at us forever
```
echo -e "chrom\tchromStart\tchromEnd\t" | cat genotypes.highcov.indels.reduced.header.bed > genotypes.highcov.indels.reduced.header.bed
```

The answer is of course awk fucked the tab delimited formatting so, need to fix that
