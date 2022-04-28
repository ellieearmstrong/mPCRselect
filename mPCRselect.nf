#!/usr/bin/env nextflow

process removeSamples {

	// Remove extraneous samples from VCF
	
	publishDir "$params.outdir/01_RemoveSamples", mode: 'copy'

	input:
	path samples from params.samples
	path raw_vcf from params.vcf
	
	output:
	path "${raw_vcf.simpleName}.smp.recode.vcf.gz" into samp_vcf_ch
	path "${raw_vcf.simpleName}.smp.log"
	
	"""
	vcftools --gzvcf $raw_vcf --out ${raw_vcf.simpleName}.smp --remove $samples --recode
	cp .command.log ${raw_vcf.simpleName}.smp.log
	gzip ${raw_vcf.simpleName}.smp.recode.vcf
	"""
	
}

process filterGQ {

	// Remove low quality sites
	
	publishDir "$params.outdir/02_GQSNPs", mode: 'copy'
	
	input:
	path samp_vcf from samp_vcf_ch
	val minGQ from params.minGQ
	
	output:
	path "${samp_vcf.simpleName}.gq.recode.vcf.gz" into gq_vcf_ch
	path "${samp_vcf.simpleName}.gq.log"
	
	"""
	vcftools --gzvcf $samp_vcf --out ${samp_vcf.simpleName}.gq --minGQ $minGQ --recode
	cp .command.log ${samp_vcf.simpleName}.gq.log
	gzip ${samp_vcf.simpleName}.gq.recode.vcf
	"""

}

process removeSingletons {

	// Remove singletons and individual-unique doubletons
	
	publishDir "$params.outdir/03_SingletonSNPs", mode: 'copy'
	
	input:
	path gq_vcf from gq_vcf_ch
	
	output:
	path "${gq_vcf.simpleName}.sng.recode.vcf.gz" into sng_vcf_ch
	path "${gq_vcf.simpleName}.singletons"
	path "${gq_vcf.simpleName}.sng.log"
	
	"""
	vcftools --gzvcf $gq_vcf --singletons --out ${gq_vcf.simpleName}
	vcftools --gzvcf $gq_vcf --exclude-positions ${gq_vcf.simpleName}.singletons --recode --out ${gq_vcf.simpleName}.sng
	cp .command.log ${gq_vcf.simpleName}.sng.log
	gzip ${gq_vcf.simpleName}.sng.recode.vcf
	"""

}

process removeMissingIndv {

	// Remove individuals with too much missing data
	
	publishDir "$params.outdir/04_MissingIndvSNPs", mode: 'copy'
	
	input:
	path sng_vcf from sng_vcf_ch
	val indvMissing from params.indvMissing
	
	output:
	path "${sng_vcf.simpleName}.indv.recode.vcf" into indv_vcf_ch
	path "${sng_vcf.simpleName}.indv.log"
	
	"""
	vcftools --gzvcf $sng_vcf --missing-indv --out ${sng_vcf.simpleName}
	${params.mPCRbindir}/indvmiss.rb ${sng_vcf.simpleName}.imiss $indvMissing > missing_samples.txt
	vcftools --gzvcf $sng_vcf --remove missing_samples.txt --recode --out ${sng_vcf.simpleName}.indv
	cp .command.log ${sng_vcf.simpleName}.indv.log
	"""

}

process cullSNPs {

	// Remove SNPs that have another SNP within the cull-distance
	
	publishDir "$params.outdir/05_CullSNPs", mode: 'copy'
	
	input:
	path indv_vcf from indv_vcf_ch
	val cull from params.cull
	
	output:
	path "${indv_vcf.simpleName}.cull.vcf.gz" into cull_vcf_ch
	path "${indv_vcf.simpleName}.cull.log"
	
	"""
	${params.mPCRbindir}/Culling.py $indv_vcf $cull > ${indv_vcf.simpleName}.cull.vcf
	vcftools --vcf ${indv_vcf.simpleName}.cull.vcf --out ${indv_vcf.simpleName}.cull
	cp .command.log ${indv_vcf.simpleName}.cull.log
	gzip ${indv_vcf.simpleName}.cull.vcf
	"""
	
}

process filterMappability {

	// Remove SNPs in regions of low mappability
	
	publishDir "$params.outdir/06_MapSNPs", mode: 'copy'
	
	input:
	path cull_vcf from cull_vcf_ch
	path map_bed from params.map_bed
	
	output:
	path "${cull_vcf.simpleName}.map.vcf.gz" into map_vcf_ch
	path "${cull_vcf.simpleName}.map.log"
	
	"""
	bedtools sort -i $cull_vcf -header > cull_sort.vcf
	bedtools sort -i $map_bed -header > map_sort.bed
	bedtools intersect -a cull_sort.vcf -b map_sort.bed -v -header -sorted > ${cull_vcf.simpleName}.map.vcf
	vcftools --vcf ${cull_vcf.simpleName}.map.vcf --out ${cull_vcf.simpleName}.map
	cp .command.log ${cull_vcf.simpleName}.map.log
	gzip *vcf
	"""
	
}

process filterSites {

	// Remove indels, non-bialleleic SNPs, sites with too much missing data
	// Remove sites without PASS flag
	
	publishDir "$params.outdir/07_SiteFilterSNPs", mode: 'copy'
	
	input:
	path map_vcf from map_vcf_ch
	val siteMissing from params.siteMissing
	
	output:
	path "${map_vcf.simpleName}.site.recode.vcf.gz" into site_vcf_ch
	path "${map_vcf.simpleName}.site.log"
	
	"""
	vcftools --gzvcf $map_vcf --out ${map_vcf.simpleName}.site --remove-filtered-geno-all --remove-indels --min-alleles 2 --max-alleles 2 --max-missing $siteMissing --recode
	cp .command.log ${map_vcf.simpleName}.site.log
	gzip ${map_vcf.simpleName}.site.recode.vcf
	"""
	
}

process filterChr {

	// Retain sites only on chr in list
	// Modified from filterChr of RatesTools 0.3 (Campana & Armstrong 2020-2021)
	
	publishDir "$params.outdir/08_ChrSNPs", mode: 'copy'
	
	input:
	path site_vcf from site_vcf_ch
	path chrs from params.chr_file
	
	output:
	path "${site_vcf.simpleName}.chr.recode.vcf.gz" into chr_vcf_ch
	path "${site_vcf.simpleName}.chr.log"
	
	"""
	chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Awkwardly make into a --chr command-list > ${site_vcf.simpleName}.site.log
	vcftools --gzvcf $site_vcf --recode --out ${site_vcf.simpleName}.chr \$chr_line
	cp .command.log ${site_vcf.simpleName}.chr.log
	gzip ${site_vcf.simpleName}.chr.recode.vcf
	"""

}

process thinSNPs {

	// Thin remaining SNPs by distance using VCFtools
	
	publishDir "$params.outdir/09_ThinSNPs", mode: 'copy'
	
	input:
	path chr_vcf from chr_vcf_ch
	val thin from params.thin
	
	output:
	path "${chr_vcf.simpleName}.thin.recode.vcf.gz" into thin_vcf_ch, thin_vcf_ch2, thin_vcf_ch3, thin_vcf_ch4
	path "${chr_vcf.simpleName}.thin.log"
	
	"""
	vcftools --gzvcf $chr_vcf --out ${chr_vcf.simpleName}.thin --recode --thin $thin
	cp .command.log ${chr_vcf.simpleName}.thin.log
	gzip ${chr_vcf.simpleName}.thin.recode.vcf
	"""
	
}

Channel
	.fromPath(params.sspecies)
	.splitCsv(header:true)
	.map { row -> tuple(row.Sample, row.Sspecies) }
	.set { sspecies_ch }

sspecies_list_ch = sspecies_ch.groupTuple(by: 1)


process splitSubspecies {

	// Split VCF by assigned subspecies
	
	publishDir "$params.outdir/10_SspeciesSNPs", mode: 'copy'
	
	input:
	path thin_vcf from thin_vcf_ch
	tuple val(samples), val(sspecies) from sspecies_list_ch
	
	output:
	path "${thin_vcf.simpleName}_${sspecies}.recode.vcf.gz" into sspec_vcf_ch
	path "${thin_vcf.simpleName}_${sspecies}.log"
	
	script:
	samplelist = ""
	for (i in samples) {
		samplelist = samplelist + " --indv " + i
	}
	"""
	vcftools --gzvcf $thin_vcf --out ${thin_vcf.simpleName}_${sspecies} --recode${samplelist}
	cp .command.log ${thin_vcf.simpleName}_${sspecies}.log
	gzip ${thin_vcf.simpleName}_${sspecies}.recode.vcf
	"""

}

process optimizePi {

	// Choose sites with the highest pi values within subspecies
	
	publishDir "$params.outdir/11_SspeciesPi", mode: 'copy'
	
	input:
	path sspec_vcf from sspec_vcf_ch
	val minPi from params.minPi
	val maxPiSNPs from params.maxPiSNPs
	val maxPiPC from params.maxPiPC
	
	output:
	path "${sspec_vcf.simpleName}.pi.recode.vcf.gz" into pi_vcf_ch
	path "${sspec_vcf.simpleName}.pi.log"
	
	"""
	vcftools --gzvcf $sspec_vcf --site-pi --out ${sspec_vcf.simpleName}
	${params.mPCRbindir}/get_pi_sites.rb ${sspec_vcf.simpleName}.sites.pi $minPi $maxPiSNPs $maxPiPC > pi_sites.txt
	vcftools --gzvcf $sspec_vcf --positions pi_sites.txt --recode --out ${sspec_vcf.simpleName}.pi 
	cp .command.log ${sspec_vcf.simpleName}.pi.log
	gzip ${sspec_vcf.simpleName}.pi.recode.vcf
	"""

}
/*
process smartPCA {

	// Choose sites that have highest SNP weightings for separating populations using EIGENSOFT smartPCA
	// EIGENSOFT has a hard limit of 99 chr and specific chr naming requirements, so first need to map chr to arbitrary IDs
	// NEEDS PCA, a way to remap the identified SNPs back
	
	publishDir "$params.outdir/PCASNPs", mode: 'copy'
	
	input:
	path thin_vcf from thin_vcf_ch2
	
	"""
	${params.mPCRbindir}/remap_chr.rb $thin_vcf 1> remapped.vcf 2> chr_maps.csv
	vcftools --vcf remapped.vcf --plink --out ${thin_vcf.simpleName}
	${params.mPCRbindir}/write_parfile.rb ${thin_vcf.simpleName} 1> convertf_parfile 2> smartpca_parfile
	convertf -p convertf_parfile
	smartpca -p smartpca_parfile
	"""

}
*/

process plinkPCA {

	// Choose sites that have highest SNP weightings for separating populations using PLINK PCA
	// Remap chr names to ensure compatibility
	
	publishDir "$params.outdir/12_PCASNPs", mode: 'copy'
	
	input:
	path thin_vcf from thin_vcf_ch2
	val maxPCASNPs from params.maxPCASNPs
	val maxPCAPC from params.maxPCAPC
	
	output:
	path "${thin_vcf.simpleName}.pca.recode.vcf.gz" into pca_vcf_ch
	path "${thin_vcf.simpleName}.pca.log"
	
	"""
	${params.mPCRbindir}/remap_chr.rb $thin_vcf 1> remapped.vcf 2> chr_maps.csv
	vcftools --vcf remapped.vcf --plink --out ${thin_vcf.simpleName}
	plink --file ${thin_vcf.simpleName} --pca var-wts --allow-extra-chr
	${params.mPCRbindir}/get_PCA_snps.rb plink.eigenvec.var chr_maps.csv $maxPCASNPs $maxPCAPC > pca_sites.txt
	vcftools --gzvcf $thin_vcf --positions pca_sites.txt --out ${thin_vcf.simpleName}.pca --recode
	cp .command.log ${thin_vcf.simpleName}.pca.log
	gzip ${thin_vcf.simpleName}.pca.recode.vcf
	"""

}

process fstSNPs {

	// Choose sites that have highest Fst values between populations using VCFtools
	// First have to convert the CSV file of subspecies assignments to individual population lists
	// Using get_pi_sites.rb because the code would be identical for Fst
	
	publishDir "$params.outdir/13_FstSNPs", mode: 'copy'
	
	input:
	path thin_vcf from thin_vcf_ch3
	path sspecies from params.sspecies
	val minFst from params.minFst
	val maxFstSNPs from params.maxFstSNPs
	val maxFstPC from params.maxFstPC
		
	output:
	path "${thin_vcf.simpleName}.fst.recode.vcf.gz" into fst_vcf_ch
	path "${thin_vcf.simpleName}.fst.log"
	
	"""
	${params.mPCRbindir}/split_pops.rb $sspecies
	popline=''
	for i in pop*.txt; do popline+=`echo '  --weir-fst-pop '\$i`; done
	vcftools --gzvcf $thin_vcf\$popline
	${params.mPCRbindir}/get_pi_sites.rb out.weir.fst $minFst $maxFstSNPs $maxFstPC > fst_snps.txt
	vcftools --gzvcf $thin_vcf --positions fst_snps.txt --out ${thin_vcf.simpleName}.fst --recode
	cp .command.log ${thin_vcf.simpleName}.fst.log
	gzip ${thin_vcf.simpleName}.fst.recode.vcf
	"""

}

// Concatenate the SNP datasets for uniquing
selected_snps_ch = pi_vcf_ch.mix(pca_vcf_ch, fst_vcf_ch)

process finalSNPs {

	// Merge datasets and find most observed SNP selections
	
	publishDir "$params.outdir/14_FinalSNPs", mode: 'copy'
	
	input:
	path sel_snps from selected_snps_ch.collect()
	val maxSNPs from params.maxSNPs
	path thin_vcf from thin_vcf_ch4
	
	output:
	path "${thin_vcf.simpleName}.fin.recode.vcf.gz" into fin_vcf_ch
	path "${thin_vcf.simpleName}.fin.log"
	
	script:
	filelist = sel_snps.join("\t")
	"""
	input=`echo $filelist`
	${params.mPCRbindir}/get_best_snps.rb $maxSNPs \$input > best_snps.txt
	vcftools --gzvcf $thin_vcf --positions best_snps.txt --out ${thin_vcf.simpleName}.fin --recode
	cp .command.log ${thin_vcf.simpleName}.fin.log
	gzip ${thin_vcf.simpleName}.fin.recode.vcf
	"""
	
}
