#!/usr/bin/env nextflow

process removeSamples {

	// Remove extraneous samples from VCF
	
	publishDir "$params.outdir/01_RemoveSamples", mode: 'copy'

	input:
	path(raw_vcf)
	path(samples)
	
	output:
	path "${raw_vcf.simpleName}.smp.vcf.gz", emit: vcf
	path "${raw_vcf.simpleName}.smp.log"
	
	"""
	vcftools --gzvcf $raw_vcf --remove $samples --recode -c | gzip > ${raw_vcf.simpleName}.smp.vcf.gz
	cp .command.log ${raw_vcf.simpleName}.smp.log
	
	"""
	
}

process filterGQ {

	// Remove low quality sites
	
	publishDir "$params.outdir/02_GQSNPs", mode: 'copy'
	
	input:
	path samp_vcf
	
	output:
	path "${samp_vcf.simpleName}.gq.vcf.gz", emit: vcf
	path "${samp_vcf.simpleName}.gq.log"
	
	"""
	vcftools --gzvcf $samp_vcf --minGQ ${params.minGQ} --recode -c | gzip > ${samp_vcf.simpleName}.gq.vcf.gz
	cp .command.log ${samp_vcf.simpleName}.gq.log
	"""

}

process removeSingletons {

	// Remove singletons and individual-unique doubletons
	
	publishDir "$params.outdir/03_SingletonSNPs", mode: 'copy'
	
	input:
	path gq_vcf
	
	output:
	path "${gq_vcf.simpleName}.sng.vcf.gz", emit: vcf
	path "${gq_vcf.simpleName}.singletons"
	path "${gq_vcf.simpleName}.sng.log"
	
	"""
	vcftools --gzvcf $gq_vcf --singletons --out ${gq_vcf.simpleName}
	vcftools --gzvcf $gq_vcf --exclude-positions ${gq_vcf.simpleName}.singletons --recode -c |  gzip > ${gq_vcf.simpleName}.sng.vcf.gz
	cp .command.log ${gq_vcf.simpleName}.sng.log
	
	"""

}

process removeMissingIndv {

	// Remove individuals with too much missing data
	
	publishDir "$params.outdir/04_MissingIndvSNPs", mode: 'copy'
	
	input:
	path sng_vcf
	
	output:
	path "${sng_vcf.simpleName}.indv.vcf.gz", emit: vcf
	path "${sng_vcf.simpleName}.indv.log"
	path "${sng_vcf.simpleName}.imiss"
	
	"""
	vcftools --gzvcf $sng_vcf --missing-indv --out ${sng_vcf.simpleName}
	indvmiss.rb ${sng_vcf.simpleName}.imiss ${params.indvMissing} > missing_samples.txt
	vcftools --gzvcf $sng_vcf --remove missing_samples.txt --recode -c | gzip > ${sng_vcf.simpleName}.indv.vcf.gz
	cp .command.log ${sng_vcf.simpleName}.indv.log
	"""

}

process cullSNPs {

	// Remove SNPs that have another SNP within the cull-distance
	
	publishDir "$params.outdir/05_CullSNPs", mode: 'copy'
	
	input:
	path indv_vcf
	
	output:
	path "${indv_vcf.simpleName}.cull.vcf.gz", emit: vcf
	path "${indv_vcf.simpleName}.cull.log"
	
	"""
	Culling.py <(gunzip -c $indv_vcf) ${params.cull} | gzip > ${indv_vcf.simpleName}.cull.vcf.gz
	vcftools --gzvcf ${indv_vcf.simpleName}.cull.vcf.gz --out ${indv_vcf.simpleName}.cull
	cp .command.log ${indv_vcf.simpleName}.cull.log
	"""
	
}

process filterMappability {

	// Remove SNPs in regions of low mappability
	
	publishDir "$params.outdir/06_MapSNPs", mode: 'copy'
	
	input:
	path(cull_vcf)
	path(map_bed)
	
	output:
	path "${cull_vcf.simpleName}.map.vcf.gz", emit: vcf
	path "${cull_vcf.simpleName}.map.log"
	
	"""
	bedtools sort -i $cull_vcf -header > cull_sort.vcf
	bedtools sort -i $map_bed -header > map_sort.bed
	bedtools intersect -a cull_sort.vcf -b map_sort.bed -v -header -sorted | gzip > ${cull_vcf.simpleName}.map.vcf.gz
	vcftools --gzvcf ${cull_vcf.simpleName}.map.vcf.gz --out ${cull_vcf.simpleName}.map
	cp .command.log ${cull_vcf.simpleName}.map.log
	"""
	
}

process filterSites {

	// Remove indels, non-bialleleic SNPs, sites with too much missing data
	// Remove sites without PASS flag
	
	publishDir "$params.outdir/07_SiteFilterSNPs", mode: 'copy'
	
	input:
	path map_vcf
	
	output:
	path "${map_vcf.simpleName}.site.recode.vcf.gz", emit: vcf
	path "${map_vcf.simpleName}.site.log"
	
	"""
	vcftools --gzvcf $map_vcf --remove-filtered-geno-all --remove-indels --min-alleles 2 --max-alleles 2 --max-missing ${params.siteMissing} --recode -c | gzip > ${map_vcf.simpleName}.site.recode.vcf.gz
	cp .command.log ${map_vcf.simpleName}.site.log
	"""
	
}

process filterChr {

	// Retain sites only on chr in list
	// Modified from filterChr of RatesTools 0.3 (Campana & Armstrong 2020-2021)
	
	publishDir "$params.outdir/08_ChrSNPs", mode: 'copy'
	
	input:
	path(site_vcf)
	path(chrs)
	
	output:
	path "${site_vcf.simpleName}.chr.recode.vcf.gz", emit: vcf
	path "${site_vcf.simpleName}.chr.log"
	
	"""
	chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Make into a --chr command-list > ${site_vcf.simpleName}.site.log
	vcftools --gzvcf $site_vcf --recode \$chr_line -c | gzip > ${site_vcf.simpleName}.chr.recode.vcf.gz
	cp .command.log ${site_vcf.simpleName}.chr.log
	"""

}

process thinSNPs {

	// Thin remaining SNPs by distance using VCFtools
	
	publishDir "$params.outdir/09_ThinSNPs", mode: 'copy'
	
	input:
	path chr_vcf
	
	output:
	path "${chr_vcf.simpleName}.thin.recode.vcf.gz", emit: vcf
	path "${chr_vcf.simpleName}.thin.log"
	
	"""
	vcftools --gzvcf $chr_vcf --recode --thin ${params.thin} -c | gzip > ${chr_vcf.simpleName}.thin.recode.vcf.gz
	cp .command.log ${chr_vcf.simpleName}.thin.log
	"""
	
}

process splitSubspecies {

	// Split VCF by assigned subspecies
	
	publishDir "$params.outdir/10_SspeciesSNPs", mode: 'copy'
	
	input:
	path(thin_vcf)
	tuple val(samples), val(sspecies)
	
	output:
	path "${thin_vcf.simpleName}_${sspecies}.recode.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}_${sspecies}.log"
	
	script:
	samplelist = ""
	for (i in samples) {
		samplelist = samplelist + " --indv " + i
	}
	"""
	vcftools --gzvcf $thin_vcf --recode${samplelist} -c | gzip > ${thin_vcf.simpleName}_${sspecies}.recode.vcf.gz
	cp .command.log ${thin_vcf.simpleName}_${sspecies}.log
	"""

}

process optimizePi {

	// Choose sites with the highest pi values within subspecies
	
	publishDir "$params.outdir/11_SspeciesPi", mode: 'copy'
	
	input:
	path sspec_vcf
	
	output:
	path "${sspec_vcf.simpleName}.pi.recode.vcf.gz", emit: vcf
	path "${sspec_vcf.simpleName}.pi.log"
	
	"""
	vcftools --gzvcf $sspec_vcf --site-pi --out ${sspec_vcf.simpleName}
	get_pi_sites.rb ${sspec_vcf.simpleName}.sites.pi ${params.minPi} ${params.maxPiSNPs} ${params.maxPiPC} > pi_sites.txt
	vcftools --gzvcf $sspec_vcf --positions pi_sites.txt --recode -c | gzip > ${sspec_vcf.simpleName}.pi.recode.vcf.gz
	cp .command.log ${sspec_vcf.simpleName}.pi.log
	"""

}

process plinkPCA {

	// Choose sites that have highest SNP weightings for separating populations using PLINK PCA
	// Remap chr names to ensure compatibility
	
	publishDir "$params.outdir/12_PCASNPs", mode: 'copy'
	
	input:
	path thin_vcf
	
	output:
	path "${thin_vcf.simpleName}.pca.recode.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.pca.log"
	
	"""
	remap_chr.rb $thin_vcf 1> remapped.vcf 2> chr_maps.csv
	plink2 --vcf remapped.vcf --pca allele-wts --allow-extra-chr --freq
	get_PCA_snps.rb plink.eigenvec.var chr_maps.csv ${params.maxPCASNPs} ${params.maxPCAPC} > pca_sites.txt
	vcftools --gzvcf $thin_vcf --positions pca_sites.txt --recode -c | gzip > ${thin_vcf.simpleName}.pca.recode.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.pca.log
	"""

}

process fstSNPs {

	// Choose sites that have highest Fst values between populations using VCFtools
	// First have to convert the CSV file of subspecies assignments to individual population lists
	// Using get_pi_sites.rb because the code would be identical for Fst
	
	publishDir "$params.outdir/13_FstSNPs", mode: 'copy'
	
	input:
	path(thin_vcf)
	path(sspecies)
		
	output:
	path "${thin_vcf.simpleName}.fst.recode.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.fst.log"
	
	"""
	split_pops.rb $sspecies
	popline=''
	for i in pop*.txt; do popline+=`echo '  --weir-fst-pop '\$i`; done
	vcftools --gzvcf $thin_vcf\$popline
	get_pi_sites.rb out.weir.fst ${params.minFst} ${params.maxFstSNPs} ${params.maxFstPC} > fst_snps.txt
	vcftools --gzvcf $thin_vcf --positions fst_snps.txt --recode -c | gzip > ${thin_vcf.simpleName}.fst.recode.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.fst.log
	"""

}

process finalSNPs {

	// Merge datasets and find most observed SNP selections
	
	publishDir "$params.outdir/14_FinalSNPs", mode: 'copy'
	
	input:
	path(sel_snps)
	path(thin_vcf)
	
	output:
	path "${thin_vcf.simpleName}.fin.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.fin.log"
	
	script:
	filelist = sel_snps.join("\t")
	"""
	input=`echo $filelist`
	get_best_snps.rb ${params.maxSNPs} \$input > best_snps.txt
	vcftools --gzvcf $thin_vcf --positions best_snps.txt --recode -c | gzip > ${thin_vcf.simpleName}.fin.recode.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.fin.log
	"""
	
}

process makePrimers {

	// Make primer sets from selected SNPs
	
	publishDir "$params.outdir/15_mPCRPrimers", mode: "copy"
	
	input:
	tuple path(fin_snps), path(refseq)
	
	output:
	
	when:
	params.makePrimers == 1
	
	"""
	"""

}

process makeBaits {

	// Make bait sets from selected SNPs using BaitsTools
	
	publishDir "$params.outdir/16_Baits", mode: "copy"
	
	input:
	tuple path(fin_snps), path(refseq)
	
	output:
	path ${fin_snps.simpleName}*
	
	when
	params.makeBaits == 1
	
	"""
	baitstools vcf2baits -i $fin_snps -r $refseq -e -o ${fin_snps.simpleName} ${params.baitsparams}
	"""

}

workflow {
	main:
		removeSamples(params.vcf, params.samples)
		filterGQ(removeSamples.out.vcf)
		removeSingletons(filterGQ.out.vcf)
		removeMissingIndv(removeSingletons.out.vcf)
		cullSNPs(removeMissingIndv.out.vcf)
		filterMappability(cullSNPs.out.vcf, params.map_bed)
		filterSites(filterMappability.out.vcf)
		filterChr(filterSites.out.vcf, params.chr_file)
		thinSNPs(filterChr.out.vcf)
		splitSubspecies(thinSNPs.out.vcf, Channel.fromPath(params.sspecies).splitCsv(header:true).map { row -> tuple(row.Sample, row.Sspecies) }.groupTuple(by: 1))
		optimizePi(splitSubspecies.out.vcf)
		plinkPCA(thinSNPs.out.vcf)
		fstSNPs(thinSNPs.out.vcf, params.sspecies)
		/* selected_snps_ch = optimizePi.out.vcf.mix(plinkPCA.out.vcf, fstSNPs.out.vcf).collect() // Concatenate the SNP datasets for uniquing
		finalSNPs(selected_snps_ch, thinSNPs.out.vcf)
		if (params.makePrimers == 1) { makePrimers(tuple finalSNPs.out.vcf, channel.fromPath(params.refseq)) }
		if (params.makeBaits == 1) { makeBaits(tuple finalSNPs.out.vcf, channel.fromPath(params.refseq)) } */
		
}
