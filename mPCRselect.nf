#!/usr/bin/env nextflow

process filterGQ {

	// Remove low quality sites
	
	publishDir "$params.outdir/01_GQSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/02_SingletonSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/03_MissingIndvSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/04_CullSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/05_MapSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/06_SiteFilterSNPs", mode: 'copy'
	
	input:
	path map_vcf
	
	output:
	path "${map_vcf.simpleName}.site.vcf.gz", emit: vcf
	path "${map_vcf.simpleName}.site.log"
	
	"""
	vcftools --gzvcf $map_vcf --remove-filtered-geno-all --remove-indels --min-alleles 2 --max-alleles 2 --max-missing ${params.siteMissing} --recode -c | gzip > ${map_vcf.simpleName}.site.vcf.gz
	cp .command.log ${map_vcf.simpleName}.site.log
	"""
	
}

process filterChr {

	// Retain sites only on chr in list
	// Modified from filterChr of RatesTools 0.3 (Campana & Armstrong 2020-2021)
	
	publishDir "$params.outdir/07_ChrSNPs", mode: 'copy'
	
	input:
	path(site_vcf)
	path(chrs)
	
	output:
	path "${site_vcf.simpleName}.chr.vcf.gz", emit: vcf
	path "${site_vcf.simpleName}.chr.log"
	
	"""
	chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Make into a --chr command-list > ${site_vcf.simpleName}.site.log
	vcftools --gzvcf $site_vcf --recode \$chr_line -c | gzip > ${site_vcf.simpleName}.chr.vcf.gz
	cp .command.log ${site_vcf.simpleName}.chr.log
	"""

}

process thinSNPs {

	// Thin remaining SNPs by distance using VCFtools
	
	publishDir "$params.outdir/08_ThinSNPs", mode: 'copy'
	
	input:
	path chr_vcf
	
	output:
	path "${chr_vcf.simpleName}.thin.vcf.gz", emit: vcf
	path "${chr_vcf.simpleName}.thin.log"
	
	"""
	vcftools --gzvcf $chr_vcf --recode --thin ${params.thin} -c | gzip > ${chr_vcf.simpleName}.thin.vcf.gz
	cp .command.log ${chr_vcf.simpleName}.thin.log
	"""
	
}

process plinkLD {

	// Remove sites in LD using PLINK2
	
	publishDir "$params.outdir/11_PlinkLD", mode: 'copy'
	
	input:
	path thin_vcf
	
	output:
	path "${thin_vcf.simpleName}.pruned.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.ld.log"
	path "${thin_vcf.simpleName}.pruned.fam"
	path "${thin_vcf.simpleName}.pruned.bim"
	path "${thin_vcf.simpleName}.pruned.bed"
	
	"""
	plink2 --vcf $thin_vcf --maf ${params.minMAF} -indep-pairwise ${params.plinkLD_indep_pairwise} --bad-ld --allow-extra-chr --set-all-var-ids '@:#' --make-bed --out tmp
	plink2 --bfile tmp --extract tmp.prune.in --make-bed --allow-extra-chr --export vcf-4.2 --out ${thin_vcf.simpleName}.pruned
	gzip ${thin_vcf.simpleName}.pruned.vcf
	cp .command.log ${thin_vcf.simpleName}.ld.log
	
	"""
	
}

process splitPopulations {

	// Split VCF by assigned populations
	
	publishDir "$params.outdir/09_PopulationsSNPs", mode: 'copy'
	
	input:
	path(thin_vcf)
	tuple val(samples), val(population)
	
	output:
	path "${thin_vcf.simpleName}_${population}.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}_${population}.raw", emit: raw
	path "${thin_vcf.simpleName}_${population}.log"
	
	script:
	samplelist = ""
	for (i in samples) {
		samplelist = samplelist + " --indv " + i
	}
	"""
	vcftools --gzvcf $thin_vcf --recode${samplelist} -c | gzip > ${thin_vcf.simpleName}_${population}.vcf.gz
	plink2 --vcf ${thin_vcf.simpleName}_${population}.vcf.gz --allow-extra-chr --export A --out ${thin_vcf.simpleName}_${population}
	cp .command.log ${thin_vcf.simpleName}_${population}.log
	"""

}

process optimizePi {

	// Choose sites with the highest pi values within subspecies/populations
	
	publishDir "$params.outdir/10_PopulationsPi", mode: 'copy'
	
	input:
	path pop_vcf
	
	output:
	path "${pop_vcf.simpleName}.pi.vcf.gz", emit: vcf
	path "${pop_vcf.simpleName}.pi.log"
	
	"""
	vcftools --gzvcf $pop_vcf --site-pi --out ${pop_vcf.simpleName}
	get_pi_sites.rb ${pop_vcf.simpleName}.sites.pi ${params.minPi} ${params.maxPiSNPs} ${params.maxPiPC} > pi_sites.txt
	vcftools --gzvcf $pop_vcf --positions pi_sites.txt --recode -c | gzip > ${pop_vcf.simpleName}.pi.vcf.gz
	cp .command.log ${pop_vcf.simpleName}.pi.log
	"""

}

process fstSNPs {

	// Choose sites that have highest Fst values between populations using VCFtools
	// First have to convert the CSV file of population assignments to individual population lists
	// Using get_pi_sites.rb because the code would be identical for Fst
	
	publishDir "$params.outdir/12_FstSNPs", mode: 'copy'
	
	input:
	path(thin_vcf)
	path(populations)
		
	output:
	path "${thin_vcf.simpleName}.fst.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.fst.log"
	
	"""
	split_pops.rb $populations
	popline=''
	for i in pop*.txt; do popline+=`echo '  --weir-fst-pop '\$i`; done
	vcftools --gzvcf $thin_vcf\$popline
	get_pi_sites.rb out.weir.fst ${params.minFst} ${params.maxFstSNPs} ${params.maxFstPC} > fst_snps.txt
	vcftools --gzvcf $thin_vcf --positions fst_snps.txt --recode -c | gzip > ${thin_vcf.simpleName}.fst.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.fst.log
	"""

}

process makeFstPlots {

	input:
	tuple path(pop1_raw), path(pop2_raw), val(rep)
	
	output:
	
	
	"""
	#!/usr/bin/env bash
	let totalind1=`wc -l $pop1_raw | cut -f2 -w`-1
	let totalind2=`wc -l $pop2_raw | cut -f2 -w`-1
	let totalsnps=`head -n1 $pop1_raw | grep -o \"\\t\" | wc -l`-6
	outstem=${pop1_raw.simpleName}_${pop2_raw.simpleName}_rep${rep}
	make_fst_plots.R \$totalind1 \$totalind2 \$totalsnps $pop1_raw $pop2_raw \$outstem ${params.Dcontrol} ${params.FstIndResamples} ${params.MaxFstEval} ${params.FstSNPInterval}
	"""

}

process finalSNPs {

	// Merge datasets and find most observed SNP selections
	
	publishDir "$params.outdir/13_FinalSNPs", mode: 'copy'
	
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
	vcftools --gzvcf $thin_vcf --positions best_snps.txt --recode -c | gzip > ${thin_vcf.simpleName}.fin.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.fin.log
	"""
	
}

process makePrimers {

	// Make primer sets from selected SNPs
	
	publishDir "$params.outdir/14_mPCRPrimers", mode: "copy"
	
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
	
	publishDir "$params.outdir/15_Baits", mode: "copy"
	
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
		filterGQ(params.vcf)
		removeSingletons(filterGQ.out.vcf)
		removeMissingIndv(removeSingletons.out.vcf)
		cullSNPs(removeMissingIndv.out.vcf)
		filterMappability(cullSNPs.out.vcf, params.map_bed)
		filterSites(filterMappability.out.vcf)
		filterChr(filterSites.out.vcf, params.chr_file)
		thinSNPs(filterChr.out.vcf)
		plinkLD(thinSNPs.out.vcf)
		splitPopulations(plinkLD.out.vcf, Channel.fromPath(params.populations).splitCsv(header:true).map { row -> tuple(row.Sample, row.Population) }.groupTuple(by: 1))
		optimizePi(splitPopulations.out.vcf)
		fstSNPs(plinkLD.out.vcf, params.populations)
		fst_ch = splitPopulations.out.vcf.combine(splitPopulations.out.vcf).filter { it[0] != it[1]}
		fst_ch.map { it -> it.sort() }
		fst_ch.collect().unique().view()
		//makeFstPlots(splitPopulations.out.vcf.collect())
		selected_snps_ch = optimizePi.out.vcf.mix(fstSNPs.out.vcf).collect() // Concatenate the SNP datasets for uniquing
		finalSNPs(selected_snps_ch, plinkLD.out.vcf)
		if (params.makePrimers == 1) { makePrimers(tuple finalSNPs.out.vcf, channel.fromPath(params.refseq)) }
		if (params.makeBaits == 1) { makeBaits(tuple finalSNPs.out.vcf, channel.fromPath(params.refseq)) }
		
}
