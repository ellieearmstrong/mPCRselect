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
	
	script:
	if (params.minGQ == "NULL")
		"""
		ln -s $samp_vcf ${samp_vcf.simpleName}.gq.vcf.gz
		vcftools --gzvcf $samp_vcf
		cp .command.log ${samp_vcf.simpleName}.gq.log
		"""
	else
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
	path "${sng_vcf.simpleName}.imiss", optional: true
	
	script:
	if (params.indvMissing == "NULL")
		"""
		ln -s $sng_vcf ${sng_vcf.simpleName}.indv.vcf.gz
		vcftools --gzvcf $sng_vcf
		cp .command.log ${sng_vcf.simpleName}.indv.log
		"""
	else
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
	
	script:
	if (params.cull == "NULL")
		"""
		ln -s $indv_vcf ${indv_vcf.simpleName}.cull.vcf.gz
		vcftools --gzvcf $indv_vcf
		cp .command.log ${indv_vcf.simpleName}.cull.log
		"""
	else
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
	
	publishDir "$params.outdir/08_ChrSNPs", mode: 'copy'
	
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
	
	publishDir "$params.outdir/09_ThinSNPs", mode: 'copy'
	
	input:
	path chr_vcf
	
	output:
	path "${chr_vcf.simpleName}.thin.vcf.gz", emit: vcf
	path "${chr_vcf.simpleName}.thin.log"
	
	script:
	if (params.thin == "NULL")
		"""
		ln -s $chr_vcf ${chr_vcf.simpleName}.thin.vcf.gz
		vcftools --gzvcf $chr_vcf
		cp .command.log ${chr_vcf.simpleName}.thin.log
		"""
	else
		"""
		vcftools --gzvcf $chr_vcf --recode --thin ${params.thin} -c | gzip > ${chr_vcf.simpleName}.thin.vcf.gz
		cp .command.log ${chr_vcf.simpleName}.thin.log
		"""
	
}

process plinkLD {

	// Remove sites in LD using PLINK2
	
	publishDir "$params.outdir/10_PlinkLD", mode: 'copy'
	
	input:
	path thin_vcf
	
	output:
	path "${thin_vcf.simpleName}.pruned.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.pruned.log"
	path "${thin_vcf.simpleName}.pruned.fam", optional: true
	path "${thin_vcf.simpleName}.pruned.bim", optional: true
	path "${thin_vcf.simpleName}.pruned.bed", optional: true
	
	script:
	if (params.plinkLD_indep_pairwise == "NULL")
		"""
		ln -s $thin_vcf ${thin_vcf.simpleName}.pruned.vcf.gz
		vcftools --gzvcf $thin_vcf
		cp .command.log ${thin_vcf.simpleName}.pruned.log
		"""
	else
	"""
	plink2 --vcf $thin_vcf --maf ${params.minMAF} -indep-pairwise ${params.plinkLD_indep_pairwise} --bad-ld --allow-extra-chr --chr-set ${params.haploidN} --set-all-var-ids '@:#' --make-bed --out tmp
	plink2 --bfile tmp --extract tmp.prune.in --make-bed --allow-extra-chr --chr-set ${params.haploidN} --export vcf-4.2 --out ${thin_vcf.simpleName}.pruned
	gzip ${thin_vcf.simpleName}.pruned.vcf
	cp .command.log ${thin_vcf.simpleName}.pruned.log
	"""
	
}

process splitPopulations {

	// Split VCF by assigned populations
	
	publishDir "$params.outdir/11_PopulationsSNPs", mode: 'copy'
	
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
	plink2 --vcf ${thin_vcf.simpleName}_${population}.vcf.gz --allow-extra-chr --chr-set ${params.haploidN} --export A --out ${thin_vcf.simpleName}_${population}
	cp .command.log ${thin_vcf.simpleName}_${population}.log
	"""

}

process optimizePi {

	// Choose sites with the highest pi values within subspecies/populations
	
	publishDir "$params.outdir/12_PopulationsPi", mode: 'copy'
	
	input:
	path pop_vcf
	
	output:
	path "${pop_vcf.simpleName}.pi.vcf.gz", emit: vcf
	path "${pop_vcf.simpleName}.pi.log"
	
	"""
	vcftools --gzvcf $pop_vcf --site-pi --out ${pop_vcf.simpleName}
	get_pi_sites.rb ${pop_vcf.simpleName}.sites.pi ${params.minPi} ${params.maxPiSNPs} ${params.maxPiPC} > pi_sites.txt
	vcftools --gzvcf $pop_vcf --positions pi_sites.txt --recode -c | bgzip > ${pop_vcf.simpleName}.pi.vcf.gz
	cp .command.log ${pop_vcf.simpleName}.pi.log
	"""

}

process fstSNPs {

	// Choose sites that have highest Fst values between populations using VCFtools
	// First have to convert the CSV file of population assignments to individual population lists
	// Using get_pi_sites.rb because the code would be identical for Fst
	
	publishDir "$params.outdir/13_FstSNPs", mode: 'copy'
	
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

	publishDir "$params.outdir/14_FstPlots", mode: 'copy'

	input:
	tuple path(pop1_raw), path(pop2_raw), val(rep)
	
	output:
	path "${pop1_raw.simpleName}_${pop2_raw.simpleName}_rep${rep}*"
	path "${pop1_raw.simpleName}_${pop2_raw.simpleName}_rep${rep}.highFst.csv", emit: fst_csv
	
	"""
	#!/usr/bin/env bash
	outstem=${pop1_raw.simpleName}_${pop2_raw.simpleName}_rep${rep}
	make_fst_plots.R $pop1_raw $pop2_raw \$outstem ${params.Dcontrol} ${params.FstIndResamples} ${params.MaxFstEval} ${params.FstSNPInterval}
	"""

}

process fstFinalSNPs {

	// Merge datasets and find most observed Fst SNP selections
	
	publishDir "$params.outdir/15_FstFinalSNPs", mode: 'copy'
	
	input:
	path(sel_snps)
	path(thin_vcf)
	
	output:
	path "${thin_vcf.simpleName}.finFst.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.finFst.log"
	
	script:
	filelist = sel_snps.join("\t")
	"""
	input=`echo $filelist`
	get_best_snps.rb ${params.maxFstSNPs} \$input > best_snps.txt
	vcftools --gzvcf $thin_vcf --positions best_snps.txt --recode -c | bgzip > ${thin_vcf.simpleName}.finFst.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.finFst.log
	"""
	
}

process piFinalSNPs {

	// Merge datasets and find most observed Pi SNP selections
	
	publishDir "$params.outdir/16_PiFinalSNPs", mode: 'copy'
	
	input:
	path(sel_snps)
	path(thin_vcf)
	
	output:
	path "${thin_vcf.simpleName}.finPi.vcf.gz", emit: vcf
	path "${thin_vcf.simpleName}.finPi.log"
	
	script:
	filelist = sel_snps.join("\t")
	"""
	input=`echo $filelist`
	get_best_snps.rb ${params.maxFstSNPs} \$input > best_snps.txt
	vcftools --gzvcf $thin_vcf --positions best_snps.txt --recode -c | bgzip > ${thin_vcf.simpleName}.finPi.vcf.gz
	cp .command.log ${thin_vcf.simpleName}.finPi.log
	"""
	
}

process concatFinalSNPs {

	// Concatenate Fst and Pi SNPs into a single dataset
	
	publishDir "$params.outdir/17_Fst_Pi_SNPs", mode: 'copy'
	
	input:
	path(fst_snps)
	path(pi_snps)
	
	output:
	path "${fst_snps.simpleName}.fst_pi.vcf.gz", emit: vcf
	path "${fst_snps.simpleName}.fst_pi.log"
	
	"""
	bcftools index $fst_snps
	bcftools index $pi_snps
	bcftools concat -ad all $fst_snps $pi_snps | gzip > ${fst_snps.simpleName}.fst_pi.vcf.gz
	vcftools --gzvcf ${fst_snps.simpleName}.fst_pi.vcf.gz
	cp .command.log ${fst_snps.simpleName}.fst_pi.log
	"""

}

process indexRef {

	// Index reference sequence for NGS-PrimerPlex
	
	input:
	path(refseq)
	
	output:
	path "$refseq.*"
	
	"""
	bwa index $refseq
	"""

}

process makePrimers {

	// Make primer sets from selected SNPs using NGS-PrimerPlex
	
	publishDir "$params.outdir/18_mPCRPrimers", mode: "copy"
	
	input:
	path(fin_snps)
	path(refseq)
	path('*')
	
	output:
	path "${fin_snps.simpleName}.npp.txt"
	
	
	"""
	vcf_2_ngsprimerplex.rb $fin_snps > ${fin_snps.simpleName}.npp.txt
	NGS_primerplex.py -regions ${fin_snps.simpleName}.npp.txt -ref $refseq -ad1 ${params.primerSeq1} -ad2 ${params.primerSeq2} -run ${fin_snps.simpleName} ${params.NPP_params}
	"""

}

process makeBaits {

	// Make bait sets from selected SNPs using BaitsTools
	
	publishDir "$params.outdir/19_Baits", mode: "copy"
	
	input:
	path(fin_snps)
	path(refseq)
	
	output:
	path "${fin_snps.simpleName}*"
	
	"""
	baitstools vcf2baits -i $fin_snps -r $refseq -e -o ${fin_snps.simpleName} ${params.baitsparams}
	"""

}

process plinkPCA {

	// Verify that chosen SNPs replicate structure of data using PCA
	
	publishDir "$params.outdir/20_PCAs", mode: 'copy'
	
	input:
	path(original_vcf)
	path(filtered_vcf)
	path(fst_vcf)
	path(pi_vcf)
	path(fst_pi_vcf)
	
	output:
	path '*eigenv*'
	path 'plink2.pca.log'
	
	"""
	#!/usr/bin/env bash
	for vcf in *vcf.gz; do
		vcftools --gzvcf \$vcf --min-alleles 2 --max-alleles 2 -c --recode | gzip > \${vcf%.vcf.gz}.biallelic.vcf.gz
		plink2 --vcf \${vcf%.vcf.gz}.biallelic.vcf.gz --pca biallelic-var-wts --allow-extra-chr --chr-set ${params.haploidN} --bad-freqs --out \${vcf%.vcf.gz}.biallelic
	done
	cp .command.log plink2.pca.log
	"""

}

workflow {
	main:
		if (params.samples == 'NULL') {
			filterGQ(params.vcf)
		} else {
			removeSamples(params.vcf, params.samples)
			filterGQ(removeSamples.out.vcf)
		}
		removeSingletons(filterGQ.out.vcf)
		removeMissingIndv(removeSingletons.out.vcf)
		cullSNPs(removeMissingIndv.out.vcf)
		if (params.map_bed == "NULL") {
			filterSites(cullSNPs.out.vcf)
		} else {
			filterMappability(cullSNPs.out.vcf, params.map_bed)
			filterSites(filterMappability.out.vcf)
		}
		if (params.chr_file == "NULL") {
			thinSNPs(filterSites.out.vcf)
		} else {
			filterChr(filterSites.out.vcf, params.chr_file)
			thinSNPs(filterChr.out.vcf)
		}
		plinkLD(thinSNPs.out.vcf)
		splitPopulations(plinkLD.out.vcf, Channel.fromPath(params.populations).splitCsv(header:true).map { row -> tuple(row.Sample, row.Population) }.groupTuple(by: 1))
		optimizePi(splitPopulations.out.vcf)
		fstSNPs(plinkLD.out.vcf, params.populations)
		fst_ch = splitPopulations.out.raw.combine(splitPopulations.out.raw).filter { it[0] != it[1]}.map { it -> it.sort() }.unique().combine(Channel.of(1..params.Fst_plot_repet))
		makeFstPlots(fst_ch)
		fst_selected_snps_ch = makeFstPlots.out.fst_csv.mix(fstSNPs.out.vcf).collect() // Concatenate the Fst SNP datasets for uniquing
		fstFinalSNPs(fst_selected_snps_ch, plinkLD.out.vcf)
		piFinalSNPs(optimizePi.out.vcf.collect(), plinkLD.out.vcf)
		concatFinalSNPs(fstFinalSNPs.out.vcf, piFinalSNPs.out.vcf)
		if (params.makePrimers == 1) {
			indexRef(params.refseq)
			makePrimers(concatFinalSNPs.out.vcf, params.refseq, indexRef.out) 
		}
		if (params.makeBaits == 1) { makeBaits(concatFinalSNPs.out.vcf, params.refseq) }
		if (params.samples == 'NULL') {
			plinkPCA(params.vcf, thinSNPs.out.vcf, fstFinalSNPs.out.vcf, piFinalSNPs.out.vcf, concatFinalSNPs.out.vcf)
		} else {
			plinkPCA(removeSamples.out.vcf, thinSNPs.out.vcf, fstFinalSNPs.out.vcf, piFinalSNPs.out.vcf, concatFinalSNPs.out.vcf)
		}
}
