#!/usr/bin/env ruby

# get_PCA_snps.rb from mPCRselect version 0.3.1
# Michael G. Campana, 2022-2024
# Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to mPCRselect;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# We politely request that this work be cited as:
# Armstrong EE, Li C, Campana MG, Ferrari T, Kelley JL, Petrov DA, Solari KA, Mooney JA.
# In prep. Recommendations for population and individual diagnostic SNP selection in non-
# model species.

# Script to get highest PC-loading sites from PLINK PCA
# Assumes plink PCA performed using the following pseudocode:
# vcftools --gzvcf <vcf.gz> --min-alleles 2 --max-alleles 2 --plink --out <stem>
# plink2 ---pedmap <stem> --pca biallelic-var-wts --allow-extra-chr --bad-freqs

# Usage: ruby get_PCA_sites.rb <plink .eigenvec.var file> <hard-maximum-number-of-sites> <soft-percent-sites-to-retain>

# Restore chr name mappings and get PC1/PC2 loadings
@pos_count = 0.0 # Total number of sites
@snps = [] # Array of SNPs
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if line[0].chr != '#' # Remove header line
			line_arr = line.split
			chr = line_arr[1].split(':')[0]
			pos = line_arr[1].split(':')[1]
			loading1 = line_arr[4].to_f.abs # Column bug fix from mPCRselect v. 0.2.0
			loading2 = line_arr[5].to_f.abs # Column bug fix from mPCRselect v. 0.2.0
			@snps.push([chr + "\t" + pos, loading1, loading2])
			@pos_count += 1.0
		end
	end
end

# Determine actual number of sites to retain
soft_limit = (@pos_count * ARGV[2].to_f).to_i
@limit = 0
ARGV[1].to_i < soft_limit ? @limit = ARGV[1].to_i : @limit = soft_limit
@snps_sort = [] # Array of all candidates sorted by loadings
if @limit > 0
	# Sort by absolute value of PC1 loadings
	@snps_sort_PC1 = @snps.sort { |a,b| b[1] <=> a[1] }
	@snps_sort_PC2 = @snps.sort { |a,b| b[2] <=> a[2] }
	# Interleave output by loadings
	for i in 0 ... @pos_count.to_i # ... to subtract one value
		@snps_sort.push(@snps_sort_PC1[i])
		@snps_sort.push(@snps_sort_PC2[i])
	end
	# Remove duplicate values
	@snps_sort.uniq!
	# Print final results
	for i in 0 ... @limit # ... to subtract one value
		puts @snps_sort[i][0]
	end
end
