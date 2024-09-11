#!/usr/bin/env ruby

# Script to get highest PC-loading sites from PLINK PCA and restore original chr names
# Usage: ruby get_PCA_sites.rb <plink file> <chr mapping file> <hard-maximum-number-of-sites> <soft-percent-sites-to-retain>

# Get chr name mappings to restore
@mappings = {} # Hash of chr mappings to restore
File.open(ARGV[1]) do |f1|
	while line = f1.gets
		line_arr = line.strip.split(',')
		@mappings[line_arr[1]] = line_arr[0]
	end
end
# Restore chr name mappings and get PC1/PC2 loadings
@pos_count = 0.0 # Total number of sites
@snps = [] # Array of SNPs
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if line[0].chr != '#' # Remove header line
			line_arr = line.split
			chr = @mappings[line_arr[0]]
			pos = line_arr[1].split(':')[1]
			loading1 = line_arr[2].to_f.abs
			puts loading1
			loading2 = line_arr[3].to_f.abs
			@snps.push([chr + "\t" + pos, loading1, loading2])
			@pos_count += 1.0
		end
	end
end

# Determine actual number of sites to retain
soft_limit = (@pos_count * ARGV[3].to_f).to_i
@limit = 0
ARGV[2].to_i < soft_limit ? @limit = ARGV[2].to_i : @limit = soft_limit
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
