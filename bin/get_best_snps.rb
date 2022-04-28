#!/usr/bin/env ruby

# Script to get the SNPs that are identified in the most datasets
# Usage: ruby get_best_snps.rb <max total SNPs> <files...>

require 'zlib'

@maxSNPs = ARGV[0].to_i
@snps = {} # Count of SNP occurrences keyed by chr and pos

for infile in ARGV[1..-1]
	Zlib::GzipReader.open(infile) do |f1|
		while line = f1.gets
			if line[0].chr != '#'
				line_arr = line.split
				site = line_arr[0..1].join("\t")
				if @snps.keys.include?(site)
					@snps[site]+=1
				else
					@snps[site]=1
				end
			end
		end
	end
end
@snps_sort = @snps.sort { |a,b| b[1] <=> a[1] }
@maxSNPs = @snps_sort.size if @snps_sort.size < @maxSNPs
# Get most observed sites
for i in 0 ... @maxSNPs
	puts @snps_sort[i][0]
end