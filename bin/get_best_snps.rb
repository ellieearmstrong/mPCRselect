#!/usr/bin/env ruby

# Script to get the SNPs that are identified in the most datasets
# Usage: ruby get_best_snps.rb <max total SNPs> <files...>

require 'zlib'

#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.8: Campana 2018
def gz_file_open(file) # Determine whether input file is gzipped or not and set method to open it
	if file[-3..-1] == ".gz"
		yield Zlib::GzipReader.open(file)
	else
		yield File.open(file)
	end
end
#-----------------------------------------------------------------------------------------------

@maxSNPs = ARGV[0].to_i
@snps = {} # Count of SNP occurrences keyed by chr and pos

for infile in ARGV[1..-1]
	gz_file_open(infile) do |f1|
		while line = f1.gets
			if line[0].chr != '#'
				line_arr = line.split
				if line_arr.size == 1 # Determine if input is VCF or site list
					site = line_arr[0].split("_")[0].gsub("\"","").gsub(":","\t") # Use PLINK snp ID and remove SNP call for site lists
				else
					site = line_arr[2].split("_")[0].gsub(":","\t") # Use PLINK snp ID and remove SNP call for VCFs
				end
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