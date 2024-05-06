#!/usr/bin/env ruby

# Script to convert a VCF to NGS-PrimerPlex input format
# Usage: ruby vcf_2_ngsprimerplex.rb <input.vcf>

require 'zlib'

Zlib::GzipReader.open(ARGV[0]) do |f1|
	while line = f1.gets
		if line[0].chr != "#"
			line_arr = line.split
			puts line_arr[0] + "\t" + (line_arr[1].to_i - 1).to_s + "\t" + (line_arr[1].to_i + 1).to_s + "\t" + line_arr[2] + "\t1\tB\tW"
		end
	end
end