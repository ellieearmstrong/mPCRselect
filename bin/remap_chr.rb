#!/usr/bin/env ruby

# Remap chromosome numbers for EIGENSOFT convertf
# Usage: ruby remap_chr.rb <vcf> 1> <remapped_vcf> 2> <remappings.csv>

require 'zlib'

@current_chr = '' # Current chromosome in VCF
@remapped_chr = 0 # Value to remap current chromosome
Zlib::GzipReader.open(ARGV[0]) do |f1|
	while line = f1.gets
		if line[0].chr == '#'
			puts line
		else
			line_arr = line.split("\t")
			chr = line_arr[0]
			if @current_chr != chr
				@current_chr = chr
				@remapped_chr += 1
				$stderr.puts @current_chr + ',' + @remapped_chr.to_s
			end
			puts @remapped_chr.to_s + "\t" + line_arr[1..-1].join("\t")
		end
	end
end