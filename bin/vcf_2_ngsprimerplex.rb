#!/usr/bin/env ruby

# vcf_2_ngsprimperplex.rb from mPCRselect version 0.3.0
# Michael G. Campana, 2024
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