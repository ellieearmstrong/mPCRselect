#!/usr/bin/env ruby

# make_wormhole.rb from mPCRselect version 0.4.0
# Michael G. Campana, 2026
# Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to mPCRselect;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# We politely request that this work be cited as:
# Armstrong EE, Li C, Campana MG, Ferrari T, Kelley JL, Petrov DA, Solari KA, Mooney JA.
# 2025. A pipeline and recommendations for population and individual diagnostic SNP 
# selection in non-model species. Mol. Ecol. Resour. 25: e14048. 
# doi: 10.1111/1755-0998.14048.

# Script to convert a tsv-formatted slopped target into a multiplex wormhole target.csv
# Usage: make_wormhole <file.tsv> <slopdist> > <out.csv>

puts 'SEQUENCE_ID,SEQUENCE_TEMPLATE,SEQUENCE_TARGET'
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		line_arr = line.strip.split
		beginning = line_arr[0].split(":")[1].split("-")[0].to_i
		ending = line_arr[0].split(":")[1].split("-")[1].to_i
		# Discard sequences that are not the full 401 bp length
		unless ending - beginning < (ARGV[1].to_i * 2) + 1
			puts line_arr.join(",") + ",\"" + (ARGV[1].to_i + 1).to_s + ",1\""
		end
	end
end