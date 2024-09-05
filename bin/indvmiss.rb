#!/usr/bin/env ruby

# indvmiss.rb from mPCRselect version 0.3.1
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

# Script to identify samples that have too much missing data from a VCFtools idepth file
# Usage is ruby indvmiss.rb <file> <missingness in float>

@header = true # Switch to ignore header line
@max_missingness = ARGV[1].to_f
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if @header
			@header = false
		else
			sample = line.split[0]
			missingness = line.strip.split[4].to_f
			if missingness > @max_missingness
				puts line
			end
		end
	end
end