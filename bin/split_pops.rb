#!/usr/bin/env ruby

# split_pops.rb from mPCRselect version 0.3.1
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

# Script to split a CSV of population assignments into population lists
# Usage: ruby split_pops.rb <pop.csv>

@pops = {} # Hash of samples keyed by pop assignment
@header = true # Switch to discard the header

File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if @header
			@header = false
		else
			line_arr = line.strip.split(',')
			sample = line_arr[0]
			pop = line_arr[1]
			if @pops.keys.include?(pop)
				@pops[pop].push(sample)
			else
				@pops[pop] = [sample]
			end
		end
	end
end
for key in @pops.keys
	File.open('pop' + key + '.txt', 'w') do |f2|
		f2.puts @pops[key].join("\n")
	end
end