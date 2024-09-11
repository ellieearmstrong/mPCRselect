#!/usr/bin/env ruby

# get_pi_sites.rb from mPCRselect version 0.3.2
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

# Script to get highest pi sites from VCFtools site-pi file
# Usage: ruby get_pi_sites.rb <file> <minimum-Pi> <hard-maximum-number-of-sites> <soft-percent-sites-to-retain>

@header = true # Switch to ignore header line
@pos_count = 0.0 # Total number of sites
@sites = [] # Array of candidate sites above minimum pi
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if @header
			@header = false
		elsif line.strip != "" # Discount blank lines
			@pos_count += 1
			pi = line.strip.split[2].to_f
			@sites.push([pi,line]) if pi >= ARGV[1].to_f
			@pos_count += 1.0
		end
	end
end
# Determine actual number of sites to retain
soft_limit = (@pos_count * ARGV[3].to_f).to_i
@limit = 0
ARGV[2].to_i < soft_limit ? @limit = ARGV[2].to_i : @limit = soft_limit
@limit = @sites.size if @limit > @sites.size
if @limit > 0
	@sorted_sites = @sites.sort { |a,b| b[0] <=> a[0] }
	for i in 0 ... @limit # ... to subtract one value
		puts @sorted_sites[i][1]
	end
end