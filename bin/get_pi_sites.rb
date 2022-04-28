#!/usr/bin/env ruby

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