#!/usr/bin/env ruby

# Script to identify samples that have too much missing data from a VCFtools idepth file
# Usage is ruby indvmiss.rb <file> <missingess in float>

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