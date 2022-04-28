#!/usr/bin/env ruby

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