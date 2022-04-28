#!/usr/bin/env ruby

# Write parfiles for EIGENSOFT convertf and smartpca
# Usage: ruby write_parfile.rb <filestub> 1> <convertf parfile> 2> <smartpca_parfile>


puts 'genotypename:    ' + ARGV[0] + '.ped'
puts 'snpname:         ' + ARGV[0] + '.map'
puts 'indivname:       ' + ARGV[0] + '.ped'
puts 'outputformat:    EIGENSTRAT'
puts 'genotypeoutname: ' + ARGV[0] + '.eigenstratgeno'
puts 'snpoutname:      ' + ARGV[0] + '.snp'
puts 'indivoutname:    ' + ARGV[0] + '.ind'
puts 'familynames:     NO'

$stderr.puts 'genotypename: ' + ARGV[0] + '.eigenstratgeno'
$stderr.puts 'snpname:    ' + ARGV[0] + '.snp'
$stderr.puts 'indivname:  ' + ARGV[0] + '.ind'
$stderr.puts 'evecoutname:        ' + ARGV[0] + '.evec.txt'
$stderr.puts 'evaloutname:        ' + ARGV[0] + '.eval.txt'
