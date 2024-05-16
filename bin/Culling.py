#!/usr/bin/env python3

# Culling.py from mPCRselect version 0.3.0
# Katherine A. Solari, Michael G. Campana, 2022-2024
# Stanford University and Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to mPCRselect;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# We politely request that this work be cited as:
# Armstrong EE, Li C, Campana MG, Ferrari T, Kelley JL, Petrov DA, Solari KA, Mooney JA.
# In prep. Recommendations for population and individual diagnostic SNP selection in non-
# model species.

# Modifications by M.G. Campana, 21 Apr 2022
# Header output rather than ignored
# Using system import to allow variation in culling distance. No need for fileinput.
# sys.argv[1] is VCF, sys.argv[2] is culling distance

import sys

cull = int(sys.argv[2])

previousline = ""
currentline = ""
nextline = ""

for line in open(sys.argv[1]):
  if "#" in line:
    print(line[:-1])

  previousline = currentline
  currentline = nextline
  nextline = line
  
  if previousline != "":
    partsPreviousLine = previousline.split("\t")
    partsCurrentLine = currentline.split("\t")
    partsNextLine = nextline.split("\t")

    if len(partsPreviousLine) < 2 or len(partsCurrentLine) < 2 or len(partsNextLine) < 2:
      continue

    isValid = True
    if partsPreviousLine[0] == partsCurrentLine[0]:
      if abs(int(partsPreviousLine[1]) - int(partsCurrentLine[1])) <= cull:
        isValid = False

    if partsNextLine[0] == partsCurrentLine[0]:
      if abs(int(partsNextLine[1]) - int(partsCurrentLine[1])) <= cull:
        isValid = False

    if isValid:
      print(line[:-1])
