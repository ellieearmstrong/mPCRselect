#!/usr/bin/env python

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
