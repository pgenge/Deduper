#!/usr/bin/env python
import re

def parsecigar_adjustpos(cigar, strand, samposition):
    """This function parses the cigar string to adjust the position of the read to the reference
    based on a these factors: Soft Clipping, Minus strand"""
    #check strand
    matched = 0
    dels = 0
    skipped = 0
    softclip = 0
    if strand == 'plus':
      if 'S' in cigar:
        S = re.findall(r'^([0-9]+)S+', cigar)
        if len(S) > 0:
          softclip = int(S[0])
        else:
          softclip = 0
        adjusted_pos = samposition - softclip
        return adjusted_pos
    else:
      if 'S' in cigar[-1]:
        S = re.findall(r'([0-9]+)S', cigar)
        if len(S) > 0:
          softclip = int(S[0])
        else:
          softclip = 0
      if 'M' in cigar:
        match = re.findall(r'([0-9]+)M', cigar)
        matches = [int(i) for i in match]
        matched = sum(matches)
      if 'D' in cigar:
        deletion = re.findall(r'([0-9]+)D', cigar)
        deletions = [int(i) for i in deletion]
        dels = sum(deletions)
      if 'N' in cigar:
        skip = re.findall(r'([0-9]+)N', cigar)
        skips = [int(i) for i in skip]
        skipped = sum(skips)
    adjusted_pos = samposition + matched + dels + skipped + softclip
    return adjusted_pos

#test cases
print(parsecigar_adjustpos("15S40M", "plus",20)) #5
print(parsecigar_adjustpos("10M4D1N5S", "minus", 10)) #30
print(parsecigar_adjustpos("10M1N5S", "minus", 10)) #26
print(parsecigar_adjustpos("10M4D5S", "minus", 10)) #29
print(parsecigar_adjustpos("10M", "minus", 10)) #20
print(parsecigar_adjustpos("10S", "minus", 10)) #20
print(parsecigar_adjustpos("5S10M", "minus", 10)) #20
print(parsecigar_adjustpos("10M5S", "minus", 10)) #25
print(parsecigar_adjustpos("5S10M1N1I", "plus", 10)) #5
print(parsecigar_adjustpos("10M5S", "plus", 10)) #10
print(parsecigar_adjustpos("10M5S", "minus", 10)) #25