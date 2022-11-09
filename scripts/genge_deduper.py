#!/usr/bin/env python

#example command: ./genge_deduper.py -u <UMI_file.txt>.txt -f <input.sam> -o <deduped_file_name.sam>
#command to run on test file
  # ./scripts/genge_deduper.py -f test.sam -u STL96.txt -o deduped.sam

#command to sort sam file via samtools documentation
  #samtools sort C1_SE_uniqAlign.sam -o C1_SE_uniqAlign.sorted.sam

#command to run on C1_SE_uniqAlign.sorted.sam
  #./scripts/genge_deduper.py -f C1_SE_uniqAlign.sorted.sam -u STL96.txt -o deduped_C1_SE_uniqAlign.sorted.sam

#import all required modules
import argparse
import re

def get_args():
    #optional argument
    parser = argparse.ArgumentParser(description='dedupe')
    parser.add_argument('-f', '--file', help = 'Specify sorted SAM file as input file -f', required= True)
    parser.add_argument('-u', '--umi', help = 'Specify known UMIs file as input', required= True)
    parser.add_argument('-o', '--outfile', help = 'Specify output SAM file for unique reads (de-duplicated) in .sam format', required= True)
    parse = parser.parse_args()
    return parse
#allows method calling
args = get_args()

#open all input files in parallel all at once
samfile = open(args.file, 'r')
umifile = open(args.umi, 'r')
deduped = open(args.outfile, 'w')
#these outputs are optional but next time wouldn't hard code, would include in argparse
invalid_umi_file = open('invalid_UMI.sam', 'w')
pcr_duplicates = open('pcr_duplicates.sam', 'w')

#make empty set of SAM file columns
unique_reads = set()

#make empty set of known UMIs
UMIs = set()
#populate UMIs set with known UMIs
for line in umifile:
  umi = line.strip('\n')
  UMIs.add(umi)

#DO NOT NEED THIS FUNCTION but I coded it so... commenting out just in case
# #grab relevant SAM file columns
# def get_cols(samfileline):
#   """This function takes each line of a SAM file splits it by column to capture relevant information 
#   for determining whether the read is a unique read or PCR duplicate(QNAME (UMI), FLAG(bitwise flag) RNAME (chr), 
#   POS (position mapped on reference), CIGAR (cigar string for parsing to determine true position) """
#   for line in samfile:
#     if '@' not in line:
#       line = samfileline.readline().split('\n').strip('\t').strip(':')
#       #grab 1st column = qname or UMI
#       qname = line[0]
#       #grab 2nd column = bitwise flag
#       bitflag = line[1]
#       #grab 3rd column = rname or chr or scaffold
#       chr_rname = line[2]
#       #grab 4th column = position on reference
#       pos = int(line[3])
#       #grab 6th column = cigar string
#       cigar = line[5]
#       return [qname, bitflag, chr_rname, pos, cigar]

#check bitwise flag for strand: either plus strand or minus strand
def whichstrand(bitflag):
    '''This function uses the bitwise flag to check which strand the read is on and returns either plus or minus'''
    if (int(bitflag) & 16 == 16):
      strand = 'minus'
    else:
      strand = 'plus'
    return strand

# CIGAR STRING INFO and assumptions
# M = match
# I = insertion to ref
# D = deletion from ref
# N = skipped from ref
# S = soft clipping


#parse cigar string left
def parsecigar_adjustpos(cigar, strand, samposition):
    """This function parses the cigar string to adjust the position of the read to the reference
    based on a these factors: Soft Clipping (left and right), Minus strand, Match, Insertions, Deletions"""
    #check strand
    #initialize integers at 0
    matched = 0
    dels = 0
    skipped = 0
    softclip = 0
    #check if it is the plus strand
    if strand == 'plus':
      #if it is a  plus strand and has an S on the left side (indicating soft clipping left)
      if 'S' in cigar:
        S = re.findall(r'^([0-9]+)S+', cigar)
        #check if the list has anything in it and if it does make the integer "softclip" = to the the value of the soft clipping
        if len(S) > 0:
          softclip = int(S[0])
        #if those above requirement are not met, set it to 0 so that the position isn't adjusted inappropriately
        else:
          softclip = 0
        #if there is a case of left soft clipping on the plus strand then adjust the 5' start position to subtract by the value of the soft clipping
        adjusted_pos = samposition - softclip
        return adjusted_pos
    #if it is the minus strand
    else:
    #check for right soft clipping, if there is right soft clipping set the integer "softclip" to that value, if not set to 0
      if 'S' in cigar[-1]:
        S = re.findall(r'([0-9]+)S', cigar)
        if len(S) > 0:
          softclip = int(S[0])
        else:
          softclip = 0
      #Check for Matches, Deletions, and Insertions, if there any of these add the value of these to their respective integer
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
    #find the sum of all of these cases if they exist (otherwise 0) and adjust the position to be the true position to include these cases
    adjusted_pos = samposition + matched + dels + skipped + softclip
    return adjusted_pos


#counters for records
not_dupe = 0
pcr_dupe = 0
invalidUMI = 0
total = 0
header = 0

#indices so we know what element of the samfile columns to call
qname = 0
bitflag = 1
rname = 2
pos = 3
cigarstr = 5
#index from samcols after splitting by ":" to get the UMI only
umi_index = 7

#chromosomes seen to track such that we can clear the unique_read set based on whether that chr has been seen or not
chr_seen = ''

#dedupe
while True:
  #get each line
  line = samfile.readline().strip('\n')
  #if no string in line stop or it is the end of the file and exit
  if line == "":
    break
  #grab the header lines and count, after that continue to non-header lines ie mapped reads
  if line.startswith('@'):
    deduped.write(f'{line}\n')
    header += 1
    continue
#split the columns of the sam file, and count the numnber of total lines
  samcols = line.split('\t')
  total += 1
  #grab the UMI of the mapped read by splitting this chunk by the colon and getting the last thing in the string which is the UMI
  umi = samcols[qname].split(':')[umi_index] 
  #check if the UMI is a known UMI that is in our set "UMIs", if not known count and write out to a file
  if umi not in UMIs:
    invalidUMI += 1
    invalid_umi_file.write(f'{line}\n')
  #grab the chromosome number of the read and check if it has not been seen, if it hasn't already been seen clear the unique_reads set to save memory
  else:
    chr = samcols[rname]
    if chr != chr_seen:
      chr_seen = chr
      unique_reads.clear()
    #if the read has a valid UMI, check the strand of the read to adjust the position of the read to the reference
    else:
      strand = whichstrand(samcols[bitflag])
      #capture the cigar string from the sam file columns
      cigar = samcols[cigarstr]
      #adjust the position based on strand (plus or minus), cigar string parsing, and left-most-mapped position as seen in the original sam file
      position = parsecigar_adjustpos(cigar, strand, int(samcols[pos]))
      #add only UNIQUE records to the unique_reads set, increment the count of unique_reads, and write out to a new file of only unique reads
      if (umi, strand, position) not in unique_reads:
        unique_reads.add((umi, strand, position))
        not_dupe += 1
        deduped.write(f'{line}\n')
      #otherwise increment the counter for pcr duplicates and write out to a file
      else:
        pcr_dupe += 1
        pcr_duplicates.write(f'{line}\n')


#close files
umifile.close()
samfile.close()
deduped.close()
invalid_umi_file.close()
pcr_duplicates.close()

#Deduper Summary
print(f'Number of Header Line: {header}') #test file expected:24
print(f'Total Records: {total}') #test file expected: 76
print(f'Number of Unique Reads: {not_dupe}') #test file expected:54
print(f'Number of Duplicates: {pcr_dupe}') #test file expected:21
print(f'Number of Invalid UMIs found: {invalidUMI}') #test file expected:1

#command to get stats per chromosome on C1_SE_uniqAlign.sorted.sam
# cat deduped_C1_SE_uniqAlign.sorted.sam | grep -v '^@' | awk '{print $3}' | uniq -c | sort -k2 -V