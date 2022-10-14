- Define the problem

We need to differentiate our reads from PCR duplication artifacts vs unique transcript reads to accurately reflect transcript abundace in our samples. 
The problem with PCR duplicates is that they will have the same alignment position on the same chromosome on the same strand with the same start/stop, and UMI. 
We need to find a way to remove these duplicates because these are identical molecules that are from the same transcript. 
These duplicates can introduce bias into our data and artificially inflate the coverage of one transcript over the other. 
Another challenge will be having to adjust the position of the aligned read to 
the reference based on the cigar string.

- Write example test cases:

- Potential test cases to include in the sam file test input

    - same UMIs duplicate on plus and minus
    - same UMIs and everything else soft slipping plus and minus and duplicate
    - reads with the different UMIs with the same pos, strand, and chromosome
    - different chr same everything else and
    - same chr but everything else different
    - different pos same everything else
    - same pos different everything else
    - different strand same everything else
    - same strand different everything else


- Develop your algorithm using pseudocode

Pseudocode

```

Read in sorted SAM File.
Loop over each line.
Save each column of interest (QNAME, RNAME, POS, FLAG) as a variable.
Use the bitwise flag mapped/unmapped function to determine if the read was mapped to the genome.
Use the cigar string parsing function and pass it the cigar string and bitwise flag strand function to determine the strand that the read is on/mapped to 
and adjust the position of the read on the reference to be the "true" position
accounting for soft clipping, skipped regions, deletions, insertions, and matches and mismatches.
For each strand that is unique add the QNAME (UMI), RNAME (chr#), Adjusted POSITION, FLAG to a tuple and add that tuple to a set.
For each strand that is not unique, we will discard the read interpreting it as a PCR duplicate.
Lastly, we'll write out the contents of our set to a new deduped sam file.

```

- Determine high level functions
    - Description
    - Function headers
    - Test examples (for individual functions)
    - Return statement

```
Bitwise Flag Functions

bitmapped(float)
'''This function check if the read is mapped to the reference'''
if bit 4 = 0 read is mapped
    move onto next step
if bit 4 = 1 read is unmapped
    go to next line

Test case:

bitstrand(float)
'''This function checks if ((flag & 16) == 16 : rev_comp = True minus strand'''

```

Counting CIGAR String, this deduper will be using inclusive numbering for the basepair start position

``` 
Parse CIGAR string function

def cigstr(str, bitstrand(float))
'''This function takes the CIGAR String and accounts for the strand from a SAM file and parses the symbols 
M, I, D, N, and S to determine the actual 5' start position on the reference of our aligned read'''

Test case:
minus strand, input "5S10M1024N" output "1034" 

```


