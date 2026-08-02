# Lab Notebook — Demulitplex

**Base Directory**

`/projects/bgmp/hvsev/bioinfo/Bi622/Demultiplex` (within Talapas)

**Environment / Versions:**
- Compute environment:
  - bgmp compute nodes 

- Software/package versions:
  
  `bash 4.4.20`

  `python 3.14.6`

**Data Source:**

Sequencing data generated from the 2017 BGMP cohort's library preps
```
DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz
```

---

### [07-24-2026] 

Wrote pseudocode to help define the problem and understand the goal. See `Assignment-the-first/pseudo.md`

### [07-30-2026]

Initial data exploration!

Biological reads = R1, R4
index reads = R2, R3

Length of reads: 
Used following bash commands to determine how long the reads in each file are

`zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc` -> 102 -> 101 w/o new line char

`zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc` -> 9 -> 8 w/o new line char

`zcat 1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc` -> 9 -> 8 w/o new line char

`zcat 1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc` -> 102 -> 101 w/o new line char

Investigated quality score lines to determine phred encoding -- found # linked to N's in the sequence line, indicative of phred-33 which includes # (equal to a qual score of 2)

### [07-31-2026]

Wrote a python script, `dist_per_n.py`, to create distribution histograms for each of the read files.

```
#!/bin/bash/env python

import matplotlib.pyplot as plt
import argparse
import bioinfo

def get_args():
    parser = argparse.ArgumentParser(description="Program to modify file name")
    parser.add_argument("-f", "--file_name", help="Desired filename", type=str)
    parser.add_argument("-r", "--read_number", help="Read number of file", type=int)
    return parser.parse_args()

args = get_args()


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    while len(lst) < 101:
        lst.append(value)
    return lst
        
dist_list: list = []
dist_list = init_list(dist_list)

def populate_list(file: str) -> tuple[list, int]:
    """This function takes a FASTQ file, counts the number of lines, and filters for the quality score line. All of the values at index i 
    in each quality score line are then summed and added to an empty list."""
    sum_list = init_list([])      #create list with 101 zeros
    with open(file, "r") as fh:
        phred_letter=0      #initialize phred letter to start at character zero in line
        num_lines=0         #initialize line counter
        for line in fh:     #iterate through the lines in the fastq file
            line = line.strip()      #strips whitespace
            num_lines+=1             #add one to line counter
            if num_lines%4 == 0:     #pull out quality score lines
                for index, phred_letter in enumerate(line):      #iterate through the character(phred letter) at index i in each qual score line
                    sum_list[index] += bioinfo.convert_phred(phred_letter)      #at index i in the list of zeros, add the converted numerical phred value
    return sum_list, num_lines       #return the list of summed phred value scores at each index and the total number of lines in the fastq file

dist_list, num_lines = populate_list(args.file_name)

for index, total in enumerate(dist_list):      #iterate through the sums(total) in the numbered my_list
    mean= total/(num_lines/4)      #define mean as the total divided by the total number of lines divided by 4(accounting for fastq file record length)
    dist_list[index] = mean          #set the newfound mean as the value at that index in my_list

x = range(101)
y = dist_list

plt.scatter(x, y, color="m")
plt.title("Mean Quality Score at Each Position")
plt.xlabel("Position Number")
plt.ylabel("Quality Score")
plt.savefig(f'dist_R{args.read_number}.png')
```

### [08-01-2026]

Wrote a bash script, `dist.sh`, to run `dist_per_n.py`. Ran into some issues with the input FASTQ files being zipped, after some research, I found I needed to `import gzip`, use `with gzip.open` instead of `with open`, and use "rt" instead of just "r" to be able to read the file as text and not bits.

Updated `dist_per_n.py`
```
#!/bin/bash/env python

import matplotlib.pyplot as plt
import argparse
import bioinfo
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Program to modify file name")
    parser.add_argument("-f", "--file_name", help="Desired filename", type=str)
    parser.add_argument("-r", "--read_number", help="Read number of file", type=int)
    return parser.parse_args()

args = get_args()


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    while len(lst) < 101:
        lst.append(value)
    return lst
        
dist_list: list = []
dist_list = init_list(dist_list)

def populate_list(file: str) -> tuple[list, int]:
    """This function takes a FASTQ file, counts the number of lines, and filters for the quality score line. All of the values at index i 
    in each quality score line are then summed and added to an empty list."""
    sum_list = init_list([])      #create list with 101 zeros
    with gzip.open(file, "rt") as fh:
        phred_letter=0      #initialize phred letter to start at character zero in line
        num_lines=0         #initialize line counter
        for line in fh:     #iterate through the lines in the fastq file
            line = line.strip()      #strips whitespace
            num_lines+=1             #add one to line counter
            if num_lines%4 == 0:     #pull out quality score lines
                for index, phred_letter in enumerate(line):      #iterate through the character(phred letter) at index i in each qual score line
                    sum_list[index] += bioinfo.convert_phred(phred_letter)      #at index i in the list of zeros, add the converted numerical phred value
    return sum_list, num_lines       #return the list of summed phred value scores at each index and the total number of lines in the fastq file

dist_list, num_lines = populate_list(args.file_name)

for index, total in enumerate(dist_list):      #iterate through the sums(total) in the numbered my_list
    mean= total/(num_lines/4)      #define mean as the total divided by the total number of lines divided by 4(accounting for fastq file record length)
    dist_list[index] = mean          #set the newfound mean as the value at that index in my_list

x = range(101)
y = dist_list

plt.scatter(x, y, color="m")
plt.title("Mean Quality Score at Each Position")
plt.xlabel("Position Number")
plt.ylabel("Quality Score")
plt.savefig(f'dist_R{args.read_number}.png')
```

`dist.sh` -- modified -f and -r flag for each file
```
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=dist_R1

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v python dist_per_n.py -f $R1 -r 1
```

**Commands run:**
```bash
$ sbatch dist.sh
```

**Job resource usage (`/usr/bin/time -v` summary from Talapas):**

`dist.sh` -- dist_R1
```
Command being timed: 
Elapsed (wall clock) time (h:mm:ss or m:ss):
Maximum resident set size (kbytes):
Percent of CPU this job got:
Exit status:
```
`dist.sh` -- dist_R2
```
Command being timed: 
Elapsed (wall clock) time (h:mm:ss or m:ss):
Maximum resident set size (kbytes):
Percent of CPU this job got:
Exit status:
```
`dist.sh` -- dist_R3
```
Command being timed: 
Elapsed (wall clock) time (h:mm:ss or m:ss):
Maximum resident set size (kbytes):
Percent of CPU this job got:
Exit status:
```
`dist.sh` -- dist_R4
```
Command being timed: 
Elapsed (wall clock) time (h:mm:ss or m:ss):
Maximum resident set size (kbytes):
Percent of CPU this job got:
Exit status:
```
