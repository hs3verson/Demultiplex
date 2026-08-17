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

plt.bar(x, y, color="cornflowerblue")
plt.title("Mean Quality Score at Each Position")
plt.xlabel("Position Number")
plt.ylabel("Quality Score")
plt.savefig(f'')
```

### [08-01-2026]

Wrote a bash script, `dist.sh`, to run `dist_per_n.py`. Ran into some issues with the input FASTQ files being zipped, after some research, I found I needed to `import gzip`, use `with gzip.open` instead of `with open`, and use "rt" instead of just "r" to be able to read the file as text and not bits. I also realized I would need separate scripts for the barcodes and the biological reads since the barcodes are only 8 nt long, I changed the plot colors for these to distinguish the plots more (cornflower blue for bio reads and light coral for barcodes). I also added an argument to argparse to be able to create better output png titles and plot titles.

To run the scripts simultaneously I ended up creating four different scripts for the four read files: `R1_dist.sh`, `R2_dist.sh`, `R3_dist.sh`, and `R4_dist.sh`
Wrapped `barcodes_dist.py` for R2 and R3 bash scripts and `bioreads_dist.py` for R1 and R4 bash scripts.

`bioreads_dist.py`
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

plt.bar(x, y, color="cornflowerblue")
plt.title(f"Mean Quality Score at Each Position - R{args.read_number}")
plt.xlabel("Position Number")
plt.ylabel("Quality Score")
plt.savefig(f'R{args.read_number}_dist.png')
```
`barcodes_dist.py`
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
    with 8 values of 0.0.'''
    while len(lst) < 8:
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

x = range(8)
y = dist_list

plt.bar(x, y, color="lightcoral")
plt.title(f"Mean Quality Score at Each Position - R{args.read_number}")
plt.xlabel("Position Number")
plt.ylabel("Quality Score")
plt.savefig(f'R{args.read_number}_dist.png')
```

**Scripts run:**

`R1_dist.sh`
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
`R2_dist.sh`
```
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=dist_R2

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v python dist_per_n.py -f $R2 -r 2
```
`R3_dist.sh`
```
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=dist_R3

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v python dist_per_n.py -f $R3 -r 3
```
`R4_dist.sh`
```
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=dist_R4

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v python dist_per_n.py -f $R4 -r 4
```

**Commands run:**
```bash
$ sbatch R1_dist.sh --output=Assignment-the-first/slurm-45933513.out
```
```bash
$ sbatch R2_dist.sh --output=Assignment-the-first/slurm-45933485.out
```
```bash
$ sbatch R3_dist.sh --output=Assignment-the-first/slurm-45933486.out
```
```bash
$ sbatch R4_dist.sh --output=Assignment-the-first/slurm-45933514.out
```

**Job resource usage (`/usr/bin/time -v` summary from Talapas):**

`R1_dist.sh`
```
Command being timed: "python bioreads_dist.py -f /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R1_001.fastq.gz -r 1"
Elapsed (wall clock) time (h:mm:ss or m:ss): 41:00.34
Maximum resident set size (kbytes): 71468
Percent of CPU this job got: 99%
Exit status: 0
```
`R2_dist.sh` 
```
Command being timed: "python barcodes_dist.py -f /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R2_001.fastq.gz -r 2"
Elapsed (wall clock) time (h:mm:ss or m:ss): 6:29.76
Maximum resident set size (kbytes): 70548
Percent of CPU this job got: 99%
Exit status: 0
```
`R3_dist.sh` 
```
Command being timed: "python barcodes_dist.py -f /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R3_001.fastq.gz -r 3"
Elapsed (wall clock) time (h:mm:ss or m:ss): 6:30.68
Maximum resident set size (kbytes): 70592
Percent of CPU this job got: 99%
Exit status: 0
```
`R4_dist.sh` 
```
Command being timed: "python bioreads_dist.py -f /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R4_001.fastq.gz -r 4"
Elapsed (wall clock) time (h:mm:ss or m:ss): 41:15.45
Maximum resident set size (kbytes): 74596
Percent of CPU this job got: 99%
Exit status: 0
```

### [08-05-2026]

Wrote a python script, `Assignment-the-third/demux.py`, to demultiplex the reads based on the logic from my pseudocode. Created a bash wrapper, `Assignment-the-third/demux.sh`, to run the python script.

**Scripts run:**
`demux.py`
```
#!/usr/bin/env python

import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to modify input four FASTQ files")
    parser.add_argument("-f1", "--R1", help="Desired R1 file", type=str)
    parser.add_argument("-f2", "--R2", help="Desired R2 file", type=str)
    parser.add_argument("-f3", "--R3", help="Desired R3 file", type=str)
    parser.add_argument("-f4", "--R4", help="Desired R4 file", type=str)
    parser.add_argument("-p", "--file_path", help="Desired file path for output files", type=str)
    return parser.parse_args()

args = get_args()

def read_record(fastq) -> list:
    '''Takes a fastq file and saves each record as a list, with each line as a separate item in the list'''
    record = []
    for i in range(4):
        record.append(fastq.readline().strip())
    return(record)

comp_dict = {"A":"t","T":"a","G":"c","C":"g", "N":"n"}
def rev_comp(seq: str) -> str:
    '''Takes a sequence string, finds the complementary sequence, and reverses it'''
    seq = seq.upper()
    new_seq = ""
    for base in seq:
        new_seq += comp_dict[base]
    new_seq = new_seq.upper()
    new_seq = new_seq[::-1]
    return new_seq

# create set of the 24 barcodes
barcodes_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
barcodes_set = set()

with open(barcodes_file, "r") as fh:
    for line in fh:
        if not line.startswith("sample"):
            column = line.strip().split("\t")
            bar_seq = column[4]
            barcodes_set.add(bar_seq)


# loop to create output files to be written to
files_dict = {}
for barcode in barcodes_set:
    files_dict[barcode] = [open(f"{args.file_path}/{barcode}_R1.fq", "w"), open(f"{args.file_path}/{barcode}_R2.fq", "w")]

unk_R1 = open(f"{args.file_path}/unk_R1.fq", "w")
unk_R2 = open(f"{args.file_path}/unk_R2.fq", "w")
hopped_R1 = open(f"{args.file_path}/hopped_R1.fq", "w")
hopped_R2 = open(f"{args.file_path}/hopped_R2.fq", "w")


def demux(fastq_R1, fastq_R2, fastq_R3, fastq_R4):
    '''Takes four fastq files: Read 1, Index 1, Index 2, and Read 2, sorts the reads to matched, hopped, and unknown files based 
    on the indexes, and prints out counts and percentages'''
    # initialize matched and hopped ditionaries to keep counts
    matched_dict = {}
    hopped_dict = {}
    # initialize counter for records sent to unknown
    unk_count = 0
    with gzip.open(fastq_R1, "rt") as fwdfile, gzip.open(fastq_R4, "rt") as revfile, gzip.open(fastq_R2, "rt") as index1, gzip.open(fastq_R3, "rt") as index2:
        record_count = 0                  # initialize record count
        while True:
            # create lists using read record fxn for each file
            R1 = read_record(fwdfile)
            R4 = read_record(revfile)
            R2 = read_record(index1)
            R3 = read_record(index2)
            fwd_bc = R2[1]                # define fwd_bc as the sequence line from the R2 record list
            rev_bc = rev_comp(R3[1])      # define rev_bc as the reverse complement of the sequence line from the R3 record list
            if R1[0] == '':
                break                     # break out of loop once R1 is empty 
            # sort barcodes not in set (ones with N's or otherwise) to unknown
            if (fwd_bc not in barcodes_set) or (rev_bc not in barcodes_set):
                unk_count += 1
                unk_R1.write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                unk_R2.write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            # sort matched barcodes to matched
            elif fwd_bc == rev_bc:
                bc_pair = f'{fwd_bc}-{rev_bc}'
                if bc_pair not in matched_dict:
                    matched_dict[bc_pair] = 1
                else:
                    matched_dict[bc_pair] += 1
                files_dict[f'{fwd_bc}'][0].write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                files_dict[f'{fwd_bc}'][1].write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            # sort hopped barcodes to hopped
            else:
                bc_pair = f'{fwd_bc}-{rev_bc}'
                if bc_pair not in hopped_dict:
                    hopped_dict[bc_pair] = 1
                else:
                    hopped_dict[bc_pair] += 1
                hopped_R1.write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                hopped_R2.write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            record_count += 1      # increment record count by 1
    # loop over matched dictionary to print counts and percentages for each of the 24 barcodes
    for bc_pair in matched_dict:
        percentage = round(((matched_dict[bc_pair])/record_count)*100, 2)
        print(f'{bc_pair}: total count = {matched_dict[bc_pair]}, percentage = {percentage}%')
    print('\n')
    # loop over hopped dictionary to print counts and percentages for each hopped combination seen
    for bc_pair in hopped_dict:
        percentage = round(((hopped_dict[bc_pair])/record_count)*100, 2)
        print(f'{bc_pair}: total count = {hopped_dict[bc_pair]}, percentage = {percentage}%')
    # print the number of reads sent to the unknown file
    print(f'\nThe number of unknown reads = {unk_count}')


# R1 = "../TEST-input_FASTQ/R1_test.fq.gz"
# R2 = "../TEST-input_FASTQ/R2_test.fq.gz"
# R3 = "../TEST-input_FASTQ/R3_test.fq.gz"
# R4 = "../TEST-input_FASTQ/R4_test.fq.gz"

demux(args.R1, args.R2, args.R3, args.R4)
```
`demux.sh`
```
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=demux

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

OFP=/scratch/bgmp/hvsev/demux

/usr/bin/time -v python demux.py -f1 $R1 -f2 $R2 -f3 $R3 -f4 $R4 -p $OFP
```

**Commands run:**
```bash
$ sbatch demux.sh --output=Assignment-the-third/slurm-46005207.out
```

**Job resource usage**

`demux.sh`
```
Command being timed: "python demux.py -f1 /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R1_001.fastq.gz -f2 /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R2_001.fastq.gz -f3 /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R3_001.fastq.gz -f4 /projects/bgmp/shared/2017_sequencing//1294_S1_L008_R4_001.fastq.gz -p /scratch/bgmp/hvsev/demux"
Elapsed (wall clock) time (h:mm:ss or m:ss): 42:12.53
Maximum resident set size (kbytes): 246488
Percent of CPU this job got: 62%
Exit status: 0
```

Based on this summary I learned I did not need to be using 8 CPU's for this job. 

Lastly, I created a markdown file, `Assignment-the-third/answers.md`, to put all the data on counts and percentages I had printed to standard out.

