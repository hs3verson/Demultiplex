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