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

def read_record(fastq):
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
    # initialize matched and hopped ditionaries to keep counts
    matched_dict = {}
    hopped_dict = {}
    # initialize counter for records sent to unknown
    unk_count = 0
    with gzip.open(fastq_R1, "rt") as fwdfile, gzip.open(fastq_R4, "rt") as revfile, gzip.open(fastq_R2, "rt") as index1, gzip.open(fastq_R3, "rt") as index2:
        record_count = 0
        while True:
            R1 = read_record(fwdfile)
            R4 = read_record(revfile)
            R2 = read_record(index1)
            R3 = read_record(index2)
            fwd_bc = R2[1]
            rev_bc = rev_comp(R3[1])
            if R1[3] == '':
                break
            if (fwd_bc not in barcodes_set) or (rev_bc not in barcodes_set):
                unk_count += 1
                unk_R1.write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                unk_R2.write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            elif fwd_bc == rev_bc:
                bc_pair = f'{fwd_bc}-{rev_bc}'
                if bc_pair not in matched_dict:
                    matched_dict[bc_pair] = 1
                else:
                    matched_dict[bc_pair] += 1
                files_dict[f'{fwd_bc}'][0].write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                files_dict[f'{fwd_bc}'][1].write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            else:
                bc_pair = f'{fwd_bc}-{rev_bc}'
                if bc_pair not in hopped_dict:
                    hopped_dict[bc_pair] = 1
                else:
                    hopped_dict[bc_pair] += 1
                hopped_R1.write(f"{R1[0]} {fwd_bc}-{rev_bc}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
                hopped_R2.write(f"{R4[0]} {fwd_bc}-{rev_bc}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
            record_count += 1
    for bc_pair in matched_dict:
        percentage = round(((matched_dict[bc_pair])/record_count)*100, 2)
        print(f'{bc_pair}: total count = {matched_dict[bc_pair]}, percentage = {percentage}%')
    print('\n')
    for bc_pair in hopped_dict:
        percentage = round(((hopped_dict[bc_pair])/record_count)*100, 2)
        print(f'{bc_pair}: total count = {hopped_dict[bc_pair]}, percentage = {percentage}%')
    print(f'\nThe number of unknown reads = {unk_count}')


# R1 = "../TEST-input_FASTQ/R1_test.fq.gz"
# R2 = "../TEST-input_FASTQ/R2_test.fq.gz"
# R3 = "../TEST-input_FASTQ/R3_test.fq.gz"
# R4 = "../TEST-input_FASTQ/R4_test.fq.gz"

demux(args.R1, args.R2, args.R3, args.R4)
