#!/usr/bin/env python

# Author: Hannah Severson hvsev@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = ['A','G','T','C']
RNA_bases = ['A','G','U','C']

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33      # converts letter to associated ASCII score and subtracts 33

def qual_score(phred_score: str) -> float:
    '''Takes a phred score line(string) and returns the average numerical
    quality score for the line.'''
    total = 0      # initialize total sum
    for char in phred_score:         # iterates over each character in phred score line
        score = convert_phred(char)  # use convert_phred fxn to find the associated value
        total += score               # add that score to the total sum
    average = total / len(phred_score)     # find the average by dividing the total by the length of the line
    return average

def validate_base_seq(seq: str, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNA_bases = set('ATGCatcg')      # define the DNA bases as A, T, C, and G --case insensitive
    RNA_bases = set('AUGCaucg')      # define the RNA bases as A, U, C, and G --case insensitive
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str, RNAflag=False):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq, RNAflag), "String contains invalid characters - are you sure you used a DNA/RNA sequence?"
    seq = seq.upper()        # converts sequence to all upper case
    Gs = seq.count("G")      # defines Gs as the number of Gs in the sequence
    Cs = seq.count("C")      # defines Cs as the number of Cs in the sequence
    return (Gs+Cs)/len(seq)  # adds the G and C counts and divdes by the length of the sequence

def calc_median(lst: list):
    '''Given a sorted list, returns the median value of the list'''
    lst.sort()      #sort list numerically
    list_len = len(lst)      #define list_len as the length of the list
    # if list length is an even number
    if list_len%2==0:      
        median = (lst[list_len//2 - 1] + lst[list_len//2]) / 2
    # if list length is an odd number
    else:
        median = lst[list_len//2]
    return median

def oneline_fasta(fasta_file, fasta_out):
    '''Given a FASTA file, returns a new FASTA file with the sequence line as one line.'''
    with open(fasta_file, "r") as f_in, open(fasta_out, "w") as f_out:      # opens input FASTA file and creates output FASTA file 
        for i, line in enumerate(f_in):      # iterates over each line in the FASTA file
            line = line.strip()              # strips line
            if i == 0:                       # for the first line in the file
                f_out.write(f'{line}\n')     # write out the line with a new line character at the end
            elif line.startswith(">"):       # for every header line after the first
                f_out.write(f'\n{line}\n')   # move to a new line, then write out the line with a new line character at the end
            else:                            # for sequence lines
                f_out.write(f'{line}')       # write out the line

def fastq_record(fastq_file):
    '''Takes a FASTQ file and stores records as a list where each item is a line'''
    with open(fastq_file, "r") as fh:
        while True:
            record_lst = [fh.readline().strip('\n') for i in range(4)]
            if record_lst[0]=="":
                break   


if __name__ == "__main__":
    # tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    # tests for qual_score
    assert qual_score("CCCFFF") == 35.5, "wrong avg qual score for 'CCCFFFF'"
    assert qual_score("()*&%$") == 6.0, "wrong avg qual score for '()*&%$'"
    assert qual_score("IJJIG") == 40.0, "wrong avg qual score for 'IJJIGJ'"
    print("Your qual_score function is working! Nice job")
    # tests for validate_base_seq
    assert validate_base_seq("AAUAGAU", True), "RNA test failed"
    assert validate_base_seq("AATAGAT"), "DNA test failed"
    assert validate_base_seq("Apples and bananas!")==False, "Non-nucleic test failed"
    print("Your validate_base_seq function is working! Nice job")
    # tests for gc_content
    assert gc_content("GGGCCC") == 1.0
    assert gc_content("UUUUU", RNAflag=True) == 0.0
    assert gc_content("GCGCGCAT") == 0.75
    print("Your gc_content function is working! Nice job")
    # tests for calc_median
    assert calc_median([3,24,47]) == 24, "calc_median function does not work for odd length list"
    assert calc_median([30,40]) == 35, "calc_median function does not work for even length list"
    assert calc_median([7,5,4,8,6]) == 6, "calc_median function does not work for on unsorted list"
    assert calc_median([12,13,14,15,16,17]) == 14.5
    print("Your calc_median function is working! Nice job")