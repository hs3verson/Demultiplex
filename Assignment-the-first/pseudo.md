```
shebang

open indexes file
    loop over it
        save indexes to a set

initialize matched dictionary
initialize hopped dictionary
open four read files (R1, R2, R3, R4)
    readline sequence line for R1, R2, R3, R4
        if R2 or R3 barcodes contain N's 
            write R1 record R1 unknown file and append index to header
            +1 to unknown counter
            write R4 record R4 unknown file and append index to header
        elif R2 and R3 barcodes match (reverse complement fxn)
            if R2 and R3 reversed are in set
                write R1 record R1 matched file and append index to header
                add to matched dictionary -- R2/R3 is key and value is +1
                write R4 record R4 matched file and append index to header
            elif R2 and R3 reversed not in set
                write R1 record R1 unknown file and append index to header
                +1 to unknown counter
                write R4 record R4 unknown file and append index to header
        elif R2 and R3 barcodes don't match
            if R2 and R3 reversed are in set
                write R1 record R1 hopped file and append index to header
                add to hopped dictionary -- R2/R3 is key and value is +1
                write R4 record R4 hopped file and append index to header
        else
            write R1 record R1 unknown file and append index to header
            +1 to unknown counter
            write R4 record R4 unknown file and append index to header

def rev_comp(seq: str) -> str:
    '''Takes a sequence string, finds the complementary sequence, and reverses it'''
    seq = seq.replace("A", "u").replace(
            "C", "g").replace("U", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq
Input: ACTGTG
Expected: CACAGT
```