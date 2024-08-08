#!/usr/bin/env python
import bioinfo
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f1", help="Name of 1st input file", required=True, type=str)
    parser.add_argument("-f2", help="Name of 2nd input file", required=True, type=str)
    parser.add_argument("-f3", help="Name of 3rd input file", required=True, type=str)
    parser.add_argument("-f4", help="Name of 4th input file", required=True, type=str)
    parser.add_argument("-q", help="Quality score threshold", required=True, type=int)
    parser.add_argument("-m", help="Matched index txt file", required=True, type=str)
    return parser.parse_args()

args = get_args()
f1=args.f1
f2=args.f2
f3=args.f3
f4=args.f4
q=args.q
m=args.m

def rev_comp(sequence: str) -> str:
    '''Takes a DNA sequence (index) and returns the reverse complement'''

    #initializing
    pair = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    reverse_comp = ""

    #reversing a given string
    reverse = sequence[::-1]
    
    #giving the complement of the newly reversed string
    for base in reverse:
        reverse_comp += pair[base]

    #return the reverse complement of a string
    return reverse_comp

#sample test!
#print(rev_comp("CAT"))
#input: CAT
#output: ATG

def edit_header(R1_header: str, R4_header: str, R2_seq: str, R3_seq: str) -> tuple[str, str]:
    '''Takes a header, adds two indexes to the end, returns a new header'''

    #header equals header plus both indexes
    R1_header += " " + R2_seq + "-" + R3_seq
    R4_header += " " + R2_seq + "-" + R3_seq

    return R1_header, R4_header

#sample test!
#print(edit_header("Header1", "Header2", "AAAA", "TTTT"))
#input: Header1, Header2, AAAA, TTTT
#output: Header1 AAAA-TTTT, Header2 AAAA-TTTT

def check_q(R2_qscores: str, R3_qscores: str, qscore_thresh: int) -> bool:
    '''Takes the quality scores of a record in reads 2 and 3, averages them, then checks them against a quality score threshold.  Returns True if the score is above the threshold, False if it is below'''
    R2 = bioinfo.qual_score(R2_qscores)
    R3 = bioinfo.qual_score(R3_qscores)
    if R2 < qscore_thresh or R3 < qscore_thresh:
        return True
    else:
        return False

#sample test!
#print(check_q("#AJAAA#","JJJJE#A", 20))
#Input: #AJAAA#, JJJJE#A, 20
#Output: True

matched_indexes = set()
unknown = 0
matched_dict: dict = {}
hopped_dict: dict = {}
file_handles: dict = {}

#convert matched index text file (known indexes) to a set
with open(m, "r") as matched:
    for line in matched:
        line = line.strip('\n')
        line = line.split('\t')
        matched_indexes.add(line[1])
        matched_indexes.add(line[3])
        matched_indexes.add(line[5])

#open all 4 input files (reads 1-4)
with gzip.open(f1, "rt") as f1, gzip.open(f2, "rt") as f2, gzip.open(f3, "rt") as f3, gzip.open(f4, "rt") as f4:
    for index in matched_indexes:
        fh1 = open(f'{index}_R1.fastq', "a")
        fh2 = open(f'{index}_R2.fastq', "a")
        file_handles[index] = fh1,fh2

    
    fh1 = open("hopped_R1.fastq", "a")
    fh2 = open("hopped_R2.fastq", "a")
    file_handles["hopped"] = fh1,fh2
    fh3 = open("unknown_R1.fastq", "a")
    fh4 = open("unknown_R2.fastq", "a")
    file_handles["unknown"] = fh3, fh4

    
    while True:
    #read 4 lines (a record) from all files
        f1header, f1seq, f1plus, f1qscore = f1.readline().strip(), f1.readline().strip(), f1.readline().strip(), f1.readline().strip()
        f2header, f2seq, f2plus, f2qscore = f2.readline().strip(), f2.readline().strip(), f2.readline().strip(), f2.readline().strip()
        f3header, f3seq, f3plus, f3qscore = f3.readline().strip(), f3.readline().strip(), f3.readline().strip(), f3.readline().strip()
        f4header, f4seq, f4plus, f4qscore = f4.readline().strip(), f4.readline().strip(), f4.readline().strip(), f4.readline().strip()

        #if the line is empty, break the loop (reached end on file)
        if f1header == "":
            break

        #reverse complement R3 to match R2
        f3seq = rev_comp(f3seq) # type: ignore

        f2seq_f3seq = f'{f2seq}-{f3seq}'

        #edit the headers to add index pair to end of headers in both R1 and R4
        f1header, f4header = edit_header(f1header, f4header, f2seq, f3seq) # type: ignore
        
        #check if index 1 AND index 2 are in the set of 24 known indexes
        if f2seq not in matched_indexes or f3seq not in matched_indexes:
            unknown += 1

            unk_R1 = file_handles["unknown"][0]#add R1 to unk_R1.fq
            unk_R4 = file_handles["unknown"][1]#add R4 to unk_R4.fq
            unk_R1.write(f'{f1header}\n{f1seq}\n{f1plus}\n{f1qscore}\n')
            unk_R4.write(f'{f4header}\n{f4seq}\n{f4plus}\n{f4qscore}\n')
            continue

        elif check_q(f2qscore, f3qscore, q): # type: ignore
            unknown += 1

            unk_R1 = file_handles["unknown"][0]#add R1 to unk_R1.fq
            unk_R4 = file_handles["unknown"][1]#add R4 to unk_R4.fq
            unk_R1.write(f'{f1header}\n{f1seq}\n{f1plus}\n{f1qscore}\n')
            unk_R4.write(f'{f4header}\n{f4seq}\n{f4plus}\n{f4qscore}\n')
            continue

        elif f2seq != f3seq:
            if f2seq_f3seq in hopped_dict:
                hopped_dict[f2seq_f3seq] += 1
            else:
                hopped_dict[f2seq_f3seq] = 1

            hopped_R1 = file_handles["hopped"][0]
            hopped_R4 = file_handles["hopped"][1]
            hopped_R1.write(f'{f1header}\n{f1seq}\n{f1plus}\n{f1qscore}\n')#add R1 to hopped_R1.fq
            hopped_R4.write(f'{f4header}\n{f4seq}\n{f4plus}\n{f4qscore}\n')#add R4 to hopped_R4.fq
            continue

        else: #if its made it this far, then it is matched, passes quality score check, and known
            if f2seq_f3seq in matched_dict:
                matched_dict[f2seq_f3seq] += 1
            else:
                matched_dict[f2seq_f3seq] = 1

            matched_R1 = file_handles[f2seq][0]
            matched_R4 = file_handles[f2seq][1]
            matched_R1.write(f'{f1header}\n{f1seq}\n{f1plus}\n{f1qscore}\n')#add R1 to matched_R1.fq
            matched_R4.write(f'{f4header}\n{f4seq}\n{f4plus}\n{f4qscore}\n')#add R4 to matched_R4.fq

    for key in file_handles:
        R1 = file_handles[key][0]
        R2 = file_handles[key][1]
        R1.close()
        R2.close()

#Percentage of reads from each sample
#Overall amount of index swapping

hopped = 0
for key in hopped_dict:
    hopped += hopped_dict[key]

matched = 0
for key in matched_dict:
    matched += matched_dict[key]

total_reads = unknown + hopped + matched

print(matched_dict)
print(hopped_dict)

for key in matched_dict:
    print(f'{key}: {matched_dict[key]}\tPercent of reads from each sample: {matched_dict[key]/total_reads * 100}')

print(f'Amount of unknown reads: {unknown}\tPercent of unknown reads: {unknown/total_reads * 100}')
print(f'Amount of hopped reads: {hopped}\tPercent of hopped reads: {hopped/total_reads * 100}')
print(f'Amount of matched reads: {matched}\tPercent of matched reads: {matched/total_reads * 100}')
