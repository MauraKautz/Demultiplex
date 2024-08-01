#!/usr/bin/env python

import bioinfo
import argparse
import matplotlib.pyplot as plt
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="Name of input file", required=True, type=str)
    parser.add_argument("-l", help="Length of the nucleotides", required=True, type=int)
    parser.add_argument("-o", help="Changing histogram ouput file", required=True, type=str)
    return parser.parse_args()

args = get_args()
f=args.f
l=args.l
o=args.o

mean_qscores: list = [0]*l

with gzip.open(f, "rt") as fh:
    i=0
    for line in fh:
        i+=1
        line = line.strip('\n') 
        if i%4 == 0:
            for x in range(0,len(line)):
                score=bioinfo.convert_phred(line[x])
                mean_qscores[x]+=score
    print(mean_qscores)

print(f'# Base Pair\tMean Quality Score')
num_records = int(i/4)
for index, sum in enumerate(mean_qscores):
    mean = (sum / (num_records))
    mean_qscores[index]= mean
    print(index, mean, sep="\t")

x=range(0,l)
y=mean_qscores

plt.bar(x,y)

plt.xlabel("Number of Bases")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Value at Each Base")
plt.savefig(o)

