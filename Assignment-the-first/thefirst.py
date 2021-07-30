#!/usr/bin/env python

import argparse 
import gzip
import numpy as np

def convert_phred(letter):
    """Converts a single character into a phred score"""
    dec = ord(letter)
    return dec - 33

def get_args():
    parse = argparse.ArgumentParser(description="A program to find the distribution of bases")
    parse.add_argument("-f", "--filename", help="the name of the file", required=True)
    parse.add_argument("-l", "--length", help="the length of the reads or indexes", type=int, required=True)
    parse.add_argument("-o", "--output", help="the name of the output file", required=True)
    return parse.parse_args()

args = get_args()

#use a numpy array for faster run time
score_list = np.zeros(args.length, dtype=int)
# open the file with gzip and read each quality line and keep a running total of each base
with gzip.open(args.filename, 'rt') as fh:
    counter = 0
    while True:
        header = fh.readline()
        if header == '':
            break
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        #add the quality score for each base to the numpy array
        for i,score in enumerate(qual.strip()):
            score_list[i] += convert_phred(score)
        counter += 1

#divide the whole list by the number of records to get the mean quality score at each base
score_list = score_list/counter

#create a list of base pair sites
bp = []
for i in range(args.length):
    bp.append(i)
    
#create a bar plot for the average quality score at each site    
import matplotlib.pyplot as plt
plt.figure(figsize=(15,5))
plt.bar(bp, score_list, edgecolor='b')

plt.title("Distribution of Average Quality Scores: %s" % args.output)
plt.xlabel("Base Pair")
plt.ylabel("Average Quality Score")

plt.savefig('%s.png' % args.output)


