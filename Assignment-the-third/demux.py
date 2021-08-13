#!/usr/bin/env python

import demux_functions as df
import argparse
import gzip

# use argparse to get user input of the file names and the index and read cutoffs
def get_args():
    parse = argparse.ArgumentParser(description="A program to find the distribution of bases")
    parse.add_argument("-r1", "--read1", help="the name of the read 1 file", required=True)
    parse.add_argument("-i1", "--index1", help="the name of the index 1 file", required=True)
    parse.add_argument("-i2", "--index2", help="the name of the index 2 file", required=True)
    parse.add_argument("-r2", "--read2", help="the name of the read 2 file", required=True)
    parse.add_argument("-ic", "--index_cutoff", help="The quality score cuttoff for the index scores", type=int, required=True)
    parse.add_argument("-sc", "--seq_cutoff", help="The quality score cuttoff for the quality scores", type=int, required=True)
    return parse.parse_args()

args = get_args()

# Create a set of all the 24 indexes (Faster look up)
index_set = {'GTAGCGTA', 'CGATCGAT', 'GATCAAGG', 'AACAGCGA', 'TAGCCATG', 'CGGTAATC', 'CTCTGGAT', 'TACCGGAT', 'CTAGCTCA', 'CACTTCAC', 'GCTACTCT', 'ACGATCAG', 'TATGGCAC', 'TGTTCCGT', 'GTCCTAAG', 'TCGACAAG', 'TCTTCGAC', 'ATCATGCG', 'ATCGTGGT', 'TCGAGAGT', 'TCGGATTC', 'GATCTTGC', 'AGAGTCCA', 'AGGATAGC'}

# create a dict with read 1 and read 2 for each index as keys so the file handles will be the values
index_dict = {}
for index in index_set:
    index_dict[str(index + '_r1')] = ''
    index_dict[str(index + '_r2')] = ''

#open output files for each index read
for index in index_dict:
    read = '{}.fastq.gz'.format(index)
    index_dict[index] = (gzip.open(read, 'wt'))

input_list = [args.read1, args.index1, args.index2, args.read2]

# create a dict with the read number as the key and file handle as the values 
input_dict = {}
input_dict['read1'] = (gzip.open(input_list[0], 'rt'))
input_dict['index1'] = (gzip.open(input_list[1], 'rt'))
input_dict['index2'] = (gzip.open(input_list[2], 'rt'))
input_dict['read2'] = (gzip.open(input_list[3], 'rt'))

#create and open bad output files in a dict
bad_dict = {}
bad_dict['bad_r1'] = gzip.open('bad_read1.fastq.gz', 'wt')
bad_dict['bad_r2'] = gzip.open('bad_read2.fastq.gz', 'wt')
bad_dict['swap_r1'] = gzip.open('swapped_read1.fastq.gz', 'wt')
bad_dict['swap_r2'] = gzip.open('swapped_read2.fastq.gz', 'wt')


record = [] #list of all four record
one_rec = [] #each record

# flag for if input files are empty
flag = True
while flag:
    record = []
    for handle in input_dict.values():
        # if line is empty, set flag to false and close all the input and output files
        fr_line = handle.readline().strip()
        if fr_line == '':
            flag = False
            for value in index_dict.values():
                value.close()
            for value in input_dict.values():
                value.close()   
            for value in bad_dict.values():
                value.close()
            break

        else: 
            # create a list for the record from each input file
            one_rec.append(fr_line)
            one_rec.append(handle.readline().strip())
            one_rec.append(handle.readline().strip())
            one_rec.append(handle.readline().strip())
            #append the record to the list of records 
            record.append(one_rec)
            one_rec = []
    
    # there are records in the list, go through the logic 
    if len(record) != 0:

        #create and add the combined indexes to the heads lines of the reads
        index1 = record[1][1]
        index2 = df.reverse_comp(record[2][1])
        comb_header =  index1 + "-" + index2
        record[0][0] = str(record[0][0] + " " + comb_header)
        record[3][0] = str(record[3][0] + " " +  comb_header)

        #if there is an 'N' in either index
        if df.contain_N(comb_header) != 0:
            # write the read1 record to the read1-bad file
            for i in record[0]:
                bad_dict['bad_r1'].write(i + "\n")
            # write the read2 record to the read2-bad file
            for i in record[3]:
                bad_dict['bad_r2'].write(i + "\n")

        #the average quality score of index < index_cutoff OR qscore of sequence < seq_cutoff
        elif df.average_qual(record[1][3]) < args.index_cutoff or df.average_qual(record[2][3]) < args.index_cutoff or df.average_qual(record[0][3]) < args.seq_cutoff or df.average_qual(record[3][3]) < args.seq_cutoff:
            # write the read1 record to the read1-bad file
            for i in record[0]:
                bad_dict['bad_r1'].write(i + "\n")
            # write the read2 record to the read2-bad file
            for i in record[3]:
                bad_dict['bad_r2'].write(i + "\n")

        # if both the indexes are real
        elif index1 in index_set and index2 in index_set:
            if index1 == index2:
                # write the read1 record to the read1-index file
                i1 = str(index1 + "_r1")
                for i in record[0]:
                    index_dict[i1].write(i + "\n")
                
                # write the read2 record to the read2-index file
                i2 = str(index2 + "_r2")
                for i in record[3]:
                    index_dict[i2].write(i + "\n")

            elif index1 != index2: #it is swapped
                # write the read1 record to the read1-swapped file
                for i in record[0]:
                    bad_dict['swap_r1'].write(i + "\n")
                # write the read2 record to the read2-swapped file
                for i in record[0]:
                    bad_dict['swap_r2'].write(i + "\n")

        else:
            #anything else
            #write the read1 record to the read1-bad file
            for i in record[0]:
                bad_dict['bad_r1'].write(i + "\n")
                
            # write the read2 record to the read2-bad file
            for i in record[3]:
                bad_dict['bad_r2'].write(i + "\n")
