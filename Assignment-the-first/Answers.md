# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. <img width="1196" alt="Index1" src="https://user-images.githubusercontent.com/83670690/127715600-65a5b341-02e3-48fe-8286-caf011fa8ba3.png">
    3. <img width="1195" alt="Index2" src="https://user-images.githubusercontent.com/83670690/127715631-ea2ba60a-2f59-4824-9b12-13e968764779.png">
    4. <img width="1198" alt="Read1" src="https://user-images.githubusercontent.com/83670690/127721430-8f5892a9-f420-49ca-8915-58dfafa1be91.png">
    5. <img width="1197" alt="Read2" src="https://user-images.githubusercontent.com/83670690/127721785-3450f960-fb58-4460-aa7b-ed1476cf1d03.png">




    b. This data will be used for differential gene expression meaning that there is already an assembled genome. For the biological reads, I would use a cutoff    quality score of 20. This means that on average 1 base for every 100 will be incorrect. Our reads are 101 base pairs long so on average only 1 base will be incorrect. Even if 2 bases are wrong, the read should still align to the genome. The worst case is that something will not align so this is not a big problem. On the other hand, the indexes need to be very accurate so I would pick a cutoff score of 30, which means 1 base out of 1000 will be incorrect on average. If the indexes are sequenced incorrectly and the read is assigned to the wrong sample, that would bias your downstream results. The consequence of the indexes being wrong is much higher so that is why they should have a higher cutoff value. 

    c. 
    ```
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | grep -F "N" | wc -l
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep -F "N" | wc -l
        3976613 indexes with an N (Index 1)
        3328051 indexes with an N (Index 2)
    ```
    
## Part 2
1. Define the problem: The problem is that all the reads are in the same file but we want them seperated by index. The indexes are not attached to the reads, they are in a different file. Additionally, some of the indexes may be swapped. In this data, all the reads are dual indexed which means that each read should have the same barcode as the paired read. If the indexes don't match, that means the indexes were swapped. Additionally, not all the reads and indexes have great quality and we want to do some preliminary data cleaning now. 

2. Describe output: The total output should have 52 files. For each read, there should 1 file with the swapped index reads and 1 file with non-existent indexes and bad quality indexes. There should 4 total files, 2 for read 1 and 2 for read 2. There will also be 24 files for each read, one for each index. These files will contain all the reads with matching, and good quality indexes. They will also have the indexes attached to the header of the record. There will be 24 files for read 1 and 24 files for read 2, for a total of 48 files. There should be 52 total files.  

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
Uploaded.

4. Pseudocode

```
Create a set of all the 24 indexes (Faster look up)
Create a set of all the input files (Read 1, Read 2, Index 1, Index 2)

For loop for the input files to read
    open the files

For loop for the indexes
    create a read1 and read 2 file for index
    open both files for writing

While there is a record in the input files
    read in read1, read2, index1 and index 2 record --> store in list (need indexes)
    combine index 1 and the reverse complement of index 2 with a dash between
    change header of index1 in list to header with added indexes
    change header of index2 in list to header with added indexes

    if there is an 'N' in either index
        write the read1 record to the read1-bad file
        write the read2 record to the read2-bad file
        continue back to the while loop

    else if the average quality score of index < index_cutoff OR qscore of sequence < seq_cutoff
        write the read1 record to the read1-bad file
        write the read2 record to the read2-bad file
        continue back to the while loop

    else if Index1 and Index2 are in set of indexes
        if index1 == index2
            write the read1 record to the read1-index file
            write the read2 record to the read2-index file
            continue back to the while loop

        else if index1 != index2 (it is swapped)
            write the read1 record to the read1-swapped file
            write the read2 record to the read2-swapped file
            continue back to the while loop

    else (anything else)
        write the read1 record to the read1-bad file
        write the read2 record to the read2-bad file

For loop to close all 52 files
```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

```
covert_phred(score: string) -> int:
    '''
    This function takes a ASCII value representing the Quality score and converts it to the decimal value - 33 (Phred 33)

    parameter: string score
    return: int decimal_score
    '''
    Input: A
    Expected output: 32

    return decimal_score
```    
```
average_qual(qual_line: string) -> float:
    '''
    This function take a quality string and returns the average quality score of that line

    parameter: string qual_line
    return: float avg_qual
    '''
    Input: AAAIII
    Expected output: 36

    return avg_qual
```
```
contain_N(index: string) -> int:
    '''
    This function takes a string of 2 indexes seperated with a '-' and returns the number of Ns it contains

    parameter: string index
    return: int num_ns
    '''
    Input: AACTTGCC-ANCTNGCC
    Expected output: 2

    return num_ns
```
```
if_match(index: string) -> bool:
    '''
    This function takes a string of 2 indexes seperated with a '-' and returns True if they match and False is they are different 

    parameter: string index
    return: bool match_flag
    '''
    Input: AACTTGCC-AACTTGCC
    Expected output: True

    return match_flag
```
```
reverse_comp(index: string) -> string
    '''
    This function takes a string of 1 index and returns a string of the reverse complement 

    parameter: string index
    return: string rev_comp
    '''
    Input: AACTTGCC
    Expected output: GGCAAGTT

    return rev_comp
```
```
if_real_index(index: string) -> bool:
    '''
    This function takes a string of 2 indexes seperated with a '-' and returns True if are both in the index set and False is they are not in the index set

    parameter: string index
    return: bool real_flag
    '''
    Input: AACTTGCC-AACTTGCC
    Expected output: False

    return real_flag
```
