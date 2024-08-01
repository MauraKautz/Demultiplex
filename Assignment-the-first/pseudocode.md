Problem, Input, and Output
Parse through four fastq files (R1:forward read, R2:index 1, R3:index 2, R4:reverse read) and move records to appropriate outputs based on indexes:
    If indexes match each other AND match one in given list of 24 known indexes, move record to files for matched indexes
    If indexes are not found in the list of 24, move record to files for unknown indexes
    If indexes are in list of 24 but not matched to each other, move record to files for hopped indexes

Start by opening all four input files (R1, R2, R3, R4)
    
    Read four lines (one record) from all files at the same time

        if the line is an empty string (end of file) then break

    Reverse complement R3 (index 2) to match index 1. Indexes need to have the same sequence.  If they have different sequences, index hopping has occurred.

    Next, edit the headers to include the sequence of the index-pair to the header of both reads in all fastq files

    Check if index 1 AND index 2 are in the set of 24 known indexes

        If one or both are NOT in the known set, create a counter for unknown indexes (unknwon += 1) and then add R1 record to unk_R1.fq and R4 record to unk_R2.fq

        If both ARE in the known set, check if the average q score of the record is above a given threshold (will be decided later)

            If the record does not have a q score above the threshold, return to the unknown path

            If the record does pass the q score threshold, check if the indexes are identical (should be easy to check bc reverse complement function as beginning will have them in the same direction)

                If the indexes are NOT identical, create a dictionary for hopped indexes and set a counter, then add R1 to hopped_R1.fq and R4 to hopped_R2.fq

                If the indexes ARE identical, create a dictionary for matched indexes and set a counter, then add R! to index_R1.fq and index_R2.fq

Unless the line is an empty string, return to beginning of loop for each new record and continue the same process for all records in all 4 files

Functions:

def rev_comp(index 2: str) -> str:
    ```Takes index 2 (R3) and reverses it to be in the same direction as index 1 (R2)```
    return reverse complement
Inuput: ATGCT
Expected Output: AGCAT

def edit_header(F1.header, F2.seq, F4.header, F3.seq: str) -> str:
    ```Adds the sequence of index 1 and index 2 to all header lines```
    return edited header
Input: @ABC123
Expected output: @ABC123-AAAAAAA-TTTTTTTTT

def check_q(F2.qscores, F3.qscores, qscore threshold: bool) -> bool:
    ```Checks if the qscores are above a sufficient threshold```
    return true or false
Input: 36, 32, 20
Expected output: True
