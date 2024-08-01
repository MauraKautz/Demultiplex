#!/usr/bin/env python

# Author: Maura Kautz mkautz@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.2"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("ATCGNatcgn")
RNA_bases = set("AUCGNaucgn")

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Calculate average quality score of Phred string'''
    sum=0   #initialize sum to be zero    
    for x in range(0,len(phred_score)):         #iterate through phred scores, incrementing each time
        sum+=convert_phred(phred_score[x])
    return sum / len(phred_score)

def validate_base_seq(seq: str,RNAflag=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq), "String contains invalid characters - are you sure you used a nucleotide sequence?"
    seq = seq.upper()
    return (seq.count("G")+seq.count("C"))/len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    list_length=len(lst)
    if list_length %2 == 0:
        i1 = list_length//2 - 1
        i2 = list_length//2
        med = (lst[i1] + lst[i2]) / 2
    else:
        i1 = list_length//2
        med = lst[i1]
    return(med)


def oneline_fasta(inputfile: str) -> str:
    '''Moves FASTA amino acids into one long line'''
    with open(f'new_{inputfile}', "w") as fh:
        with open(inputfile, "r") as fh1:
            linenum=0
            while True:
                line=fh1.readline().strip()
                if line=="":
                    break

                if line.startswith(">") and linenum == 0:
                    fh.write(f'{line}\n')
                    linenum += 1
                elif line.startswith(">") and linenum != 0:
                    fh.write(f'\n{line}\n')
                    linenum += 1
                else:
                    fh.write(line)
                    linenum += 1
    return f'new_{inputfile}'

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
