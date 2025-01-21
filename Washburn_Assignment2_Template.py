# This script was created by Tucker J Lancaster
import sys
import textwrap


#finished
def valid_DNA_sequence(DNA):
    # This function takes in a string containing possible DNA sequence and returns True if all nucleotides are valid (A,a,C,c,G,g,T,t) and False if any nucleotide is invalid
    #DNA = DNA.strip()
    valid = "AGTCagtc"
    for letter in DNA: 
        if letter not in valid: 
            return False
    return True

def print_DNA_sequence(DNA, mode):
    # This function takes in a string containing valid DNA sequence and a mode that specifies the output format. Prints info to a screen, maximum characters per line
    for i in range(0,3):
        #This loop is used to create the three possible frames on the forward strand
        DNA_sequence = DNA[i:]# Create DNA sequence for current frame. Ensure it is divisible by 3
        if len(DNA_sequence)%3 == 1:
            DNA_sequence = DNA_sequence[:-1]
        elif len(DNA_sequence)%3 == 2:
            DNA_sequence = DNA_sequence[:-2]
            #print(DNA_sequence)
        else:
            DNA_sequence = DNA_sequence
            #print(DNA_sequence)
        
        translated_sequence = translate(DNA_sequence, mode)
    
        # #Print out direction (5' to 3') and frame to screen
        frame = list(range(3))
        
                
        if mode == "DNA":
            translated_sequence = textwrap.wrap(translated_sequence, width = 60, subsequent_indent = " ")
            DNA_sequence = textwrap.wrap(DNA_sequence, width = 60)
            new_dict = dict(zip(DNA_sequence, translated_sequence))
            print("5' to 3' Frame: {}".format(i))
            for key in new_dict:
                #print("5' to 3' Frame: {}".format(i))
                print(key)
                print(new_dict[key])        
        else:
            print("5' to 3' Frame: {}".format(i))
            print(translated_sequence)
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
          
    # #this is for the reverse strand     
    rev_DNA_sequence = reverse_complement(DNA)
    for i in range(0,3):
        #This loop is used to create the three possible frames on the forward strand
        DNA_sequence = rev_DNA_sequence[i:] # Create DNA sequence for current frame. Ensure it is divisible by 3
        if len(DNA_sequence)%3 == 1:
            DNA_sequence = DNA_sequence[:-1]
           
        elif len(DNA_sequence)%3 == 2:
            DNA_sequence = DNA_sequence[:-2]
           
        else:
            DNA_sequence = DNA_sequence

        translated_sequence = translate(DNA_sequence, mode)
        # Print out direction (5' to 3') and frame to screen
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
        if mode == "DNA":
            translated_sequence = textwrap.wrap(translated_sequence, width = 60, subsequent_indent= " ")
            DNA_sequence = textwrap.wrap(DNA_sequence, width = 60)
            new_dict = dict(zip(DNA_sequence, translated_sequence))
            print("3' to 5' Frame: {}".format(i))
            for key in new_dict:
                #new_dict[key] = " " + new_dict[key]
                print(key)
                print(new_dict[key])
        else:
            print("3' to 5' Frame: {}".format(i))
            print(translated_sequence)
           
    
        
        

#this works now  
if len(sys.argv) == 2:
    mode = sys.argv[1]
    mode = mode.upper()

    

def translate(DNA_sequence, mode):
    # Create dictionaries that translate codons into amino acids of appropriate format
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
    'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    #If VERBOSE, modify the start and stop codons and modify all codons to add space after
    if mode == 'VERBOSE':
        codontable.update({'ATG': 'Met', 'TAA': 'Stop', 'TGA':'Stop', 'TAG':'Stop'})
        for values in codontable:
            codontable[values] = codontable[values] + " "
        
        
    #If DNA, modify all codons to add space before and after
    if mode == 'DNA':
        for values in codontable:
            codontable[values] = " " + codontable[values] + " " 
    
    # Loop through DNA sequence codons
    out_seq = ''
    if len(DNA_sequence)%3 == 0:
        for i in range(0, len(DNA_sequence), 3):
            codon = DNA_sequence[i:i + 3]
            codon = codon.upper()
            out_seq += codontable[codon]
        return out_seq
    
    
     
#finished
def reverse_complement(DNA_sequence):
    # Returns string containing reverse complement of DNA sequence
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    #check about .upper() to DNA_sequence
    seq_list = list(DNA_sequence.upper())
    for i in range(len(seq_list)):
        nuc = seq_list[i]
        pair = complement[nuc]
        seq_list[i] = pair

    seq_list.reverse()
    
    return "".join(seq_list)



# Part I: Determine if the user has entered the appropriate number of argments when they called the script (one). Determine if the user entered one of three valid options for the mode. If there is an error in either of these, print out informative error messages indicating which error was made, what the three valid options are, and then quit the program.

# Part II: Create loop to query user for DNA sequence

# Part III: Determine if user wants to exit program

# Part IV: Determine if user input is valid DNA sequence. If DNA sequence is not valid, print error message and allow user to enter new DNA sequence. 

# Part V: Print out 6 translated frames to the screen in appropriate format

mode_list = ["VERBOSE", "DNA", "COMPACT"]
while True:

    if len(sys.argv) != 2:
        print("Invalid number of options.\nUsage: python3 assignment2_template.py <mode> \nMode can be one of the following options: \n\t VERBOSE \n\t COMPACT \n\t DNA")
        sys.exit()
    elif mode not in mode_list:
        print(f"{mode} Is not a valid option\nUsage: python3 assignment2_template.py <mode> \nMode can be one of the following options: \n\t VERBOSE \n\t COMPACT \n\t DNA")
        sys.exit()
    else:
        break
    
# # Part II: Create loop to query user for DNA sequence
while True:
    DNA = input("Enter DNA sequence (or exit to quit the program): " )
    if DNA == 'Exit' or DNA == 'exit':
        sys.exit()
    elif DNA == "":
        continue
    else:
        #valid_DNA_sequence(DNA)
        if valid_DNA_sequence(DNA) == True:
            #maybe change this to the print(print_DNA_sequence(DNA, mode)) function 
            print_DNA_sequence(DNA, mode)
        else:
            print("Invalid DNA sequence. Characters must be one of A, a, T, t, G, g, C, c")
    


    














