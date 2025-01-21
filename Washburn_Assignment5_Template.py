"""
Assignment 5: Automatic Primer Identification

#python version: Python 3.9.18
#used git clone to install primer3 and placed primer3_core into a directory in my $PATH

In this assignment, you will practice using python's subprocess library, the pyfasta library, and the primer3
command-line-interface to automatically identify primers. This template provides a potential route to a completed
script, but the implementation details are entirely up to you. The checker script, as you will see, does little more
than confirm that your script prints the correct primers to the console for a handful of argument combinations.
Author: Tucker J Lancaster
"""
#python version: Python 3.9.18
#used git clone to install primer3 and placed primer3_core into a directory in my $PATH


import argparse
import subprocess
from pyfasta import Fasta 

#make sure to use primer3 environment
"""
If you have not already, install pyfasta and primer3 into the current environment using the following commands:
conda install -c bioconda primer3
conda install -c bioconda pyfastaz
"""

"""
Create your parser and arguments. Your script should expect three positional arguments: First, the name of a
fasta file, second, the chromosome you want to amplify, and third, the position you want amplified. Enforce that
the first two arguments are of type str, and the third is of type int.
"""


parser = argparse.ArgumentParser(description = 'This program runs primer3 to identify the best two primers to amplify a given sequence.', usage = "Assignment5_Template.py [-h] fasta_file chromosome position") #fix this
parser.add_argument('fasta_file', help='Enter a genome file containing DNA sequence', type=str)
parser.add_argument('chromosome', help = 'Enter the chromosome you want amplified', type = str) 
parser.add_argument('position', help = 'Enter the position you want amplified', type = int)
args = parser.parse_args()

"""
Open the fasta file and read in the entirety of the sequence. You can use whatever data structure you prefer.
"""

#https://github.com/brentp/pyfasta
fasta_file = Fasta(args.fasta_file)

"""
Identify the sequence 500 bp upstream and downstream of the requested position and store it as a string.
"""

chrom = args.chromosome
position = args.position

seq = fasta_file[chrom][position-500:position+501]


"""
Add braces to the DNA sequence 100 bp upstream and downstream of the requested position. Alternatively, you can use 
an argument in the input file to specify this
"""

#seq = ''.join((seq[:400],'[',seq[400:]))
#seq = ''.join((seq[:601],']',seq[601:]))


"""
Create an input file containing the DNA sequence and specify a product length of 600 - 800 base pairs
"""

#https://stackoverflow.com/questions/48959098/how-to-create-a-new-text-file-using-python
#check target position sequence 
with open('assignment5_primer3.txt', 'w') as file: 
    file.write('SEQUENCE_ID=assignment5_primer3\n' + 'SEQUENCE_TEMPLATE=' + seq + '\n' + 'PRIMER_PICK_LEFT_PRIMER=1\n' + 'PRIMER_PICK_RIGHT_PRIMER=1\n' + 'PRIMER_OPT_SIZE=20\n' + 'PRIMER_MIN_SIZE=18\n' + 'PRIMER_MAX_SIZE=27\n' + 'PRIMER_PRODUCT_SIZE_RANGE=600-800\n' +  'SEQUENCE_TARGET=400,201\n' + '=')

"""   
Execute the primer3_core command on the input file you created.
"""

info = subprocess.run(['primer3_core', 'assignment5_primer3.txt'], capture_output = True)

""" 
Read in and parse the output to identify the sequence of best two primers. Print these out to the user
"""
primer_out = info.stdout.decode('utf-8')
primer_out = primer_out.split('\n')
left = primer_out[16]
right = primer_out[17]

#https://stackoverflow.com/questions/12572362/how-to-get-a-string-after-a-specific-substring
print(left.split("=",1)[1])
print(right.split("=",1)[1])