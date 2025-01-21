"""
Assignment 6: Protein Structure Analysis

In this template, the ret_transmembrane() function, as well as the __name__=='__main__' block, have been written
for you and should not be modified. The remaining 11 functions are missing portions that you will need to provide.
While this may seem like a lot, most of the functions require only a few lines of code to complete (i.e., the task
has been broken into bitesize pieces for your convenience). For functions with return values, these values are initially
to "None" in the template as a placeholder. For functions without return values, you will instead see a "pass"
statement. In both cases, please replace these lines with meaningful code when you write your functions.
Assignment Author: Tucker J Lancaster
"""

import pdb
import Bio.PDB
from Bio.Align import substitution_matrices, MultipleSeqAlignment
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import argparse


def ret_transmembrane(loc):
    """
    This function is already implemented. This function takes in the coordinate of a residue from
    the b2 inactive receptor and returns the transmembrane receptor the residue belongs to or
    False if the residue is part of the intercellular or extracellular loop
    NOTE: You can only use this function on residue coordinates from the b2 inactive receptor

    :param loc: residue coordinate
    :type loc: int
    :return: transmembrane receptor name, or False if the residue is part of a loop
    :rtype: str or bool
    """

    if loc >= 31 and loc <= 56:
        return 'TM1'
    if loc >= 70 and loc <= 92:
        return 'TM2'
    if loc >= 104 and loc <= 129:
        return 'TM3'
    if loc >= 150 and loc <= 171:
        return 'TM4'
    if loc >= 198 and loc <= 220:
        return 'TM5'
    if loc >= 275 and loc <= 298:
        return 'TM6'
    if loc >= 306 and loc <= 327:
        return 'TM7'
    return False


def ret_ligand(pdb, ligand_name):
    """
    This function is already implemented. It takes in a pdb structure object and ligand name, and returns the
    associated ligand residue

    :param pdb: pdb object
    :type pdb: Bio.PDB.Structure.Structure
    :param ligand_name: name of the ligand of interest
    :type ligand_name: str
    :return: ligand residue object if ligand name is valid, else 'None'
    :rtype: Bio.PDB.Residue.Residue or Nonetype
    """

    
    for residue in pdb.get_residues():
        if residue.resname == ligand_name:
            return residue
    return None


def open_pdb_files(ref_pdb_file, mov_pdb_file):
    """
    Reads in the ref (activated) and mov (inactivated) pdb files, and returns the structures they encode as
    Bio.PDB.Structure.Structure objects
    Hint: use the Bio.PDB.PDBParser class, and the associated get_structure method

    :param ref_pdb_file: filename for the inactivated pdb
    :type ref_pdb_file: str
    :param mov_pdb_file: filename for the activated pdb
    :type mov_pdb_file: str
    :return ref_pdb: inactivated pdb object
    :rtype ref_pdb: Bio.PDB.Structure.Structure
    :return mov_pdb: activated pdb object
    :rtype mov_pdb: Bio.PDB.Structure.Structure
    """
    
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_pdb = pdb_parser.get_structure('inactive', ref_pdb_file)
    mov_pdb = pdb_parser.get_structure('active', mov_pdb_file)
    return ref_pdb, mov_pdb


def create_polypeptides(ref_pdb, mov_pdb):
    """
    take in ref_pdb and mov_pdb objects (the objects returned by open_pdb_files), convert each to a polypeptide object,
    and return the two polypeptide objects
    Hint: Use the Bio.PDB.PPBuilder object, and the associated build_peptides method

    :param ref_pdb: inactivated pdb object
    :type ref_pdb: Bio.PDB.Structure.Structure
    :param mov_pdb: activated pdb object
    :type mov_pdb: Bio.PDB.Structure.Structure
    :return ref_peptides: inactivated peptide object
    :rtype ref_peptides: list
    :return mov_peptides: activated peptide object
    :rtype mov_peptides: list
    """
    
    pep=Bio.PDB.PPBuilder()
    ref_peptides = pep.build_peptides(ref_pdb)
    mov_peptides = pep.build_peptides(mov_pdb)
    return ref_peptides, mov_peptides


def create_sequences(ref_peptides, mov_peptides):
    """
    take in ref and mov peptide objects of the sort returned by create_polypeptides, and return them as ref and mov
    sequence objects. Note that one of the polypeptide objects has a gap. You will need to append the two sequences
    together
    Hint: use the get_sequence method of ref_peptides and mov_peptides

    :param ref_peptides: inactivated peptide object
    :type ref_peptides: list
    :param mov_peptides: activated peptide object
    :type mov_peptides: list
    :return ref_seq: activated sequence
    :rtype ref_seq: Bio.Seq.Seq
    :return mov_seq: inactivated sequence
    :rtype mov_seq: Bio.Seq.Seq
    """


    ref_seq = ref_peptides[0].get_sequence()
    
    #https://www.biostars.org/p/92665/
    mov_list= [i.get_sequence() for i in mov_peptides]
    mov1 = mov_list[0]
    mov2 = mov_list[1]
    mov_seq= mov1 + mov2
    return ref_seq, mov_seq

   


def create_pairwise_alignment(ref_seq, mov_seq):
    """
    take in ref and mov sequence objects (of the type returned by create_sequence) and return a pairwise alignment of
    the two sequences. Use the BLOSUM62 matrix with a gap open penalty of -10 and a gap extend penalty of -0.5
    Hint: use substitution_matrices.load to get the substitution matrix. Use pairwise2.align.globalds to get the
    pairwise alignment.

    :param ref_seq: inactivated sequence
    :type ref_seq: Bio.Seq.Seq
    :param mov_seq: activated sequence
    :type mov_seq: Bio.Seq.Seq
    :return t_alignment: pairwise alignment object
    :rtype: Bio.pairwise2.Alignment
    """
    
    matrix = substitution_matrices.load(name = 'BLOSUM62')
    gap_open = -10
    gap_extend = -0.5
    align = pairwise2.align.globalds(ref_seq, mov_seq, matrix, gap_open, gap_extend)
    t_alignment = align[0]
    return t_alignment


def create_structure_alignment(t_alignment):
    """
    take in a paired alignment (as returned by create_pairwise_alignment) and return a structure alignment object.
    Hint: make a MultipleSeqAlignment object out of the pairwise alignment, then the Bio.PDB.StructureAlignment method
    to generate a structure_alignment.

    :param t_alignment: pairwise alignment object
    :type t_alignment: Bio.pairwise2.Alignment
    :return structure_alignment: Structure alignment object
    :rtype structure_alignment: Bio.PDB.StructureAlignment.StructureAlignment
    """
   
    alignment_paired = MultipleSeqAlignment([Bio.SeqRecord.SeqRecord(t_alignment[0]), Bio.SeqRecord.SeqRecord(t_alignment[1])])
    structure_alignment = Bio.PDB.StructureAlignment(alignment_paired, ref_pdb, mov_pdb)
    return structure_alignment


def create_atom_list(structure_alignment):
    """
    Create a list of atoms to use for the structural alignment. Only use residues that are in transmembrane regions,
    using the ret_transmembrane function I have implemented for you .Also, use the four atoms in each amino acid that
    are present in all 20 amino acids (hint two C, 1 N, and 1 O atom). The structure_alignment.duos attribute may be
    useful

    :param structure_alignment: structure alignment object
    :type structure_alignment: Bio.PDB.StructureAlignment.StructureAlignment
    :return ref_atoms: list of atoms in the ref to use for alignment
    :rtype ref_atoms: list
    :return mov_atoms: list of atoms in the mov to use for alignment
    :rtype mov_atoms: list
    """
    

    ref_atom_list = []
    mov_atom_list = []
    for duo in structure_alignment.duos:
        res1 = duo[0]
        res2 = duo[1]
        if res1 and res2:
            location = res1.id[1]
            if ret_transmembrane(location):
                ref_atom_list.append(res1.child_list[0:4])
                mov_atom_list.append(res2.child_list[0:4])
    
    #https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    ref_atoms = [atom for sub_atom in ref_atom_list for atom in sub_atom]
    mov_atoms = [atom for sub_atom in mov_atom_list for atom in sub_atom]
    return ref_atoms, mov_atoms

   


def create_superimposer(ref_atoms, mov_atoms, mov_pdb):
    """
    take in the ref_atoms list, mov_atoms list, and mov_pdb pdb structure object, and return the associated aligned
    Superimposer object

    :param ref_atoms: list of atoms in the ref to use for alignment
    :type ref_atoms: list
    :param mov_atoms: list of atoms in the mov to use for alignment
    :type mov_atoms: list
    :param mov_pdb: inactivated pdb object
    :type mov_pdb: Bio.PDB.Structure.Structure
    :return super_imposer: aligned superimposer object
    :rtype super_imposer: Bio.PDB.Superimposer.Superimposer
    """
   
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, mov_atoms)
    super_imposer.apply(mov_pdb.get_atoms())
    return super_imposer


def save_aligned_pdb(mov_pdb, out_pdb_file):
    """
    save the aligned mov_pdb object to out_pdb_file
    Hint: create a Bio.PDB.PDBIO() object, then the set_structure method of that object with mov_pdb as the argument,
    and finally the save method of the Bio.PDB.PDBIO() object to save it.

    :param mov_pdb: aligned pdb structure object 
    :type mov_pdb: io.PDB.Structure.Structure
    :param out_pdb_file: destination file for results
    :type out_pdb_file: str
    """

   
    io = Bio.PDB.PDBIO()
    io.set_structure(mov_pdb)
    io.save(out_pdb_file)
    


def avg_residue_distance(res1, res2):
    """
    This function takes in two residues and returns the average distance between the four common atoms found in all
    amino acids. This method is used to compare how far a residue moves when activated by binding to an agonist ligand

    :param res1: first residue
    :type res1: Bio.PDB.Residue.Residue
    :param res2: second residue
    :type res2: Bio.PDB.Residue.Residue
    :return: average distance
    :rtype: float
    """

    res1_list = res1.child_list[0:4] 
    res2_list = res2.child_list[0:4]  
    distance_list = []
    for x in range(4):
        distance_list.append(res1_list[x] - res2_list[x])


    avg_dist = ((sum(distance_list))/4)
    return avg_dist


def min_residue_distance(res1, res2):
    """
    This method takes in two residues and returns the minimum distance between the two residues This method is used to
    identify residues that are close to the ligand

    :param res1: first residue
    :type res1: Bio.PDB.Residue.Residue
    :param res2: second residue
    :type res2: Bio.PDB.Residue.Residue
    :return: minimum distance
    :rtype: float
    """
    
    atom_dis = []
    for i in res1:
        for j in res2:
            atom_dis.append(j - i)
    min_dist = min(atom_dis)
    return min_dist


if __name__ == '__main__':
    """
    This section of the code is implemented for you, but reading through (and understanding) it will help you with the
    rest of the assignment. Do not make any modifications in this section
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--Inactivated_pdb', type=str, help='Name of a pdb file containing reference pdb',
                        default='B2R_Inactive.pdb')
    parser.add_argument('--Activated_pdb', type=str, help='Name of a pdb file containing activated pdb',
                        default='B2R_Active.pdb')
    parser.add_argument('--Moved_pdb', type=str,
                        help='Name of a pdb file to store the structurally aligned activated pdb',
                        default='B2R_Active_moved.pdb')
    args = parser.parse_args()

    # save the arguments to different names for convenience
    ref_pdb_file = args.Inactivated_pdb
    mov_pdb_file = args.Activated_pdb
    out_pdb_file = args.Moved_pdb

    # use the user-defined function open_pdb_files to read in the ref and mov pdbs
    ref_pdb, mov_pdb = open_pdb_files(ref_pdb_file, mov_pdb_file)

    # use the user-defined function create_polypeptides to convert the pdb objects to polypeptide objects
    ref_peptides, mov_peptides = create_polypeptides(ref_pdb, mov_pdb)

    # use the user-defined function create_sequences to create sequences from the two polypeptide objects
    ref_seq, mov_seq = create_sequences(ref_peptides, mov_peptides)

    # use the user-defined function create_pairwise_alignment to create a pairwise alignment from the sequence objects
    t_alignment = create_pairwise_alignment(ref_seq, mov_seq)

    # use the user-defined function create_structure_alignment to create a StructureAlignment object
    structure_alignment = create_structure_alignment(t_alignment)

    # use the user-defined function create_atom_list to create a list of atoms to use for the structural alignment.
    ref_atoms, mov_atoms = create_atom_list(structure_alignment)

    # use the user-defined function create_superimposer to create an aligned Superimposer object
    super_imposer = create_superimposer(ref_atoms, mov_atoms, mov_pdb)

    # use the user-defined function save_aligned_pdb to save the aligned pdb object to out_pdb_file
    save_aligned_pdb(mov_pdb, out_pdb_file)

    # Determine how far each transmembrane residue moves when activated. Y have already identified the duos of
    # residues that are homologous. Use these along with the avg_residue_distance function to determine how far each
    # residue moves. Plot this information using  the residue position on the x axis and the distance it moved on the
    # y axis

    distances = {}
    # The keys in the distances dict are the transmembrane domains (TM1 - TM7). The value are tuples with the first
    # element the residue position and second element the distance moved
    for duo in structure_alignment.duos:
        resR = duo[0]
        resM = duo[1]
        if resR and resM is not None:
            loc = resR.get_id()[1]
            if ret_transmembrane(loc):
                try:
                    distances[ret_transmembrane(loc)].append((loc, avg_residue_distance(resR, resM)))
                except KeyError:
                    distances[ret_transmembrane(loc)] = [(loc, avg_residue_distance(resR, resM))]

    # Plot distances between the aligned pdb files (Don't change anything)
    for key in distances:
        a, = plt.plot([x[0] for x in distances[key]], [x[1] for x in distances[key]], label='Active/Inactive',
                      color='black')

    plt.xlabel('Position (amino acids)')
    plt.ylabel('Distance from reference (Angstroms)')
    plt.text(35, 6, 'tm1')
    plt.text(75, 6, 'tm2')
    plt.text(110, 6, 'tm3')
    plt.text(155, 6, 'tm4')
    plt.text(205, 6, 'tm5')
    plt.text(280, 6, 'tm6')
    plt.text(310, 6, 'tm7')

    plt.show()

    # Create a command file for ChimeraX to color each residue that is within 4 Angstroms red. You will also need to
    # implement the min_residue_distance for this to work

    ligand_residue = ret_ligand(ref_pdb, 'LIG')
    with open('Ligand_residues.cxc', 'w') as f:
        for res in ref_pdb.get_residues():
            if min_residue_distance(res, ligand_residue) < 4:
                print('color :' + str(res.get_id()[1]) + ' red', file=f)
