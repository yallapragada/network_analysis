'''
script to remove duplicates after concatenating all 10 proteins and create new fasta files with deduplicated sequences
usage: python duplicate_remover.py [folder containing aligned files with .afasta extension
'''


import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
import graph_analysis_util as util

strains_by_year = {}


def concat_sequences(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10):

    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')
    sequences3 = AlignIO.read(file3, 'fasta')
    sequences4 = AlignIO.read(file4, 'fasta')
    sequences5 = AlignIO.read(file5, 'fasta')
    sequences6 = AlignIO.read(file6, 'fasta')
    sequences7 = AlignIO.read(file7, 'fasta')
    sequences8 = AlignIO.read(file8, 'fasta')
    sequences9 = AlignIO.read(file9, 'fasta')
    sequences10 = AlignIO.read(file10, 'fasta')

    complete_sequences = []

    for sequence1 in sequences1:
        strain_name = util.get_strain_name(sequence1)
        sequence2 = util.get_matching_sequence(sequences2, strain_name=strain_name)
        sequence3 = util.get_matching_sequence(sequences3, strain_name=strain_name)
        sequence4 = util.get_matching_sequence(sequences4, strain_name=strain_name)
        sequence5 = util.get_matching_sequence(sequences5, strain_name=strain_name)
        sequence6 = util.get_matching_sequence(sequences6, strain_name=strain_name)
        sequence7 = util.get_matching_sequence(sequences7, strain_name=strain_name)
        sequence8 = util.get_matching_sequence(sequences8, strain_name=strain_name)
        sequence9 = util.get_matching_sequence(sequences9, strain_name=strain_name)
        sequence10 = util.get_matching_sequence(sequences10, strain_name=strain_name)

        if (sequence2 and sequence3 and sequence4 and sequence5 and sequence6 and sequence7 and sequence8 and sequence9 and sequence10):
            update_year_counts(sequence1)
            complete_sequence=[]
            complete_sequence.append(util.get_strain_name(sequence1))
            complete_sequence.append(sequence1.seq+sequence2.seq+sequence3.seq+sequence4.seq+sequence5.seq+sequence6.seq+sequence7.seq+sequence8.seq+sequence9.seq+sequence10.seq)
            complete_sequences.append(complete_sequence)

    print_year_counts()
    return complete_sequences


def update_year_counts(record):

    strain_name = util.get_strain_name(record)
    year        = strain_name.split('/')[-1]
    if year in strains_by_year:
        strains_by_year[year].append(strain_name)
    else:
        strains_by_year[year]   = [strain_name]

def print_year_counts():
    # what are the counts per year?
    for key, value in strains_by_year.items():
        print(key, len(value))

def sequence_cleaner(complete_sequences):
    # Create our hash table to add the sequences
    unique_sequences={}

    for complete_sequence in complete_sequences:
        # Take the current sequence
        sequence = complete_sequence[1].upper()
        strain = complete_sequence[0]

        if sequence not in unique_sequences and strain not in unique_sequences.values():
            unique_sequences[sequence] = strain

    return unique_sequences


def create_unique_sequences(fasta_file, unique_strains_file):
    sequences = AlignIO.read(fasta_file, 'fasta')
    strains = list(util.read_file(unique_strains_file))[0]
    new_fasta_sequences = []

    for sequence in sequences:
        strain_name = util.get_strain_name(sequence)
        if strain_name in strains:
            new_fasta_sequences.append(sequence)
            strains.remove(strain_name)

    return new_fasta_sequences


def run_dup_remover():

    folder = sys.argv[1]
    file1 = folder + os.sep + 'ha.afasta'
    file2 = folder + os.sep + 'na.afasta'
    file3 = folder + os.sep + 'm1.afasta'
    file4 = folder + os.sep + 'm2.afasta'
    file5 = folder + os.sep + 'np.afasta'
    file6 = folder + os.sep + 'pb1.afasta'
    file7 = folder + os.sep + 'pb2.afasta'
    file8 = folder + os.sep + 'pa.afasta'
    file9 = folder + os.sep + 'ns1.afasta'
    file10 = folder + os.sep + 'ns2.afasta'
    complete_sequences = concat_sequences(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10)
    print('number of complete sequences ', len(complete_sequences))
    print('length of complete sequence', len(complete_sequences[0][1]))

    unique_sequences_dict = sequence_cleaner(complete_sequences)
    strains_100 = [strain for strain in unique_sequences_dict.values()]
    util.write_strains_to_csv(strains_100, folder + os.sep + 'unique_strains_100')
    print('number of sequences_100 ', len(strains_100))


def run_unique_fasta_creator():

    folder = sys.argv[1]
    unique_strains_file = folder + os.sep + 'unique_strains_100.csv'
    proteins = ['ha', 'na', 'm1', 'm2', 'np', 'pb1', 'pb2', 'pa', 'ns1', 'ns2']

    for protein in proteins:
        fasta_file = folder + os.sep + protein + '.afasta'
        unique_sequences = create_unique_sequences(fasta_file, unique_strains_file)
        print('number of unique sequences for ', protein, ' ', len(unique_sequences))
        unique_file = folder + os.sep + 'u_' + protein + '.afasta'
        msa = MultipleSeqAlignment(unique_sequences)
        AlignIO.write(msa, handle=unique_file, format='fasta')


if __name__ == '__main__':
    run_dup_remover()
    run_unique_fasta_creator()
