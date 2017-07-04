'''
when you download strain data with all proteins, you get a single fasta file
that contains sequences of all 10 proteins for a given strain
this script splits a big fasta file that contains all 10 proteins of influenza
to 10 seperate fasta files per protein
'''

import graph_analysis_util as gau
import sys
from Bio import SeqIO
import random

def regular_split():
    strains_file = sys.argv[1]
    gau.split_fasta(strains_file)

#get list of strain names; randomly pick k strains per year if
#no. of strains > k for that year
#includes only strains < max_year
def get_strain_names(fasta_file, k, max_year):

    strain_names    = []
    strains_by_year = {}

    for record in SeqIO.parse(fasta_file, "fasta"):

        description_list = record.description.split('|')

        if (len(description_list)<4):
            print("protein description too short ", record.description)
            continue

        if 'HA' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            year        = strain_name.split('/')[-1]
            if year in strains_by_year:
                strains_by_year[year].append(strain_name)
            else:
                strains_by_year[year]   = [strain_name]

    # what are the counts per year?
    for key, value in strains_by_year.items():
        print(key, len(value))

    for key, value in strains_by_year.items():
        try:
            year = int(key)
        except Exception:
            year = 0
        if year < max_year:
            print(key + ' is less than ' + str(max_year))
            if len(value) > k:
                small_list = random.sample(value, k=k)
                for strain_name in small_list:
                    strain_names.append(strain_name)
            else:
                for strain_name in value:
                    strain_names.append(strain_name)
        else:
            print(key + ' is greater than ' + str(max_year))

    return strain_names

#splits big fasta file to protein fasta files
def custom_split():

    fasta_file = sys.argv[1]
    strain_names = get_strain_names(fasta_file, k=100, max_year=2020)

    ha_sequences = []
    na_sequences = []
    np_sequences = []
    m1_sequences = []
    m2_sequences = []
    pb1_sequences = []
    pb2_sequences = []
    ns1_sequences = []
    ns2_sequences = []
    pa_sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        description_list = record.description.split('|')

        if (len(description_list) < 4):
            print("protein description too short ", record.description)
            continue

        if 'HA' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                ha_sequences.append(record)

        if 'NA' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                na_sequences.append(record)

        if 'NP' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                np_sequences.append(record)

        if 'M1' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                m1_sequences.append(record)

        if 'M2' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                m2_sequences.append(record)

        if 'PB1 Polymerase (basic) protein 1' in description_list[5] or 'PB1 Polymerase (basic) protein 1' in \
                description_list[4]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                pb1_sequences.append(record)

        if 'PB2' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                pb2_sequences.append(record)

        if 'NS1 Non-structural protein 1' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                ns1_sequences.append(record)

        if 'NS2 Non-structural protein 2' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                ns2_sequences.append(record)

        if 'PA Polymerase (acidic) protein' in description_list[5]:
            strain_name = gau.get_strain_name(record)
            if strain_name in strain_names:
                pa_sequences.append(record)

    output_handle = open("ha.fasta", "w")
    print("length of ha sequences ", len(ha_sequences))
    SeqIO.write(ha_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("na.fasta", "w")
    print("length of na sequences ", len(na_sequences))
    SeqIO.write(na_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("m1.fasta", "w")
    print("length of m1 sequences ", len(m1_sequences))
    SeqIO.write(m1_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("m2.fasta", "w")
    print("length of m2 sequences ", len(m2_sequences))
    SeqIO.write(m2_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("np.fasta", "w")
    print("length of np sequences ", len(np_sequences))
    SeqIO.write(np_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("pb1.fasta", "w")
    print("length of pb1 sequences ", len(pb1_sequences))
    SeqIO.write(pb1_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("pb2.fasta", "w")
    print("length of pb2 sequences ", len(pb2_sequences))
    SeqIO.write(pb2_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("ns1.fasta", "w")
    print("length of ns1 sequences ", len(ns1_sequences))
    SeqIO.write(ns1_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("ns2.fasta", "w")
    print("length of ns2 sequences ", len(ns2_sequences))
    SeqIO.write(ns2_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("pa.fasta", "w")
    print("length of pa sequences ", len(pa_sequences))
    SeqIO.write(pa_sequences, output_handle, "fasta")
    output_handle.close()


if __name__ == '__main__':
    custom_split()