from Bio import SeqIO, AlignIO
import sys
from numpy import transpose, array
from pandas import DataFrame
from minepy import MINE
from scipy.stats import mode
import networkx as nx
import csv
import datetime
import matplotlib.pyplot as plt

import operator
from bokeh.plotting import figure, show, output_file
import math


# split a fasta file that contains several different proteins into multiple fasta files
# one for protein
def split_fasta(fasta_file):
    ha_sequences = []
    na_sequences = []
    np_sequences = []
    m1_sequences = []
    m2_sequences = []
    pb1_sequences = []
    pb2_sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        description_list = record.description.split('|')

        if (len(description_list)<4):
            print(record.description)
            continue

        if 'HA' in description_list[5]:
            ha_sequences.append(record)

        if 'NA' in description_list[5]:
            na_sequences.append(record)

        if 'NP' in description_list[5]:
            np_sequences.append(record)

        if 'M1' in description_list[5]:
            m1_sequences.append(record)

        if 'M2' in description_list[5]:
            m2_sequences.append(record)

        if 'PB1 Polymerase (basic) protein 1' in description_list[5]:
            pb1_sequences.append(record)

        if 'PB2' in description_list[5]:
            pb2_sequences.append(record)

    output_handle = open("ha.fasta", "w")
    SeqIO.write(ha_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("na.fasta", "w")
    SeqIO.write(na_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("m1.fasta", "w")
    SeqIO.write(m1_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("m2.fasta", "w")
    SeqIO.write(m2_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("np.fasta", "w")
    SeqIO.write(np_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("pb1.fasta", "w")
    SeqIO.write(pb1_sequences, output_handle, "fasta")
    output_handle.close()

    output_handle = open("pb2.fasta", "w")
    SeqIO.write(pb2_sequences, output_handle, "fasta")
    output_handle.close()


def generate_binary_sequence(fasta_list):
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in fasta_list]))]
    return array([[1 if x==mod[i] else 0 for i,x in enumerate(y)] for y in fasta_list])


def get_consensus_sequence(sequences):
    raw_sequences = [x.seq for x in sequences]
    mod = [mode(y)[0][0] for y in transpose(array([list(z) for z in raw_sequences]))]
    consensus_sequence = ''.join(mod)
    return consensus_sequence


def convert_01(sequence, consensus_sequence):
    return [1 if x==consensus_sequence[i] else 0 for i,x in enumerate(sequence)]


def convert_all_to_123(sequences):
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in sequences]))]
    consensus_sequence = ''.join(mod)
    sequences_123 = [convert_to_123(sequence, consensus_sequence) for sequence in sequences]
    # sequences_123 = [list(sequence_123) for sequence_123 in sequences_123_array]
    return sequences_123


def convert_to_123(sequence, consensus_sequence):
    new_sequence = [change_to_123(sequence[i], consensus_sequence[i]) for i in range(len(sequence))]
    return new_sequence


def change_to_123(residue, consensus_residue):
    if (residue == consensus_residue):
        return 0

    # initialize Transformed residue values to unknown to begin with
    residueT = 'U'
    consensus_residueT = 'U'

    types = {'X': ['R', 'K', 'D', 'E', 'P', 'N'],  ## surface
             'Y': ['Q', 'H', 'S', 'T', 'G'],  ## neutral
             'Z': ['A', 'I', 'L', 'F', 'V', 'Y', 'C', 'M', 'W']}  ## buried

    for t in types.keys():
        if residue in types[t]:
            residueT = t
        if consensus_residue in types[t]:
            consensus_residueT = t

    if (residueT == 'U'):
        return 0

    if (residueT == consensus_residueT):
        return 1
    else:
        if (residueT == 'Y' or consensus_residueT == 'Y'):
            return 2
        else:
            return 3


def get_central_node(G):
    centrality_dict = nx.degree_centrality(G)
    degree_centrality_nodes   = sorted(centrality_dict.items(), key=operator.itemgetter(1), reverse=True)
    print('nodes with highest degree centrality')
    if (len(degree_centrality_nodes)>200):
        central_nodesX = [node[0] for node in degree_centrality_nodes[:200]]
    else:
        central_nodesX = [node[0] for node in degree_centrality_nodes]
    for node in central_nodesX:
        print(node)
    return central_nodesX

def get_eigen_node(G):
    eigen_dict = nx.eigenvector_centrality(G)
    eigen_nodes = sorted(eigen_dict.items(), key=operator.itemgetter(1), reverse=True)
    print('nodes with highest eigen centrality')
    if (len(eigen_nodes)>50):
        eigen_nodesX = [node[0] for node in eigen_nodes[:50]]
    else:
        eigen_nodesX = [node[0] for node in eigen_nodes]
    for node in eigen_nodesX:
        print(node)
    return eigen_nodesX


def add_attribute(graph):
    protein_name_dict = {}
    for node in nx.nodes(graph):
        protein_name_dict[node] = node.split('_')[0]
    nx.set_node_attributes(graph, 'protein', protein_name_dict)

def create_graph(dataframe, filename):
    graph = nx.from_pandas_dataframe(dataframe, 'x', 'y', 'weight')
    add_attribute(graph)
    nx.write_graphml(graph, filename+'.graphml')

def get_sequences_from_alignment(alignment):
    sequences = [x.seq for x in alignment]
    return sequences


def perform_transpose(binary_sequences):
    transposed_list = transpose(array(binary_sequences)).tolist()
    return transposed_list


def create_DF_from_dict(dictionary):
    return DataFrame(dictionary)


def find(list_of_tuples, value):
    try:
        return next(x for x in list_of_tuples if value in x)
    except:
        print("***")
        return None


# Given a SeqRecord, return strain_name
def get_strain_name(record):
    strain_full_str = record.description.split('|')[3]
    if (strain_full_str.startswith('Strain')):
        strain = strain_full_str.split(':')[1]
    else:
        strain = record.description.split('|')[4].split(':')[1]
    return strain


def get_matching_sequence(sequences, strain_name):
    for sequence in sequences:
        if (get_strain_name(sequence) == strain_name):
            return sequence
    return None


def create_01_sequences(file1, file2):
    print(file1, file2)
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    p1_sequences_01 = []  # list of p1 sequences
    p2_sequences_01 = []  # list of p2 sequences

    consensus_sequence1 = get_consensus_sequence(sequences1)
    consensus_sequence2 = get_consensus_sequence(sequences2)

    for sequence1 in sequences1:
        strain_name = get_strain_name(sequence1)
        sequence2 = get_matching_sequence(sequences2, strain_name=strain_name)
        if (sequence2):
            p1_sequences_01.append(convert_01(sequence1.seq, consensus_sequence1))
            p2_sequences_01.append(convert_01(sequence2.seq, consensus_sequence2))

    return p1_sequences_01, p2_sequences_01


# create sequences with 0,1,2,3 values for protein1 and protein2 sequences
def create_0123_sequences(file1, file2):
    print(file1, file2)
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    p1_sequences_0123 = []  # list of p1 sequences
    p2_sequences_0123 = []  # list of p2 sequences

    consensus_sequence1 = get_consensus_sequence(sequences1)
    consensus_sequence2 = get_consensus_sequence(sequences2)

    for sequence1 in sequences1:
        strain_name = get_strain_name(sequence1)
        sequence2 = get_matching_sequence(sequences2, strain_name=strain_name)
        if (sequence2):
            p1_sequences_0123.append(convert_to_123(sequence1.seq, consensus_sequence1))
            p2_sequences_0123.append(convert_to_123(sequence2.seq, consensus_sequence2))

    return p1_sequences_0123, p2_sequences_0123


def create_df_from_dict(dictionary):
    return DataFrame(dictionary)


def perform_mic_2p(p1_sequences, p2_sequences, p1, p2, cutoff=0.5):
    mic_scores = []
    p1_sequences_t = transpose(array([list(z) for z in p1_sequences])).tolist()
    p2_sequences_t = transpose(array([list(z) for z in p2_sequences])).tolist()

    for idx1, record1 in enumerate(p1_sequences_t):
        for idx2, record2 in enumerate(p2_sequences_t):
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(record1, record2)
            if (mine.mic() > float(cutoff)):
                mic_score = {}
                mic_score['x'] = p1+'_'+str(idx1+1)
                mic_score['y'] = p2+'_'+str(idx2+1)
                mic_score['weight'] = mine.mic()
                mic_scores.append(mic_score)

    print('done computing mics for ', p1, p2)
    return mic_scores

def generate_output_filename(output_dir, basename):
    suffix = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = "_".join([basename, suffix])
    output_filename=output_dir+"\\"+filename+".txt"
    return output_filename

def entropy(text):
    frequencies={}
    shannon=0.00
    for each in text:
        try:
            frequencies[each]+=1
        except:
            frequencies[each]=1
    textlen=len(text)
    for k,v in frequencies.items():
        if (k!='-'):
            freq  =  1.0*v/textlen
            shannon+=freq*math.log2(freq)
    shannon*=-1
    return shannon


def entropy_all_positions(sequences):
    sequencesT = transpose(sequences)
    scores = [entropy(x) for x in sequencesT]
    for score in scores:
        if (score>1):
            print("SCORE ", score)
    return scores


def plotEntropies(entropies, residues, title):
    plot_simple_bokeh_line(y=entropies, x=list(residues), filename=title+".html", width=400, height=300, xaxis="residues", yaxis="entropies", title=title)


def plot_simple_bokeh_line(x,y, filename, height, width, xaxis, yaxis, title):
    output_file(filename)
    p=figure(plot_width=width, plot_height=height,title=title)
    p.line(x,y, color="green")
    p.xaxis.axis_label=xaxis
    p.yaxis.axis_label=yaxis
    show(p)

def read_file(filename):
    f = open(filename)
    csv_f = csv.reader(f)
    return csv_f

def run_entropies():
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    p1 = file1.split('.')[0]
    p2 = file2.split('.')[0]

    alignment1 = AlignIO.read(file1, 'fasta')
    alignment2 = AlignIO.read(file2, 'fasta')

    sequences1 = [x.seq for x in alignment1]
    sequences2 = [x.seq for x in alignment2]

    residues1 = range(0, len(sequences1[0]))
    residues2 = range(0,len(sequences2[0]))

    entropies1 = entropy_all_positions(sequences1)
    entropies2 = entropy_all_positions(sequences2)

    plotEntropies(entropies1, residues1, p1)
    plotEntropies(entropies2, residues2, p2)

#print(convert_01('ABCDEFGHIJ', 'AACDEEGHIJ'))
#split_fasta('strains_2859_2013to15.fasta')
#split_fasta('strains_230_2015.fasta')
#split_fasta('strains_865_2013.fasta')
#split_fasta('strains_310_90_to_04.fasta')
#split_fasta('strains_432_05_to_07.fasta')