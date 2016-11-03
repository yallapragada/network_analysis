import csv
import datetime
import math
import operator
import sys

import networkx as nx
from Bio import SeqIO, AlignIO
from minepy import MINE
from numpy import transpose, array
from pandas import DataFrame
from scipy.stats import mode

from plot_util import plotEntropies


def load_blosum_matrix(matrix_filename):

    with open(matrix_filename) as matrix_file:
        matrix_content = matrix_file.read()
    lines = matrix_content.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    blosum_matrix = {}

    for row in lines:
        entries = row.split()
        row_name = entries.pop(0)
        blosum_matrix[row_name] = {}

        if len(entries) != len(columns):
            print('Improper entry number in row')

        for column_name in columns:
            blosum_matrix[row_name][column_name] = entries.pop(0)

    return blosum_matrix


def lookup_score(blosum_matrix, a, b):
    a = a.upper()
    b = b.upper()

    if a not in blosum_matrix or b not in blosum_matrix[a]:
        print(a, ' or ', b, ' not in matrix ')

    x = blosum_matrix[a][b]
    return x

# this function is partitioning the input list into n smaller chunks
def partition(list, n):
    q, r = divmod(len(list), n)
    indices = [q*i + min(i, r) for i in range(n+1)]
    return [list[indices[i]:indices[i + 1]] for i in range(n)]


# this function is performing a simple split of the complete file
# does not split by protein
def split_fasta_file_to_smaller(fasta_file, num_chunks):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append(record)
    random.shuffle(records)
    chunks = partition(records, num_chunks)
    for i in range(len(chunks)):
        output_handle = open("chunk"+str(i)+".fasta", "w")
        SeqIO.write(chunks[i], output_handle, "fasta")
        output_handle.close()


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
    ns1_sequences = []
    ns2_sequences = []
    pa_sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        description_list = record.description.split('|')

        if (len(description_list)<4):
            print("protein description too short ", record.description)
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

        if 'PB1 Polymerase (basic) protein 1' in description_list[5] or 'PB1 Polymerase (basic) protein 1' in description_list[4]:
            pb1_sequences.append(record)

        if 'PB2' in description_list[5]:
            pb2_sequences.append(record)

        if 'NS1 Non-structural protein 1' in description_list[5]:
            ns1_sequences.append(record)

        if 'NS2 Non-structural protein 2' in description_list[5]:
            ns2_sequences.append(record)

        if 'PA Polymerase (acidic) protein' in description_list[5]:
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


# blosum62 based encoding >=0 implies 0; else 1
def convert_01_B62b(sequence, consensus_sequence, b62_matrix):
    new_sequence = [change_to_01_B62b(sequence[i], consensus_sequence[i], b62_matrix) for i in range(len(sequence))]
    return new_sequence


def change_to_01_B62b(residue, consensus_residue, b62_matrix):
    b62_value = lookup_score(b62_matrix, residue, consensus_residue)
    if int(b62_value) > 0:
        return 0
    else:
        return 1


# blosum62 based encoding >=0 implies 0; else 1
def convert_01_B62a(sequence, consensus_sequence, b62_matrix):
    new_sequence = [change_to_01_B62a(sequence[i], consensus_sequence[i], b62_matrix) for i in range(len(sequence))]
    return new_sequence


def change_to_01_B62a(residue, consensus_residue, b62_matrix):
    b62_value = lookup_score(b62_matrix, residue, consensus_residue)
    if int(b62_value) >= 0:
        return 0
    else:
        return 1

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


# get top n nodes with highest degree
def get_highest_degree(G, n):
    centrality_dict = nx.degree_centrality(G)
    degree_centrality_nodes   = sorted(centrality_dict.items(), key=operator.itemgetter(1), reverse=True)
    if (len(degree_centrality_nodes)> int(n)):
        central_nodesX = [node[0] for node in degree_centrality_nodes[:int(n)]]
    else:
        central_nodesX = [node[0] for node in degree_centrality_nodes]
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

def write_strains_to_csv(list_of_strains, title):

    csv.register_dialect(
        'mydialect',
        delimiter=',',
        quotechar='"',
        doublequote=True,
        skipinitialspace=True,
        lineterminator='\n',
        quoting=csv.QUOTE_MINIMAL)

    f = open(title + ".csv", 'wt')

    try:
        writer = csv.writer(f, dialect='mydialect')
        writer.writerow(list_of_strains)
    finally:
        f.close()


# Given a SeqRecord, return strain_name
def get_strain_name(record):
    strain_full_str = record.description.split('|')[3]
    if (strain_full_str.startswith('Strain')):
        strain = strain_full_str.split(':')[1]
    else:
        strain = record.description.split('|')[4].split(':')[1]
    return strain


def get_matching_sequence(records, strain_name):
    for record in records:
        if (get_strain_name(record) == strain_name):
            return record
    return None


def create_01_sequences(file1, file2):
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


# create 0,1 sequences based on blosum 62
def create_B01a_sequences(file1, file2):
    blosum_matrix = load_blosum_matrix('blosum62.txt')
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    p1_sequences_B01a = []  # list of p1 sequences
    p2_sequences_B01a = []  # list of p2 sequences

    consensus_sequence1 = get_consensus_sequence(sequences1)
    consensus_sequence2 = get_consensus_sequence(sequences2)

    for sequence1 in sequences1:
        strain_name = get_strain_name(sequence1)
        sequence2 = get_matching_sequence(sequences2, strain_name=strain_name)
        if (sequence2):
            p1_sequences_B01a.append(convert_01_B62a(sequence1.seq, consensus_sequence1, blosum_matrix))
            p2_sequences_B01a.append(convert_01_B62a(sequence2.seq, consensus_sequence2, blosum_matrix))

    return p1_sequences_B01a, p2_sequences_B01a


# create 0,1 sequences based on blosum 62
def create_B01b_sequences(file1, file2):
    blosum_matrix = load_blosum_matrix('blosum62.txt')
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    p1_sequences_B01b = []  # list of p1 sequences
    p2_sequences_B01b = []  # list of p2 sequences

    consensus_sequence1 = get_consensus_sequence(sequences1)
    consensus_sequence2 = get_consensus_sequence(sequences2)

    for sequence1 in sequences1:
        strain_name = get_strain_name(sequence1)
        sequence2 = get_matching_sequence(sequences2, strain_name=strain_name)
        if (sequence2):
            p1_sequences_B01b.append(convert_01_B62b(sequence1.seq, consensus_sequence1, blosum_matrix))
            p2_sequences_B01b.append(convert_01_B62b(sequence2.seq, consensus_sequence2, blosum_matrix))

    return p1_sequences_B01b, p2_sequences_B01b


def create_df_from_dict(dictionary):
    return DataFrame(dictionary)

#how many nodes per protein
#how many edges between proteins
def create_report(df):
    proteins = ['ha', 'na', 'np', 'm1', 'm2', 'pb1', 'pb2', 'ns1', 'ns2', 'pa']
    for protein in proteins:
        df1 = list(df['x'][(df['p1'] == protein)].values.ravel())
        df2 = list(df['y'][(df['p2'] == protein)].values.ravel())
        df3=list(set(df1+df2))
        print('# unique residues of ', protein, ' = ', len(df3))

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
                mic_score['p1'] = p1
                mic_score['p2'] = p2
                mic_score['weight'] = mine.mic()
                mic_scores.append(mic_score)

    #print('computed ', len(mic_scores), ' mics for ', p1, p2, 'for cutoff ', cutoff)
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
    return scores


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
#split_fasta_file_to_smaller('strains_800_name_A.fasta', 2)
#split_fasta('strains_random_400.fasta')
#split_fasta('strains_800_name_A_subset.fasta')
#split_fasta('strains_800_random.fasta')
#split_fasta('all_strains_comp_gen.fasta')

def test_dict():
    dict1 = {'x':'ha_10', 'y':'na_11', 'p1':'ha', 'p2':'na', 'weight':'0.8'}
    dict2 = {'x':'ha_11', 'y':'na_11', 'p1':'ha', 'p2':'na', 'weight':'0.8'}
    dict3 = {'x': 'ha_10', 'y': 'na_13', 'p1': 'ha', 'p2': 'na', 'weight': '0.8'}
    dict4 = {'x': 'ha_14', 'y': 'na_1', 'p1': 'ha', 'p2': 'na', 'weight': '0.8'}
    dict5 = {'x': 'ha_15', 'y': 'na_2', 'p1': 'ha', 'p2': 'na', 'weight': '0.8'}
    dict6 = {'x': 'ha_21', 'y': 'na_13', 'p1': 'ha', 'p2': 'na', 'weight': '0.8'}
    dict7 = {'x': 'ha_21', 'y': 'na_2', 'p1': 'ha', 'p2': 'na', 'weight': '0.8'}
    dict8 = {'x': 'np_3', 'y': 'na_2', 'p1': 'np', 'p2': 'na', 'weight': '0.8'}
    dict9 = {'x': 'np_4', 'y': 'na_2', 'p1': 'np', 'p2': 'na', 'weight': '0.8'}
    dict10 = {'x': 'np_5', 'y': 'ha_11', 'p1': 'np', 'p2': 'ha', 'weight': '0.8'}

    dictionaries = []
    dictionaries.append(dict1)
    dictionaries.append(dict2)
    dictionaries.append(dict3)
    dictionaries.append(dict4)
    dictionaries.append(dict5)
    dictionaries.append(dict6)
    dictionaries.append(dict7)
    dictionaries.append(dict8)
    dictionaries.append(dict9)
    dictionaries.append(dict10)

    df=create_df_from_dict(dictionaries)
    create_report(df)

def read_csv_strains(file):
    strains = list(read_file(file))[0]


#read_csv('unique_strains_100.csv')

#test_dict()

#def test_convert_B62a():
    sequence = 'EGRMNYYWTLVEPGDKITFEA'
    consensus_sequence = 'EGRMNYYWTGGGPGDKITFEA'
    blosum_matrix = load_blosum_matrix('blosum62.txt')
    print(blosum_matrix)
    print(convert_01_B62a(sequence, consensus_sequence, blosum_matrix))

#test_convert_B62a()


def read_graphml(infilename):
    graph = nx.read_graphml(infilename)
    return graph


def write_graphml(outgraph, infilename):
    nx.write_graphml(outgraph, 'p_'+infilename)