from Bio import SeqIO, AlignIO
from sys import argv
from numpy import transpose, array
from pandas import DataFrame
from minepy import MINE
from scipy.stats import mode
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import operator

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

        if ('HA' in description_list[5]):
            ha_sequences.append(record)

        if ('NA' in description_list[5]):
            na_sequences.append(record)

        if ('NP' in description_list[5]):
            np_sequences.append(record)

        if ('M1' in description_list[5]):
            m1_sequences.append(record)

        if ('M2' in description_list[5]):
            m2_sequences.append(record)

        if ('PB1' in description_list[5]):
            pb1_sequences.append(record)

        if ('PB2' in description_list[5]):
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

def get_consensus_sequence(sequences):
    raw_sequences = [x.seq for x in sequences]
    mod = [mode(y)[0][0] for y in transpose(array([list(z) for z in raw_sequences]))]
    consensus_sequence = ''.join(mod)
    return consensus_sequence


def convert_all_to_123(sequences):
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in sequences]))]
    consensus_sequence = ''.join(mod)
    sequences_123 = [convert_to_123(sequence, consensus_sequence) for sequence in sequences]
    # sequences_123 = [list(sequence_123) for sequence_123 in sequences_123_array]
    return sequences_123

def convert_to_123(sequence, consensus_sequence):
    new_sequence = [change_to_123(sequence[i], consensus_sequence[i]) for i in range(len(sequence))]
    return new_sequence

def get_central_node(G):
    centrality_dict = nx.degree_centrality(G)
    print('centrality_dict ', centrality_dict)
    central_nodes   = sorted(centrality_dict.items(), key=operator.itemgetter(1), reverse=True)
    print('sorted centrality_dict ', central_nodes)
    for node in central_nodes:
        print(node)
    central_nodesX  = [node[0] for node in central_nodes]
    return central_nodesX

def get_between_node(G):
    between_dict = nx.betweenness_centrality(G)
    print('between_dict', between_dict)
    between_nodes = sorted(between_dict.items(), key=operator.itemgetter(1), reverse=True)[:20]
    between_nodesX = [node[0] for node in between_nodes]
    return between_nodesX

def get_eigen_node(G):
    eigen_dict = nx.eigenvector_centrality(G)
    eigen_nodes = sorted(eigen_dict.items(), key=operator.itemgetter(1), reverse=True)[:20]
    eigen_nodesX = [node[0] for node in eigen_nodes]
    return eigen_nodesX

def createGraph(
                dataframe,
                labels=None, graph_layout='shell', node_size=1600, node_color='blue', node_alpha=0.3,
                node_text_size=12,
                edge_color='blue', edge_alpha=0.3, edge_tickness=1,
                edge_text_pos=0.3,
                text_font='sans-serif'):
    G = nx.from_pandas_dataframe(dataframe, 'x', 'y', 'mic')
    node_names={}

    print("graph info")
    print(nx.info(G))

    graph_pos=nx.shell_layout(G)

    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size,alpha=node_alpha,node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,font_family=text_font)

    plt.show()

    highest_degree = get_central_node(G)
    sub = G.subgraph(highest_degree)

    print("centrality sub-graph info")
    print(nx.info(sub))

    # draw graph
    nx.draw_networkx_nodes(sub,graph_pos,node_size=node_size,alpha=node_alpha,node_color=node_color)
    nx.draw_networkx_edges(sub,graph_pos,width=edge_tickness,alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(sub, graph_pos,font_size=node_text_size,font_family=text_font)

    plt.show()
    #----------------

    between_node = get_between_node(G)
    sub = G.subgraph(between_node)

    print("betweenness sub-graph info")
    print(nx.info(sub))

    # draw graph
    nx.draw_networkx_nodes(sub,graph_pos,node_size=node_size,alpha=node_alpha,node_color=node_color)
    nx.draw_networkx_edges(sub,graph_pos,width=edge_tickness,alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(sub, graph_pos,font_size=node_text_size,font_family=text_font)

    plt.show()

    '''
    eigen_node = get_eigen_node(G)
    sub = G.subgraph(eigen_node)

    print("eigen sub-graph info")
    print(nx.info(sub))

    graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(sub,graph_pos,node_size=node_size,alpha=node_alpha,node_color=node_color)
    nx.draw_networkx_edges(sub,graph_pos,width=edge_tickness,alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(sub,graph_pos,font_size=node_text_size,font_family=text_font)

    plt.show()

    return sub
    '''

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

# create sequences with 0,1,2,3 values for protein1 and protein2 sequences
def create_concat_sequence(file1, file2):
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

def createDFFromDict(dictionary):
    return DataFrame(dictionary)

def performMIC2p(p1_sequences_0123, p2_sequences_0123,p1,p2):
    mic_scores = []
    p1_sequences_0123T = transpose(array([list(z) for z in p1_sequences_0123])).tolist()
    p2_sequences_0123T = transpose(array([list(z) for z in p2_sequences_0123])).tolist()

    print(len(p1_sequences_0123T), len(p1_sequences_0123T[0]))
    print(len(p2_sequences_0123T), len(p2_sequences_0123T[0]))

    for idx1, record1 in enumerate(p1_sequences_0123T):
        for idx2, record2 in enumerate(p2_sequences_0123T):
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(record1, record2)
            if (mine.mic() > 0.2):
                mic_score = {}
                mic_score['x'] = p1+str(idx1)
                mic_score['y'] = p2+str(idx2)
                mic_score['mic'] = mine.mic()
                mic_scores.append(mic_score)
    return mic_scores

def run():
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    file1 = argv[1]
    file2 = argv[2]
    p1 = file1.split('.')[0]
    p2 = file2.split('.')[0]
    p1_sequences_0123, p2_sequences_0123 = create_concat_sequence(file1, file2)
    mic_scores = performMIC2p(p1_sequences_0123, p2_sequences_0123, p1, p2)
    micDF = createDFFromDict(mic_scores)
    createGraph(micDF)

run()
#split_fasta('strains_2859_2013to15.fasta')
