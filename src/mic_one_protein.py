import sys
from Bio import AlignIO
from Bio import SeqIO
from minepy import MINE
from numpy import array, transpose, var
from odo.convert import higher_precision_freqs
from scipy.stats import mode
from pandas import DataFrame
import networkx as nx
import matplotlib.pyplot as plt
import operator

def generate_binary_sequence(sequences):
    """ Generates a binary sequence out of a fasta MSA
    """
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in sequences]))]
    binary_sequence_array = array([[1 if x==mod[i] else 0 for i,x in enumerate(y)] for y in sequences])
    binary_sequences = [list(binary_sequence) for binary_sequence in binary_sequence_array]
    return binary_sequences

def conver_all_to_123(sequences):
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in sequences]))]
    consensus_sequence = ''.join(mod)
    sequences_123 = [convert_to_123(sequence, consensus_sequence) for sequence in sequences]
    #sequences_123 = [list(sequence_123) for sequence_123 in sequences_123_array]
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

def remove_zeros(data):
    dataT = transpose(data)
    print(dataT.shape)
    indices = [index for index, x in enumerate(dataT) if var(x) != 0]
    for i, index in enumerate(indices):
        index_map[i] = index
    data_trimmed = transpose(array([x for x in dataT if var(x) != 0]))
    return data_trimmed

#use biopython SeqIO to parse fasta file
#returns a list of sequences
def get_sequences(fasta_file):
    sequences = [x.seq for x in SeqIO.parse(fasta_file, "fasta")]
    return sequences

def get_sequences_from_alignment(alignment):
    sequences = [x.seq for x in alignment]
    return sequences

def perform_transpose(binary_sequences):
    transposed_list = transpose(array(binary_sequences)).tolist()
    return transposed_list

def performMIC(transposed_list):
    mic_scores=[]
    for counter1 in range(0, len(transposed_list)-1):
        for counter2 in range(counter1+1, len(transposed_list)):
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(transposed_list[counter1], transposed_list[counter2])
            if (mine.mic()>0.6):
                mic_score={}
                mic_score['x']=counter1
                mic_score['y']=counter2
                mic_score['mic']=mine.mic()
                mic_scores.append(mic_score)
    return mic_scores

def createDFFromDict(dictionary):
    df = DataFrame(dictionary)
    return DataFrame(dictionary)

def centrality_sort(centrality_dict):
    return sorted(centrality_dict.items(), key=operator.itemgetter(1), reverse=True)

def createGraph(dataframe):
    G = nx.from_pandas_dataframe(dataframe, 'x', 'y', 'mic')
    node_names={}
    #convert 1.0, 2.0.. to 1, 2 for node names
    for i, node in enumerate(G.nodes()):
        node_names[node]=str(int(node))
    G = nx.relabel_nodes(G, node_names)
    print("graph info")
    print(nx.info(G))
    degree_sorted = centrality_sort(nx.degree_centrality(G))
    highest_degree = [node[0] for node in degree_sorted[-20:]]
    sub = G.subgraph(highest_degree)
    print("sub-graph info")
    print(nx.info(sub))
    return sub

index_map = {}

def run_mic():
    alignment = AlignIO.read(sys.argv[1], 'fasta')
    sequences = get_sequences_from_alignment(alignment)
    binary_sequences = generate_binary_sequence(sequences)
    print(len(binary_sequences), len(binary_sequences[0]))
    binary_sequences = remove_zeros(binary_sequences)
    print(len(binary_sequences), len(binary_sequences[0]))
    transposed_list = perform_transpose(binary_sequences)
    print(len(transposed_list), len(transposed_list[0]))
    mic_scores=performMIC(transposed_list)
    micDF=createDFFromDict(mic_scores)
    print(micDF.head(5))
    G=createGraph(micDF)
    nx.draw_networkx(G)
    plt.axis('off')
    plt.show()

def run_mic_123():
    alignment = AlignIO.read(sys.argv[1], 'fasta')
    sequences = get_sequences_from_alignment(alignment)
    sequences123 = conver_all_to_123(sequences)
    print(len(sequences123), len(sequences123[0]))
    sequences123 = remove_zeros(sequences123)
    print(len(sequences123), len(sequences123[0]))
    transposed_list = transpose(sequences123)
    print(len(transposed_list), len(transposed_list[0]))
    mic_scores=performMIC(transposed_list)
    micDF=createDFFromDict(mic_scores)
    print(micDF.head(5))
    G=createGraph(micDF)
    nx.draw_networkx(G)
    plt.axis('off')
    plt.show()

def unit_test():
    print(change_to_123('Y', 'Y'))
    print(change_to_123('W', 'A'))
    print(change_to_123('M', 'H'))
    print(change_to_123('C', 'P'))
    print(convert_to_123('AWCIMMAAYGSVQ', 'ADNTMVAILGSVY'))

run_mic_123()
