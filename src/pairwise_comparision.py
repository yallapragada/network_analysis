# perform pairwise comparision of residues at two positions with relatively high mic
# in p1, p2 to get an understanding of pair wise metrics;
# example: postion 10 in p1 and position 32 in p2
# W,C = 1320; S, C = 4550; V, N = 8; rest = 0 total sequences = 1320+4550+8

import graph_analysis_util as util
import os, sys
from Bio import AlignIO
from numpy import transpose, array


def perform_residue_analysis(file1, file2, position1, position2):
    residue_combination_count = {}
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    sequences1_t = transpose(array([list(z) for z in sequences1])).tolist()
    sequences2_t = transpose(array([list(z) for z in sequences2])).tolist()

    zipped_list = list(zip(sequences1_t[position1], sequences2_t[position2]))

    for element in zipped_list:
        residue_combination = str(element[0]+element[1])
        if residue_combination in residue_combination_count.keys():
            residue_combination_count[residue_combination] = int(residue_combination_count[residue_combination])+1
        else:
            residue_combination_count[residue_combination] = 1

    print(residue_combination_count)
    print('--------------------------')


def get_best_edges_from_graph(g, n):
    highest_degree_nodes = util.get_highest_degree(g, n)

    best_edges = []
    for node in highest_degree_nodes:
        best_edges_this_node = [(u, v, d['weight']) for (u, v, d) in g.edges(node, data=True) if d['weight'] >= 0.9]
        for edge in best_edges_this_node:
            best_edges.append(edge)

    best_edges.sort(key=lambda x:x[2], reverse=True)

    return best_edges[:10]


def run_pairwise_comparision():
    graphml_file = sys.argv[1]
    folder = sys.argv[2]

    g = util.read_graphml(graphml_file)
    top_ten_edges = get_best_edges_from_graph(g, 10)

    for edge in top_ten_edges:
        residue1 = int(edge[0].split('_')[1])-1
        protein1 = edge[0].split('_')[0]
        residue2 = int(edge[1].split('_')[1])-1
        protein2 = edge[1].split('_')[0]
        print(protein1, protein2, residue1, residue2)
        file1 = folder + os.sep + protein1 + '.afasta'
        file2 = folder + os.sep + protein2 + '.afasta'
        perform_residue_analysis(file1, file2, residue1, residue2)


run_pairwise_comparision()
