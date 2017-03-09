# perform pairwise comparision of residues at two positions
# in p1, p2 to get an understanding of pair wise metrics;
# example: postion 10 in p1 and position 32 in p2
# W,C = 1320; S, C = 4550; V, N = 8; rest = 0 total sequences = 1320+4550+8
from sympy.functions.elementary.complexes import sign

import graph_analysis_util as util
import os, sys
from Bio import AlignIO
from numpy import transpose, array
import plot_util as pu


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

    return residue_combination_count


def get_best_edges_from_graph(g, n):
    highest_degree_nodes = util.get_highest_degree(g, n)

    best_edges = []
    for node in highest_degree_nodes:
        best_edges_this_node = [(u, v, d['weight']) for (u, v, d) in g.edges(node, data=True) if d['weight'] >= 0.9]
        for edge in best_edges_this_node:
            best_edges.append(edge)

    best_edges.sort(key=lambda x:x[2], reverse=True)

    return best_edges[:10]


def run_pairwise_counts_manual_edges():
    folder = sys.argv[2]
    significant_residues_file = 'pairs.txt'
    lines = util.read_file(significant_residues_file)
    logo_folder = 'C:\\uday\\gmu\\correlations\\results\\10proteins\\logos'
    for line in lines:
        print(line)
        protein1 = line[0]
        residue1 = int(line[1])-1
        protein2 = line[2]
        residue2 = int(line[3])-1
        file1 = folder + os.sep + protein1 + '.afasta'
        file2 = folder + os.sep + protein2 + '.afasta'
        residue_comparision_count = perform_residue_analysis(file1, file2, residue1, residue2)
        residue1_str = protein1+str(residue1+1)
        residue2_str = protein2+str(residue2+1)
        filename = residue1_str + '_' + residue2_str + '.bat'
        pu.create_image_magick_script(residue_comparision_count, residue1_str, residue2_str, logo_folder, filename)


def run_pairwise_counts_edges(start, end, reverse):
    rcc_list = []
    graphml_file = sys.argv[1]
    data_folder = sys.argv[2]

    g = util.read_graphml(graphml_file)
    top_n_edges = util.get_best_edges(g, start=start, end=end, reverse=reverse)

    #logo_folder = 'C:\\uday\\gmu\\correlations\\results\\10proteins\\logos'

    for edge in top_n_edges:
        residue1 = int(edge[0].split('_')[1])-1
        protein1 = edge[0].split('_')[0]
        residue2 = int(edge[1].split('_')[1])-1
        protein2 = edge[1].split('_')[0]
        file1 = data_folder + os.sep + protein1 + '.afasta'
        file2 = data_folder + os.sep + protein2 + '.afasta'
        rcc = perform_residue_analysis(file1, file2, residue1, residue2)
        print(rcc)
        print(protein1, protein2, residue1+1, residue2+1)
        rcc_list.append(rcc)
        #residue1_str = protein1 + str(residue1+1)
        #residue2_str = protein2 + str(residue2+1)
        #filename = residue1_str + '_' + residue2_str + '.bat'
        #pu.create_image_magick_script(rcc, residue1_str, residue2_str, logo_folder, filename)
    return rcc_list


def create_pairwise_count_plots_top_edges():
    out_folder = sys.argv[3]
    rcc_list = run_pairwise_counts_edges(start=0, end=50, reverse=True)
    pu.plot_residue_counts(rcc_list=rcc_list, title="residue combination count plot for top 50 edges (mic ~0.99)", folder=out_folder, mic='0.99')
    rcc_list = run_pairwise_counts_edges(start=100000, end=100050, reverse=True)
    pu.plot_residue_counts(rcc_list=rcc_list, title="residue combination count plot for edges with mic ~0.8", folder=out_folder, mic='0.8')
    rcc_list = run_pairwise_counts_edges(start=200000, end=200050, reverse=True)
    pu.plot_residue_counts(rcc_list=rcc_list, title="residue combination count plot for edges with mic ~0.51", folder=out_folder, mic='0.51')
    rcc_list = run_pairwise_counts_edges(start=240000, end=240050, reverse=True)
    pu.plot_residue_counts(rcc_list=rcc_list, title="residue combination count plot for edges with mic ~0.32", folder=out_folder, mic='0.32')


create_pairwise_count_plots_top_edges()
#run_pairwise_counts_manual_edges()