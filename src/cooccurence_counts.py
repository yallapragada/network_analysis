# perform pairwise comparision of residues at two positions
# in p1, p2 to get an understanding of pair wise metrics;
# example: postion 10 in p1 and position 32 in p2
# W,C = 1320; S, C = 4550; V, N = 8; rest = 0 total sequences = 1320+4550+8


import graph_analysis_util as util
import os, sys
from Bio import AlignIO
from numpy import transpose, array
import plot_util as pu
import operator


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


#cooccurence counts for top 50 and bottom 50 edges in a graphml
def get_cooccurences(graphml_file, fasta_folder):
    top_50_edges_cooccurence_counts     = []
    bottom_50_edges_cooccurence_counts  = []

    g = util.read_graphml(graphml_file)
    top_n_edges = util.get_edges(g, start=0, end=50, reverse=True)
    bottom_n_edges = util.get_edges(g, start=0, end=50, reverse=False)

    for edge in top_n_edges:
        residue1 = int(edge[0].split('_')[1])-1
        protein1 = edge[0].split('_')[0]
        residue2 = int(edge[1].split('_')[1])-1
        protein2 = edge[1].split('_')[0]
        file1 = fasta_folder + os.sep + protein1 + '.afasta'
        file2 = fasta_folder + os.sep + protein2 + '.afasta'
        rcc = perform_residue_analysis(file1, file2, residue1, residue2)
        print(edge)
        print(sorted(rcc.items(), key=operator.itemgetter(1), reverse=True))
        top_50_edges_cooccurence_counts.append(rcc)

    print('=====================')
    for edge in bottom_n_edges:
        residue1 = int(edge[0].split('_')[1])-1
        protein1 = edge[0].split('_')[0]
        residue2 = int(edge[1].split('_')[1])-1
        protein2 = edge[1].split('_')[0]
        file1 = fasta_folder + os.sep + protein1 + '.afasta'
        file2 = fasta_folder + os.sep + protein2 + '.afasta'
        rcc = perform_residue_analysis(file1, file2, residue1, residue2)
        print(sorted(rcc.items(), key=operator.itemgetter(1), reverse=True))
        bottom_50_edges_cooccurence_counts.append(rcc)

    return top_50_edges_cooccurence_counts, bottom_50_edges_cooccurence_counts

#cooccurences for top 50 and bottom 50 in a graphml
def plot_cooccurences(graphml_file, fasta_folder, dataset):
    out_folder = 'C:\\uday\\gmu\\correlations\\results\\10proteins\\cooccurence'
    top_50, bottom_50 = get_cooccurences(graphml_file=graphml_file, fasta_folder=fasta_folder)

    title = "cooccurence plot for top 50 edges in " + dataset
    filename = dataset + 'T'
    pu.plot_residue_counts(rcc_list=top_50, title=title, folder=out_folder, filename=filename)

    title = "cooccurence plot for bottom 50 edges in " + dataset
    filename = dataset + 'B'
    pu.plot_residue_counts(rcc_list=bottom_50, title=title, folder=out_folder, filename=filename)


def run_cooccurence_all_datasets():
    #input_file = 'location of cooccurence_input.txt'
    input_file = sys.argv[1]
    lines = util.read_file(input_file)
    for line in lines:
        fasta_folder = line[0]
        graphml_file = line[1]
        dataset = line[2]
        plot_cooccurences(graphml_file=graphml_file, fasta_folder=fasta_folder, dataset=dataset)
        print('************************')

if __name__ == '__main__':
    run_cooccurence_all_datasets()
