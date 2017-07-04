import matplotlib.pyplot as plt
import matplotlib
import networkx as nx
from graph_analysis_util import read_graphml
import os
import sys

# script to create degree and clustering plots from a graphml
# usage: python create_degree_clustering_plots.py [output_folder_location] [path to graphml] [title for plot]

#degree distribution
def create_clustering_plot(g, folder, title):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    clustering_sequence = sorted(nx.clustering(g).values(), reverse=True)  # degree sequence
    plt.hist(clustering_sequence, color="#3F5D7D")
    plt.xlim([min(clustering_sequence), max(clustering_sequence)])
    plt.title("Clustering histogram for " + title, fontsize=18)
    plt.ylabel("frequency", fontsize=14)
    plt.xlabel("clustering coefficient", fontsize=14)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plot_file = folder + os.sep + 'CL_' + title + '.png'
    plt.savefig(plot_file,format='png',bbox_inches="tight")
    plt.close()


#degree distribution
def create_degree_plot(g, folder, title):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    degree_sequence = sorted(nx.degree(g).values(), reverse=True)  # degree sequence
    print(degree_sequence)
    plt.hist(degree_sequence, color="#3F5D7D")
    plt.xlim([min(degree_sequence), max(degree_sequence)])
    plt.title("Degree histogram for " + title, fontsize=18)
    plt.ylabel("frequency", fontsize=14)
    plt.xlabel("degree", fontsize=14)
    plt.tick_params(axis='both', labelsize=12)
    plt.tight_layout()
    plot_file = folder + os.sep + title + '.png'
    plt.savefig(plot_file,format='png',bbox_inches="tight")
    plt.close()


def run(folder, infilename, title):
    ingraph = read_graphml(folder + '\\' + infilename)
    create_clustering_plot(ingraph, folder=folder, title=title)
    create_degree_plot(ingraph, folder=folder, title=title)


if __name__ == '__main__':
    folder      = sys.argv[1]
    infilename  = sys.argv[2]
    title       = sys.argv[3]
    run(folder=folder, infilename=infilename, title=title)
