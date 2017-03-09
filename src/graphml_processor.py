from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import networkx as nx
from graph_analysis_util import read_graphml, write_macro_graphml, get_top_node, write_cliques_to_csv
from plot_util import plot_mic_histogram
import sys
import os


#degree distribution
def create_degree_plot(g):
    degree_sequence = sorted(nx.degree(g).values(), reverse=True)  # degree sequence

    plt.loglog(degree_sequence, 'b-', marker='o')
    plt.title("Degree rank plot")
    plt.ylabel("degree")
    plt.xlabel("rank")

    #plt.savefig("degree_histogram.png")
    plt.show()


#macro plot
def create_plot(g, plot_file, title):
    colors = ['r', 'b', 'g', 'c', 'm', 'r', 'b', 'g', 'c', 'm']
    graph_pos=nx.shell_layout(sorted(g.nodes()))

    labels = [node+'_'+str(g.node[node]['count']) for node in g.nodes_iter()]
    labelsD = {node:node+'('+str(g.node[node]['count'])+')' for node in g.nodes_iter()}

    #g2=nx.relabel_nodes(g,labelsD)

    weights = [d['weight'] if 'weight' in d else 0 for u, v, d in g.edges(data=True)]
    norm_weights = [10.0*(float(i)-min(weights))/max(weights) for i in weights]

    # draw graph
    nx.draw_networkx_nodes(g,graph_pos, node_size=2400, node_color=colors)
    nx.draw_networkx_edges(g,graph_pos,width=norm_weights)
    nx.draw_networkx_labels(g,graph_pos, labelsD)

    plt.axis('off')
    plt.title(title)
    plt.savefig(plot_file, bbox_inches='tight')


def create_protein_graph(ingraph):
    pnodes = {}
    pnodes = defaultdict(lambda: 0, pnodes)

    pedges = {}
    pedges = defaultdict(lambda: 0, pedges)

    outgraph = nx.Graph()

    for node in ingraph.nodes_iter():
        pnodes[ingraph.node[node]['protein']]+=1

    for u,v,d in ingraph.edges(data=True):
        key=(ingraph.node[u]['protein'], ingraph.node[v]['protein'])
        pedges[key]+=1

    for key in pnodes.keys():
        outgraph.add_node(key)

    edges = combinations(pnodes.keys(), 2)

    outgraph.add_nodes_from(pnodes.keys())
    outgraph.add_edges_from(edges)
    nx.set_node_attributes(outgraph, 'count', pnodes)
    nx.set_edge_attributes(outgraph, 'weight', pedges)

    return outgraph

def draw_histogram_big_file(big_infile, out_folder, dataset):
    mics_all=[]
    with open(big_infile, 'r') as inputfile:
        for line in inputfile:
            if '<data key="d1">' in line:
                mic=line[line.find('>') + 1: line.find('</')]
                mic_float=float(mic)
                mics_all.append(mic_float)
    plot_file = big_infile.split('.')[0]
    plot_mic_histogram(mics=mics_all, folder=out_folder, file=plot_file, dataset=dataset, suffix='01', bins=40)


def run(folder, infilename, title):
    ingraph = read_graphml(folder + os.sep + infilename)
    try:
        #avg_cluster_coeff = nx.average_clustering(ingraph)
        #print('average clustering for ' + title + " = " + str(avg_cluster_coeff))

        avg_deg_coeff = nx.average_degree_connectivity(ingraph)
        print('average degree for ' + title + ' = ' + str(avg_deg_coeff))
    except Exception as e:
        print (e.__str__())


    '''
    outgraph = create_protein_graph(ingraph)
    write_macro_graphml(outgraph, folder, infilename)
    plot_file = infilename.split('.')[0] + '.png'
    create_plot(outgraph, folder, plot_file, title)
    '''

if __name__ == '__main__':
    folder     = sys.argv[1]
    infilename = sys.argv[2]
    title      = sys.argv[3]
    run(folder=folder, infilename=infilename, title=title)

#draw_histogram_big_file('uh7n9_01_01.graphml', 'C:\\uday\\gmu\\correlations\\results\\10proteins\\histograms')
#draw_histogram_big_file('400a_001_01.graphml', 'C:\\Users\\uday\\pycharm_projects\\network_analysis\\data')
#run('u800d_05_B01b.graphml')
#compare_two_graphs('all_05_01.graphml', '800d_05_01.graphml')
#get_max_clique(sys.argv[1])
#create_pajek(sys.argv[1])
#create_clique_top_node(sys.argv[1])