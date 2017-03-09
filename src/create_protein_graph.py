#input: graphml of correlated mutations
#output: plot of protein level graph depicted nodes/edges
#for each of the 10 proteins

from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import networkx as nx
from graph_analysis_util import read_graphml, write_macro_graphml, get_top_node, write_cliques_to_csv
import matplotlib
import sys, os


#macro plot
def create_plot(g, folder, plot_file, title):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

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
    plt.title(title, fontsize=18)
    plot_file = folder + os.sep + plot_file
    plt.tight_layout()
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()

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


def run(folder, infilename, title):
    ingraph = read_graphml(folder + os.sep + infilename)
    outgraph = create_protein_graph(ingraph)
    write_macro_graphml(outgraph, folder, infilename)
    plot_file = infilename.split('.')[0] + '.png'
    create_plot(outgraph, folder, plot_file, title)

if __name__ == '__main__':
    folder     = sys.argv[1]
    infilename = sys.argv[2]
    title      = sys.argv[3]
    run(folder=folder, infilename=infilename, title=title)