from collections import defaultdict
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx

from graph_analysis_util import read_graphml, write_graphml

def create_plot(g, plot_file):
    graph_pos=nx.shell_layout(g)

    labels = [node+'_'+str(g.node[node]['count']) for node in g.nodes_iter()]
    print(labels)

    labelsD = {node:node+'('+str(g.node[node]['count'])+')' for node in g.nodes_iter()}
    print(labelsD)

    #g2=nx.relabel_nodes(g,labelsD)

    weights = [g[u][v]['weight'] for u, v in g.edges()]
    norm_weights = [10.0*float(i)/max(weights) for i in weights]

    # draw graph
    nx.draw_networkx_nodes(g,graph_pos, node_size=2400)
    nx.draw_networkx_edges(g,graph_pos,width=norm_weights)
    nx.draw_networkx_labels(g,graph_pos, labelsD)

    plt.axis('off')
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

    print(pedges)
    for key in pnodes.keys():
        outgraph.add_node(key)

    edges = combinations(pnodes.keys(), 2)
    outgraph.add_nodes_from(pnodes.keys())
    outgraph.add_edges_from(edges)
    nx.set_node_attributes(outgraph, 'count', pnodes)
    nx.set_edge_attributes(outgraph, 'weight', pedges)

    return outgraph

#compare nodes in two graphs
def compare_two_graphs(file1, file2):
    g1 = read_graphml(file1)
    g2 = read_graphml(file2)

    nodes1 = g1.nodes()
    nodes2 = g2.nodes()

    matches = []
    in1only = []
    in2only = []
    for node in nodes1:
        if node in nodes2:
            matches.append(node)
        else:
            in1only.append(node)

    for node in nodes2:
        if node not in nodes1:
            in2only.append(node)

    print(len(nodes1), len(nodes2), len(matches), len(in1only), len(in2only))

def run(infilename):
    ingraph = read_graphml(infilename)
    outgraph = create_protein_graph(ingraph)
    write_graphml(outgraph, infilename)
    plot_file = infilename.split('.')[0] + '.png'
    create_plot(outgraph, plot_file)

run('all_05_01.graphml')
#compare_two_graphs('all_05_01.graphml', '800d_05_01.graphml')