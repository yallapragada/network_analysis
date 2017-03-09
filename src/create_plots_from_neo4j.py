#crete node/edge count vs cut-off plots
#generates bokeh html that is included in web application

from py2neo import Graph
import matplotlib
import matplotlib.pyplot as plt
from bokeh.io import output_file, show, gridplot, save
from bokeh.plotting import figure
from graph_analysis_util import read_datasets_csv

# get handle to local neo4j instance
def get_connection_to_neo4j():
    graph = Graph(password='flu')
    return graph


TEMPLATE_FOLDER = 'C:\\Users\\uday\\pycharm_projects\\network_analysis_web\\corrmut\\templates\\'
PNG_FOLDER      = 'C:\\uday\\gmu\\correlations\\results\\10proteins\\node_edge_counts\\'

graph = get_connection_to_neo4j()
datasets = read_datasets_csv()

#get number of edges in the network for different cutoff values
def get_edge_counts_for_cutoffs(dataset):
    cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98]
    edge_counts = []
    for cutoff in cutoffs:
        query_part1 = 'MATCH (n:' + dataset
        query_part2 = ''')-[r]->(m)
        where toFloat(r.mic)>{cutoff}
        return count(r) as edges
        '''
        query = query_part1 + query_part2
        results = graph.data(query, cutoff=cutoff, dataset=dataset)
        edge_counts.append(results[0]['edges'])
    return cutoffs, edge_counts

# get number of nodes in the network for different cutoff values
def get_node_counts_for_cutoffs(dataset):
    cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98]
    node_counts = []
    for cutoff in cutoffs:
        query_part1 = 'MATCH (n:' + dataset
        query_part2 = ''')-[r]-(m)
        where toFloat(r.mic)>{cutoff}
        return count(distinct(n)) as nodes
        '''
        query = query_part1 + query_part2
        results = graph.data(query, cutoff=cutoff, dataset=dataset)
        node_counts.append(results[0]['nodes'])
    return cutoffs, node_counts


def plot_edgecounts_matplotlib_sub(datasets, colors, filename):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    i=0
    fig, ax = plt.subplots(1)
    for dataset in datasets:
        ax.set_xlim([0.1, 1])
        x,y=get_edge_counts_for_cutoffs(dataset=dataset)
        ax.plot(x, y, label=dataset, color=colors[i], linewidth=2.0)
        i=i+1
        ax.set_title('Edge counts for ' + filename, fontsize=18)
        ax.set_xlabel("mic cut-off", fontsize=14)
        ax.set_ylabel("# edges", fontsize=14)
        ax.legend(loc='best', fontsize=14)
    plt.tight_layout()
    plt.savefig(PNG_FOLDER + filename + '_edges.png', format='png', bbox_inches="tight")
    plt.close()

def plot_edgecounts_matplotlibe_main():

    colors = ["#3F5D7D", "#AA5D66", "#CD5D93"]

    datasets_csv = read_datasets_csv()
    datasets = sorted([dataset['name'] for dataset in datasets_csv])

    for dataset in datasets:
        plot_edgecounts_matplotlib_sub([dataset], colors, dataset)

    datasets = []
    filename = 'HUMAN_H3N2'
    for dataset in datasets_csv:
        if dataset['name'].startswith('HUMAN_H3N2'):
            datasets.append(dataset['name'])
    plot_edgecounts_matplotlib_sub(datasets, colors, filename='HUMAN_H3N2')

    datasets = []
    filename = 'HUMAN_H1N1'
    for dataset in datasets_csv:
        if dataset['name'].startswith('HUMAN_H1N1'):
            datasets.append(dataset['name'])
    plot_edgecounts_matplotlib_sub(datasets, colors, filename='HUMAN_H1N1')

    datasets = []
    filename = 'SWINE'
    for dataset in datasets_csv:
        if dataset['name'].startswith('SWINE'):
            datasets.append(dataset['name'])
    plot_edgecounts_matplotlib_sub(datasets, colors, filename='SWINE')

def plot_nodecounts_matplotlib_sub(datasets, colors, filename):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    i=0
    fig, ax = plt.subplots(1)
    for dataset in datasets:
        ax.set_xlim([0.1, 1])
        x,y=get_node_counts_for_cutoffs(dataset=dataset)
        ax.plot(x, y, label=dataset, color=colors[i], linewidth=2.0)
        i=i+1
        ax.set_title('Node counts for ' + filename, fontsize=18)
        ax.set_xlabel("mic cut-off", fontsize=14)
        ax.set_ylabel("# nodes", fontsize=14)
        ax.legend(loc='best', fontsize=14)
    plt.tight_layout()
    plt.savefig(PNG_FOLDER + filename + '_nodes.png', format='png', bbox_inches="tight")
    plt.close()

def plot_nodecounts_matplotlibe_main():

    colors = ["#3F5D7D", "#AA5D66", "#CD5D93"]

    datasets_csv = read_datasets_csv()
    datasets = sorted([dataset['name'] for dataset in datasets_csv])

    for dataset in datasets:
        plot_nodecounts_matplotlib_sub([dataset], colors, dataset)

    datasets = []
    filename = 'HUMAN_H3N2'
    for dataset in datasets_csv:
        if dataset['name'].startswith('HUMAN_H3N2'):
            datasets.append(dataset['name'])
    plot_nodecounts_matplotlib_sub(datasets, colors, filename='HUMAN_H3N2')

    datasets = []
    filename = 'HUMAN_H1N1'
    for dataset in datasets_csv:
        if dataset['name'].startswith('HUMAN_H1N1'):
            datasets.append(dataset['name'])
    plot_nodecounts_matplotlib_sub(datasets, colors, filename='HUMAN_H1N1')

    datasets = []
    filename = 'SWINE'
    for dataset in datasets_csv:
        if dataset['name'].startswith('SWINE'):
            datasets.append(dataset['name'])
    plot_nodecounts_matplotlib_sub(datasets, colors, filename='SWINE')

#simple line plot for plotting cutoff vs #edges
def plot_edgecounts_bokeh_line():
    datasets_csv = read_datasets_csv()
    datasets = sorted([dataset['name'] for dataset in datasets_csv])
    plots = []
    for dataset in datasets:
        output_file(PNG_FOLDER + dataset + '_edges.html')
        x,y=get_edge_counts_for_cutoffs(dataset=dataset)
        p = figure(title='edge counts for ' + dataset, tools="pan,wheel_zoom,box_zoom,previewsave,reset")
        p.left[0].formatter.use_scientific=False
        p.xaxis.axis_label='cut-off'
        p.yaxis.axis_label='# edges'
        p.line(x,y, color="green")
        save(p)
        plots.append(p)
    output_file(TEMPLATE_FOLDER + 'bokeh_edgecount.html')
    gp = gridplot(plots, ncols=2, tools="pan,wheel_zoom,box_zoom,previewsave,reset")
    save(gp)

#simple line plot of (x,y) using bokeh
def plot_nodecounts_bokeh_line():
    datasets_csv = read_datasets_csv()
    datasets = sorted([dataset['name'] for dataset in datasets_csv])
    plots = []
    for dataset in datasets:
        output_file(PNG_FOLDER + dataset + '_nodes.html')
        x,y=get_node_counts_for_cutoffs(dataset=dataset)
        p = figure(title='node counts for ' + dataset, tools="pan,wheel_zoom,box_zoom,previewsave,reset")
        p.xaxis.axis_label='cut-off'
        p.yaxis.axis_label='# nodes'
        p.line(x,y, color="green")
        save(p)
        plots.append(p)
    output_file(TEMPLATE_FOLDER + 'bokeh_nodecount.html')
    gp = gridplot(plots, ncols=2, tools="pan,wheel_zoom,box_zoom,save,reset")
    save(gp)

if __name__ == '__main__':
    plot_edgecounts_matplotlibe_main()
    plot_nodecounts_matplotlibe_main()
    #plot_edgecounts_bokeh_line()
    #plot_nodecounts_bokeh_line()