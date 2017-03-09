import graph_analysis_util as util
import os, sys
from Bio import AlignIO
from dssp_parse import DSSPData
import plot_util
from scipy import stats

def get_edges_from_graph(g, protein, residue):
    node = protein + '_' + residue
    edges = g.edges(node)
    return edges


def get_residues_from_graph(g):
    in_residues = {}
    proteins = ['ha', 'na', 'm1', 'm2', 'np', 'pb1', 'pb2', 'pa', 'ns1', 'ns2']

    for protein in proteins:
        in_residues.setdefault(protein,[])

    for node in g.nodes_iter():
        residue_no = str(node).split('_')[1]
        protein = str(node).split('_')[0]
        in_residues[protein].append(residue_no)

    return in_residues


def in_out_entropies(entropies, in_residues_one_p, folder, protein):
    tot_entropy=0.0
    in_network_entropy=0.0
    out_network_entropy=0.0
    in_network_y    =[]
    in_network_x    =[]
    colors = []
    length = len(entropies)

    for id, entropy in enumerate(entropies):
        tot_entropy+=entropy
        if str(id+1) in in_residues_one_p:
            in_network_entropy+=entropy
            in_network_y.append(entropy)
            in_network_x.append(str(id+1))
            colors.append('r')
        else:
            out_network_entropy+=entropy
            in_network_y.append(entropy)
            in_network_x.append(str(id + 1))
            colors.append('g')

    len_out_network = len(entropies) - len(in_residues_one_p)
    if in_network_entropy > 0:
        in_network_entropy_avg = in_network_entropy/len(in_residues_one_p)
    else:
        in_network_entropy_avg = 0.0
    out_network_entropy_avg = out_network_entropy/len_out_network
    entropy_avg = tot_entropy/length

    #plot_util.plot_in_out_entropies(in_network_x, in_network_y, folder, protein, colors)
    return entropy_avg, in_network_entropy_avg, out_network_entropy_avg


def run_in_out_entropy_comparision_all_datasets():
    avg_entropies = []
    input_file = 'C:\\Users\\uday\\pycharm_projects\\network_analysis\\data\\entropy_input.txt'
    lines = util.read_file(input_file)
    for line in lines:
        folder = line[0]
        graphml_file = line[1]
        dataset = line[2]
        avg_entropy_dict=run_in_out_entropy_comparision(dataset=dataset, graphml_file=graphml_file, folder=folder)
        avg_entropies.append(avg_entropy_dict)
    util.write_avg_entropies_to_csv(avg_entropies)

#folder = folder that contains aligned fasta files
#graphml = contains in-network residues
#dataset = name of dataset
def run_in_out_entropy_comparision(dataset, graphml_file, folder):
    g = util.read_graphml(graphml_file)
    in_residues = get_residues_from_graph(g)

    avg_entropy_dict = {}
    avg_entropy_dict['dataset'] = dataset
    proteins = ['ha', 'na', 'm1', 'm2', 'np', 'pb1', 'pb2', 'pa', 'ns1', 'ns2']
    for protein in proteins:
        file = folder + os.sep + protein + '.afasta'
        alignment = AlignIO.read(file, 'fasta')
        sequences = [x.seq for x in alignment]
        entropies = util.entropy_all_positions(sequences)
        entropy_avg, in_network_entropy_avg, out_network_entropy_avg = in_out_entropies(entropies, in_residues[protein], folder, protein)
        avg_entropy_dict[protein] = entropy_avg
        print('averages for ', dataset, ' ', protein, ' ', in_network_entropy_avg, out_network_entropy_avg)

    return avg_entropy_dict


def run_in_out_acc_comparision():
    graphml_file = sys.argv[1]
    folder       = sys.argv[2]

    g = util.read_graphml(graphml_file)
    in_residues = get_residues_from_graph(g)
    dssp = DSSPData()
    run_in_out_acc_ha(in_residues, dssp, g, folder)
    run_in_out_acc_m1(in_residues, dssp, g, folder)
    run_in_out_acc_na(in_residues, dssp, g, folder)
    run_in_out_acc_np(in_residues, dssp, g, folder)
    run_in_out_acc_ns1(in_residues, dssp, g, folder)


# ha pdb/uniprot numbering delta={5,18 - ChainA} {501,344 - ChainB}
def run_in_out_acc_ha(in_residues, dssp, g, folder):
    dssp.parseDSSP('HA_1RU7.dssp')
    residue_acc_a = dssp.getAllResAccForChain('A', delta=13)
    residue_acc_b = dssp.getAllResAccForChain('B', delta=-157)

    residue_acc = {}
    residue_acc.update(residue_acc_a)
    residue_acc.update(residue_acc_b)

    sav_in_network=[]
    sav_out_network=[]
    num_edges = []

    for key, value in residue_acc.items():
        if key in in_residues['ha']:
            sav_in_network.append(int(value))
            num_edges.append(len(get_edges_from_graph(g, 'ha', key)))
        else:
            sav_out_network.append(int(value))

    plot_util.plot_in_edges_acc(num_edges, sav_in_network, folder, 'ha')
    ttest_result = stats.ttest_ind(sav_in_network, sav_out_network, equal_var=False)
    print('ttest for ha ', ttest_result)
    print('average in-network ACC for HA ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for HA ', sum(sav_out_network)/len(sav_out_network))

# m1 same numbering for pdb & uniprot
def run_in_out_acc_m1(in_residues, dssp, g, folder):
    dssp.parseDSSP('M1_1EA3.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]
    num_edges = []

    for key, value in residue_acc.items():
        if key in in_residues['m1']:
            sav_in_network.append(int(value))
            num_edges.append(len(get_edges_from_graph(g, 'm1', key)))
        else:
            sav_out_network.append(int(value))

    plot_util.plot_in_edges_acc(num_edges, sav_in_network, folder, 'm1')
    ttest_result = stats.ttest_ind(sav_in_network, sav_out_network, equal_var=False)
    print('ttest for m1 ', ttest_result)

    print('average in-network ACC for M1 ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for M1 ', sum(sav_out_network)/len(sav_out_network))

# na same numbering for pdb & uniprot
def run_in_out_acc_na(in_residues, dssp, g, folder):
    dssp.parseDSSP('NA_3BEQ.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network = []
    sav_out_network = []
    num_edges = []

    for key, value in residue_acc.items():
        if key in in_residues['na']:
            sav_in_network.append(int(value))
            num_edges.append(len(get_edges_from_graph(g, 'na', key)))
        else:
            sav_out_network.append(int(value))

    plot_util.plot_in_edges_acc(num_edges, sav_in_network, folder, 'na')
    ttest_result = stats.ttest_ind(sav_in_network, sav_out_network, equal_var=False)
    print('ttest for na ', ttest_result)

    print('average in-network ACC for NA ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NA ', sum(sav_out_network)/len(sav_out_network))


# ns1 same numbering for pdb & uniprot
def run_in_out_acc_ns1(in_residues, dssp, g, folder):
    dssp.parseDSSP('NS1_2GX9.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]
    num_edges = []

    for key, value in residue_acc.items():
        if key in in_residues['ns1']:
            sav_in_network.append(int(value))
            num_edges.append(len(get_edges_from_graph(g, 'ns1', key)))
        else:
            sav_out_network.append(int(value))

    plot_util.plot_in_edges_acc(num_edges, sav_in_network, folder, 'ns1')
    ttest_result = stats.ttest_ind(sav_in_network, sav_out_network, equal_var=False)
    print('ttest for ns1 ', ttest_result)

    print('average in-network ACC for NS1 ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NS1 ', sum(sav_out_network)/len(sav_out_network))


# np same numbering for pdb & uniprot
def run_in_out_acc_np(in_residues, dssp, g, folder):
    dssp.parseDSSP('NP_2IQH.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]

    num_edges = []

    for key, value in residue_acc.items():
        if key in in_residues['np']:
            sav_in_network.append(int(value))
            num_edges.append(len(get_edges_from_graph(g, 'np', key)))
        else:
            sav_out_network.append(int(value))

    plot_util.plot_in_edges_acc(num_edges, sav_in_network, folder, 'np')

    ttest_result = stats.ttest_ind(sav_in_network, sav_out_network, equal_var=False)
    print('ttest for np ', ttest_result)

    print('average in-network ACC for NP ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NP ', sum(sav_out_network)/len(sav_out_network))

run_in_out_entropy_comparision_all_datasets()
#run_in_out_acc_comparision()
