import graph_analysis_util as util
import os, sys
from Bio import AlignIO
from dssp_parse import DSSPData
from prody.proteins import dssp, parsePDB

def get_residues_from_graph(graphml_file):
    g = util.read_graphml(graphml_file)
    in_residues = {}
    proteins = ['ha', 'na', 'm1', 'm2', 'np', 'pb1', 'pb2', 'pa', 'ns1', 'ns2']

    for protein in proteins:
        in_residues.setdefault(protein,[])

    for node in g.nodes_iter():
        residue_no = str(node).split('_')[1]
        protein = str(node).split('_')[0]
        in_residues[protein].append(residue_no)

    return in_residues


def in_out_entropies(entropies, in_residues_one_p):
    in_network_entropy=0.0
    out_network_entropy=0.0

    for id, entropy in enumerate(entropies):
        if str(id+1) in in_residues_one_p:
            in_network_entropy+=entropy
        else:
            out_network_entropy+=entropy

    len_out_network = len(entropies) - len(in_residues_one_p)
    in_network_entropy_avg = in_network_entropy/len(in_residues_one_p)
    out_network_entropy_avg = out_network_entropy/len_out_network

    return in_network_entropy_avg, out_network_entropy_avg


def run_in_out_entropy_comparision():
    folder       = sys.argv[1]
    graphml_file = sys.argv[2]
    in_residues = get_residues_from_graph(graphml_file)

    proteins = ['ha', 'na', 'm1', 'm2', 'np', 'pb1', 'pb2', 'pa', 'ns1', 'ns2']
    for protein in proteins:
        file = folder + os.sep + 'u_' + protein + '.afasta'
        alignment = AlignIO.read(file, 'fasta')
        sequences = [x.seq for x in alignment]
        entropies = util.entropy_all_positions(sequences)
        in_network_entropy_avg, out_network_entropy_avg = in_out_entropies(entropies, in_residues[protein])
        print('averages for ', protein, ' ', in_network_entropy_avg, out_network_entropy_avg)


def run_in_out_acc_comparision():
    graphml_file = sys.argv[1]
    in_residues = get_residues_from_graph(graphml_file)
    print(in_residues['m1'])
    dssp = DSSPData()
    run_in_out_acc_ha(in_residues, dssp)
    run_in_out_acc_m1(in_residues, dssp)
    run_in_out_acc_na(in_residues, dssp)
    run_in_out_acc_np(in_residues, dssp)
    run_in_out_acc_ns1(in_residues, dssp)


# ha pdb/uniprot numbering delta={5,18 - ChainA} {501,344 - ChainB}
def run_in_out_acc_ha(in_residues, dssp):
    dssp.parseDSSP('HA_1RU7.dssp')
    residue_acc_a = dssp.getAllResAccForChain('A', delta=13)
    residue_acc_b = dssp.getAllResAccForChain('B', delta=-157)

    residue_acc = {}
    residue_acc.update(residue_acc_a)
    residue_acc.update(residue_acc_b)

    sav_in_network=[]
    sav_out_network=[]

    for key, value in residue_acc.items():
        if key in in_residues['ha']:
            sav_in_network.append(int(value))
        else:
            sav_out_network.append(int(value))

    print('average in-network ACC for HA ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for HA ', sum(sav_out_network)/len(sav_out_network))

# m1 same numbering for pdb & uniprot
def run_in_out_acc_m1(in_residues, dssp):
    dssp.parseDSSP('M1_1EA3.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]

    for key, value in residue_acc.items():
        if key in in_residues['m1']:
            sav_in_network.append(int(value))
        else:
            sav_out_network.append(int(value))

    print('average in-network ACC for M1 ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for M1 ', sum(sav_out_network)/len(sav_out_network))

# na same numbering for pdb & uniprot
def run_in_out_acc_na(in_residues, dssp):
    dssp.parseDSSP('NA_3BEQ.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]

    for key, value in residue_acc.items():
        if key in in_residues['na']:
            sav_in_network.append(int(value))
        else:
            sav_out_network.append(int(value))

    print('average in-network ACC for NA ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NA ', sum(sav_out_network)/len(sav_out_network))

# ns1 same numbering for pdb & uniprot
def run_in_out_acc_ns1(in_residues, dssp):
    dssp.parseDSSP('NS1_2GX9.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]

    for key, value in residue_acc.items():
        if key in in_residues['ns1']:
            sav_in_network.append(int(value))
        else:
            sav_out_network.append(int(value))

    print('average in-network ACC for NS1 ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NS1 ', sum(sav_out_network)/len(sav_out_network))


# np same numbering for pdb & uniprot
def run_in_out_acc_np(in_residues, dssp):
    dssp.parseDSSP('NP_2IQH.dssp')
    residue_acc = dssp.getAllResAccForChain('A', delta=0)

    sav_in_network=[]
    sav_out_network=[]

    for key, value in residue_acc.items():
        if key in in_residues['np']:
            sav_in_network.append(int(value))
        else:
            sav_out_network.append(int(value))

    print('average in-network ACC for NP ', sum(sav_in_network)/len(sav_in_network))
    print('average out-network ACC for NP ', sum(sav_out_network)/len(sav_out_network))

#run_in_out_entropy_comparision()
run_in_out_acc_comparision()
