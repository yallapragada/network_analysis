from py2neo import Graph
from py2neo import Node, Relationship
import pandas as pd
import sys
import graph_analysis_util as util
import glob, os
consensus_sequences = {}


def get_consensus_sequences(fasta_folder):
    os.chdir(fasta_folder)
    for filename in glob.glob('*.afasta'):
        protein = filename.split('.')[0].upper()
        filename = fasta_folder + os.sep + filename
        consensus_sequences[protein] = util.get_consensus_sequence_fasta(filename).upper()


# get handle to local neo4j instance
def get_connection_to_neo4j():
    graph = Graph(password='flu')
    return graph


# create pandas dataframe from csv
def create_pd_df_from_csv(filename):
    csv_df = pd.read_csv(filename)
    return csv_df


# given node type and residue number, return node
def find_node(graph, type, unique_id):
    node = graph.find_one(type, property_key='unique_id', property_value=unique_id)
    return node


#create database unique id
def create_unique_id(dataset, protein, residue_number):
    return dataset + '_' + protein + '_' + residue_number


# create nodes for each unique source & target residues in csv
def create_nodes(csv_df, graph, dataset):

    source_nodes = csv_df.x.unique()
    target_nodes = csv_df.y.unique()

    for node in source_nodes:
        tokens = str(node).split('_')
        protein = tokens[0].upper()
        residue_number = tokens[1]

        unique_id = create_unique_id(dataset=dataset, protein=protein, residue_number=residue_number)
        # don't create a node if one exists already
        if find_node(graph=graph, type=protein, unique_id=unique_id) is None:
            aa = consensus_sequences[protein][int(residue_number)-1]
            aa3 = util.get_3letter_aa(aa)
            print(protein, residue_number, aa3)
            node = Node(protein, dataset, residue_number=residue_number, aa=aa3, protein=protein, dataset=dataset, unique_id=unique_id)
            graph.create(node)
        else:
            print('node exists for %s %s' % (protein, residue_number))

    for node in target_nodes:
        tokens = str(node).split('_')
        protein = tokens[0].upper()
        residue_number = tokens[1]
        unique_id = create_unique_id(dataset=dataset, protein=protein, residue_number=residue_number)

        if find_node(graph=graph, type=protein, unique_id=unique_id) is None:
            aa = consensus_sequences[protein][int(residue_number)-1]
            aa3 = util.get_3letter_aa(aa)
            print(protein, residue_number, aa3)
            node = Node(protein, dataset, residue_number=residue_number, aa=aa3, protein=protein, dataset=dataset, unique_id=unique_id)
            graph.create(node)
        else:
            print('node exists for %s %s' % (protein, residue_number))


def create_edges(csv_df, graph, dataset):

    for index, row in csv_df.iterrows():
        p1 = str(row['p1']).upper()
        p2 = str(row['p2']).upper()

        p1_residue_number = str(row['x']).split('_')[1]
        p2_residue_number = str(row['y']).split('_')[1]

        unique_id = create_unique_id(dataset=dataset, protein=p1, residue_number=p1_residue_number)
        n1 = find_node(graph=graph, type=p1, unique_id=unique_id)

        unique_id = create_unique_id(dataset=dataset, protein=p2, residue_number=p2_residue_number)
        n2 = find_node(graph=graph, type=p2, unique_id=unique_id)

        mic = format(row['weight'], '.3f')

        edge = Relationship(n1, 'MUTATES_WITH', n2, mic=mic)
        graph.create(edge)


def load(filename, dataset):
    csv_df = create_pd_df_from_csv(filename)
    create_nodes(csv_df=csv_df, dataset=dataset, graph=graph)
    create_edges(csv_df=csv_df, graph=graph, dataset=dataset)


if __name__ == '__main__':
    graph     = get_connection_to_neo4j()
    csv_folder = sys.argv[1]
    dataset    = sys.argv[2]
    fasta_folder = sys.argv[3]
    get_consensus_sequences(fasta_folder=fasta_folder)
    os.chdir(csv_folder)
    for filename in glob.glob('*.csv'):
        print(filename)
        filename = csv_folder + os.sep + filename
        load(filename=filename, dataset=dataset)