# network_analysis
network analysis of correlated mutations in proteins

# code status
still in a "development phase"

# functionality
performs network analysis (based on "maximal information coefficient") of
correlated mutations in proteins.

# details
* download strains from fludb.org
* split the strain file to individual (protein) fasta files
* align using MUSCLE/MAFFT
* run deduplication script - duplicate_remover.py
* run graph_analysis_2p.py or graph_analysis_1p.py for computing MICs and creating CSV files with node and edge data for a given dataset
* run load_graph_db_from_csv.py to load neo4j from CSV files
* run create_networkx_graph_from_csv.py to create a graphml from CSV files
* import graphml into Gephi to create visualizations and perform analysis
* run create_protein_graph.py to create macro/protein level graphs and visualizations from graphml
* run create_degree_clustering_plots.py to create degree/clustering plots from graphml
* run create_plots_from_neo4j.py to create node and edge count plots from data loaded in neo4j
* run cooccurence_counts.py to perform residue coocurrence analysis
* run create_mic_histogram.py to create a histogram of MICs