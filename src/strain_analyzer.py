from Bio import AlignIO, SeqIO
import graph_analysis_util as util
import sys, os
import operator
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def compare_string(string1, string2):
    same = sum(map(operator.eq, string1,string2))
    match_ratio = same/len(string1)
    return match_ratio


# remove >(ratio)% identical sequences
def remove_similar_sequences(sequences, ratio, folder):
    unique_sequences = []
    unique_sequences.append(sequences[0])

    for sequence in sequences[1:]:
        similar = False
        for unique_sequence in unique_sequences:
            similarity = compare_string(sequence[1], unique_sequence[1])
            if similarity>=float(ratio):
                similar = True
                break
        if similar is False:
            unique_sequences.append(sequence)

    #convert 0.999 to 999 for appending to filename
    file_append_ratio = ratio.split('.')[1]
    strains_ratio = [unique_sequence[0] for unique_sequence in unique_sequences]
    print('number of ', ratio, ' strains ', len(strains_ratio))
    util.write_strains_to_csv(strains_ratio, folder + os.sep + 'unique_strains_' + file_append_ratio)


def total_nodes_edges(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10):

    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')
    sequences3 = AlignIO.read(file3, 'fasta')
    sequences4 = AlignIO.read(file4, 'fasta')
    sequences5 = AlignIO.read(file5, 'fasta')
    sequences6 = AlignIO.read(file6, 'fasta')
    sequences7 = AlignIO.read(file7, 'fasta')
    sequences8 = AlignIO.read(file8, 'fasta')
    sequences9 = AlignIO.read(file9, 'fasta')
    sequences10 = AlignIO.read(file10, 'fasta')

    lengths = []

    lengths.append(len(sequences1[0]))
    lengths.append(len(sequences2[0]))
    lengths.append(len(sequences3[0]))
    lengths.append(len(sequences4[0]))
    lengths.append(len(sequences5[0]))
    lengths.append(len(sequences6[0]))
    lengths.append(len(sequences7[0]))
    lengths.append(len(sequences8[0]))
    lengths.append(len(sequences9[0]))
    lengths.append(len(sequences10[0]))

    total_nodes = 0
    total_edges = 0
    for i in range(0, 10):
        print(lengths[i])
        total_nodes = total_nodes + lengths[i]
        for j in range(i+1, 10):
            total_edges = total_edges + (lengths[i] * lengths[j])

    print('total nodes ', total_nodes)
    print('total edges ', total_edges)


def concat_sequences(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10):

    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')
    sequences3 = AlignIO.read(file3, 'fasta')
    sequences4 = AlignIO.read(file4, 'fasta')
    sequences5 = AlignIO.read(file5, 'fasta')
    sequences6 = AlignIO.read(file6, 'fasta')
    sequences7 = AlignIO.read(file7, 'fasta')
    sequences8 = AlignIO.read(file8, 'fasta')
    sequences9 = AlignIO.read(file9, 'fasta')
    sequences10 = AlignIO.read(file10, 'fasta')

    complete_sequences = []

    for sequence1 in sequences1:
        strain_name = util.get_strain_name(sequence1)
        sequence2 = util.get_matching_sequence(sequences2, strain_name=strain_name)
        sequence3 = util.get_matching_sequence(sequences3, strain_name=strain_name)
        sequence4 = util.get_matching_sequence(sequences4, strain_name=strain_name)
        sequence5 = util.get_matching_sequence(sequences5, strain_name=strain_name)
        sequence6 = util.get_matching_sequence(sequences6, strain_name=strain_name)
        sequence7 = util.get_matching_sequence(sequences7, strain_name=strain_name)
        sequence8 = util.get_matching_sequence(sequences8, strain_name=strain_name)
        sequence9 = util.get_matching_sequence(sequences9, strain_name=strain_name)
        sequence10 = util.get_matching_sequence(sequences10, strain_name=strain_name)

        if (sequence2 and sequence3 and sequence4 and sequence5 and sequence6 and sequence7 and sequence8 and sequence9 and sequence10):
            complete_sequence=[]
            complete_sequence.append(util.get_strain_name(sequence1))
            complete_sequence.append(sequence1.seq+sequence2.seq+sequence3.seq+sequence4.seq+sequence5.seq+sequence6.seq+sequence7.seq+sequence8.seq+sequence9.seq+sequence10.seq)
            complete_sequences.append(complete_sequence)

    return complete_sequences

def combine_2_fastas(file1, file2):

    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    records = []
    seq_id = 1

    for sequence1 in sequences1:
        strain_name = util.get_strain_name(sequence1)
        sequence2 = util.get_matching_sequence(sequences2, strain_name=strain_name)

        if sequence2 is not None:
            records.append(SeqRecord(Seq(str(sequence1.seq+sequence2.seq)), id=str(seq_id), name='HA_NA', description='HA_NA'))
            seq_id += 1

    print(len(records))
    with open("ha_na.fasta", "w") as handle:
        SeqIO.write(records, handle, "fasta")

    handle.close()

def run(folder, ratios_list):

    file1 = folder + os.sep + 'ha.afasta'
    file2 = folder + os.sep + 'na.afasta'
    file3 = folder + os.sep + 'm1.afasta'
    file4 = folder + os.sep + 'm2.afasta'
    file5 = folder + os.sep + 'np.afasta'
    file6 = folder + os.sep + 'pb1.afasta'
    file7 = folder + os.sep + 'pb2.afasta'
    file8 = folder + os.sep + 'pa.afasta'
    file9 = folder + os.sep + 'ns1.afasta'
    file10 = folder + os.sep + 'ns2.afasta'
    complete_sequences = concat_sequences(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10)
    print('number of complete strains ', len(complete_sequences))

    processes = []
    for ratio in ratios_list:
        p = mp.Process(target=remove_similar_sequences, args=(complete_sequences, ratio, folder))
        processes.append(p)

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    print('after join')

def number_of_max_nodes_edges(folder):
    file1 = folder + os.sep + 'ha.afasta'
    file2 = folder + os.sep + 'na.afasta'
    file3 = folder + os.sep + 'm1.afasta'
    file4 = folder + os.sep + 'm2.afasta'
    file5 = folder + os.sep + 'np.afasta'
    file6 = folder + os.sep + 'pb1.afasta'
    file7 = folder + os.sep + 'pb2.afasta'
    file8 = folder + os.sep + 'pa.afasta'
    file9 = folder + os.sep + 'ns1.afasta'
    file10 = folder + os.sep + 'ns2.afasta'
    total_nodes_edges(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10)

#one time use; move it to a util file
def clean_up_swine_h1n1(filename):
    sequences = list(SeqIO.parse(filename, "fasta"))
    nfilename = 'nha.fasta'
    new_sequences = []
    for sequence in sequences:
        if len(sequence.seq) > 550:
            new_sequences.append(sequence)
    SeqIO.write(new_sequences, nfilename, 'fasta')

if __name__ == '__main__':
    folder = sys.argv[1]
    number_of_max_nodes_edges(folder)
    #file1 = 'ha.fasta'
    #clean_up_swine_h1n1(file1)
    #file2 = folder + os.sep + 'na.afasta'
    #combine_2_fastas(file1=file1, file2=file2)
    #ratios = sys.argv[2]
    #ratios_list= ratios.split(',')
    #run(folder, ratios_list)
