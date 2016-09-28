from Bio import AlignIO
import graph_analysis_util as util
import sys, os
import operator


def compare_string(string1, string2):
    same = sum(map(operator.eq, string1,string2))
    match_ratio = same/len(string1)
    return match_ratio


# remove >(ratio)% identical sequences
def remove_similar_sequences(sequences, ratio):
    unique_sequences = []
    unique_sequences.append(sequences[0])

    for sequence in sequences[1:]:
        similar = False
        for unique_sequence in unique_sequences:
            similarity = compare_string(sequence[1], unique_sequence[1])
            if similarity>=ratio:
                similar = True
                break
        if similar is False:
            unique_sequences.append(sequence)
    return unique_sequences


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

def run():

    folder = sys.argv[1]
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
    print('number of complete sequences ', len(complete_sequences))

    '''
    sequences_100 = remove_similar_sequences(complete_sequences, 1.0)
    strains_100 = [unique_sequence[0] for unique_sequence in sequences_100]
    util.write_strains_to_csv(strains_100, folder + os.sep + 'unique_strains_100')
    print('number of sequences_100 ', len(sequences_100))
    '''

    sequences_999 = remove_similar_sequences(complete_sequences, 0.999)
    strains_999 = [unique_sequence[0] for unique_sequence in sequences_999]
    util.write_strains_to_csv(strains_999, folder + os.sep + 'unique_strains_999')
    print('number of sequences_999 ', len(sequences_999))

    sequences_998 = remove_similar_sequences(complete_sequences, 0.998)
    strains_998 = [unique_sequence[0] for unique_sequence in sequences_998]
    util.write_strains_to_csv(strains_998, folder + os.sep + 'unique_strains_998')
    print('number of sequences_998 ', len(sequences_998))


run()
#print(compare_string('VCAGREF', 'VC-GREF'))