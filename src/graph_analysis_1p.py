import sys
import graph_analysis_util as util
import multiprocessing as mp
import os

'''
python script for computing intra-protein MICs for a given protein
this script can compute MICs for multiple proteins in several datasets
usage: python graph_analysis_1p.py [path to args_1p.txt]
this script processes multiple datasets based on config in args_1p.txt
#format for args_1p.txt:
input_1p.txt,uh3n2_swine_all,0.1,C:\\uday\\gmu\\correlations\\results\\10proteins\\mic_csv\\uh3n2_swine_all\\intra
input_1p.txt,[unique_identifier_for_dataset],[mic_threshold],[path to folder that contains uh3n2_human_all aligned fasta files]
where input_1p.txt contains names of proteins for which MIC should be computed
and 0.1 = MIC threshold
'''


def run_01(file, cutoff, output, out_folder):
    p = file.split('.')[0].split(os.sep)[1]
    p_sequences_01 = util.create_01_sequences_gaps_single(file)
    mic_scores = util.perform_mic_1p(p_sequences=p_sequences_01, p=p, cutoff=cutoff, out_folder=out_folder)
    output.put(mic_scores)
    print('done with run_01 for ', file, cutoff)


def run(input_file, folder, cutoff, out_folder):
    manager = mp.Manager()
    output = manager.Queue()

    all_mic_scores = []

    lines = util.read_file(input_file)

    processes = []
    for line in lines:
        print("start processing ", line[0])
        file = folder+os.sep+line[0]
        p = mp.Process(target=run_01, args=(file, cutoff, output, out_folder))
        processes.append(p)

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    print('after join')
    results = [output.get() for p in processes]
    print('no. of results ', len(results))

    for list_of_mic_scores in results:
        if len(list_of_mic_scores) > 0:
            for mic_score in list_of_mic_scores:
                all_mic_scores.append(mic_score)

    print('no. of mic_scores ', len(all_mic_scores))


if __name__ == '__main__':
    args_file = sys.argv[1]
    arglines = util.read_file(args_file)
    for argline in arglines:
        print("start processing ", argline[0], argline[1], argline[2], argline[3])
        run(argline[0], argline[1], argline[2], argline[3])
        print("done processing ", argline[0], argline[1], argline[2], argline[3])