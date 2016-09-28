import graph_analysis_util as util
import matplotlib
import matplotlib.pyplot as plt
import sys
import multiprocessing as mp


def run_0123(file1, file2, cutoff, output):
    p1 = file1.split('.')[0].split("\\")[1]
    p2 = file2.split('.')[0].split("\\")[1]

    print(p1, p2)

    p1_sequences_0123, p2_sequences_0123 = util.create_0123_sequences(file1, file2)
    mic_scores = util.perform_mic_2p(p1_sequences_0123, p2_sequences_0123, p1, p2, cutoff=cutoff)

    output.put(mic_scores)
    print('done with run_0123')


def run_01(file1, file2, cutoff, output):
    p1 = file1.split('.')[0].split("\\")[1]
    p2 = file2.split('.')[0].split("\\")[1]

    p1_sequences_01, p2_sequences_01 = util.create_01_sequences(file1, file2)
    mic_scores = util.perform_mic_2p(p1_sequences_01, p2_sequences_01, p1, p2, cutoff=cutoff)

    output.put(mic_scores)
    print('done with run_01 for ', file1, file2, cutoff)


def run(input_file, folder, cutoff, target):
    manager = mp.Manager()
    output = manager.Queue()

    all_mic_scores = []

    lines = util.read_file(input_file)
    out_filename = folder + '_' + cutoff.replace(".", "") + '_' + target

    processes = []
    for line in lines:
        print("start processing ", line[0], line[1])
        file1 = folder+"\\"+line[0]
        file2 = folder+"\\"+line[1]
        if target == '0123':
            p = mp.Process(target=run_0123, args=(file1, file2, cutoff, output))
        else:
            p = mp.Process(target=run_01, args=(file1, file2, cutoff, output))
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
    mic_df = util.create_df_from_dict(all_mic_scores)
    util.create_graph(mic_df, out_filename)


if __name__ == '__main__':
    args_file = sys.argv[1]
    arglines = util.read_file(args_file)
    for argline in arglines:
        print("start processing ", argline[0], argline[1], argline[2], argline[3])
        run(argline[0], argline[1], argline[2], argline[3])
        print("done processing ", argline[0], argline[1], argline[2], argline[3])