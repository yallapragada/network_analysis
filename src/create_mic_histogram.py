from plot_util import plot_mic_histogram
import sys
import os


def draw_mic_histogram(folder, infile, dataset):
    mics_all = []
    plot_file = 'hist_' + infile.split('.')[0]
    infile = folder + os.sep + infile
    with open(infile, 'r') as inputfile:
        for line in inputfile:
            if '<data key="d1">' in line:
                mic=line[line.find('>') + 1: line.find('</')]
                mic_float=float(mic)
                mics_all.append(mic_float)
    plot_file = folder + os.sep + plot_file + '.png'
    plot_mic_histogram(mics=mics_all, plot_file=plot_file, dataset=dataset, bins=40)


def run(folder, infile, dataset):
    draw_mic_histogram(folder=folder, infile=infile, dataset=dataset)


if __name__ == '__main__':
    folder     = sys.argv[1]
    infile     = sys.argv[2]
    dataset    = sys.argv[3]
    run(folder=folder, infile=infile, dataset=dataset)