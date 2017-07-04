# this script contains several utility plotting functions

from bokeh.io import output_file, show
from bokeh.plotting import figure
import matplotlib.pyplot as plt
import matplotlib, os
import matplotlib.patches as mpatches
import operator


image_magick_colors = ['DarkRed', 'DarkBlue', 'DarkOrange', 'DarkGreen', 'DeepPink', 'indigo', 'LightGreen',
                       'LightGoldenrod', 'magenta', 'olive', 'PaleGreen', 'yellow', 'turquoise', 'tan', 'thistle',
                       'sienna', 'seashell', 'PeachPuff', 'peru', 'MistyRose', 'MintCream', 'PaleGoldenrod', 'wheat',
                       'teal', 'SpringGreen3', 'tomato4', 'yellow4', 'BlanchedAlmond', 'CadetBlue2', 'sienna', 'seashell',
                       'PeachPuff', 'peru', 'MistyRose', 'MintCream', 'PaleGoldenrod', 'wheat',
                       'teal', 'SpringGreen3', 'tomato4', 'yellow4', 'BlanchedAlmond', 'CadetBlue2']

#logo for residue counts
def create_image_magick_script(residue_comparision_count, residue1, residue2, folder, filename):
    logo_file = filename.replace('bat', 'png')
    filename = folder + os.sep + filename
    target = open(filename, 'w')
    total = sum(residue_comparision_count.values())
    target.write('magick -background white -bordercolor white -fill green -gravity center -size 160*' + str(total) + ' ^')
    target.write('\n')
    i=0

    for key, value in sorted(residue_comparision_count.items(), key=operator.itemgetter(1), reverse=True):
        print(key, value)
        i=i+1
        ht_str = 'ht' + str(i)
        if value < 20:
            value = value + 20
        row = '( -size 160x' + str(value) + ' -fill ' + image_magick_colors[i-1] + ' -font Consolas caption:\"' + key + '\" -trim -set option:' + ht_str + ' \"%%h\" -resize 160x%%[' + ht_str + ']! ) ^'
        target.write(row)
        target.write('\n')
    target.write('-size 160x160 -fill DarkViolet -font Consolas caption:"' + residue1 + ' ' + residue2 + '"  -trim -border 0x3 -append ' + logo_file)
    target.close()

def plot_mic_histogram(mics, plot_file, dataset, bins):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")
    plt.xlim([0.1, 1])
    plt.hist(mics, bins=bins,color="#3F5D7D")
    plt.title("Histogram of mic values for " + dataset, fontsize=18)
    plt.xlabel("mic", fontsize=14)
    plt.ylabel("counts", fontsize=14)
    plt.tick_params(axis='both', labelsize=12)
    plt.tight_layout()
    plt.savefig(plot_file, format='png',bbox_inches="tight")
    plt.close()


#plot entropies of residues in protein using bokeh
def plotEntropies(entropies, residues, title):
    plot_simple_bokeh_line(y=entropies, x=list(residues), filename=title+".html", width=400, height=300, xaxis="residues", yaxis="entropies", title=title)


#simple line plot of (x,y) using bokeh
def plot_simple_bokeh_line(x,y, filename, height, width, xaxis, yaxis, title):
    output_file(filename)
    p=figure(plot_width=width, plot_height=height,title=title)
    p.line(x,y, color="green")
    p.xaxis.axis_label=xaxis
    p.yaxis.axis_label=yaxis
    show(p)


def plot_residue_counts(rcc_list, title, folder, filename):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

    for rcc in rcc_list:
        values = [x[1] for x in sorted(rcc.items(), key=operator.itemgetter(1), reverse=True)]
        plt.plot(values)
    plt.title(title)
    plt.tight_layout()
    savepath = get_save_path(folder, filename, 'png')
    plt.savefig(savepath, format='png')
    plt.close()


#matplotlib line plots of in & out entropies of residues in a network
def plot_in_out_entropies(in_network_x, in_network_y, folder, protein, colors):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

    plt.scatter(in_network_x, in_network_y, c=colors, s=5)
    plt.title(protein.upper() + ':entropy for in vs out-of network residues', fontsize=16)
    plt.xlabel('residue#', fontsize=12)
    plt.ylabel('entropy', fontsize=12)

    red_patch = mpatches.Patch(color='red', label='in-network')
    blue_patch = mpatches.Patch(color='green', label='out-of-network')
    plt.legend(handles=[red_patch, blue_patch])
    plt.tight_layout()
    savepath = get_save_path(folder, protein+'_ent', 'png')
    plt.savefig(savepath, format='png')
    plt.close()


def plot_in_edges_acc(num_edges, in_network_acc, folder, protein):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

    plt.scatter(num_edges, in_network_acc, s=5)
    plt.title(protein.upper() + ':acc vs num-edges for in-network residues', fontsize=16)
    plt.xlabel('num of edges', fontsize=12)
    plt.ylabel('acc', fontsize=12)
    plt.tight_layout()
    savepath = get_save_path(folder, protein+'_acc', 'png')
    plt.savefig(savepath, format='png')
    plt.close()


def get_save_path(directory, filename, format):
    filename = "%s.%s" % (filename, format)
    if not directory:
        savepath = filename
    else:
        savepath = os.path.join(directory, filename)
    return savepath
