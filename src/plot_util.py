from bokeh.io import output_file, show
from bokeh.plotting import figure
import matplotlib.pyplot as plt
import matplotlib, os

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


#matplotlib line plots of in & out entropies of residues in a network
def plot_in_out_entropies(in_network_x, in_network_y, folder, protein, colors):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

    #plt.gca().set_color_cycle(['red', 'green'])
    plt.scatter(in_network_x, in_network_y, c=colors, s=5)
    plt.title(protein.upper() + ':entropy for in vs out-of network residues', fontsize=16)
    plt.xlabel('residue#', fontsize=12)
    plt.ylabel('entropy', fontsize=12)
    plt.tight_layout()
    savepath = get_save_path(folder, protein+'_ent', 'png')
    plt.savefig(savepath, format='png')
    plt.close()


def plot_in_edges_acc(num_edges, in_network_acc, folder, protein):
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="Times New Roman")

    #plt.gca().set_color_cycle(['red', 'green'])
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
