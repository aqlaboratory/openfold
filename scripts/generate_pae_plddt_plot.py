# taken from: https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb#scrollTo=jQUP8Ab3RN7s
import argparse
import sys
import pickle as pkl
#from zipfile import Path
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt, colors as cols, cm as cm, rcParams, font_manager
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.table import table
from matplotlib.gridspec import GridSpec
import json
from sys import exit
import os
from Bio import PDB as pdb
from Bio import SeqIO
import io
import json
from json import encoder

encoder.FLOAT_REPR = lambda o: format(o, '.2f')

# plot size, in inches.
plot_size = 16

plot_increment = "200"  # @param[10,25,50,100,250,500]
plot_increment = int(plot_increment)


# Define class for AlphaFold metadata file and class methods
class AlphaFoldMetaData(object):
    def __init__(self, name, PathToFile, FastaSequence=None, ranking=None):
        # Define attributes
        self.name = name
        self.PathToFile = PathToFile
        self.FastaSequence = FastaSequence
        self.saving_filename = name
        self.saving_pathname = self.PathToFile.split(self.saving_filename)[0]
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)


class AlphaFoldPickle(AlphaFoldMetaData):

    def __init__(self, name, PathToFile, FastaSequence=None, ranking=None):
        super().__init__(name, PathToFile, FastaSequence, ranking)  # Define attributes
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)
        self.data = []
        self.PAE = None

        # Extract pickled data
        with (open("{}".format(self.PathToFile), "rb")) as openfile:
            while True:
                try:
                    self.data.append(pkl.load(openfile))
                except EOFError:
                    break

        # Try statement accounts for data run using non-pTM models, with no PAE output
        try:
            self.PAE = self.data[0]['predicted_aligned_error'].round(2)
        except:
            print("PAE model data not present. To access this performance metric, run AlphaFold"
                  "using pTM-enabled models.")

        # Define pLDDT
        self.pLDDT = self.data[0]['plddt'].round(2)
        self.max_pae = self.data[0]['max_predicted_aligned_error']
        self.ptm = self.data[0]['ptm_score']
        self.iptm = self.data[0]['iptm_score']

    def save_to_json(self):
        # save pkl to json format as colabfold
        colab_data = {}
        colab_data['plddt'] = list(np.around(np.array(self.pLDDT.tolist()), 2))
        colab_data['pae'] = list(np.around(np.array(self.PAE.tolist()), 2))
        colab_data['max_pae'] = self.max_pae
        colab_data['ptm'] = self.ptm
        colab_data['iptm'] = self.iptm

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        with open('{}/{}.json'.format(self.saving_pathname, self.saving_filename), "w") as outfile:
            outfile.write(json.dumps(colab_data, cls=NumpyEncoder))


def plot_pLDDT(outdir, name, model1, model2, model3, prot1len, size_in_inches=3.5, axis_label_increment=100):
    m1_x = list(range(0, len(model1.pLDDT), 1))
    m1_y = list(model1.pLDDT)
    m2_x = list(range(0, len(model2.pLDDT), 1))
    m2_y = list(model2.pLDDT)
    m3_x = list(range(0, len(model3.pLDDT), 1))
    m3_y = list(model3.pLDDT)

    plt.figure(figsize=(size_in_inches, (size_in_inches / 2)))
    ticks = np.arange(0, len(model1.pLDDT), axis_label_increment)
    plt.xticks(ticks)
    plt.yticks()
    plt.title(name, size=20, fontweight="bold")
    plt.xlabel("Residue index", size=16, fontweight="bold")
    plt.ylabel("Predicted LDDT", size=16, fontweight="bold")
    plt.plot(m1_x, m1_y, '-b', label='model1')
    plt.plot(m2_x, m2_y, '-m', label='model2')
    plt.plot(m3_x, m3_y, '-g', label='model3')

    plt.vlines(x=prot1len, ymin=0, ymax=100, colors='k', linestyles='--')

    plt.legend(loc='lower right')
    plt.savefig('{}/{}_pLDDT.png'.format(outdir, name), dpi=300)


def plot_paE(outdir, name, model1, model2, model3, prot1len, interface_df, size_in_inches=3.5, axis_label_increment=200):

    # data = [
    #     [0.742, 376, 64, 83, 92, 2, 4, 8, 6.0],
    #     [0.742, 348, 69, 86, 92, 2, 3, 6, 6.0],
    #     [0.018, 3, 54, 58, 63, 14, 15, 15, 5.7]
    # ]
    #
    # columns = (
    # 'pdockq', 'ncontacts', 'plddt_min', 'plddt_avg', 'plddt_max', 'pae_min', 'pae_avg', 'pae_max', 'distance_avg')
    #
    # df = pd.DataFrame(
    #     data,
    #     columns=list(columns)
    # )


    def draw_subplot(name, ax, model, prot1len, display_scale=False):
        ticks = np.arange(0, model.PAE[1].size, axis_label_increment)
        img_ax = ax.imshow(model.PAE, cmap="bwr")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_title(name, size=20, fontweight="bold")
        ax.set_xlabel("Residue index", size=16, fontweight="bold")
        ax.set_ylabel("Residue index", size=16, fontweight="bold")
        ax.axvline(x=prot1len, color='k', linewidth=4)
        ax.axhline(y=prot1len, color='k', linewidth=4)
        return img_ax

    nrows = 1
    height_ratios = [1]
    if interface_df is not None:
        nrows = 2
        height_ratios = [1, 2]

    fig = plt.figure(figsize=(12, 10), layout="constrained")
    gs1 = GridSpec(nrows, 4, figure=fig, width_ratios=[1,1,1,0.1], height_ratios=height_ratios)

    models = [model1, model2, model3]

    ax1 = fig.add_subplot(gs1[0, 0])
    im1 = draw_subplot(f'model1', ax1, models[0], prot1len)

    ax2 = fig.add_subplot(gs1[0, 1])
    im2 = draw_subplot(f'model2', ax2, models[1], prot1len)

    ax3 = fig.add_subplot(gs1[0, 2])
    im3 = draw_subplot(f'model3', ax3, models[2], prot1len)

    ax4 = fig.add_subplot(gs1[0, 3])
    mesh = ax4.pcolormesh(models[2].PAE, cmap="bwr")
    scale = fig.colorbar(mesh, ax4, label="Predicted error (Å)")
    scale.set_label(label="Predicted error (Å)", size=14, fontweight="bold")

    if interface_df is not None:
        ax5 = fig.add_subplot(gs1[1, :])
        ax5.axis('off')
        ax5.axis('tight')
        rows = ['model %d' % x for x in (1, 2, 3)]
        tbl = ax5.table(
            cellText=interface_df.values[:,2:],
            rowLabels=rows,
            colLabels=list(interface_df.columns)[2:],
            loc="upper center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(14)
        tbl.auto_set_column_width([0, 1, 2, 3, 4, 5, 6, 7, 8])

    # Save plot
    plt.savefig('{}/{}_PAE.png'.format(outdir, name), dpi=300)


def generate_plots(fasta, pkl1, pkl2, pkl3, outdir, name, interface):
    model1_results = AlphaFoldPickle(name, pkl1)
    model2_results = AlphaFoldPickle(name, pkl2)
    model3_results = AlphaFoldPickle(name, pkl3)

    def get_multimer_prot1_len(f):
        with open(f) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return len(record.seq)

    prot1len = get_multimer_prot1_len(fasta)

    print("Generating plddt plot")
    plot_pLDDT(outdir, name, model1_results, model2_results, model3_results, prot1len, size_in_inches=plot_size,
               axis_label_increment=plot_increment)

    print("Generating PAE plot")
    df = None
    if interface is None:
        print("No interface file provided, will not output interface table")
    elif os.path.exists(interface):
        df = pd.read_csv(interface, sep=",")
    else:
        print(f"Unable to create pandas dataframe with provided interface file {interface}")

    plot_paE(outdir, name, model1_results, model2_results, model3_results, prot1len, df, size_in_inches=plot_size,
             axis_label_increment=plot_increment)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', dest='fasta', required=True)
    parser.add_argument('--model1_pkl', dest='model1_pkl', required=True)
    parser.add_argument('--model2_pkl', dest='model2_pkl', required=True)
    parser.add_argument('--model3_pkl', dest='model3_pkl', required=True)
    parser.add_argument('--output_dir', dest='output_dir', required=True)
    parser.add_argument('--basename', dest='basename', required=True)
    parser.add_argument('--interface', dest='interface', required=False)
    args = parser.parse_args()

    generate_plots(args.fasta, args.model1_pkl, args.model2_pkl, args.model3_pkl, args.output_dir,
                            args.basename, args.interface)
    # generate_plddt_plot(args.fasta, args.model1_pkl, args.model2_pkl, args.model3_pkl, args.output_dir, args.basename)
    # generate_pae_plot(args.fasta, args.model1_pkl, args.model2_pkl, args.model3_pkl, args.output_dir, args.basename)
