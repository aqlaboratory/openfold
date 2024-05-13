# taken from: https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb#scrollTo=jQUP8Ab3RN7s
import argparse
import sys
import pickle as pkl
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt, colors as cols, cm as cm, rcParams, font_manager
from mpl_toolkits.axes_grid1 import ImageGrid
import json
from sys import exit
import os
from Bio import SeqIO
import io
import json
from json import encoder

encoder.FLOAT_REPR = lambda o: format(o, '.2f')

# plot size, in inches.
plot_size = 16

# @markdown Input value to increment plot axes by (this may need finetuning based on output)
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


def plot_paE(outdir, name, model1, model2, model3, prot1len, size_in_inches=3.5, axis_label_increment=200):
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

    fig = plt.figure(figsize=(size_in_inches, size_in_inches))

    grid = ImageGrid(fig, 111,  # as in plt.subplot(111)
                     nrows_ncols=(1, 3),
                     axes_pad=0.15,
                     share_all=False,
                     cbar_location="right",
                     cbar_mode="single",
                     cbar_size="7%",
                     cbar_pad=0.15,
                     )

    models = [model1, model2, model3]

    cnt = 1
    for ax, model in zip(grid, models):
        im = draw_subplot(f'model{cnt}', ax, model, prot1len)
        cnt += 1

    scale = ax.cax.colorbar(im, label="Predicted error (Å)")
    scale.set_label(label="Predicted error (Å)", size=14, fontweight="bold")
    # Save plot
    plt.savefig('{}/{}_PAE.png'.format(outdir, name), dpi=300)


def generate_pae_plot(fasta, pkl1, pkl2, pkl3, outdir, name):
    model1_results = AlphaFoldPickle(name, pkl1)
    model1_results.saving_pathname = outdir
    model1_results.saving_filename = name

    model2_results = AlphaFoldPickle(name, pkl2)
    model2_results.saving_pathname = outdir
    model2_results.saving_filename = name

    model3_results = AlphaFoldPickle(name, pkl3)
    model3_results.saving_pathname = outdir
    model3_results.saving_filename = name

    def get_multimer_prot1_len(f):
        with open(f) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return len(record.seq)

    prot1len = get_multimer_prot1_len(fasta)
    # results.write_pLDDT_file()
    print("Plotting pLDDT for {}".format(name))
    plot_paE(outdir, name, model1_results, model2_results, model3_results, prot1len, size_in_inches=plot_size,
             axis_label_increment=plot_increment)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', dest='fasta', required=True)
    parser.add_argument('--model1_pkl', dest='model1_pkl', required=True)
    parser.add_argument('--model2_pkl', dest='model2_pkl', required=True)
    parser.add_argument('--model3_pkl', dest='model3_pkl', required=True)
    parser.add_argument('--output_dir', dest='output_dir', required=True)
    parser.add_argument('--basename', dest='basename', required=True)
    args = parser.parse_args()

    generate_pae_plot(args.fasta, args.model1_pkl, args.model2_pkl, args.model3_pkl, args.output_dir, args.basename)
