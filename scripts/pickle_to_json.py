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


def generate_json(pkl1, outdir, name, model_nbr):
    model1_results = AlphaFoldPickle(name, pkl1)
    model1_results.saving_pathname = outdir
    # "${NAME}_model_${model}_multimer_v3_relaxed"
    model1_results.saving_filename = f"{name}_model_{model_nbr}_multimer_v3_relaxed"
    print(f"Saving model{model_nbr} in json format")
    model1_results.save_to_json()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_pkl', dest='model_pkl', required=True)
    parser.add_argument('--output_dir', dest='output_dir', required=True)
    parser.add_argument('--basename', dest='basename', required=True)
    parser.add_argument('--model_nbr', dest='model_nbr', required=True)
    args = parser.parse_args()

    generate_json(args.model_pkl, args.output_dir, args.basename, args.model_nbr)
