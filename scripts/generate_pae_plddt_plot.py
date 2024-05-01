# taken from: https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb#scrollTo=jQUP8Ab3RN7s
import argparse
import sys
import pickle as pkl
#from zipfile import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, colors as cols, cm as cm, rcParams, font_manager
import json
from sys import exit
import os
from Bio import PDB as pdb
import io



# plot size, in inches.
plot_size = 16

# @markdown Input value to increment plot axes by (this may need finetuning based on output)
plot_increment = "50"  # @param[10,25,50,100,250,500]
plot_increment = int(plot_increment)


# Define class for AlphaFold metadata file and class methods
class AlphaFoldMetaData(object):
    def __init__(self, name, PathToFile, FastaSequence=None, ranking=None):
        # Define attributes
        self.name = name
        self.PathToFile = PathToFile
        self.FastaSequence = FastaSequence
        self.saving_filename = self.PathToFile.split("/")[-1].split(".")[0]
        self.saving_pathname = self.PathToFile.split(self.saving_filename)[0]
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)

    # Generate a plot of pLDDT value
    def plot_pLDDT(self, size_in_inches=3.5, axis_label_increment=100):
        x = list(range(0, len(self.pLDDT), 1))
        y = list(self.pLDDT)

        # Use standard AlphaFold colors
        cmap = cols.LinearSegmentedColormap.from_list("", ["red", "orange", "yellow", "cornflowerblue", "blue"])

        plt.figure(figsize=(size_in_inches, (size_in_inches / 2)))
        ticks = np.arange(0, len(self.pLDDT), axis_label_increment)
        plt.xticks(ticks)
        plt.yticks()
        plt.title(self.name, size=20, fontweight="bold")
        plt.xlabel("Residue index", size=16, fontweight="bold")
        plt.ylabel("Predicted LDDT", size=16, fontweight="bold")
        plt.scatter(x, y, c=y, cmap=cmap, s=5)
        plt.clim(0, 100)
        scale = plt.colorbar(shrink=0.5)
        scale.set_label(label="Predicted LDDT", size=12, fontweight="bold")
        # Save to directory with pickle file in
        plt.savefig('{}/{}_pLDDT.png'.format(self.saving_pathname, self.saving_filename), dpi=300)

        # Generate a plot from PAE measurements

    def plot_PAE(self, size_in_inches=3.5, axis_label_increment=100):
        ticks = np.arange(0, self.PAE[1].size, axis_label_increment)
        plt.figure(figsize=(size_in_inches, size_in_inches))
        PAE = plt.imshow(self.PAE, cmap="bwr")
        plt.xticks(ticks)
        plt.yticks(ticks)
        plt.title(self.name, size=20, fontweight="bold")
        plt.xlabel("Residue index", size=16, fontweight="bold")
        plt.ylabel("Residue index", size=16, fontweight="bold")
        scale = plt.colorbar(PAE, shrink=0.5)
        scale.set_label(label="Predicted error (Å)", size=14, fontweight="bold")

        # Save plot
        plt.savefig('{}/{}_PAE.png'.format(self.saving_pathname, self.saving_filename), dpi=300)

        # Generate dataframe from PAE data and save to csv
        pd_PAE = pd.DataFrame(self.PAE)
        pd_PAE.to_csv('{}/{}_PAE.csv'.format(self.saving_pathname, self.saving_filename))
        pd_PAE.to_json('{}/{}_PAE.json'.format(self.saving_pathname, self.saving_filename))

class AlphaFoldPickle(AlphaFoldMetaData):

    def __init__(self, PathToFile, FastaSequence=None, ranking=None):
        super().__init__(PathToFile, FastaSequence, ranking)  # Define attributes
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
            self.PAE = self.data[0]['predicted_aligned_error']
        except:
            print("PAE model data not present. To access this performance metric, run AlphaFold"
                  "using pTM-enabled models.")

        # Define pLDDT
        self.pLDDT = self.data[0]['plddt']

    # Generate a ChimeraX attribute file from pLDDT measurements
    def write_pLDDT_file(self):
        seqMismatch = False
        pd_lDDT = pd.DataFrame(self.pLDDT)
        # Name dataframe column
        pd_lDDT.columns = ["pLDDT"]

        # If the fasta file was provided:
        if self.FastaSequence != None:

            # Open the fasta file in read mode
            with (open("{}".format(self.FastaSequence), "r")) as openfile:
                fasta = openfile.read()

            # Delete header line and remove line-breaks
            sequence = fasta.split("\n", 1)[1].replace("\n", "")

            # Check that the lengths of the two sequences match
            if len(sequence) != len(pd_lDDT):

                # If not, ignore the fasta file
                print(
                    "Length of sequence in fasta file provided ({}) does not match length of sequence used in AlphaFold prediction ({}). Ignoring fasta file.".format(
                        len(sequence), len(pd_lDDT)))
                seqMismatch = True
            # If they do,
            else:
                # Convert the fasta sequence into a residue list
                list_sequence = []
                for item in sequence:
                    list_sequence.append(item)

                # Convert the list into a pandas series
                pd_sequence = pd.Series(list_sequence)

                # Insert the series into the dataframe at column 1 to act as labels for the data
                pd_lDDT.insert(0, "Residue", pd_sequence)

        # Otherwise, remind user to check that they have used corret input files
        else:
            print("Number of residues for which pLDDT is provided: ", len(pd_lDDT),
                  "If this does not match the length of your sequence, please double check the input file.")

        # Tell python not to elide middle rows of dataframe when printing to std.out
        pd.set_option("display.max_rows", None, "display.max_columns", None)

        # Save dataframe to ./outputfiles with same name as original pickle and .csv extension
        pd_lDDT.to_csv('{}/{}_pLDDT.csv'.format(self.saving_pathname, self.saving_filename))
        # Delete residue ID
        if self.FastaSequence != None and seqMismatch == False:
            lDDT_table = pd_lDDT.drop('Residue', axis=1)
        else:
            lDDT_table = pd_lDDT

        # Initialise list to store Chimera-style residue identifiers (":x", where x = residue number)
        residue_list = []

        # Populate this list
        for residue in range(0, len(lDDT_table)):
            residue_list.append(":{}".format(residue + 1))

        # Save to pandas format
        chimerax_numbering = pd.Series(residue_list)

        # Insert in the first column of the dataframe, to satisfy ChimeraX formatting
        lDDT_table.insert(0, 'Numbering', chimerax_numbering)

        # Tidy indices so the first label is 1 not 0
        pd_lDDT.index += 1

        # Create a file to save the Chimera attribute output into

        with (open('{}/{}_lDDT.txt'.format(self.saving_pathname, self.saving_filename), 'w+')) as openfile:

            # Write file header in correct format
            openfile.write('attribute: pLDDTvalue\nmatch mode: 1-to-1\nrecipient: residues\n')

            # Iterate over rows of dataframe, writing residue ID and lDDT value to file with correct formatting
            for i, row in lDDT_table.iterrows():
                openfile.write("\t{}\t{}\n".format(row['Numbering'], row['pLDDT']))

        return pd_lDDT


class AlphaFoldJson:
    def __init__(self, PathToDirectory):
        self.PathToDirectory = PathToDirectory
        self.RankingDebug = []
        try:
            with open("{}/ranking_debug.json".format(self.PathToDirectory)) as jsonfile:
                self.RankingDebugRaw = json.load(jsonfile)
            for index in enumerate(self.RankingDebugRaw['order']):
                self.RankingDebug.append(index)
        except:
            exit(
                "To use batch processing, please ensure that the ranking_debug.json file and the result_model_n.pkl files are present in the directory issued in the command. Exiting AlphaPickle now...")


class AlphaFoldPDB(AlphaFoldMetaData):
    def loadCleanStructure(self, id, filePath):
        standardResidues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                            "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

        parser = pdb.PDBParser()
        parsedStructure = parser.get_structure(id, filePath)
        for chain in parsedStructure.get_chains():
            removeResidues = list()
            for i, residue in enumerate(chain.get_residues()):
                if residue.resname not in standardResidues:
                    removeResidues.append(residue.id)
                    print(residue.id)
            [chain.detach_child(id) for id in removeResidues]

        return parsedStructure

    def extractPLDDT(self, PDBobject):
        pLDDT = []
        for residue in PDBobject.get_residues():
            i = 0
            for atom in residue.get_atoms():
                while i < 1:
                    pLDDT.append(atom.bfactor)
                    i += 1
        pLDDT_series = pd.Series(pLDDT)
        return pLDDT_series

    def __init__(self, PathToFile, FastaSequence=None, ranking=None):
        super().__init__(PathToFile, FastaSequence, ranking)
        # Define attributes
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)
        self.structure = self.loadCleanStructure("test", PathToFile)
        self.pLDDT = self.extractPLDDT(self.structure)
        self.data = []
        self.PAE = None

    def PDB_write_pLDDT(self):
        residueNumbers = pd.Series(range(1, len(self.pLDDT) + 1))
        if len(residueNumbers) != len(self.pLDDT):
            print("Oops")
        else:
            pd_lDDT = pd.DataFrame(self.pLDDT)
            pd_lDDT.columns = ["pLDDT"]
            pd_lDDT.insert(0, "Residue", residueNumbers)
            pd_lDDT.to_csv('{}/{}_pLDDT.csv'.format(self.saving_pathname, self.saving_filename))


class AlphaFoldPAEJson(AlphaFoldMetaData):
    def extractPAEfromJson(self, PathToFile):

        with open(PathToFile, 'r') as file:
            jsonstring = json.load(file)
            if 'predicted_aligned_error' in jsonstring[0]:
                pae = jsonstring[0]['predicted_aligned_error']
            else:
                residue1 = jsonstring[0]['residue1']
                residue2 = jsonstring[0]['residue2']
                pae = jsonstring[0]['distance']

        if 'predicted_aligned_error' in jsonstring[0]:
            paeArray = np.array(pae)
        else:
            paeArray = np.ones((max(residue1), (max(residue2))))
            for i, j, n in zip(residue1, residue2, pae):
                paeArray[int(i - 1), int(j - 1)] = n

        return paeArray

    def __init__(self, PathToFile, FastaSequence=None, ranking=None):
        super().__init__(PathToFile, FastaSequence, ranking)
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)

        self.PAE = self.extractPAEfromJson(PathToFile)
        self.pLDDT = None


# use_files_from_google_drive = False        #@param {type:"boolean"}
#
# if not use_files_from_google_drive:
#     print("Select PAE files for upload")
#     PAEfiles =  list(google.colab.files.upload().keys())
#     print("Select pLDDT files for upload")
#     pLDDTfiles =  list(google.colab.files.upload().keys())
# else:
#     #print("Select PAE files for upload")
#     path_to_PAE_file_in_drive = ""        #@param {type:"string"}
#     if ":" in path_to_PAE_file_in_drive:
#         path_to_PAE_file_in_drive = path_to_PAE_file_in_drive.split(":")
#     else:
#         path_to_PAE_file_in_drive = [path_to_PAE_file_in_drive]
#     path_to_pLDDT_file_in_drive = ""        #@param {type:"string"}
#     if ":" in path_to_pLDDT_file_in_drive:
#         path_to_pLDDT_file_in_drive = path_to_pLDDT_file_in_drive.split(":")
#     else:
#         path_to_pLDDT_file_in_drive = [path_to_pLDDT_file_in_drive]

def generate_plots(pkl, outdir, name):
    results = AlphaFoldPickle(name,pkl, None)
    results.saving_pathname = outdir
    results.saving_filename = name
    if type(results.PAE) == np.ndarray:
        print("Plotting PAE for {} and saving to csv".format(pkl))
        results.plot_PAE(size_in_inches=plot_size, axis_label_increment=plot_increment)

    results = AlphaFoldPickle(name,pkl, None)
    results.saving_filename = name
    results.saving_pathname = outdir
    results.write_pLDDT_file()
    print("Plotting pLDDT for {} and saving to csv".format(pkl))
    results.plot_pLDDT(size_in_inches=plot_size, axis_label_increment=plot_increment)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pkl', dest='input_pkl', required=True)
    parser.add_argument('--output_dir', dest='output_dir', required=True)
    parser.add_argument('--basename', dest='basename', required=True)
    args = parser.parse_args()

    generate_plots(args.input_pkl, args.output_dir, args.basename)
