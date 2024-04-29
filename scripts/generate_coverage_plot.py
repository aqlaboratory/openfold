# taken from: https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
# and https://github.com/sokrypton/ColabFold/blob/main/colabfold/batch.py
import argparse
import os
from pathlib import Path
import pickle as pkl
import numpy as np
from matplotlib import pyplot as plt


def plot_msa_v2(feature_dict, sort_lines=True, dpi=100):
    seq = feature_dict["msa"][0]
    if "asym_id" in feature_dict:
        Ls = [0]
        k = feature_dict["asym_id"][0]
        for i in feature_dict["asym_id"]:
            if i == k:
                Ls[-1] += 1
            else:
                Ls.append(1)
            k = i
    else:
        Ls = [len(seq)]
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"]

    msa = feature_dict["msa"][:N]
    gap = msa != 21
    qid = msa == seq
    gapid = np.stack([gap[:, Ln[i]:Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:, Ln[i]:Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)
    plt.figure(figsize=(8, 5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(lines,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
               extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    return plt


def generate_coverage(fd_pkl, output_dir, name, dpi=500):
    feature_dict_pkl = []
    with (open("{}".format(fd_pkl), "rb")) as openfile:
        while True:
            try:
                feature_dict_pkl.append(pkl.load(openfile))
            except EOFError:
                break
    feature_dict = feature_dict_pkl[0]
    msa_plot = plot_msa_v2(feature_dict, sort_lines=True, dpi=dpi)
    coverage_png = os.path.join(output_dir,f"{name}_coverage.png")
    msa_plot.savefig(str(coverage_png), bbox_inches='tight')
    msa_plot.close()
    # result_files.append(coverage_png)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pkl', dest='feature_dict_pkl', required=True)
    parser.add_argument('--output_dir', dest='output_dir', required=True)
    parser.add_argument('--basename', dest='basename', required=True)
    args = parser.parse_args()

    # generate_coverage(args.input_pkl, args.output_dir, args.basename)
    print("gen coverage plot")
    generate_coverage(args.feature_dict_pkl, args.output_dir, args.basename)
