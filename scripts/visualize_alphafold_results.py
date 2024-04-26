# taken from: https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py
import glob
import math
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle


def get_pae_plddt(model_names):
    out = {}
    for i, name in enumerate(model_names):
        d = pickle.load(open(name, 'rb'))
        basename = os.path.basename(name)
        basename = basename[basename.index('model'):]
        out[f'{basename}'] = {'plddt': d['plddt'], 'pae': d['predicted_aligned_error']}
    return out


def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    plt.figure(figsize=(14, 4), dpi=100)
    ##################################################################
    plt.subplot(1, 2, 1)
    plt.title("Sequence coverage")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    ##################################################################
    plt.subplot(1, 2, 2)
    plt.title("Predicted LDDT per position")
    for model_name, value in pae_plddt_per_model.items():
        plt.plot(value["plddt"], label=model_name)
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}/{name + ('_' if name else '')}coverage_LDDT.png")
    ##################################################################

    ##################################################################
    num_models = 5  # columns
    num_runs_per_model = math.ceil(len(model_names) / num_models)
    fig = plt.figure(figsize=(3 * num_models, 2 * num_runs_per_model), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(num_runs_per_model, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    fig.tight_layout()
    plt.savefig(f"{out_dir}/{name + ('_' if name else '')}PAE.png")
    ##################################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', dest='input_dir', required=True)
    parser.add_argument('--name', dest='name')
    parser.set_defaults(name='')
    parser.add_argument('--output_dir', dest='output_dir')
    parser.set_defaults(output_dir='')
    args = parser.parse_args()

    feature_dict = pickle.load(open(f'{args.input_dir}/DTX1_1-DTX2_1_model_1_multimer_v3_output_dict.pkl', 'rb'))
    # is_multimer = ('result_model_1_multimer.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
    # is_ptm = ('result_model_1_ptm.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
    # model_names = [f'{args.input_dir}/result_model_{f}{"_multimer" if is_multimer else "_ptm" if is_ptm else ""}.pkl' for f in range(1,6)]
    model_names = sorted(glob.glob(f'{args.input_dir}/DTX1_1-DTX2_1_model_1_multimer_v3_output_dict.pkl'))

    pae_plddt_per_model = get_pae_plddt(model_names)
    generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name,
                           pae_plddt_per_model)
