import os
import json
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
import py3Dmol
from openfold.np import protein
from openfold.np.relax.utils import overwrite_b_factors
from .utils import get_config_preset_list_for_model, get_run_folder_by_id

# Color bands for visualizing plddt
PLDDT_BANDS = [
    (0, 50, '#FF7D45'),
    (50, 70, '#FFDB13'),
    (70, 90, '#65CBF3'),
    (90, 100, '#0053D6')
]

plddts = {}
pae_outputs = {}
weighted_ptms = {}


def load_prediction_results(pkl_dir, fasta_id, model_name):
    with open(f'{pkl_dir}/{fasta_id}_{model_name}_output_dict.pkl', 'rb') as f:
        prediction_result = pickle.load(f)
    return prediction_result

def get_metrics_and_visualizations(run_id, weight_set, model_name, fasta_id, relax_prediction=False):
    best_model_name, to_visualize_pdb_path = get_best_prediction_by_plddt(run_id, weight_set, model_name, fasta_id, relax_prediction)
    plot_3d_structure(to_visualize_pdb_path)
    get_plddt_pae(f"{get_run_folder_by_id(run_id)}/output/selected_predictions", best_model_name, fasta_id)


def get_best_prediction_by_plddt(run_id, weight_set, model_name, fasta_id, relax_prediction=False):
    global plddts, pae_outputs, weighted_ptms, best_unrelaxed_prot
    
    run_folder = get_run_folder_by_id(run_id)
    
    pkl_dir = os.path.join(f"{run_folder}/output", 'predictions')
    output_dir = os.path.join(f"{run_folder}/output", 'selected_predictions')
    os.makedirs(output_dir, exist_ok=True)
    
    config_preset_list = get_config_preset_list_for_model(weight_set, model_name)

    for config_preset_item in config_preset_list:
        prediction_result = load_prediction_results(pkl_dir, fasta_id, config_preset_item)
        mean_plddt = prediction_result['plddt'].mean()

        plddts[config_preset_item] = prediction_result['plddt']
        
        if model_name == 'multimer':
            pae_outputs[config_preset_item] = (prediction_result['predicted_aligned_error'], prediction_result['max_predicted_aligned_error'])
            weighted_ptms[config_preset_item] = prediction_result['weighted_ptm_score']

        final_atom_mask = prediction_result['final_atom_mask']

    # Find the best model according to the mean pLDDT.
    if model_name == 'monomer' or model_name == 'monomer_ptm':
        best_model_name = max(plddts.keys(), key=lambda x: plddts[x].mean())
    elif model_name =='multimer':
        best_model_name = max(weighted_ptms.keys(), key=lambda x: weighted_ptms[x])
        
    print(f"Best model is: {best_model_name}, relaxed: {relax_prediction}")
    
    best_model_plddts = plddts[best_model_name].mean()
    
    print(f"Mean PLDDT: {best_model_plddts} ")
    
    # Save the mean pLDDT 
    pred_output_path = os.path.join(output_dir, f'{fasta_id}_mean_plddt.txt')
    with open(pred_output_path, 'w') as f:
        f.write(str(best_model_plddts))
    
    unrelaxed_file_name = f'{pkl_dir}/{fasta_id}_{best_model_name}_unrelaxed.pdb'
    
    if relax_prediction:
        pdb_file_name = f'{pkl_dir}/{fasta_id}_{best_model_name}_relaxed.pdb'
    else:
        pdb_file_name = unrelaxed_file_name
    
    with open(pdb_file_name, 'r') as file:
        best_pdb = file.read()
        
    with open(unrelaxed_file_name, 'r') as file:
        best_unrelaxed_pdb_str = file.read()
        
    best_unrelaxed_prot = protein.from_pdb_string(best_unrelaxed_pdb_str)
    
    pred_output_path = os.path.join(output_dir, f'{fasta_id}_selected_prediction.pdb')
    with open(pred_output_path, 'w') as f:
        f.write(best_pdb)

    banded_b_factors = []
    for plddt in plddts[best_model_name]:
        for idx, (min_val, max_val, _) in enumerate(PLDDT_BANDS):
            if plddt >= min_val and plddt <= max_val:
                banded_b_factors.append(idx)
                break
    banded_b_factors = np.array(banded_b_factors)[:, None] * final_atom_mask
    to_visualize_pdb = overwrite_b_factors(best_pdb, banded_b_factors)
    
    visualize_output_path = os.path.join(output_dir, f'{fasta_id}_selected_prediction_visualize.pdb')
    with open(visualize_output_path, 'w') as f:
        f.write(to_visualize_pdb)

    pae_output_path = os.path.join(output_dir, f'{fasta_id}_predicted_aligned_error.json')
    if pae_outputs:
        rounded_errors = np.round(pae_outputs[best_model_name][0].astype(np.float64), decimals=1)
        indices = np.indices((len(rounded_errors), len(rounded_errors))) + 1
        indices_1 = indices[0].flatten().tolist()
        indices_2 = indices[1].flatten().tolist()
        pae_data = json.dumps([{
            'residue1': indices_1,
            'residue2': indices_2,
            'distance': rounded_errors.flatten().tolist(),
            'max_predicted_aligned_error': pae_outputs[best_model_name][1].item()
        }],
                            indent=None,
                            separators=(',', ':'))
        with open(pae_output_path, 'w') as f:
            f.write(pae_data)

    return best_model_name, visualize_output_path


def get_plddt_pae(output_dir, best_model_name, fasta_id):
    if pae_outputs:
        num_plots = 2
    else:
        num_plots = 1

    plt.figure(figsize=[8 * num_plots, 6])
    plt.subplot(1, num_plots, 1)
    plt.plot(plddts[best_model_name])
    plt.title('Predicted LDDT')
    plt.xlabel('Residue')
    plt.ylabel('pLDDT')

    if num_plots == 2:
        plt.subplot(1, 2, 2)
        pae, max_pae = list(pae_outputs.values())[0]
        plt.imshow(pae, vmin=0., vmax=max_pae, cmap='Greens_r')
        plt.colorbar(fraction=0.046, pad=0.04)
        
        total_num_res = best_unrelaxed_prot.residue_index.shape[-1]
        chain_ids = best_unrelaxed_prot.chain_index
        for chain_boundary in np.nonzero(chain_ids[:-1] - chain_ids[1:]):
            if chain_boundary.size:
                plt.plot([0, total_num_res], [chain_boundary, chain_boundary], color='red')
                plt.plot([chain_boundary, chain_boundary], [0, total_num_res], color='red')
        plt.title('Predicted Aligned Error')
        plt.xlabel('Scored residue')
        plt.ylabel('Aligned residue')

    # Save the pLDDT and predicted aligned error plots as PNG
    plt.savefig(os.path.join(output_dir, f'{fasta_id}_plddt_pae.png'))
    print(f"Saved pLDDT and predicted aligned error plots as PNG to {output_dir}")
    return plt

def plot_plddt_legend():
    thresh = [
        'Very low (pLDDT < 50)',
        'Low (70 > pLDDT > 50)',
        'Confident (90 > pLDDT > 70)',
        'Very high (pLDDT > 90)']

    colors = [x[2] for x in PLDDT_BANDS]

    plt.figure(figsize=(2, 2))
    for c in colors:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False, loc='center', fontsize=20)
    plt.xticks([])
    plt.yticks([])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.title('Model Confidence', fontsize=20, pad=20)
    plt.show()
    return plt

def plot_3d_structure(pdb_file_path):
    """Plots the 3D structure for use in a Jupyter notebook."""
    show_sidechains = True

    with open(pdb_file_path, 'r') as f:
        pdb_content = f.read()
        
    color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}
    view = py3Dmol.view(width=800, height=600)
    view.addModelsAsFrames(pdb_content)
    style = {'cartoon': {
        'colorscheme': {
            'prop': 'b',
            'map': color_map}
            }}
    if show_sidechains:
        style['stick'] = {}
    view.setStyle({'model': -1}, style)
    view.zoomTo()
    view.show()