import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from openfold.data import parsers
from .utils import get_run_folder_by_id, read_fasta_file


def get_msa_plot(run_id, fasta_id=None):
    
    run_folder = get_run_folder_by_id(run_id)
    input_msa_folder = os.path.join(f"{run_folder}/output", 'alignments')
    output_folder = os.path.join(f"{run_folder}/output", 'msa_plots')
    os.makedirs(output_folder, exist_ok=True)
    
    
    if fasta_id:
        fasta_file_path = os.path.join(f"{run_folder}/fasta_dir/tmp", f"{fasta_id}.fasta")
        sequence = read_fasta_file(fasta_file_path)
        output_file_plot = os.path.join(output_folder, f"{fasta_id}_msa_plot.png")
        create_msa_plot(fasta_id, sequence, output_file_plot, input_msa_folder)
        
    else:
        # Iterate over each subfolder in output_msa_openfold
        for subfolder in os.listdir(input_msa_folder):
            
            fasta_id = subfolder
            fasta_file_path = os.path.join(f"{run_folder}/fasta_dir/tmp", f"{fasta_id}.fasta")
            
            if os.path.exists(fasta_file_path):
                sequence = read_fasta_file(fasta_file_path)
                output_file_plot = os.path.join(output_folder, f"{fasta_id}_msa_plot.png")
                create_msa_plot(fasta_id, sequence, output_file_plot, input_msa_folder)

def create_msa_plot(fasta_id, original_sequence, output_file, input_msa_folder):
    # Path to the search results files
    search_results_files = {
        original_sequence: {
            'uniref90': os.path.join(input_msa_folder, f'{fasta_id}/uniref90_hits.sto'),
            'mgnify': os.path.join(input_msa_folder, f'{fasta_id}/mgnify_hits.sto'),
            'smallbfd': os.path.join(input_msa_folder, f'{fasta_id}/bfd_uniclust_hits.a3m'),
            # Add other databases if needed
        },
    }

    MAX_HITS_BY_DB = {
        'uniref90': 10000,
        'mgnify': 501,
        'smallbfd': 5000,
    }

    msas_by_seq_by_db = {seq: {} for seq in search_results_files.keys()}
    full_msa_by_seq = {seq: [] for seq in search_results_files.keys()}

    # Function to parse the MSA files
    def parse_msa_file(file_path, file_type):
        if file_type == 'sto':
            with open(file_path, 'r') as f:
                sto_content = f.read()
            msa_obj = parsers.parse_stockholm(sto_content)
        elif file_type == 'a3m':
            with open(file_path, 'r') as f:
                a3m_content = f.read()
            msa_obj = parsers.parse_a3m(a3m_content)
        return msa_obj

    # Load the search results from files
    for seq_name, db_files in search_results_files.items():
        print(f'Loading results for sequence: {fasta_id}:{seq_name}')
        for db_name, result_file in db_files.items():
            print(f'  Loading database: {db_name}')
            file_type = result_file.split('.')[-1]
            msa_obj = parse_msa_file(result_file, file_type)

            msas, del_matrix, targets = msa_obj.sequences, msa_obj.deletion_matrix, msa_obj.descriptions
            db_msas = parsers.Msa(msas, del_matrix, targets)
            if db_msas:
                if db_name in MAX_HITS_BY_DB:
                    db_msas.truncate(MAX_HITS_BY_DB[db_name])
                msas_by_seq_by_db[seq_name][db_name] = db_msas
                full_msa_by_seq[seq_name].extend(msas)
                msa_size = len(set(msas))
                print(f'{msa_size} Sequences Found in {db_name}')

    # Deduplicate full MSA and calculate total MSA size
    for seq_name in full_msa_by_seq.keys():
        full_msa_by_seq[seq_name] = list(dict.fromkeys(full_msa_by_seq[seq_name]))
        total_msa_size = len(full_msa_by_seq[seq_name])
        print(f'\n{total_msa_size} Sequences Found in Total for {seq_name}\n')

    # Visualize the results
    fig = plt.figure(figsize=(12, 3))
    max_num_alignments = 0

    for seq_idx, seq_name in enumerate(search_results_files.keys()):
        full_msas = full_msa_by_seq[seq_name]
        deduped_full_msa = list(dict.fromkeys(full_msas))
        total_msa_size = len(deduped_full_msa)

        aa_map = {restype: i for i, restype in enumerate('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')}
        msa_arr = np.array([[aa_map[aa] for aa in seq] for seq in deduped_full_msa])
        num_alignments, num_res = msa_arr.shape
        plt.plot(np.sum(msa_arr != aa_map['-'], axis=0), label=f'Chain {seq_idx}')
        max_num_alignments = max(num_alignments, max_num_alignments)

    plt.title('Per-Residue Count of Non-Gap Amino Acids in the MSA')
    plt.ylabel('Non-Gap Count')
    plt.yticks(range(0, max_num_alignments + 1, max(1, int(max_num_alignments / 3))))
    plt.legend()

    # Save the plot to a file
    plt.savefig(output_file)
    print(f'MSA plot saved to: {output_file}')
    return fig