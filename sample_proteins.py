import json
import os
import random
import matplotlib.pyplot as plt


def select_proteins_from_PDB():
    with open('/data/sk844/pdb_seqres.txt', 'r') as f:
        lines = f.readlines()
        max_idx = len(lines)
        print(max_idx)
        fastas = []
        for r in range(1):
            label_idx = random.randrange(0, max_idx, 2)
            print(label_idx)
            fasta = lines[label_idx] + lines[label_idx + 1]
            print(fasta)
            fastas.append(fasta)
        print(len(fastas))
    for fasta in fastas:
        with open('output_dir/sample_proteins/' + fasta[1:7] + '.fasta', 'w') as f:
            f.write(fasta)
    print('File start')
    print(repr(lines[0]))
    print(repr(lines[1]))
    print(repr(lines[2]))

def measure_time():
    timings_all = []
    rootdir = "../alphafold/output_dir/sample_proteins"
    for dir in sorted(os.listdir(rootdir)):
        print(dir)
        if os.path.isfile(dir):
            continue
        with open(os.path.join(rootdir, dir, 'msas', 'timings_msa.json'), 'r') as f:
            timings = json.load(f)
            timings_all.append((dir, timings))

    for protein_name, protein_dict in timings_all:
        with open(os.path.join('output_dir','sample_proteins', protein_name + '.fasta'), 'r') as f:
            fasta_lines = f.readlines()
            protein_length = len(fasta_lines[1].strip())
            protein_dict['length'] = protein_length
    print(timings_all)
    timings_msas = [prot_dict['msas']/60 for _, prot_dict in timings_all if prot_dict['msas']]
    timings_templates_featurization = [prot_dict['templates']/60 for _, prot_dict in timings_all if prot_dict['templates']]
    prot_lengths = [prot_dict['length'] for _, prot_dict in timings_all if prot_dict['msas']]
    # Plot alignment times
    plt.scatter(prot_lengths, timings_msas)
    plt.title("Alignment times")
    plt.xlabel("Protein size")
    plt.ylabel("Execution times (minutes)")
    plt.show()
    # Plot template featurization times
    plt.scatter(prot_lengths, timings_templates_featurization)
    plt.title("Template Featurization times")
    plt.xlabel("Protein size")
    plt.ylabel("Execution times (minutes)")
    plt.show()
    print("Average alignment times: ", sum(timings_msas)/len(timings_msas))
    print("Template featurization times: ", sum(timings_templates_featurization)/len(timings_templates_featurization))

    protein_lengths = []
    for file in sorted(os.listdir(os.path.join('output_dir','sample_proteins'))):
        with open(os.path.join('output_dir', 'sample_proteins', file), 'r') as f:
            fasta_lines = f.readlines()
            protein_lengths.append(len(fasta_lines[1].strip()))
    print('Protein lenghts')
    print(protein_lengths)

    #Now dealing with the execution times on 1-core
    rootdir2 = "../alphafold/output_dir/sample_proteins_1core"
    timings_granular = {}
    for dir in sorted(os.listdir(rootdir2)):
        print(dir)
        if os.path.isfile(dir):
            continue
        with open(os.path.join(rootdir2, dir, 'msas', 'timings_msa.json'), 'r') as f:
            timings = json.load(f)
            timings['msas_1core'] = timings.pop('msas')
            timings['templates_1core'] = timings.pop('templates')
            timings_granular[dir] = timings

    parallel_to_1core_ratio_msa = []
    parallel_to_1core_ratio_template_featurization = []
    protein_lengths_1_core = []
    for protein_name, protein_dict in timings_all:
        if protein_name not in timings_granular:
            continue
        protein_dict.update(timings_granular[protein_name])
        parallel_to_1core_ratio_msa.append(protein_dict['msas_1core'] / protein_dict['msas'])
        parallel_to_1core_ratio_template_featurization.append(protein_dict['templates_1core'] / protein_dict['templates'])
        protein_lengths_1_core.append(protein_dict['length'])
        print(protein_name, protein_dict['msas_1core'] / protein_dict['msas'], '\t', protein_dict)
    print('Average 1-core multiplier for MSAs: ', sum(parallel_to_1core_ratio_msa)/len(parallel_to_1core_ratio_msa))
    print('Average 1-core multiplier for template featurization: ',
          sum(parallel_to_1core_ratio_template_featurization) / len(parallel_to_1core_ratio_template_featurization))

    # Seeings alignment times wrt protein lengths
    plt.scatter(protein_lengths_1_core, parallel_to_1core_ratio_msa)
    plt.show()

    # Now, let's deal with template search and featurization times
    hhsearch_times_1_core = []
    for protein_name, protein_dict in timings_granular.items():
        hhsearch_times_1_core.append(protein_dict['hhsearch'])
    print('Average hhsearch times on 1-core (seconds): ', sum(hhsearch_times_1_core) / len(hhsearch_times_1_core))


if __name__ == '__main__':
    measure_time()

