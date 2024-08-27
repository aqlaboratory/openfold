import random
import string
import os


def generate_random_run_id(length=6):
    run_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))
    print(f"Run ID: {run_id}")
    return run_id

def generate_random_sequence_name(length=8):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

def write_sequences_to_fasta(sequences, fasta_dir_root_path):
    # os.makedirs(f'inference/{fasta_dir}', exist_ok=True)
    fasta_file = os.path.join(fasta_dir_root_path, "sequences.fasta")
    
    with open(fasta_file, 'w') as f:
        for seq in sequences:
            sequence_name = f"sequence_{generate_random_sequence_name()}"
            f.write(f">{sequence_name}\n")
            f.write(f"{seq}\n")
    
    return fasta_file

def validate_sequence(input_sequence, weight_set):
    # Remove all whitespaces, tabs, and end lines; convert to upper-case
    input_sequence = input_sequence.translate(str.maketrans('', '', ' \n\t')).upper()
    aatypes = set('ACDEFGHIKLMNPQRSTVWY')  # 20 standard amino acids
    allowed_chars = aatypes.union({':'})
    
    if not set(input_sequence).issubset(allowed_chars):
        raise Exception(f'Input sequence contains non-amino acid letters: {set(input_sequence) - allowed_chars}. OpenFold only supports 20 standard amino acids as inputs.')
    
    if ':' in input_sequence and weight_set != 'AlphaFold':
        raise ValueError('Input sequence is a multimer, must select Alphafold weight set')
    
    return input_sequence

def generate_individual_sequence_files(output_path):
    with open(f"{output_path}/fasta_dir/sequences.fasta", "r") as infile:
        sequence_id = None
        sequence_lines = []
        
        output_path = f"{output_path}/fasta_dir/tmp"
        os.makedirs(output_path, exist_ok=True)
        
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id is not None:
                    # Save the previous sequence to a file
                    output_file = os.path.join(output_path, f"{sequence_id}.fasta")
                    with open(output_file, "w") as outfile:
                        outfile.write(f">{sequence_id}\n")
                        outfile.write("\n".join(sequence_lines) + "\n")
                    print(f"Saved {sequence_id} to {output_file}")

                # Start a new sequence
                sequence_id = line[1:].split('.')[0]  # Remove '>' and split by '.' to remove the suffix
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if sequence_id is not None:
            output_file = os.path.join(output_path, f"{sequence_id}.fasta")
            with open(output_file, "w") as outfile:
                outfile.write(f">{sequence_id}\n")
                outfile.write("\n".join(sequence_lines) + "\n")
            print(f"Saved {sequence_id} to {output_file}")
