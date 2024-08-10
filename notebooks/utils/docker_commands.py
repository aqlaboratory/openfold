from .model import get_config_preset_list_for_model
import docker
import random
import string
import os

def generate_random_directory_name(prefix="fasta_dir_", length=6):
    return prefix + ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

def generate_random_sequence_name(length=8):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

def write_sequences_to_fasta(sequences, fasta_dir):
    os.makedirs(f'inference/{fasta_dir}', exist_ok=True)
    fasta_file = os.path.join(f'inference/{fasta_dir}', "sequences.fasta")
    
    with open(fasta_file, 'w') as f:
        for seq in sequences:
            sequence_name = f"sequence_{generate_random_sequence_name()}"
            f.write(f">{sequence_name}\n")
            f.write(f"{seq}\n")
    
    return fasta_dir

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


def run_inference(
    inference_input,
    database_dir,
    weight_set,
    model_name, 
    use_precomputed_alignments=False):

    default_fasta_dir_name = ''

    # Determine if the inference_input is a directory or a sequence string
    if os.path.isdir(inference_input):
        fasta_dir = inference_input
        default_fasta_dir_name = 'fasta_dir'
        print(f"Input is a directory: {fasta_dir}/{default_fasta_dir_name}")
    else:
        # Treat inference_input as a sequence string
        sequence_string = inference_input
        print(f"Input is a sequence string: {sequence_string}")
        
        # Validate and clean the sequence string
        sequence_string = validate_sequence(inference_input, weight_set)
        
        os.makedirs('inference', exist_ok=True)
        
        fasta_dir = f'{os.getcwd()}/inference'

        # Generate a random directory name for storing the fasta file
        default_fasta_dir_name  = generate_random_directory_name()
        print(f"Generated directory: {default_fasta_dir_name}")
        
        inference_dir = write_sequences_to_fasta(sequence_string.split(':'), default_fasta_dir_name)
        print(f"Sequences written to FASTA file in directory: {inference_dir}")
        
       
    config_preset_list = get_config_preset_list_for_model(weight_set, model_name)
    
    for config_preset in config_preset_list:
        run_inference_for_model(
            inference_input, 
            database_dir, 
            fasta_dir, 
            default_fasta_dir_name, 
            config_preset, 
            use_precomputed_alignments
        )
        
    
def run_inference_for_model(
    inference_input, 
    database_dir,
    fasta_dir,
    default_fasta_dir_name,
    config_preset, 
    use_precomputed_alignments,
    ):
    
    client = docker.from_env()

    # Set up the volumes dictionary
    volumes = {
        fasta_dir: {'bind': '/inference', 'mode': 'rw'},
    }
    
    volumes[database_dir] = {'bind': '/database', 'mode': 'rw'}

    command = [
        "python3", "/opt/openfold/run_pretrained_openfold.py",
        f"/inference/{default_fasta_dir_name}",
        "/database/pdb_mmcif/mmcif_files/",
        # Databases
        "--uniref90_database_path", "/database/uniref90/uniref90.fasta",
        "--mgnify_database_path", "/database/mgnify/mgy_clusters_2022_05.fa",
        "--pdb70_database_path", "/database/pdb70/pdb70",
        "--uniclust30_database_path", "/database/uniclust30/uniclust30_2018_08/uniclust30_2018_08",
        "--bfd_database_path", "/database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        # Search MSA
        "--jackhmmer_binary_path", "/opt/conda/bin/jackhmmer",
        "--hhblits_binary_path", "/opt/conda/bin/hhblits",
        "--hhsearch_binary_path", "/opt/conda/bin/hhsearch",
        "--kalign_binary_path", "/opt/conda/bin/kalign",
        # Inference settings
        "--model_device", "cuda:0",
        "--config_preset", config_preset,
        "--save_outputs",
        "--output_dir", f"/inference/results/{config_preset}",
    ]

    if use_precomputed_alignments:
        command.extend(["--use_precomputed_alignments", "/inference/alignments"])

    try:
        print("Running Docker container for inference...")
        container = client.containers.run(
            image="openfold:latest",
            command=command,
            volumes=volumes,
            runtime="nvidia",  # This is required for GPU usage
            remove=True,
            detach=True,  # Running in the background to allow log streaming
            stdout=True,
            stderr=True
        )
        
        # Stream logs in real-time
        for log in container.logs(stream=True):
            print(log.decode('utf-8'), end='')
        
    except docker.errors.ContainerError as e:
        print(f"ContainerError: {e}")
    except docker.errors.ImageNotFound as e:
        print(f"ImageNotFound: {e}")
    except docker.errors.APIError as e:
        print(f"APIError: {e}")