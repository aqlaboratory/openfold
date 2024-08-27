import os
import shutil
from .docker_commands import run_inference_for_model, run_msa_alignment
from .utils import validate_sequence, generate_random_run_id, write_sequences_to_fasta, generate_individual_sequence_files

default_fasta_dir_name = 'fasta_dir'

def run_inference(
    inference_input,
    database_dir,
    weight_set,
    model_name, 
    use_precomputed_alignments=False,
    run_id=None):
    
    os.makedirs('data', exist_ok=True)
    
    if use_precomputed_alignments:
        if not run_id:
            raise ValueError(f"run_id is required when using pre computed aligments.")
    
    if not run_id:
        
        # Generate a random directory name for storing the fasta file
        run_id  = generate_random_run_id()
                
        run_path = f"{os.getcwd()}/data/run_{run_id}"
        fasta_dir_root_path = f'{run_path}/{default_fasta_dir_name}'
            
        os.makedirs(fasta_dir_root_path, exist_ok=True)
            
        print(f"Fasta root directory: {fasta_dir_root_path}")
        
        # Determine if the inference_input is a directory or a sequence string
        if os.path.isfile(inference_input):
            shutil.copy2(inference_input, f"{fasta_dir_root_path}/sequences.fasta")
            
            print(f"Using input file : {inference_input}")
        else:
            # Treat inference_input as a sequence string
            sequence_string = inference_input
            print(f"Input is a sequence string: {sequence_string}")
            
            # Validate and clean the sequence string
            sequence_string = validate_sequence(inference_input, weight_set)
            
            fasta_file = write_sequences_to_fasta(sequence_string.split(':'), fasta_dir_root_path)
            print(f"Sequences written to FASTA file: {fasta_file}")
            
            
        generate_individual_sequence_files(run_path)
            
        run_msa_alignment(run_path, database_dir, default_fasta_dir_name)
            
        run_model_with(run_path, database_dir, default_fasta_dir_name, weight_set, model_name)
    
    elif run_id and use_precomputed_alignments:
        # Run inference with given Run ID
        run_path = f"{os.getcwd()}/data/run_{run_id}"
        
        if os.path.isdir(run_path):
            run_model_with(run_path, database_dir, default_fasta_dir_name, weight_set, model_name)
        else: 
            raise ValueError(f"Provided Run ID does not exists")
            
def run_model_with(run_path, database_dir, default_fasta_dir_name, weight_set, model_name):
    config_preset_list = get_config_preset_list_for_model(weight_set, model_name)
        
    for config_preset in config_preset_list:
        run_inference_for_model(
            run_path,
            database_dir, 
            default_fasta_dir_name, 
            config_preset
        )     

def get_config_preset_list_for_model(weight_set, model_name):
    if weight_set == "OpenFold" and model_name == "monomer":
        model_names = [
            'finetuning_3.pt',
            'finetuning_4.pt',
            'finetuning_5.pt',
            'finetuning_ptm_2.pt',
            'finetuning_no_templ_ptm_1.pt'
        ]
    elif weight_set == "AlphaFold" and model_name == "multimer":
        # Generate 'model_1_multimer_v3' to 'model_5_multimer_v3' using a loop
        model_names = [f'model_{i}_multimer_v3' for i in range(1, 6)]
    elif weight_set == "AlphaFold" and model_name == "monomer":
        # Generate 'model_1' to 'model_5' using a loop
        model_names = [f'model_{i}' for i in range(1, 6)]
    else:
        raise ValueError(f"Invalid combination of weight_set '{weight_set}' and model_name '{model_name}'")
    
    return model_names