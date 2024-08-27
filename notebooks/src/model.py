import os
import shutil
from .docker_commands import run_inference_for_model, run_msa_alignment
from .utils import validate_sequence, generate_random_run_id, write_sequences_to_fasta, generate_individual_sequence_files

class InferenceClientOpenFold:
    default_fasta_dir_name = 'fasta_dir'

    def __init__(self, database_dir):
        self.database_dir = database_dir

    def run_inference(self, weight_set, model_name, inference_input=None, use_precomputed_alignments=False, run_id=None):
        os.makedirs('data', exist_ok=True)

        if use_precomputed_alignments:
            if not run_id:
                raise ValueError("run_id is required when using pre-computed alignments.")
        else:
            if not inference_input:
                raise ValueError("inference_input is required to compute alignments.")

        if not run_id:
            run_id = generate_random_run_id()
            run_path = f"{os.getcwd()}/data/run_{run_id}"
            fasta_dir_root_path = f'{run_path}/{self.default_fasta_dir_name}'

            os.makedirs(fasta_dir_root_path, exist_ok=True)
            print(f"Fasta root directory: {fasta_dir_root_path}")

            if os.path.isfile(inference_input):
                shutil.copy2(inference_input, f"{fasta_dir_root_path}/sequences.fasta")
                print(f"Using input file: {inference_input}")
            else:
                sequence_string = inference_input
                print(f"Input is a sequence string: {sequence_string}")

                sequence_string = validate_sequence(inference_input, weight_set)
                fasta_file = write_sequences_to_fasta(sequence_string.split(':'), fasta_dir_root_path)
                print(f"Sequences written to FASTA file: {fasta_file}")

            generate_individual_sequence_files(run_path)
           
            self.run_msa_alignment(run_path)
            self.run_model_with(run_path, weight_set, model_name)

        elif run_id and use_precomputed_alignments:
            run_path = f"{os.getcwd()}/data/run_{run_id}"

            if os.path.isdir(run_path):
                self.run_model_with(run_path, weight_set, model_name)
            else:
                raise ValueError("Provided Run ID does not exist")

        return run_id

    def run_model_with(self, run_path, weight_set, model_name):
        config_preset_list = self.get_config_preset_list_for_model(weight_set, model_name)

        for config_preset in config_preset_list:
            run_inference_for_model(
                run_path,
                self.database_dir,
                self.default_fasta_dir_name,
                config_preset
            )
            
    def run_msa_alignment(self, run_path):
        run_msa_alignment(run_path, self.database_dir, self.default_fasta_dir_name)

    def get_config_preset_list_for_model(self, weight_set, model_name):
        if weight_set == "OpenFold" and model_name == "monomer":
            model_names = [
                'finetuning_3.pt',
                'finetuning_4.pt',
                'finetuning_5.pt',
                'finetuning_ptm_2.pt',
                'finetuning_no_templ_ptm_1.pt'
            ]
        elif weight_set == "AlphaFold" and model_name == "multimer":
            model_names = [f'model_{i}_multimer_v3' for i in range(1, 6)]
        elif weight_set == "AlphaFold" and model_name == "monomer":
            model_names = [f'model_{i}' for i in range(1, 6)]
        else:
            raise ValueError(f"Invalid combination of weight_set '{weight_set}' and model_name '{model_name}'")

        return model_names