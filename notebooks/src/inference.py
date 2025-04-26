import os
import shutil
from .docker_runner import DockerRunner
from datetime import datetime
from .utils import validate_sequence, generate_random_run_id, write_sequences_to_fasta, generate_individual_sequence_files, get_run_folder_by_id, get_config_preset_list_for_model

class InferenceClientOpenFold:
    default_fasta_dir_name = 'fasta_dir'

    def __init__(self, database_dir, docker_client):
        self.database_dir = database_dir
        self.docker_client = docker_client
        self.docker_runner = None

    def run_inference(self, weight_set, model_name, inference_input=None, use_precomputed_alignments=False, run_id=None, gpu="cuda:0"):
        os.makedirs('data', exist_ok=True)

        if use_precomputed_alignments:
            if not run_id:
                raise ValueError("run_id is required when using pre-computed alignments.")
            return self._run_with_precomputed_alignments(run_id, weight_set, model_name, gpu)
        
        if not inference_input:
            raise ValueError("inference_input is required to compute alignments.")
        
        if not run_id:
            run_id = generate_random_run_id()
        
        return self._run_new_inference(run_id, weight_set, model_name, inference_input, gpu)

    def _run_new_inference(self, run_id, weight_set, model_name, inference_input, gpu):
        # Get the current date and time in dd_MM_yy_hh_mm_ss format
        current_datetime = datetime.now().strftime('%d_%m_%y_%H_%M_%S')
        run_path = os.path.join(os.getcwd(), 'data', f'run_{current_datetime}_{run_id}')
        # run_path = os.path.join(os.getcwd(), 'data', f'run_{run_id}')
        fasta_dir_root_path = os.path.join(run_path, self.default_fasta_dir_name)
        self._initialize_docker_runner(run_path)

        self._prepare_fasta_directory(fasta_dir_root_path, inference_input, weight_set)
        generate_individual_sequence_files(run_path)

        self.run_msa_alignment()
        self.run_model_with_preset(run_path, weight_set, model_name, gpu)
        
        return run_id

    def _run_with_precomputed_alignments(self, run_id, weight_set, model_name, gpu):
        run_path = get_run_folder_by_id(run_id)
        print(f"Using pre-computed alignments from: {run_path}")
        if not os.path.isdir(run_path):
            raise ValueError(f"Provided Run ID '{run_id}' does not exist.")

        self._initialize_docker_runner(run_path)
        self.run_model_with_preset(run_path, weight_set, model_name, gpu)
        
        return run_id

    def _initialize_docker_runner(self, run_path):
        self.docker_runner = DockerRunner(self.docker_client, run_path, self.database_dir, self.default_fasta_dir_name)

    def _prepare_fasta_directory(self, fasta_dir_root_path, inference_input, weight_set):
        os.makedirs(fasta_dir_root_path, exist_ok=True)
        print(f"Fasta root directory: {fasta_dir_root_path}")

        if os.path.isfile(inference_input):
            self._copy_input_file(inference_input, fasta_dir_root_path)
        else:
            self._write_sequence_to_fasta(inference_input, fasta_dir_root_path, weight_set)

    def _copy_input_file(self, input_file, fasta_dir_root_path):
        destination = os.path.join(fasta_dir_root_path, 'sequences.fasta')
        shutil.copy2(input_file, destination)
        print(f"Using input file: {input_file}")

    def _write_sequence_to_fasta(self, sequence_string, fasta_dir_root_path, weight_set):
        validated_sequence = validate_sequence(sequence_string, weight_set)
        fasta_file = write_sequences_to_fasta(validated_sequence.split(':'), fasta_dir_root_path)
        print(f"Sequences written to FASTA file: {fasta_file}")

    def run_model_with_preset(self, run_path, weight_set, model_name, gpu):
        config_preset_list = get_config_preset_list_for_model(weight_set, model_name)
        for config_preset in config_preset_list:
            self.docker_runner.run_inference_for_model(weight_set, config_preset, gpu)
            
    def run_msa_alignment(self, cpus_per_task=32, no_tasks=1):
        self.docker_runner.run_msa_alignment(cpus_per_task=cpus_per_task, no_tasks=no_tasks)