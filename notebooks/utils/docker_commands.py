import docker
import os

def run_inference_for_model(
    run_path,
    database_dir,
    default_fasta_dir_name,
    config_preset, 
    use_precomputed_alignments=True,
    ):
    
    client = docker.from_env()
    # Set up the volumes dictionary
    volumes = {
        run_path: {'bind': '/run_path', 'mode': 'rw'},
    }
    
    volumes[database_dir] = {'bind': '/database', 'mode': 'rw'}
    
    fasta_dir = f"/run_path/{default_fasta_dir_name}/tmp"
    outpur_dir = "/run_path/output"
    precomputed_alignments_dir = "/run_path/output/alignments"
    
    command = [
        "python3", "/opt/openfold/run_pretrained_openfold.py",
        fasta_dir,
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
        "--pdb_seqres_database_path", "/database/pdb_seqres/pdb_seqres.txt",
        "--uniref30_database_path", "/database/uniref30/UniRef30_2021_03",
        "--uniprot_database_path", "/database/uniprot/uniprot.fasta",
        "--hmmsearch_binary_path", "/opt/conda/bin/hmmsearch",
        "--hmmbuild_binary_path", "/opt/conda/bin/hmmbuild",
        # Inference settings
        "--model_device", "cuda:0",
        "--config_preset", config_preset,
        "--save_outputs",
         "--output_dir", outpur_dir
    ]

    if use_precomputed_alignments:
        command.extend(["--use_precomputed_alignments", precomputed_alignments_dir])

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
        
def run_msa_alignment(
    run_path,
    database_dir,
    default_fasta_dir_name,
    cpus_per_task=32,
    no_tasks=1,
):
    client = docker.from_env()
    
    precomputed_alignments_dir = f"{run_path}/output/alignments"
    os.makedirs(precomputed_alignments_dir, exist_ok=True)

    volumes = {
        run_path: {'bind': '/run_path', 'mode': 'rw'},
    }
    
    volumes[database_dir] = {'bind': '/database', 'mode': 'rw'}
    
    fasta_dir = f"/run_path/{default_fasta_dir_name}/tmp"
    precomputed_alignments_dir_docker = "/run_path/output/alignments"
   
    # Define the Docker command
    command = [
        "python3", "/opt/openfold/scripts/precompute_alignments.py",
        fasta_dir,
        precomputed_alignments_dir_docker,
        "--uniref90_database_path", "/database/uniref90/uniref90.fasta",
        "--mgnify_database_path", "/database/mgnify/mgy_clusters_2022_05.fa",
        "--pdb70_database_path", "/database/pdb70/pdb70",
        "--uniclust30_database_path", "/database/uniclust30/uniclust30_2018_08/uniclust30_2018_08",
        "--bfd_database_path", "/database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--cpus_per_task", str(cpus_per_task),
        "--no_tasks", str(no_tasks),
        "--jackhmmer_binary_path", "/opt/conda/bin/jackhmmer",
        "--hhblits_binary_path", "/opt/conda/bin/hhblits",
        "--hhsearch_binary_path", "/opt/conda/bin/hhsearch",
        "--kalign_binary_path", "/opt/conda/bin/kalign"
    ]

    try:
        print("Running Docker container for MSA alignment...")
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
        raise e
    except docker.errors.ImageNotFound as e:
        print(f"ImageNotFound: {e}")
        raise e
    except docker.errors.APIError as e:
        print(f"APIError: {e}")
        raise e