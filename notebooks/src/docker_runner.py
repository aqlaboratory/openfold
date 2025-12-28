import os

class DockerRunner:
    def __init__(self, client, run_path, database_dir, default_fasta_dir_name):
        self.client = client
        self.run_path = run_path
        self.database_dir = database_dir
        self.default_fasta_dir_name = default_fasta_dir_name

    def _setup_volumes(self):
        return {
            self.run_path: {'bind': '/run_path', 'mode': 'rw'},
            self.database_dir: {'bind': '/database', 'mode': 'rw'}
        }

    def _stream_logs(self, container):
        for log in container.logs(stream=True):
            print(log.decode('utf-8'), end='')

    def run_inference_for_model(self, weight_set, config_preset, gpu, use_precomputed_alignments=True):
        command = self._build_inference_command(weight_set, config_preset, gpu, use_precomputed_alignments)
        self._run_container(command)

    def run_msa_alignment(self, cpus_per_task=32, no_tasks=1):
        precomputed_alignments_dir = f"{self.run_path}/output/alignments"
        os.makedirs(precomputed_alignments_dir, exist_ok=True)
        command = self._build_msa_alignment_command(cpus_per_task, no_tasks)
        self._run_container(command)

    def _build_inference_command(self, weight_set, config_preset, gpu, use_precomputed_alignments):
        fasta_dir = f"/run_path/{self.default_fasta_dir_name}/tmp"
        output_dir = "/run_path/output"
        precomputed_alignments_dir = "/run_path/output/alignments"
        
        command = [
            "python3", "/opt/openfold/run_pretrained_openfold.py",
            fasta_dir,
            "/database/pdb_mmcif/mmcif_files/",
            "--uniref90_database_path", "/database/uniref90/uniref90.fasta",
            "--mgnify_database_path", "/database/mgnify/mgy_clusters_2022_05.fa",
            "--pdb70_database_path", "/database/pdb70/pdb70",
            "--uniclust30_database_path", "/database/uniclust30/uniclust30_2018_08/uniclust30_2018_08",
            "--bfd_database_path", "/database/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
            "--jackhmmer_binary_path", "/opt/conda/bin/jackhmmer",
            "--hhblits_binary_path", "/opt/conda/bin/hhblits",
            "--hhsearch_binary_path", "/opt/conda/bin/hhsearch",
            "--kalign_binary_path", "/opt/conda/bin/kalign",
            "--pdb_seqres_database_path", "/database/pdb_seqres/pdb_seqres.txt",
            "--uniref30_database_path", "/database/uniref30/UniRef30_2021_03",
            "--uniprot_database_path", "/database/uniprot/uniprot.fasta",
            "--hmmsearch_binary_path", "/opt/conda/bin/hmmsearch",
            "--hmmbuild_binary_path", "/opt/conda/bin/hmmbuild",
            "--model_device", gpu,
            "--save_outputs",
            "--output_dir", output_dir
        ]
        
        if weight_set == "AlphaFold":
            command.extend(["--config_preset", config_preset])
            
        if weight_set == "OpenFold":
            command.extend(["--openfold_checkpoint_path", f"/database/openfold_params/{config_preset}"])

        if use_precomputed_alignments:
            command.extend(["--use_precomputed_alignments", precomputed_alignments_dir])

        return command

    def _build_msa_alignment_command(self, cpus_per_task, no_tasks):
        fasta_dir = f"/run_path/{self.default_fasta_dir_name}/tmp"
        precomputed_alignments_dir_docker = "/run_path/output/alignments"

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

        return command

    def _run_container(self, command):
        volumes = self._setup_volumes()

        try:
            print("Running Docker container...")
            container = self.client.containers.run(
                image="openfold:latest",
                command=command,
                volumes=volumes,
                runtime="nvidia",
                remove=True,
                detach=True,
                stdout=True,
                stderr=True
            )
            
            self._stream_logs(container)
        
        except docker.errors.ContainerError as e:
            print(f"ContainerError: {e}")
            raise e
        except docker.errors.ImageNotFound as e:
            print(f"ImageNotFound: {e}")
            raise e
        except docker.errors.APIError as e:
            print(f"APIError: {e}")
            raise e