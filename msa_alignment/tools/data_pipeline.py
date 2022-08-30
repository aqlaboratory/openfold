import os
import time
from typing import Optional
from multiprocessing import cpu_count, Process
import jackhmmer, hhblits, hhsearch, parsers


class AlignmentRunner:
    """Runs alignment tools and saves the results"""

    def __init__(
        self,
        jackhmmer_binary_path: Optional[str] = None,
        hhblits_binary_path: Optional[str] = None,
        hhsearch_binary_path: Optional[str] = None,
        uniref90_database_path: Optional[str] = None,
        mgnify_database_path: Optional[str] = None,
        bfd_database_path: Optional[str] = None,
        uniclust30_database_path: Optional[str] = None,
        pdb70_database_path: Optional[str] = None,
        use_small_bfd: Optional[bool] = None,
        no_cpus: Optional[int] = None,
        uniref_max_hits: int = 10000,
        mgnify_max_hits: int = 5000,
    ):
        """
        Args:
            jackhmmer_binary_path:
                Path to jackhmmer binary
            hhblits_binary_path:
                Path to hhblits binary
            hhsearch_binary_path:
                Path to hhsearch binary
            uniref90_database_path:
                Path to uniref90 database. If provided, jackhmmer_binary_path
                must also be provided
            mgnify_database_path:
                Path to mgnify database. If provided, jackhmmer_binary_path
                must also be provided
            bfd_database_path:
                Path to BFD database. Depending on the value of use_small_bfd,
                one of hhblits_binary_path or jackhmmer_binary_path must be 
                provided.
            uniclust30_database_path:
                Path to uniclust30. Searched alongside BFD if use_small_bfd is 
                false.
            pdb70_database_path:
                Path to pdb70 database.
            use_small_bfd:
                Whether to search the BFD database alone with jackhmmer or 
                in conjunction with uniclust30 with hhblits.
            no_cpus:
                The number of CPUs available for alignment. By default, all
                CPUs are used.
            uniref_max_hits:
                Max number of uniref hits
            mgnify_max_hits:
                Max number of mgnify hits
        """
        node_id = int(os.environ.get("SLURM_NODEID", 999))
        db_map = {
            "jackhmmer": {
                "binary": jackhmmer_binary_path,
                "dbs": [
                    uniref90_database_path,
                    mgnify_database_path,
                    bfd_database_path if use_small_bfd else None,
                ],
            },
            "hhblits": {
                "binary": hhblits_binary_path,
                "dbs": [
                    bfd_database_path if not use_small_bfd else None,
                ],
            },
            "hhsearch": {
                "binary": hhsearch_binary_path,
                "dbs": [
                    pdb70_database_path,
                ],
            },
        }

        for name, dic in db_map.items():
            binary, dbs = dic["binary"], dic["dbs"]
            if (binary is None and not all([x is None for x in dbs])):
                raise ValueError(
                    f"{name} DBs provided but {name} binary is None"
                )

        if (not all([x is None for x in db_map["hhsearch"]["dbs"]])
                and uniref90_database_path is None):
            raise ValueError(
                """uniref90_database_path must be specified in order to perform
                   template search"""
            )

        self.uniref_max_hits = uniref_max_hits
        self.mgnify_max_hits = mgnify_max_hits
        self.use_small_bfd = use_small_bfd

        if (no_cpus is None):
            no_cpus = cpu_count()

        self.jackhmmer_uniref90_runner = None
        if (jackhmmer_binary_path is not None and
                    uniref90_database_path is not None
                ):
            self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniref90_database_path,
                n_cpu=no_cpus,
            )
            print(f"Node{node_id}:uniref90 is prepared")

        self.jackhmmer_small_bfd_runner = None
        self.hhblits_bfd_uniclust_runner = None
        if (bfd_database_path is not None):
            if use_small_bfd:
                self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                    binary_path=jackhmmer_binary_path,
                    database_path=bfd_database_path,
                    n_cpu=no_cpus,
                )
            else:
                dbs = [bfd_database_path]
                if (uniclust30_database_path is not None):
                    dbs.append(uniclust30_database_path)
                self.hhblits_bfd_uniclust_runner = hhblits.HHBlits(
                    binary_path=hhblits_binary_path,
                    databases=dbs,
                    n_cpu=no_cpus,
                )
                print(f"Node{node_id}:bfd and uniclust_30 is prepared")

        self.jackhmmer_mgnify_runner = None
        if (mgnify_database_path is not None):
            self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=mgnify_database_path,
                n_cpu=no_cpus,
            )
            print(f"Node{node_id}:magnify is prepared")

        self.hhsearch_pdb70_runner = None
        if (pdb70_database_path is not None):
            self.hhsearch_pdb70_runner = hhsearch.HHSearch(
                binary_path=hhsearch_binary_path,
                databases=[pdb70_database_path],
                n_cpu=no_cpus,
            )
            print(f"Node{node_id}:pdb70 is prepared")

    def uniref90_pdb70_search_engine(self,
                                     file_id: str,
                                     fasta_path: str, output_dir: str):
        uniref90_start_time = time.time()
        print(f"{file_id}_uniref90 is searching")
        jackhmmer_uniref90_result = self.jackhmmer_uniref90_runner.query(
            fasta_path
        )[0]
        uniref90_end_time = time.time()
        print(f"{file_id} uniref90 is completed, time passed: {round((uniref90_end_time - uniref90_start_time) / 60, 2)} mins")
        uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(
            jackhmmer_uniref90_result["sto"],
            max_sequences=self.uniref_max_hits
        )
        uniref90_out_path = os.path.join(output_dir, "uniref90_hits.a3m")
        with open(uniref90_out_path, "w") as f:
            f.write(uniref90_msa_as_a3m)

        if (self.hhsearch_pdb70_runner is not None):
            print(f"{file_id}_pdb70 is searching")
            pdb70_start_time = time.time()
            hhsearch_result = self.hhsearch_pdb70_runner.query(
                uniref90_msa_as_a3m
            )
            pdb70_end_time = time.time()
            print(
                f"{file_id} pdb70 is completed, time passed: {round((pdb70_end_time - pdb70_start_time) / 60, 2)} mins")
            pdb70_out_path = os.path.join(output_dir, "pdb70_hits.hhr")
            with open(pdb70_out_path, "w") as f:
                f.write(hhsearch_result)

    def mgnify_search_engine(self, file_id, fasta_path, output_dir):
        print(f"{file_id}_mgnify is searching")
        mgnify_start_time = time.time()
        jackhmmer_mgnify_result = self.jackhmmer_mgnify_runner.query(
            fasta_path
        )[0]
        mgnify_end_time = time.time()
        print(
            f"{file_id} magnify is completed, time passed: {round((mgnify_end_time - mgnify_start_time) / 60, 2)} mins")
        mgnify_msa_as_a3m = parsers.convert_stockholm_to_a3m(
            jackhmmer_mgnify_result["sto"],
            max_sequences=self.mgnify_max_hits
        )
        mgnify_out_path = os.path.join(output_dir, "mgnify_hits.a3m")
        with open(mgnify_out_path, "w") as f:
            f.write(mgnify_msa_as_a3m)

    def bfd_uniclust_search_engine(self, file_id, fasta_path, output_dir):
        print(f"{file_id}_bfd is searching")
        bfd_start_time = time.time()
        hhblits_bfd_uniclust_result = (
            self.hhblits_bfd_uniclust_runner.query(fasta_path)
        )
        bfd_end_time = time.time()
        print(
            f"{file_id} bfd is completed, time passed: {round((bfd_end_time - bfd_start_time) / 60, 2)} mins")
        if output_dir is not None:
            bfd_out_path = os.path.join(
                output_dir, "bfd_uniclust_hits.a3m")
            with open(bfd_out_path, "w") as f:
                f.write(hhblits_bfd_uniclust_result["a3m"])

    def run(
        self,
        fasta_path: str,
        output_dir: str,
    ):
        """Runs alignment tools on a sequence"""
        file_id = output_dir.split("/")[-1]
        funcs = []
        process_list = []
        if (self.jackhmmer_uniref90_runner is not None):
            funcs.append(self.uniref90_pdb70_search_engine)
        if (self.jackhmmer_mgnify_runner is not None):
            funcs.append(self.mgnify_search_engine)
        if (self.use_small_bfd and self.jackhmmer_small_bfd_runner is not None):
            jackhmmer_small_bfd_result = self.jackhmmer_small_bfd_runner.query(
                fasta_path
            )[0]
            bfd_out_path = os.path.join(output_dir, "small_bfd_hits.sto")
            with open(bfd_out_path, "w") as f:
                f.write(jackhmmer_small_bfd_result["sto"])
        elif (self.hhblits_bfd_uniclust_runner is not None):
            funcs.append(self.bfd_uniclust_search_engine)
        
        for i in range(len(funcs)):
            p = Process(target=funcs[i], args=(file_id, fasta_path, output_dir))
            p.start()
            process_list.append(p)
        
        for i in process_list:
            p.join()
        