# Copyright 2023 RIKEN & Fujitsu Limited
# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A Python wrapper for OpenFold."""
import os
import signal
import subprocess
import re
import resource
from typing import Sequence

from absl import logging

from openfold.data.tools import utils


class OpenFoldInference:
    """Python wrapper of the OpenFold inference."""

    def __init__(self, script_path: str, exec_path: str = None):
        """Initializes the Python OpenFold wrapper.

        Args:
          script_path: The path to the OpenFold script.
          exec_path: The path for executing OpenFold script

        Raises:
          RuntimeError: If OpenFold script not found within the path.
        """
        self.script_path = script_path
        self.exec_path = exec_path

    def run(
            self,
            fasta_dir: str,
            template_mmcif_dir: str,
            subdir: str,
            args: object,
            timeout: float=None):
        """Run inference.

        Args:
          fasta_dir:
          template_mmcif_dir:
          timeout:

        Returns:
          A string with the alignment in a3m format.

        Raises:
          RuntimeError: If OpenFold fails.
          ValueError: If any of the sequences is less than 6 residues long.
        """
        logging.info("Run inference for {}".format(fasta_dir))

        cmd = [
            "python3",
            self.script_path,
            fasta_dir,
            template_mmcif_dir,
            "--model_device",
            args.model_device,
            "--jax_param_path",
            args.jax_param_path,
            "--output_dir",
            args.output_dir,
            "--use_precomputed_alignments",
            args.use_precomputed_alignments,
            "--max_template_date",
            args.max_template_date,
            "--kalign_binary_path",
            args.kalign_binary_path,
            "--config_preset",
            args.config_preset,
        ]
        if args.obsolete_pdbs_path is not None:
            cmd.append("--obsolete_pdbs_path")
            cmd.append(args.obsolete_pdbs_path)
        if args.release_dates_path is not None:
            cmd.append("--release_dates_path")
            cmd.append(args.release_dates_path)
        if args.data_random_seed is not None:
            cmd.append("--data_random_seed")
            cmd.append(args.data_random_seed)
        if subdir is not None:
            cmd.append("--sub_directory")
            cmd.append(subdir)

        def preexec_fn():
            os.setpgrp()
            if args.max_memory is not None:
                resource.setrlimit(
                    resource.RLIMIT_AS,
                    (args.max_memory, resource.RLIM_INFINITY))

        logging.info('Launching subprocess "%s"', " ".join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.exec_path,
            preexec_fn=preexec_fn,
        )

        with utils.timing("OpenFold inference query"):
            try:
                stdout, stderr = process.communicate(timeout=timeout)
                retcode = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired as e:
                pgid = os.getpgid(process.pid)
                logging.info(f'Terminating the whole process group (pid:{process.pid}, pgid={pgid})...')
                os.killpg(pgid, signal.SIGTERM)
                raise e

            stdout_dec = stdout.decode("utf-8")
            stderr_dec = stderr.decode("utf-8")
            logging.info(
                "OpenFold inference stdout:\n%s\n\nstderr:\n%s\n",
                stdout_dec,
                stderr_dec,
            )

        if retcode:
            raise RuntimeError(
                "OpenFold inference failed\nstdout:\n%s\n\nstderr:\n%s\n"
                % (stdout_dec, stderr_dec)
            )

        ret = dict()
        ret['inference_time'] = float(re.findall('Inference time: *(.*)', stderr_dec)[0])
        ret['relaxation_time'] = float(re.findall('Relaxation time: *(.*)', stderr_dec)[0])

        return ret
