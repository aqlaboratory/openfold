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

"""Common utilities for data pipeline tools."""
import contextlib
import datetime
import logging
import shutil
import tempfile
import time
import pickle
import os
import lz4
from typing import Optional

from openfold.data import mmcif_parsing

@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@contextlib.contextmanager
def timing(msg: str):
    logging.info("Started %s", msg)
    tic = time.perf_counter()
    yield
    toc = time.perf_counter()
    logging.info("Finished %s in %.3f seconds", msg, toc - tic)


def to_date(s: str):
    return datetime.datetime(
        year=int(s[:4]), month=int(s[5:7]), day=int(s[8:10])
    )

def load_cif(path: str, hit_pdb_code: str):
    filename, ext = os.path.splitext(os.path.basename(path))

    if ext == ".cif" and os.path.isfile(path):
        with open(path, "r") as cif_file:
            cif_string = cif_file.read()

            parsing_result = mmcif_parsing.parse(
                file_id=hit_pdb_code, mmcif_string=cif_string
            )
    else:
        pkl_path = os.path.join(os.path.dirname(path), filename + ".pkl")

        ### pickle
        # with open(pkl_path, "rb") as pkl_file:
        #     parsing_result = pickle.load(pkl_file)

        ### pickle + lz4
        with lz4.frame.LZ4FrameFile(pkl_path, 'rb') as pkl_file:
            pkl = pkl_file.read()
            parsing_result = pickle.loads(pkl)

    return parsing_result
