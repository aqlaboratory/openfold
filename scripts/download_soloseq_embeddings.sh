#!/bin/bash
#
# Copyright 2024 AlQuraishi Laboratory
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
#
# Downloads ESM-1b embeddings used to train OpenFold SoloSeq single-seq model.
#
# Usage: bash download_soloseq_embeddings.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aws &> /dev/null ; then
    echo "Error: aws could not be found. Please install aws."
    exit 1
fi

DOWNLOAD_DIR="${1}/soloseq_embeddings"
mkdir -p "${DOWNLOAD_DIR}"
aws s3 cp --no-sign-request --region us-east-1 s3://openfold/soloseq_embeddings/ "${DOWNLOAD_DIR}" --recursive
