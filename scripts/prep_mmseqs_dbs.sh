#!/bin/bash
#
# Copyright 2021 AlQuraishi Laboratory 
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
# Downloads and unzips all required data for AlphaFold.
#
# Usage: bash download_all_data.sh /path/to/download/directory
set -e

DOWNLOAD_DIR="$1"

for f in $(ls ${DOWNLOAD_DIR}/*.tar.gz)
do
  tar --extract --verbose --file="${DOWNLOAD_DIR}/${f}" \
      --directory="${DOWNLOAD_DIR}/mmseqs_dbs"
  rm "${f}"
  BASENAME="$(basename {f%%.*})"
  DB_NAME="${BASENAME}_db"
  OLD_PWD=$(pwd)
  cd "${DOWNLOAD_DIR}/mmseqs_dbs" 
  mmseqs tar2exprofiledb "${BASENAME}" "${DB_NAME}"
  mmseqs createindex "${DB_NAME}" "${DOWNLOAD_DIR}/tmp/"
  cd "${OLD_PWD}"
done


