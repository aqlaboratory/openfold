# Copyright 2023 RIKEN & Fujitsu Limited
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

import os
from openfold.data import mmcif_parsing
import pickle
import argparse
import lz4

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_list', type=str)
    parser.add_argument('output_dir', type=str)
    parser.add_argument('num_procs', type=int)
    args = parser.parse_args()

    print('settings', args)

    rank = int(os.environ.get('PMIX_RANK', '0'))
    print(f'rank = {rank} / {args.num_procs}', flush=True)

    with open(args.input_list) as f:
        lines = f.readlines()

    lines = [line.rstrip('\n') for line in lines]

    lines = lines[rank::args.num_procs]

    for path in lines:
        file_id = os.path.splitext(os.path.basename(path))[0]
        outpathtmp = os.path.join(args.output_dir, f'{file_id}.tmp')
        outpath = os.path.join(args.output_dir, f'{file_id}.pkl')

        if os.path.isfile(outpath):
            continue

        with open(path, 'r') as f:
            mmcif_string = f.read()
        mmcif = mmcif_parsing.parse(
            file_id=file_id, mmcif_string=mmcif_string,
            with_raw_string=False,
            with_structure=False,
        )

        ### pickle
        # with open(outpathtmp, 'wb') as f:
        #     pickle.dump(mmcif, f)

        ### pickle + lz4
        pkl = pickle.dumps(mmcif)
        with lz4.frame.LZ4FrameFile(outpathtmp, 'wb', compression_level=1) as f:
            f.write(pkl)


        os.rename(outpathtmp, outpath)
