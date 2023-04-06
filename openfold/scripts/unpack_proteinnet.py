import argparse
import os
from pathlib import Path


def _write_file(args, file_in_progress):
    file_id = file_in_progress[1]
    fname = file_id.upper() + ".core"
    fpath = os.path.join(args.output_dir, fname)
    with open(fpath, "w") as fp:
        fp.write('\n'.join(file_in_progress))


def main(args):
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)    

    with open(args.proteinnet_file, "r") as fp:
        proteinnet_string = fp.readlines()

    file_in_progress = []
    for line in proteinnet_string:
        if(line == "[ID]\n"):
            if(len(file_in_progress) > 0):
                _write_file(args, file_in_progress)
                file_in_progress = []

        file_in_progress.append(line.strip())

    if(len(file_in_progress) > 0):
        _write_file(args, file_in_progress)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "proteinnet_file", type=str,
        help="Path to ProteinNet file to unpack"
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Path to directory in which to output .core files"
    )

    args = parser.parse_args()

    main(args)
