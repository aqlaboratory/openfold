import argparse
import json
import os


""" Unifies databases created with create_alignment_db.py """


def main(args):
    super_index = {}
    for f in os.listdir(args.alignment_db_dir):
        if(not os.path.splitext(f)[-1] == ".index"):
            continue
        
        with open(os.path.join(args.alignment_db_dir, f), "r") as fp:
            index = json.load(fp)

        db_name = f"{os.path.splitext(f)[0]}.db"
        
        for k in index:
            super_index[k] = {
                "db": db_name,
                "files": index[k],
            }

    with open(os.path.join(args.output_dir, "super.index"), "w") as fp:
        json.dump(super_index, fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_db_dir", type=str, help="Path to directory containing alignment_dbs")
    parser.add_argument("output_dir", type=str, help="Path in which to output super index")

    args = parser.parse_args()

    main(args)
