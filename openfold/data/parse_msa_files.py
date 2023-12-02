import os, multiprocessing, argparse, pickle, tempfile
import multiprocessing.pool # Need to import multiprocessing.pool first otherwise multiprocessing.pool.Pool cannot be called
from openfold.data import parsers
import contextlib
def parse_stockholm_file(alignment_dir: str, stockholm_file: str):
    path = os.path.join(alignment_dir, stockholm_file)
    file_name,_ = os.path.splitext(stockholm_file)
    with open(path, "r") as infile:
        msa = parsers.parse_stockholm(infile.read())
        infile.close()
    # queue.put()
    return {file_name: msa}

def parse_a3m_file(alignment_dir: str, a3m_file: str):
    path = os.path.join(alignment_dir, a3m_file)
    file_name,_ = os.path.splitext(a3m_file)
    with open(path, "r") as infile:
        msa = parsers.parse_a3m(infile.read())
        infile.close()
    # queue.put({file_name: msa})
    return {file_name: msa}

def run_parse_all_msa_files_multiprocessing(stockholm_files: list, a3m_files: list, alignment_dir:str):
    msa_results={}
    processes = []
    a3m_tasks = [(alignment_dir, f) for f in a3m_files]
    sto_tasks = [(alignment_dir, f) for f in stockholm_files]
    with multiprocessing.pool.Pool(len(a3m_tasks) + len(sto_tasks)) as pool:
        a3m_results = pool.starmap_async(parse_a3m_file, a3m_tasks).get()
        sto_results = pool.starmap_async(parse_stockholm_file, sto_tasks).get()
        for res in [*a3m_results, *sto_results]:
            msa_results.update(res)
        return msa_results

def main():  
    parser = argparse.ArgumentParser(description='Process msa files in parallel')
    parser.add_argument('--alignment_dir', type=str, help='path to alignment dir')
    args = parser.parse_args()
    alignment_dir = args.alignment_dir
    stockholm_files = [i for i in os.listdir(alignment_dir) if (i.endswith('.sto') and ("hmm_output" not in i))]
    a3m_files = [i for i in os.listdir(alignment_dir) if i.endswith('.a3m')]
    msa_data = run_parse_all_msa_files_multiprocessing(stockholm_files, a3m_files, alignment_dir)
    with tempfile.NamedTemporaryFile('wb', suffix='.pkl', delete=False) as outfile:
        pickle.dump(msa_data, outfile)
        print(outfile.name)

if __name__ == "__main__":
    main()