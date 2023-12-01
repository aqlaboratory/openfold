import os, multiprocessing, argparse, json
from openfold.data import parsers
def parse_stockholm_file(alignment_dir: str, stockholm_file: str, queue: multiprocessing.Queue):
                 path = os.path.join(alignment_dir, stockholm_file)
                 file_name,_ = os.path.splitext(stockholm_file)
                 with open(path, "r") as infile:
                     msa = parsers.parse_stockholm(infile.read())
                     infile.close()
                 queue.put({file_name: msa})

def parse_a3m_file(alignment_dir: str, a3m_file: str,queue: multiprocessing.Queue):
    path = os.path.join(alignment_dir, a3m_file)
    file_name,_ = os.path.splitext(a3m_file)
    with open(path, "r") as infile:
        msa = parsers.parse_a3m(infile.read())
        infile.close()
    queue.put({file_name: msa})

def run_parse_all_msa_files_multiprocessing(stockholm_files: list, a3m_files: list, alignment_dir:str):
    print(f"#### line 764 start running in multiprocessing way")
    msa_results={}
    processes = []
    queue = multiprocessing.Queue()
    for f in stockholm_files: 
        process = multiprocessing.Process(target = parse_stockholm_file, args=(alignment_dir, f, queue))
        process.deamon = False
        processes.append(process)
        process.start()
    for f in a3m_files:
        process = multiprocessing.Process(target = parse_a3m_file, args=(alignment_dir, f, queue))
        process.daemon = False  
        processes.append(process)
        process.start()
    
    for p in processes:
        res = queue.get()
        msa_results.update(res)
        p.join()

def main(alignment_dir):  
    parser = argparse.ArgumentParser(description='Process msa files in parallel')
    parser.add_argument('alignment_dir', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
    args = parser.parse_args()
    stockholm_files = [i for i in os.listdir(alignment_dir) if (i.endswith('.sto') and ("hmm_output" not in i))]
    a3m_files = [i for i in os.listdir(alignment_dir) if i.endswith('.a3m')]
    msa_data = run_parse_all_msa_files_multiprocessing(stockholm_files, a3m_files, alignment_dir)
    return json.dumps(msa_data)

if __name__ == "__main__":
    main()