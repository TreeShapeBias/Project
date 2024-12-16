from multiprocessing import Pool, cpu_count
import subprocess
from glob import glob
import tqdm
from os.path import dirname

def work(path: str):
    out_path = dirname(path) + "/estim_tree.tre"
    out_path_log = dirname(path) + "/estim_tree.tre.log"
    with open(out_path_log, "w") as outfile:
        _ = subprocess.call(["./astral4", "-i", path, "-o", out_path, "-t", "2"], stdout=outfile) 


if __name__ == "__main__":
    paths = glob("gene_trees_concat/**/all_g_trees.tre", recursive=True)
    pool = Pool(processes=(cpu_count() - 1))

    for _ in tqdm.tqdm(pool.imap_unordered(work, paths), total=len(paths)):
        pass

    pool.close()
    pool.join()

