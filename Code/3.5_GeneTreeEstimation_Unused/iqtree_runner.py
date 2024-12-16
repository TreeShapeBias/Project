from multiprocessing import Pool, cpu_count
import subprocess
from glob import glob
import tqdm

def work(path: str):
    _ = subprocess.call(["iqtree2", "-quiet", "-s", path, "-m", "GTR", "-T", "2"]) 


if __name__ == "__main__":
    paths = glob("gene_trees/**/*.phy", recursive=True)
    pool = Pool(processes=(cpu_count() - 1))

    for _ in tqdm.tqdm(pool.imap_unordered(work, paths), total=len(paths)):
        pass

    pool.close()
    pool.join()

