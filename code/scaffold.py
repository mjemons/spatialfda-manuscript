import os
import sys
import numpy as np
import subprocess
import hashlib
import random

def str_to_int(string):
    return(int(hashlib.sha256(string.encode('utf-8')).hexdigest(), 16) % 10**7)

seed = (int(snakemake.wildcards['rep']) - 1) * 40 + str_to_int(snakemake.wildcards['typ']) * 4 + str_to_int(snakemake.wildcards['sim'])
 
np.random.seed(seed)
random.seed(seed)

# draw both rmin and rmax from a normal distribution to obtain slightly differening scaffolds

rmin = str(5)
rmax = str(20)

outdir = "-".join(snakemake.output[1].split("-")[:-1])
seed = str(seed) 
code = "software/PowerAnalysisForSpatialOmics/spatialpower/tissue_generation/random_circle_packing.py"

subprocess.run(["python", code, "-x", "1000", "-y", "1000", "--rmin", rmin, "--rmax", rmax, "-s", seed, "-o", outdir])
