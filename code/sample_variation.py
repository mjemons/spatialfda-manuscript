import numpy as np
import pandas as pd
import json
import hashlib
import random


def str_to_int(string):
        return(int(hashlib.sha256(string.encode('utf-8')).hexdigest(), 16) % 10**7)

seed = (int(snakemake.wildcards['rep']) - 1) * 100 + str_to_int(snakemake.wildcards['typ']) * 4 + 9 * str_to_int(snakemake.wildcards['sim'])

np.random.seed(seed)
random.seed(seed)

outdir = "/".join(str(snakemake.output).split("/")[:-1])

sample = np.random.uniform(-0.2,0.2)

#write to txt
with open(str(snakemake.output), "w") as f:
    f.write(str(sample))
