# code adapted from "Power analysis for spatial omics" Nature Methods 2023
# Authors: Ethan A. G. Baker, Denis Schapiro, Bianca Dumitrascu, Sanja Vickovic & Aviv Regev 

# BSD 3-Clause License

# Copyright (c) 2022, Broad Institute
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import pandas as pd
import sys
import json
import hashlib
import random

def str_to_int(string):
        return(int(hashlib.sha256(string.encode('utf-8')).hexdigest(), 16) % 10**7)

seed = (int(snakemake.wildcards['rep']) - 1) * 100 + str_to_int(snakemake.wildcards['typ']) * 4 + 9 * str_to_int(snakemake.wildcards['sim']) + str_to_int(snakemake.wildcards['fov'])

np.random.seed(seed)
random.seed(seed)

path = "software/PowerAnalysisForSpatialOmics"
sys.path.append(str(path))


outdir = "/".join(str(snakemake.output).split("/")[:-1])

import spatialpower.tissue_generation.assign_labels as assign_labels
import spatialpower.tissue_generation.visualization as viz
import spatialpower.neighborhoods.permutationtest as perm_test
import networkx as nx 

# load adjacaceny matrix of scaffold
A = np.load(str(snakemake.input['adj']))
# load centroids
C = np.load(str(snakemake.input['cent']))
STD = {"patient1": 0.25, "patient2": 0.25, "patient3": 0.25, "patient4": 0.25, "patient5": 0.25}
std = float(snakemake.wildcards["std"]) * STD[snakemake.wildcards['sim']]

with open(snakemake.input['var'], "r") as f:
    shift = float(f.read())
print("shift is " + str(shift))
prop = float(snakemake.wildcards['prop'])

def truncated_normal_sample(loc, scale):
    # make truncation symmetric
    delta = min(abs(0-loc), abs(1-loc))
    while True:
        sample = np.random.normal(loc=loc, scale=scale)
        if (loc - delta) <= sample <= (loc + delta):
            return sample

                    

prop = truncated_normal_sample(loc = prop, scale = std)

# define cell type probabilities
cell_type_probabilities = np.array([prop, 0.1, 0.1, 0.1])

#row normalise in order to have the changed proportion as fixed
cell_type_probabilities[1:] = (1 - cell_type_probabilities[0]) / len(cell_type_probabilities[1:])
print(cell_type_probabilities)

loc = float(snakemake.wildcards['prob']) + float(shift)

# obtain the simulation specific interaction probability

p_perturbed = truncated_normal_sample(loc = loc, scale = std)
print(p_perturbed)
# define neighbourhood probabilities
neighborhood_probabilities = np.array(([p_perturbed, 0.25, 0.25, 0.25],
                                        [0.25, 0.25, 0.25, 0.25],
                                        [0.25, 0.25, 0.25, 0.25],
                                        [0.25, 0.25, 0.25, 0.25]))

# row and column normalise
neighborhood_probabilities[0,1:] = (1 - neighborhood_probabilities[0,0]) / len(neighborhood_probabilities[1,1:])
# set row normalised first row to first column 
neighborhood_probabilities[:,0] = neighborhood_probabilities[0,:]
# normalise the rest of the matrix with 1-the first entry
neighborhood_probabilities[1:,1:] = (1 - neighborhood_probabilities[0,1]) / len(neighborhood_probabilities[1,1:])

print(neighborhood_probabilities)
def build_assignment_matrix(attribute_dict, n_cell_types):
    data = list(attribute_dict.items())
    data = np.array(data) # Assignment matrix
    
    B = np.zeros((data.shape[0],n_cell_types)) # Empty matrix
    
    for i in range(0, data.shape[0]):
        t = data[i,1]
        B[i,t] = 1
    
    return B 

position_dict = dict()
for i in range(0, C.shape[0]):
    position_dict[i] = C[i, :]

# create graph from the numpy array of cell locations of the scaffold
graph = nx.from_numpy_array(A)
# assign cell types
cell_assignments = assign_labels.heuristic_assignment(graph, cell_type_probabilities, neighborhood_probabilities, mode='graph', dim=1000, position_dict=position_dict, grid_size=50, revision_iters=100, n_swaps=25)

B = build_assignment_matrix(cell_assignments, n_cell_types=4)

position_dict
df_pos = np.transpose(pd.DataFrame(position_dict))
df_pos.index

df_pos["Labels"] = cell_assignments.values()

df_pos = df_pos.rename(columns={0: 'x', 1: 'y'})

df_pos["sample_id"] = snakemake.wildcards['typ'] + "_" + snakemake.wildcards['sim']
df_pos["interaction_probability"] = snakemake.wildcards['prob']
df_pos["condition"] = snakemake.wildcards['typ']

#write to csv
df_pos.to_csv(str(snakemake.output), header = True, index = True)

