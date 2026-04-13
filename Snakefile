import numpy as np 
import re
import json
np.random.seed(1234)

TYP = ["ctrl", "pert1", "pert2", "pert3"]
COMP = ["pert1", "pert2", "pert3"]
SIM = ["patient1", "patient2", "patient3", "patient4", "patient5"]
FILE = ["centroids", "adjacency"]
FOV = list(range(1,11))
SUBFOV = 1
PROB = {"ctrl": 0.5, "pert1": 0.5, "pert2": 0.4, "pert3": 0.3} 
METHOD = ["spatialFDAG", "spatialFDAL", "spicyRMM", "spicyRLM", "spaceANOVAUni", "spaceANOVAMulti", "smoppix", "intensityMM", "mxfdaFM", "mxfdaMM"]
REP = list(range(1,501))
STD = [1]
PROP = [0.15, 0.25, 0.6]

res_scaffold = list()
res_celltyping = list()
res_dftospe = list()
res_randomfov = list()

output_celltypes = list()
input_merge = list()

for val in TYP:
  input_merge += expand("outs/{rep}/{typ}/{sim}/dataframe-fov-{prob}-{std}-{prop}-{sim}-{rep}-{fov}-{subfov}.csv", typ = val, sim = SIM, fov = FOV, prob = PROB[val], subfov = SUBFOV, allow_missing=True)

container: "docker://condaforge/mambaforge:24.9.2-0"

rule all:
    input: 
            res_pvalues = expand("outs/pValues/pValues-{method}-{std}-{prop}-{comp}.rds", method = METHOD, std = STD, prop = PROP, comp = COMP),
            res_powerCurve = expand("outs/powerCurve.pdf"),
            res_iCOBRA = ["outs/TPRFDP.pdf", "outs/TPRFPR.pdf"],
            res_plotSpicyRSim = "outs/manyRROC.pdf",
            res_marginalIntensities = "outs/marginal_simulated_intensities.pdf",
            res_diabetesExample = ["outs/intensityBoxplot.pdf", "outs/heatmapComb.pdf", "outs/residualPlot.pdf", "outs/qqdeltaTh.pdf"],
            res_runtime = ["outs/runtimes.pdf", "outs/runtimes.rds"]

rule clone_simulation_repo:
    output:
        directory("software/PowerAnalysisForSpatialOmics")
    shell:
        """
        git clone https://github.com/mjemons/PowerAnalysisForSpatialOmics {output}
        """

rule scaffold_generation:
    input: software = "software/PowerAnalysisForSpatialOmics"
    output: "outs/{rep}/{typ}/{sim}/scaffold-{prob}-{sim}-centroids.npy",
            "outs/{rep}/{typ}/{sim}/scaffold-{prob}-{sim}-adjacency.npy"
    conda: "envs/ist_sim.yml"
    script: 
            "code/scaffold.py"

rule sample_variation:
    output: "outs/{rep}/{typ}/{sim}/sample-variation-{std}-{prop}.txt"
    conda: "envs/ist_sim.yml"
    script: 
            "code/sample_variation.py"

rule celltype_assignment:
    input:
        adj =  "outs/{rep}/{typ}/{sim}/scaffold-{prob}-{sim}-adjacency.npy",
        cent = "outs/{rep}/{typ}/{sim}/scaffold-{prob}-{sim}-centroids.npy",
        var = "outs/{rep}/{typ}/{sim}/sample-variation-{std}-{prob}.txt"
    output: "outs/{rep}/{typ}/{sim}/dataframe-celltype-{prob}-{std}-{prop}-{sim}-{rep}-{fov}.csv"
    conda: "envs/ist_sim.yml"
    script: 
            "code/celltyping.py"

rule random_fov_selection:
    input: file =  "outs/{rep}/{typ}/{sim}/dataframe-celltype-{prob}-{std}-{prop}-{sim}-{rep}-{fov}.csv"
    output: file = "outs/{rep}/{typ}/{sim}/dataframe-fov-{prob}-{std}-{prop}-{sim}-{rep}-{fov}-{subfov}.csv"
    conda: "envs/fov_selection.yml"
    script:
            "code/fov_selection.R"

rule marginalIntensitiesPlot:
    input:  ls = expand("outs/{rep}/df_total-{std}-{prop}.rds", rep = REP, std = STD, prop = PROP)
    output: plt = "outs/marginal_simulated_intensities.pdf"
    conda: "envs/diabetesExample.yml"
    script: 
           "code/simulated_marginal_intensities.R"
           
rule df_merge:
    input:  ls = input_merge
    output: rds = "outs/{rep}/df_total-{std}-{prop}.rds"
    conda: "envs/df_merge.yml"
    script:
            "code/df_merge.R"  
            
rule df_to_spe:
    input:  rds = "outs/{rep}/df_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    conda: "envs/df_to_spe.yml"
    script:
            "code/df_to_spe.R"   

rule install_spatialFDA:
    output: touch = "outs/spatialFDA_installed"  
    conda: "envs/spatialFDA.yml"
    script:
            "code/install_spatialFDA.R"
            
rule spatialFDAG:
    input:  software = "outs/spatialFDA_installed",
            rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/spatialFDAG-{std}-{prop}-{comp}.rds"
    conda: "envs/spatialFDA.yml"
    script:
            "code/spatialFDAG.R"   

rule spatialFDAL:
    input:  software = "outs/spatialFDA_installed",
            rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/spatialFDAL-{std}-{prop}-{comp}.rds"
    conda: "envs/spatialFDA.yml"
    script:
            "code/spatialFDAL.R"

rule install_spaceANOVA:
    output: touch = "outs/spaceANOVA_installed"
    conda: "envs/spaceANOVA.yml"
    script:
	    "code/install_spaceANOVA.R"
           
rule spaceANOVA:
    input: softare = "outs/spaceANOVA_installed",
           rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rdsUni = "outs/{rep}/spaceANOVAUni-{std}-{prop}-{comp}.rds",
            rdsMulti = "outs/{rep}/spaceANOVAMulti-{std}-{prop}-{comp}.rds" 
    conda: "envs/spaceANOVA.yml"
    threads: 1
    script:
            "code/spaceANOVA.R"

rule spicyRMM:
    input: rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/spicyRMM-{std}-{prop}-{comp}.rds",
    conda: "envs/spicyR.yml"
    threads: 1
    script:
            "code/spicyRMM.R"

rule spicyRLM:
    input:  rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/spicyRLM-{std}-{prop}-{comp}.rds",
    conda: "envs/spicyR.yml"
    threads: 1
    script:
            "code/spicyRLM.R"

rule install_smoppix:
    output: touch = "outs/smoppix_installed"
    conda: "envs/smoppix.yml"
    script:
            "code/install_smoppix.R"

rule smoppix:
    input:  software = "outs/smoppix_installed",
            rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/smoppix-{std}-{prop}-{comp}.rds",
    conda: "envs/smoppix.yml"
    threads: 1
    script:
            "code/smoppix.R"

rule intensityMM:
    input:  rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/intensityMM-{std}-{prop}-{comp}.rds",
    conda: "envs/intensity.yml"
    threads: 1
    script:
            "code/intensitytMM.R"

rule install_mxfda:
    output: touch = "outs/mxfda_installed"
    conda: "envs/mxfda.yml"
    script:
            "code/install_mxfda.R"

rule mxfdaFM:
    input: software = "outs/mxfda_installed",
           rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/mxfdaFM-{std}-{prop}-{comp}.rds",
    conda: "envs/mxfda.yml"
    threads: 1
    script:
            "code/mxfdaFM.R"

rule mxfdaMM:
    input: software = "outs/mxfda_installed",
           rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/mxfdaMM-{std}-{prop}-{comp}.rds",
    conda: "envs/mxfda.yml"
    threads: 1
    script:
            "code/mxfdaMM.R"

rule comparePvalues:
    input: ls = expand("outs/{rep}/{method}-{std}-{prop}-{comp}.rds", rep = REP, allow_missing = True),
    output: rds = "outs/pValues/pValues-{method}-{std}-{prop}-{comp}.rds"
    conda: "envs/comparePvalues.yml"
    script:
            "code/comparePvalues.R"
               
rule powerCurve:
    input: ls = expand("outs/pValues/pValues-{method}-{std}-{prop}-{comp}.rds", method = METHOD, comp = COMP, std = STD, prop = PROP)
    output: plt = "outs/powerCurve.pdf",
            supplement = "outs/supplementSimRes.pdf",
            combined = "outs/combinedSimRes.pdf"
    conda: "envs/powerCurve.yml"
    script: 
 	    "code/powerCurve.R"

rule iCOBRA:
    input: ls = expand("outs/pValues/pValues-{method}-{std}-{prop}-{comp}.rds", method = METHOD, comp = COMP, std = STD, prop = PROP) 
    output: plt = "outs/TPRFDP.pdf",
            roc = "outs/TPRFPR.pdf"
    conda: "envs/powerCurve.yml"
    script:
            "code/iCOBRA.R"

rule spicyRsim:
    input: software = "outs/spatialFDA_installed"
    output: rds = "outs/manyRsim.rds"
    conda: "envs/spatialFDA.yml"
    threads: 10
    script:
           "code/spicyRsim.R"

rule plotSpicyRSim:
     input: rds = "outs/manyRsim.rds"
     output: plt = "outs/manyRROC.pdf"
     conda: "envs/spatialFDA.yml"
     script:
            "code/plotSpicyRSim.R"

rule runtimeComparison:
     output: rds = "outs/runtimes.rds"
     conda: "envs/runtime.yml"
     script:
            "code/runtimeComparison.R"

rule runtimeComparisonPlot:
     input:  rds = "outs/runtimes.rds"
     output: plt = "outs/runtimes.pdf",
     conda: "envs/runtime.yml"
     script:
            "code/runtimeComparisonPlot.R"

rule diabetesExample:
     output: intensityBoxplot = "outs/intensityBoxplot.pdf",
             rds = "outs/diabetesExample.rds"
     conda: "envs/diabetesExample.yml"
     threads: 5
     script:
            "code/diabetesExample.R"

rule diabetesExamplePlot:
     input: rds = "outs/diabetesExample.rds"
     output: heatmap = "outs/heatmapComb.pdf",
             qcPlot = "outs/residualPlot.pdf",
             qcPlotDeltaTh = "outs/qqdeltaTh.pdf"
     conda: "envs/diabetesExample.yml"
     script:
            "code/diabetesExamplePlot.R"
