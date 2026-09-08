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
METHOD = ["spatialFDAG", "spatialFDAGAdj", "spatialFDAGNSW",
 "spatialFDAL", "spatialFDALAdj", "spatialFDALNSW", 
 "spicyRMM", "spicyRLM", "spaceANOVAUni", "spaceANOVAMulti", "smoppix", 
 "intensityMM", "mxfdaFM", "mxfdaMM"]
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

container: "docker://continuumio/miniconda3:26.5.3"

rule all:
    input: 
            res_pvalues = expand("outs/pValues/pValues-{method}-{std}-{prop}-{comp}.rds", method = METHOD, std = STD, prop = PROP, comp = COMP),
            res_powerCurve = expand("outs/powerCurve.pdf"),
            res_iCOBRA = ["outs/TPRFDP.pdf", "outs/TPRFPR.pdf",
            "outs/TPRFDPSupp.pdf", "outs/TPRFPRSupp.pdf"],
            res_plotSpicyRSim = "outs/manyRROC.pdf",
            res_marginalIntensities = "outs/marginal_simulated_intensities.pdf",
            res_diabetesExample = ["outs/intensityBoxplot.pdf", "outs/heatmapComb.pdf", 
            "outs/residualPlot.pdf", "outs/qqdeltaTh.pdf", "outs/heatmapCombSuppA.pdf",
            "outs/heatmapCombSuppB.pdf", "outs/exampleFOVs.pdf"],
            res_runtime = ["outs/runtimes.pdf", "outs/runtimes.rds"],
            res_plotSimData = ["outs/plotSimData.pdf", "outs/plotCurves.pdf"],
            res_diabetesExampleLOO = ["outs/diabetesExampleLOO.rds", 
            "outs/diabetesExampleLOOPlot.pdf"],
            res_cordsExample = ["outs/plotCords.pdf", "outs/fbplotTF.pdf", 
            "outs/fbplotVF.pdf",]

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

rule plotSimData:
    input: rds = "outs/15/spe_total-1-0.25.rds"
    output: pltSims = "outs/plotSimData.pdf",
            pltCurves = "outs/plotCurves.pdf"
    conda: "envs/spatialFDA.yml"
    script:
            "code/plotSimData.R"
            
rule spatialFDAG:
    input:  rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rdsG = "outs/{rep}/spatialFDAG-{std}-{prop}-{comp}.rds",
            rdsGNSW = "outs/{rep}/spatialFDAGNSW-{std}-{prop}-{comp}.rds",
            rdsGAdj = "outs/{rep}/spatialFDAGAdj-{std}-{prop}-{comp}.rds"
    conda: "envs/spatialFDA.yml"
    script:
            "code/spatialFDAG.R"   

rule spatialFDAL:
    input:  rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rdsL = "outs/{rep}/spatialFDAL-{std}-{prop}-{comp}.rds",
            rdsLNSW = "outs/{rep}/spatialFDALNSW-{std}-{prop}-{comp}.rds",
            rdsLAdj = "outs/{rep}/spatialFDALAdj-{std}-{prop}-{comp}.rds"
    conda: "envs/spatialFDA.yml"
    script:
            "code/spatialFDAL.R"
           
rule spaceANOVA:
    input: rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
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

rule smoppix:
    input:  rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
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

rule mxfdaFM:
    input: rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
    output: rds = "outs/{rep}/mxfdaFM-{std}-{prop}-{comp}.rds",
    conda: "envs/mxfda.yml"
    threads: 1
    script:
            "code/mxfdaFM.R"

rule mxfdaMM:
    input: rds = "outs/{rep}/spe_total-{std}-{prop}.rds"
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
            roc = "outs/TPRFPR.pdf",
            pltSupp = "outs/TPRFDPSupp.pdf",
            rocSupp = "outs/TPRFPRSupp.pdf"
    conda: "envs/powerCurve.yml"
    script:
            "code/iCOBRA.R"

rule spicyRsim:
    output: rds = "outs/manyRsim.rds"
    conda: "envs/spatialFDA.yml"
    threads: 10
    script:
           "code/spicyRsim.R"

rule plotSpicyRSim:
     input: rds = "outs/manyRsim.rds"
     output: roc = "outs/manyRROC.pdf",
             fdp = "outs/manyRFDP.pdf"
     conda: "envs/powerCurve.yml"
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

rule diabetesExampleFOVs:
     output: plt = "outs/exampleFOVs.pdf"
     conda: "envs/diabetesExample.yml"
     script:
            "code/diabetesExampleFOVs.R"

rule diabetesExampleLOO:
     output: rds = "outs/diabetesExampleLOO.rds"
     conda: "envs/diabetesExample.yml"
     threads: 15
     script:
            "code/diabetesExampleLOO.R"

rule diabetesExampleLOOPlot:
     input: rds = "outs/diabetesExampleLOO.rds"
     output: plt = "outs/diabetesExampleLOOPlot.pdf"
     conda: "envs/diabetesExample.yml"
     script:
            "code/diabetesExampleLOOPlot.R"

rule diabetesExamplePlot:
     input: rds = "outs/diabetesExample.rds"
     output: heatmap = "outs/heatmapComb.pdf",
             heatmapSuppA = "outs/heatmapCombSuppA.pdf",
             heatmapSuppB = "outs/heatmapCombSuppB.pdf",
             qcPlot = "outs/residualPlot.pdf",
             qcPlotDeltaTh = "outs/qqdeltaTh.pdf"
     conda: "envs/diabetesExample.yml"
     script:
            "code/diabetesExamplePlot.R"

rule downloadCordsData:
    output:
        rds="data/Cords/SingleCellExperiment Objects/sce_all_annotated.rds"
    params:
        zipfile="data/raw/cords2024/SingleCellExperiment_Objects.zip"
    shell:
        """
        mkdir -p "$(dirname "{params.zipfile}")"
        wget -c 'https://zenodo.org/records/7961844/files/SingleCellExperiment%20Objects.zip?download=1' \
            -O "{params.zipfile}"
        mkdir -p "data/Cords"

        python -m zipfile -e \
            "{params.zipfile}" \
            "data/Cords"
        """

rule CordsExample:
     input: rds = "data/Cords/SingleCellExperiment Objects/sce_all_annotated.rds"
     output: rds = "outs/cords_results.rds"
     conda: "envs/diabetesExample.yml"
     threads: 10
     script:
            "code/cordsExample.R"

rule CordsExamplePlot:
     input: rds = "outs/cords_results.rds",
            sce = "data/Cords/SingleCellExperiment Objects/sce_all_annotated.rds"
     output: pCords = "outs/plotCords.pdf",
             pCordsSupp = "outs/pCordsSupp.pdf",
             fbplotTF = "outs/fbplotTF.pdf",
             fbplotVF = "outs/fbplotVF.pdf",
     conda: "envs/diabetesExample.yml"
     script:
            "code/cordsExamplePlot.R"