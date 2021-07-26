""""
Seeker
Manuscript: https://academic.oup.com/nar/article/48/21/e121/5921300
Software: https://github.com/gussow/seeker
"""

import os
import sys

#CONFIG
outDirName = "seeker" #tool dir in global output directory
seekerBuild = os.path.join(workflow.basedir, "../build/seeker")

#GENERIC
include: os.path.join(workflow.basedir, "../rules/preflight.smk")

#TARGET
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_seeker_tptn.tsv"), genome=GENOMES)

#BUILD
rule build_seeker:
    output:
        os.path.join(seekerBuild, 'setup.done')
    conda:
        "../conda_environments/seeker.yaml"
    shell:
        """
        predict-metagenome
        """

#RECIPES
rule run_seeker:
    input:
        fa = os.path.join(outputdir, "{genome}.fna"),
    output:
        bed = os.path.join(outputdir, "{genome}.bed")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_seeker.txt")
    conda:
        "../conda_environments/seeker.yaml"
    params:
        os.path.join(workflow.basedir,'../'),
    resources:
        mem_mb = 8000
    shell:
        """
        export PYTHONPATH={params}
        python {scripts}/run_seeker.py -i {input.fa} -o {output.bed}
        """

rule seeker_2_tbl:
    input:
        os.path.join(outputdir, "{genome}.bed")
    output:
        os.path.join(outputdir, "{genome}", "_seeker_locs.tsv")
    run:
        infh = open(input[0],'r')
        outfh = open(output[0], 'w')
        for line in infh:
            if not line.startswith("#"):
                data = line.split("\t")
                outfh.write(f"{data[0].strip()}\t{data[1].strip()}\t{data[2].strip()}\n")
        outfh.close()

rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}", "_seeker_locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_seeker_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """