""""
ProphET
Manuscript: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0223364
Software: https://github.com/jaumlrc/ProphET
"""

import os
import sys

#CONFIG
outDirName = 'prophet' #tool dir in global output directory
prophetBuild = os.path.join(workflow.basedir, "../build/")
prophetHome = os.path.join(prophetBuild, 'ProphET')
prophetRun = os.path.join(prophetHome, 'ProphET_standalone.pl')
prophetGit = 'https://github.com/jaumlrc/ProphET.git'

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/prophET_rules.smk")

#TARGET
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_prophet_tptn.tsv"), genome=GENOMES)

#BUILD
rule build_prophet:
    output:
        prophetRun
    conda:
        "../conda_environments/prophet.yaml"
    shell:
        """
        cd {prophetBuild}
        git clone {prophetGit}
        cd {prophetHome}
        ./INSTALL.pl
        """

#RECIPES
rule run_prophet:
    input:
        fa = os.path.join(outputdir, "{genome}.prophET.fna"),
        gff = os.path.join(outputdir, "{genome}.prophET.gff"),
        req = prophetRun
    output:
        os.path.join(outputdir, '{genome}', 'phages_coords')
    conda:
        "../conda_environments/prophet.yaml"
    params:
        dir = directory(os.path.join(outputdir, '{genome}'))
    benchmark:
        os.path.join(outputdir,"benchmarks","{genome}_prophet.txt")
    resources:
        mem_mb = 8000
    shell:
        """
        perl {prophetRun} --fasta_in {input.fa} --gff_in {input.gff} --outdir {params.dir}
        """

rule prophet_2_tbl:
    input:
        os.path.join(outputdir, '{genome}/phages_coords')
    output:
        os.path.join(outputdir, "{genome}", "locs.tsv")
    run:
        infh = open(input[0],'r')
        outfh = open(output[0], 'w')
        for line in infh:
            if not line.startswith("#"):
                data = line.split("\t")
                outfh.write(f"{data[0].strip().replace('_dot_', '.')}\t{data[2].strip()}\t{data[3].strip()}\n") # reinsert the dots that were remove to appease ptrophET
        outfh.close()

rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_prophet_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """