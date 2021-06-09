"""
Virsorter2

Manuscript: https://doi.org/10.1186/s40168-020-00990-y
Software: https://github.com/jiarong/VirSorter2

"""

import os
import sys


# CONFIG
vs2Build = os.path.join(workflow.basedir, "../build/vs2")
if not os.path.exists(vs2Build):
    os.makedirs(vs2Build)
vs2DbUrl = 'https://osf.io/v46sc/download'
outDirName = 'virsorter2'

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_virsorter2_tptn.tsv"), genome=GENOMES)


# RECIPES
rule virsorter_setup:
    output:
        os.path.join(vs2Build, 'setup.done')
    conda:
        "../conda_environments/virsorter2.yaml"
    shell:
        """
        cd {vs2Build};
        wget -O db.tgz {vs2DbUrl};
        tar xvf db.tgz;
        virsorter config --init-source --db-dir=./db;
        touch {output}
        """


rule run_virsorter2:
    input:
        fna = os.path.join(outputdir, "{genome}.fna"),
        req = os.path.join(vs2Build, 'setup.done')
    output:
        os.path.join(outputdir, "{genome}.out", 'final-viral-boundary.tsv'),
    params:
        os.path.join(outputdir, "{genome}.out")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_virsorter2.txt")
    conda:
        "../conda_environments/virsorter2.yaml"
    resources:
        mem_mb = 16000
    shell:
        """
        virsorter run --use-conda-off --db-dir {vs2Build}/db -w {params} -i {input.fna} -j 1 all
        """


rule virsorter2_to_tbl:
    input:
        os.path.join(outputdir, "{genome}.out", 'final-viral-boundary.tsv'),
    output:
        temp(os.path.join(outputdir, "{genome}.out", "locs.tsv.tmp"))
    run:
        outFH = open(output[0], 'w')
        inFH = open(input[0], 'r')
        for line in inFH:
            if not line.startswith('seqname'):
                l = line.split('\t')
                outFH.write(f'{l[0]}\t{l[16]}\t{l[17]}\n')
        outFH.close()
        inFH.close()


rule virsorter2_tbl_merge:
    input:
        os.path.join(outputdir,"{genome}.out","locs.tsv.tmp")
    output:
        os.path.join(outputdir,"{genome}.out","locs.tsv")
    shell:
        "bedtools merge -i {input} > {output}"


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}.out", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_virsorter2_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
