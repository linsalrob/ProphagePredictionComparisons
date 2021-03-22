"""
PhageBoost

Manuscript: https://www.biorxiv.org/content/10.1101/2020.08.09.243022v1.full.pdf
Software: https://github.com/ku-cbd/PhageBoost

Notes:
    PhageBoost can be invoked from the command line (PhageBoost -h), however
    by default it accepts a fasta format file and performs its own gene predictions.
    The phageboost_genbank.py script is required to use the annotations provided by
    the genbank files.
"""

# CONFIG
outDirName = "phageboost"
dataDir = os.path.join(workflow.basedir, "../data")

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../scripts/preflight.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phageboost_tptn.tsv"), genome=GENOMES)


# RECIPES
rule run_phageboost:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        tsv = os.path.join(outputdir, "{genome}_phageboost.tsv")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_phageboost.txt")
    conda:
        "../conda_environments/phageboost.yaml"
    params:
        pypath = os.path.join(workflow.basedir,'../'),
        datdir = dataDir
    shell:
        """
        export PYTHONPATH={params.pypath};
        python3 scripts/phageboost_genbank.py -g {input.gen} -o {output.tsv} \
            -m {params.datdir}/model_delta_std_hacked.pickled.silent.gz
        """


rule phageboost_to_tbl:
    input:
        tsv = os.path.join(outputdir, "{genome}_phageboost.tsv")
    output:
        os.path.join(outputdir, "{genome}_phageboost_locs.tsv")
    shell:
        """
        if [ $(stat -c %s {input}) -lt 50 ]; then
            touch {output}
        else
            grep -v probability {input.tsv} | cut -f 3,4,5 > {output}
        fi
        """


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}_phageboost_locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_phageboost_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
