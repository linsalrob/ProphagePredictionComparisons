"""
Phigaro

Manuscript: https://academic.oup.com/bioinformatics/article/36/12/3882/5822875
Software: https://github.com/bobeobibo/phigaro

"""

import os
import sys

# CONFIG
outDirName = "phigaro"
pgBuild = os.path.join(workflow.basedir, "../build/phigaro")

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phigaro_tptn.tsv"), genome=GENOMES)


# BUILD PHIGARO
rule phigaro_setup:
    """
    phigaro needs a one-time setup to download the databases.
    """
    output:
        os.path.join(pgBuild, 'setup.done')
    conda:
        "../conda_environments/phigaro.yaml"
    shell:
        "printf '1\n1\nN\n' | phigaro-setup --no-updatedb && touch {output}"


# RECIPES
rule run_phigaro:
    input:
        fna = os.path.join(outputdir, "{genome}.fna"),
        req = os.path.join(pgBuild, 'setup.done')
    output:
        tsv = os.path.join(outputdir, "{genome}_phigaro", "{genome}.phigaro.tsv")
    params:
        tsv = "{genome}_phigaro" # phigaro adds the .tsv extension
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_phigaro.txt")
    conda:
        "../conda_environments/phigaro.yaml"
    shell:
        """
        cd {outputdir};
        phigaro -f {input.fna} -e tsv -o {params.tsv} --delete-shorts
        """


rule phigaro_to_tbl:
    input:
        tsv = os.path.join(outputdir, "{genome}_phigaro", "{genome}.phigaro.tsv")
    output:
        os.path.join(outputdir, "{genome}_phigaro_locs.tsv")
    shell:
        """
        if [ $(stat -c %s {input}) -lt 50 ]; then
            touch {output}
        else
            grep -v scaffold {input.tsv} | cut -f 1,2,3 > {output}
        fi
        """


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}_phigaro_locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_phigaro_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
