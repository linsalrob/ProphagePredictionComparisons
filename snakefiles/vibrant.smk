"""
Vibrant

Manuscript: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0
Software: https://github.com/AnantharamanLab/VIBRANT

"""

import os
import sys


# CONFIG
outDirName = 'vibrant'
vBuild = os.path.join(workflow.basedir, "../build/vibrant/")


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")



# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, f"{{genome}}_{outDirName}_tptn.tsv"), genome=GENOMES)


# RECIPES
rule vibrant_db_download:
    output:
        os.path.join(vBuild, 'setup.done')
    conda:
        "../conda_environments/vibrant.yaml"
    shell:
        """
        download-db.sh;
        touch {output}
        """


rule run_vibrant:
    input:
        f = os.path.join(outputdir,"{genome}.fna"),
        req = os.path.join(vBuild, 'setup.done')
    output:
        os.path.join(outputdir, '{genome}/VIBRANT_{genome}/VIBRANT_results_{genome}/VIBRANT_integrated_prophage_coordinates_{genome}.tsv')
    conda:
        "../conda_environments/vibrant.yaml"
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_vibrant.txt")
    params:
        os.path.join(outputdir, '{genome}')
    shell:
        """
        mkdir -p {params};
        cd {params};
        VIBRANT_run.py -i {input.f} -no_plot -t 1;
        """


rule vibrant_to_tbl:
    input:
        os.path.join(outputdir, '{genome}/VIBRANT_{genome}/VIBRANT_results_{genome}/VIBRANT_integrated_prophage_coordinates_{genome}.tsv')
    output:
        os.path.join(outputdir, "{genome}.out", "locs.tsv")
    run:
        outFH = open(output[0], 'w')
        inFH = open(input[0], 'r')
        for line in inFH:
            l = line.split('\t')
            if l[0] != 'scaffold':
                outFH.write(f'{l[0]}\t{l[5]}\t{l[6]}\n')
        outFH.close()
        inFH.close()


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}.out", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_vibrant_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """