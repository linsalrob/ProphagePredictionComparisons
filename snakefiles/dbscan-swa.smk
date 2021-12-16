"""
DBSCAN-SWA

Manuscript: https://www.biorxiv.org/content/10.1101/2020.07.12.199018v1
Software: https://github.com/HIT-ImmunologyLab/DBSCAN-SWA

"""

import os
import sys

# CONFIG
outDirName = 'dbscan-swa'
dbsBuild = os.path.join(workflow.basedir, "../build/")
dbsHome = os.path.join(dbsBuild, 'DBSCAN')
dbsRun = os.path.join(dbsHome, 'bin/dbscan-swa.py')
dlUrl = 'https://github.com/beardymcjohnface/DBSCAN-SWA-1.git'



# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_dbscan-swa_tptn.tsv"), genome=GENOMES)


rule build_dbscan_swa:
    """
    clone git, make executable, run test to do the first-time database installation, delete test
    """
    output:
        dbsRun
    conda:
        "../conda_environments/dbscan-swa.yaml"
    shell:
        """
        cd {dbsBuild}
        git clone {dlUrl}
        cd DBSCAN-SWA-1/
        cd {dbsBuild}
        rm -rf DBSCAN/ && mkdir DBSCAN
        mv DBSCAN-SWA-1/* DBSCAN/
        rm -rf DBSCAN-SWA-1/
        cd DBSCAN/
        wget -O db.tar.gz https://cloudstor.aarnet.edu.au/plus/s/2nfXwDm9YuWNy5C/download
        tar xvf db.tar.gz
        cd {dbsHome}
        chmod u+x -R bin/
        chmod u+x -R software/
        cd {dbsHome}/test
        # sed 's/num_threads = 20/num_threads = thread_num/' {dbsRun} \
        #     | sed 's/thread_num <- 10//' \
        #     | sed "s/.*0 prophage region was detected in the query bacterial genome.*/\\t\\topen(save_prophage_summary_file,'w').close()\\n\\t\\tsys.exit(0)/" \
        #     > tmp.py && mv tmp.py {dbsRun}
        python {dbsRun} --input NC_007054.fasta --output yeet --thread_num 1
        rm -r yeet
        """


rule run_dbscan_swa:
    input:
        fa = os.path.join(outputdir,"{genome}.fna"),
        req = dbsRun
    output:
        os.path.join(outputdir, '{genome}/bac_DBSCAN-SWA_prophage_summary.txt')
    conda:
        "../conda_environments/dbscan-swa.yaml"
    params:
        os.path.join(outputdir, '{genome}')
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_dbscan-swa.txt")
    log:
        os.path.join(outputdir, '{genome}/dbscan.log')
    resources:
        mem_mb = 8000
    shell:
        "python {dbsRun} --input {input.fa} --output {params} --thread_num 1 &> {log}"


rule dbscan_swa_2_tbl:
    input:
        os.path.join(outputdir,'{genome}/bac_DBSCAN-SWA_prophage_summary.txt')
    output:
        os.path.join(outputdir,"{genome}","locs.tsv")
    run:
        infh = open(input[0],'r')
        outfh = open(output[0], 'w')
        line = infh.readline()
        for line in infh:
            l = line.split('\t')
            id = l[1].split()
            outfh.write(f'{id[0]}\t{l[3]}\t{l[4]}\n')
        outfh.close()


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_dbscan-swa_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
