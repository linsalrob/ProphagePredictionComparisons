"""
Virsorter

Manuscript: https://peerj.com/articles/985/
Software: https://github.com/simroux/VirSorter

"""

import os
import sys

# CONFIG
vs1Build = os.path.join(workflow.basedir, "../build/vs1")
if not os.path.exists(vs1Build):
    os.makedirs(vs1Build)
vs1DbUrl = 'https://cloudstor.aarnet.edu.au/plus/s/m55PsF0siDDWI7o/download'
vs1DbTar = 'virsorter-data-v2.tar.gz'
vs1Db = os.path.join(vs1Build, 'virsorter-data')

outDirName = "virsorter"

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_virsorter_tptn.tsv"), genome=GENOMES)


# RECIPES
rule virsorter_db:
    output:
        os.path.join(vs1Db,'VirSorter_Readme.txt')
    shell:
        """
        cd {vs1Build};
        curl -o {vs1DbTar} {vs1DbUrl};
        tar xvf {vs1DbTar};
        rm {vs1DbTar};
        """


rule run_virsorter:
    input:
        fna = os.path.join(outputdir, "{genome}.fna"),
        req = os.path.join(vs1Db,'VirSorter_Readme.txt')
    output:
        c1 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-1.gb"),
        c2 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-2.gb"),
        c3 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-3.gb"),
        c4 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_prophages_cat-4.gb"),
        c5 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_prophages_cat-5.gb"),
        c6 = os.path.join(outputdir,"{genome}_virsorter","Predicted_viral_sequences/VIRSorter_prophages_cat-6.gb")
    params:
        odir = os.path.join(outputdir, "{genome}_virsorter")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_virsorter.txt")
    conda:
        "../conda_environments/virsorter.yaml"
    resources:
        mem_mb = 8000
    shell:
        """
        wrapper_phage_contigs_sorter_iPlant.pl --ncpu 1 -f {input.fna} --db 1 --wdir {params.odir} --data-dir {vs1Db}
        """


rule virsorter_to_tbl:
    input:
        c1 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-1.gb"),
        c2 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-2.gb"),
        c3 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_cat-3.gb"),
        c4 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_prophages_cat-4.gb"),
        c5 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_prophages_cat-5.gb"),
        c6 = os.path.join(outputdir, "{genome}_virsorter", "Predicted_viral_sequences/VIRSorter_prophages_cat-6.gb")
    output:
        os.path.join(outputdir, "{genome}_virsorter", "locs.tsv")
    run:
        import re
        out = open(output[0], 'w')
        for f in [input.c1, input.c2, input.c3]:
            infh = open(f, 'r')
            for line in infh:
                l = line.split()
                if l[0]=='LOCUS':
                    out.write(f'{l[1]}\t1\t{l[2]}\n')
            infh.close()
        for f in [input.c4, input.c5, input.c6]:
            infh = open(f, 'r')
            for line in infh:
                if line.startswith('LOCUS'):
                    l = line.split()
                    l[1] = re.sub('VIRSorter_|-cat_.','',l[1])
                    l[1] = re.sub('_gene.*gene_\d*-|-','\t',l[1])
                    out.write(f'{l[1]}\n')
            infh.close()
        out.close()


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}_virsorter", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_virsorter_tptn.tsv")
    params:
        os.path.join(workflow.basedir, '../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
