"""
Phispy (v2)

Manuscript: (v1) https://academic.oup.com/nar/article/40/16/e126/1027055
Software: (v2) https://github.com/linsalrob/PhiSpy

"""

# CONFIG
outputdir = "phispy_tptn"

# GENERIC CONFIG/RECIPES
include: "../scripts/preflight.smk"


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phispy_tptn.tsv"), genome=GENOMES)


# RECIPES
rule run_phispy:
    input:
        g = os.path.join(test_genomes, "{genome}.gb.gz")
    params:
        o = os.path.join(outputdir, "{genome}.phispy")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}.benchmarks.txt")
    output:
        temporary(os.path.join(outputdir, "{genome}.phispy", "bacteria.fasta")),
        temporary(os.path.join(outputdir, "{genome}.phispy", "bacteria.gbk")),
        temporary(os.path.join(outputdir, "{genome}.phispy", "phage.fasta")),
        os.path.join(outputdir, "{genome}.phispy", "phage.gbk"),
        os.path.join(outputdir, "{genome}.phispy", "phispy.log"),
    conda:
        "../conda_environments/phispy.yaml"
    resources:
        mem_mb = 8000
    shell:
        """
        PhiSpy.py -o {params.o} --output_choice 4 {input.g}
        """


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        phg = os.path.join(outputdir, "{genome}.phispy", "phage.gbk")
    output:
        tp = os.path.join(outputdir, "{genome}_phispy_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
