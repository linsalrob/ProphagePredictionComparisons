

outputdir = "phispy_tptn"


include: "../scripts/preflight.smk"


rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phispy_tptn.tsv"), genome=GENOMES)

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
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
