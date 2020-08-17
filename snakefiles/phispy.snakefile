

outputdir = "phispy_tptn"

phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))

rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phispy_tptn.tsv"), genome=GENOMES)

rule run_phispy:
    input:
        g = os.path.join(phispydir, "{genome}.gb.gz")
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
    shell:
        """
        PhiSpy.py -o {params.o} --output_choice 4 {input.g}
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        phg = os.path.join(outputdir, "{genome}.phispy", "phage.gbk")
    output:
        tp = os.path.join(outputdir, "{genome}_phispy_tptn.tsv")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
