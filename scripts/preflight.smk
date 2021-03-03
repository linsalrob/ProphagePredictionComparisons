
# config for all tests
test_genomes = os.path.join(workflow.basedir, "../genbank")
EdwardsLab = os.path.join(workflow.basedir, "../EdwardsLab")
scripts = os.path.join(workflow.basedir, "../scripts")
GENOMES, = glob_wildcards(os.path.join(test_genomes, '{genome}.gb.gz'))



rule convert_gb_to_fna:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        fna = os.path.join(outputdir, "{genome}.fna")
    shell:
        """
        export PYTHONPATH={EdwardsLab};
        python3 {EdwardsLab}/bin/genbank2sequences.py -g {input.gen} -n {output.fna}
        """
