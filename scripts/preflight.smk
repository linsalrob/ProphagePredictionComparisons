
# config for all tests
test_genomes = os.path.join(workflow.basedir, "../genbank")
scripts = os.path.join(workflow.basedir, "../scripts")
GENOMES, = glob_wildcards(os.path.join(test_genomes, '{genome}.gb.gz'))



rule convert_gb_to_fna:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        fna = os.path.join(outputdir, "{genome}.fna")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 {scripts}/genbank2sequences.py -g {input.gen} -n {output.fna}
        """