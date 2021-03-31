
# config for all tests
test_genomes = os.path.join(workflow.basedir, "../genbank")
scripts = os.path.join(workflow.basedir, "../scripts")
GENOMES, = glob_wildcards(os.path.join(test_genomes, '{genome}.gb.gz'))

testDir = os.path.join(workflow.basedir, '../tests')
outputdir = os.path.join(testDir, outDirName)
if not os.path.exists(testDir):
    os.mkdir(testDir)


rule convert_gb_to_fna:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        fna = os.path.join(outputdir, "{genome}.fna")
    params:
        env = os.path.join(workflow.basedir,'../'),
        out = os.path.join(outputdir, "{genome}")
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params.env};
        python3 {scripts}/genbank2sequences.py -g {input.gen} -n {params.out}
        """
