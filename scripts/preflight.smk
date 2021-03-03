
# config for all tests
test_genomes = os.path.join(workflow.basedir, "../genbank")
EdwardsLab = os.path.join(workflow.basedir, "../EdwardsLab")
scripts = os.path.join(workflow.basedir, "../scripts")
GENOMES, = glob_wildcards(os.path.join(test_genomes, '{genome}.gb.gz'))

# Virsorter2
virsorter2Build = os.path.join(workflow.basedir, "../build/vs2")
vs2DbUrl = 'https://osf.io/v46sc/download'



