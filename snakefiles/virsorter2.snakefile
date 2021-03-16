
# CONFIG
vs2Build = os.path.join(workflow.basedir, "../build/vs2")
if not os.path.exists(vs2Build):
    os.makedirs(vs2Build)
vs2DbUrl = 'https://osf.io/v46sc/download'
outputdir = os.path.join(workflow.basedir, "../virsorter2_tests")

# GENERIC CONFIG/RECIPES
include: "../scripts/preflight.smk"


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_virsorter2_tptn.tsv"), genome=GENOMES)


# RECIPES
rule virsorter_setup:
    output:
        os.path.join(vs2Build, 'setup.done')
    conda:
        "../conda_environments/virsorter2.yaml"
    shell:
        """
        cd {vs2Build};
        wget -O db.tgz {vs2DbUrl};
        tar xvf db.tgz;
        virsorter config --init-source --db-dir=./db;
        touch {output}
        """


rule run_virsorter2:
    input:
        fna = os.path.join(outputdir, "{genome}.fna"),
        req = os.path.join(vs2Build, 'setup.done')
    output:
        os.path.join(outputdir, "{genome}.out", 'final-viral-boundary.tsv'),
    params:
        os.path.join(outputdir, "{genome}.out")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_virsorter2.txt")
    conda:
        "../conda_environments/virsorter2.yaml"
    resources:
        mem_mb = 16000,
        cpus = 4,
        time_min = 480
    shell:
        """
        virsorter run --use-conda-off  -w {params} -i {input.fna} -j 4 all
        """


rule virsorter2_to_tbl:
    input:
        os.path.join(outputdir, "{genome}.out", 'final-viral-boundary.tsv'),
    output:
        os.path.join(outputdir, "{genome}.out", "locs.tsv")
    run:
        outFH = open(output[0], 'w')
        inFH = open(input[0], 'r')
        for line in inFH:
            l = line.split()
            if l[13]=='0':
                outFH.write(f'{l[0]}\t{l[16]}\t{l[17]}\n')
        outFH.close()
        inFH.close()


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}.out", "locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_virsorter2_tptn.tsv")
    shell:
        """
        python3 {scripts}/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
