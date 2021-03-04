

# CONFIG
outputdir = "phage_finder_tests"
pfBuild = os.path.join(workflow.basedir, "../build/phage_finder")
pfRun = os.path.join(pfBuild, 'phage_finder_v2.1/bin/phage_finder_v2.1.sh')
dlUrl = 'https://cloudstor.aarnet.edu.au/plus/s/LZAWr3htZbZc1uF/download'
dlTar = 'phage_finder_v2.1.tar.gz'

# GENERIC CONFIG/RECIPES
include: "../scripts/preflight.smk"


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, "{genome}_phage_finder", "{genome}_phage_finder_tptn.tsv"), genome=GENOMES)


# RECIPES
rule build_phage_finder:
    output:
        pfRun
    shell:
        """
        cd {pfBuild};
        curl -o {dlTar} {dlUrl};
        tar xvf {dlTar};
        rm {dlTar};
        """


rule run_phage_finder:
    input:
        fna = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.fna"),
        faa = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.faa"),
        pfi = os.path.join(outputdir, "{genome}_phage_finder", "phage_finder_info.txt"),
    output:
        os.path.join(outputdir, "{genome}_phage_finder", "PFPR_tab.txt")
    params:
        os.path.join(outputdir, "{genome}_phage_finder")
    benchmark:
        os.path.join(outputdir, "benchmarks", "{genome}_phage_finder.txt")
    conda:
        "../conda_environments/phage_finder.yaml"
    shell:
        """
        cd {params} && touch error.log formatdb.log && {pfRun} {wildcards.genome}
        """


rule phage_finder_to_tbl:
    input:
        tsv = os.path.join(outputdir, "{genome}_phage_finder", "PFPR_tab.txt")
    output:
        os.path.join(outputdir, "{genome}_phage_finder", "{genome}_phage_finder_locs.tsv")
    shell:
        """
        set +e
        G=$(grep -v ^# {input.tsv})
        exitcode=$?
        if [ $exitcode == 0 ]; then
            IFS=$"\n"
            for X in $(echo $G); do 
                echo $X | cut -f 1,4,5 > {output}
            done
        elif [ $exitcode == 1 ]; then
            touch {output}
        else
            exit $exitcode
        fi
        """


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        tbl = os.path.join(outputdir, "{genome}_phage_finder", "{genome}_phage_finder_locs.tsv")
    output:
        tp = os.path.join(outputdir, "{genome}_phage_finder", "{genome}_phage_finder_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
