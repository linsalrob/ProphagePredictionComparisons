"""
Phage_Finder

Manuscript: https://academic.oup.com/nar/article/34/20/5839/3100473
Software: http://phage-finder.sourceforge.net/

"""

# CONFIG
outputdir = "phage_finder_tests"
pfBuild = os.path.join(workflow.basedir, "../build/phage_finder")
pfHome = os.path.join(pfBuild, 'phage_finder_v2.1')
pfRun = os.path.join(pfBuild, 'phage_finder_v2.1/bin/Phage_Finder_v2.1.pl')
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


rule convert_gb_to_fna_faa:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        fna = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.fna"),
        faa = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.faa"),
        pfi = os.path.join(outputdir, "{genome}_phage_finder", "phage_finder_info.txt"),
    params:
        dir = os.path.join(workflow.basedir,'../'),
        prefix = os.path.join(outputdir, "{genome}_phage_finder", "{genome}")
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params.dir};
        python3 {scripts}/genbank2sequences.py -g {input.gen} -n {params.prefix} -a {params.prefix} --phage_finder {output.pfi}
        """


rule run_phage_finder:
    input:
        fna = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.fna"),
        faa = os.path.join(outputdir, "{genome}_phage_finder", "{genome}.faa"),
        pfi = os.path.join(outputdir, "{genome}_phage_finder", "phage_finder_info.txt"),
        req = pfRun
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
        # hmm searches
        for i in `cat {pfHome}/hmm3.lst`; do
            hmmsearch {pfHome}/PHAGE_HMM3s_dir/$i.HMM {input.faa};
        done > combined.out
        
        # blast
        blastall -p blastp -d {pfHome}/DB/phage_10_02_07_release.db -m 8 -e 0.001 -i {input.faa} \
            -o ncbi.out -v 4 -b 4 -a 2 -F F
        
        # tRNA scan
        tRNAscan-SE -B -o tRNAscan.out {input.fna}
        
        # aragorn
        aragorn -m -o tmRNA_aragorn.out {input.fna}
        
        # phage_finder
        {pfRun} -t ncbi.out -i {input.pfi} -r tRNAscan.out -n tmRNA_aragorn.out -A {input.fna} -S
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
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -r {input.tbl} > {output.tp}
        """
