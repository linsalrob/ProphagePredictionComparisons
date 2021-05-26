# RECIPES
rule dl_pvogs:
    output:
        os.path.join(psBuild, 'pVOGs.hmm')
    shell:
        """
        cd {psBuild};
        wget http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz;
        tar -zxvf AllvogHMMprofiles.tar.gz;
        cat AllvogHMMprofiles/* > pVOGs.hmm;
        rm -r AllvogHMMprofiles;
        rm AllvogHMMprofiles.tar.gz;
        """


rule run_phispy:
    input:
        g = os.path.join(test_genomes, "{genome}.gb.gz"),
        req = os.path.join(psBuild, 'pVOGs.hmm')
    params:
        o = os.path.join(outputdir, f"{{genome}}.{outDirName}"),
        t = training_set,
        v = use_pvogs
    benchmark:
        os.path.join(outputdir, "benchmarks", f"{{genome}}_{outDirName}.txt")
    output:
        temporary(os.path.join(outputdir, f"{{genome}}.{outDirName}", "bacteria.fasta")),
        temporary(os.path.join(outputdir, f"{{genome}}.{outDirName}", "bacteria.gbk")),
        temporary(os.path.join(outputdir, f"{{genome}}.{outDirName}", "phage.fasta")),
        os.path.join(outputdir, f"{{genome}}.{outDirName}", "phage.gbk"),
        os.path.join(outputdir, f"{{genome}}.{outDirName}", "phispy.log"),
    conda:
        "../conda_environments/phispy.yaml"
    resources:
        mem_mb = 8000
    shell:
        """
        PhiSpy.py --threads 1 -o {params.o} {params.t} {params.v} --output_choice 4 {input.g}
        """


rule count_tp_tn:
    input:
        gen = os.path.join(test_genomes, "{genome}.gb.gz"),
        phg = os.path.join(outputdir, f"{{genome}}.{outDirName}", "phage.gbk")
    output:
        tp = os.path.join(outputdir, f"{{genome}}_{outDirName}_tptn.tsv")
    params:
        os.path.join(workflow.basedir,'../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params};
        python3 scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """