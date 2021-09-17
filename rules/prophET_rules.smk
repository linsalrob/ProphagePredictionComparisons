# config for all tests
test_genomes = os.path.join(workflow.basedir, "../genbank")
scripts = os.path.join(workflow.basedir, "../scripts")
GENOMES, = glob_wildcards(os.path.join(test_genomes, '{genome}.gb.gz'))

testDir = os.path.join(workflow.basedir, '../tests')
outputdir = os.path.join(testDir, outDirName)
if not os.path.exists(testDir):
    os.mkdir(testDir)


#INPUT PREPROCESS
rule unpack_gbk:
    input:
        gbk_gz = os.path.join(test_genomes, "{genome}.gb.gz")
    output:
        tmp_gbk = temp(os.path.join(outputdir, "{genome}.prophET.tmp.gbk"))
    shell:
        """
        gzip -dkc {input.gbk_gz} > {output.tmp_gbk}
        """

rule sanitize_accesions:  # if accessions include dots they must be replaced (_dot_) to avoid problem during prophET run
    input:
        tmp_gbk = os.path.join(outputdir, "{genome}.prophET.tmp.gbk")
    output:
        fixed_tmp_gbk = temp(os.path.join(outputdir, "{genome}.prophET.fixed.tmp.gbk"))
    params:
        env = os.path.join(workflow.basedir, '../')
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        export PYTHONPATH={params.env};
        python3 {scripts}/sanitize_genbank.py --in {input.tmp_gbk} --out {output.fixed_tmp_gbk}
        """

rule convert_gb_to_fna:
    input:
        tmp_gbk = os.path.join(outputdir, "{genome}.prophET.fixed.tmp.gbk")
    output:
        fna = os.path.join(outputdir, "{genome}.prophET.fna")
    params:
        env = os.path.join(workflow.basedir, '../'),
        out = os.path.join(outputdir, "{genome}.prophET")
    conda:
        "../conda_environments/roblib.yaml"
    shell:
        """
        echo {workflow.basedir}
        export PYTHONPATH={params.env};
        python3 {scripts}/genbank2sequences.py -g {input.tmp_gbk} -n {params.out}
        """

rule get_prophet_gff:
    input:
        tmp_gbk = os.path.join(outputdir, "{genome}.prophET.fixed.tmp.gbk")
    output:
        tmp_gff = temp(os.path.join(outputdir,"{genome}.prophET.tmp.gff")),
        gff = os.path.join(outputdir, "{genome}.prophET.gff"),
    params:
        perl5lib = os.path.join(scripts,'GFFLib'),
        standardisation_script = os.path.join(scripts,'GFFLib/gff_rewrite.pl'),
        logfile = os.path.join(outputdir,"gbk2gff.conversion.log")
    conda:
        "../conda_environments/prophet.yaml"
    shell:
        """
        export PERL5LIB={params.perl5lib}
        echo {input.tmp_gbk} 2>> {params.logfile}
        bp_genbank2gff3.pl {input.tmp_gbk} --outdir stdout > {output.tmp_gff} 2>> {params.logfile}
        perl {params.standardisation_script} --input {output.tmp_gff} --output {output.gff} --add_missing_features 2>> {params.logfile}
        """