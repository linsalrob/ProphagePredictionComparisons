
# CONFIG
TOOLS = ['phage_finder','phageboost','phigaro','phispy','virsorter','virsorter2']
outDirName = ''

# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../scripts/preflight.smk")


rule all:
    input:
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_benchmarks.tsv"),
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_tptn.tsv")


rule combine_tptn:
    input:
        expand(os.path.join(testDir, '{tool}/{genome}_{tool}_tptn.tsv'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_tptn.tsv")
    params:
        script = os.path.join(scripts, 'parse_prophage_predictions.pl'),
        test = testDir
    run:
        out = open(output[0],'w')
        out.write(f'Prophage Caller\tGenome\tTP\tTN\tFP\tFN\tAccuracy\tPrecision\tRecall\tSpecificity\tf1 score\n')
        for tool in TOOLS:
            pipeJob = f'perl {params.script} -d {params.test}/{tool} -c {tool}'
            for line in shell(pipeJob,iterable=True):
                out.write(line)
        out.close()


rule combine_benchmarks:
    input:
        expand(os.path.join(testDir, '{tool}/benchmarks/{genome}.benchmarks.txt'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_benchmarks.tsv")
    params:
        script = os.path.join(scripts, 'parse_benchmarks.pl'),
        test = testDir
    run:
        out = open(output[0], 'w')
        out.write(f'Prophage\tCaller Genome\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time')
        for tool in TOOLS:
            pipeJob = f'perl {params.script} -d {params.test}/{tool}/benchmarks -c {tool}'
            for line in shell(pipeJob,iterable=True):
                out.write(line)
        out.close()

