
# CONFIG
outDirName = ''
TOOLS = {
    'dbscan-swa' : 'DBSCAN-SWA',
    'phage_finder' : 'Phage Finder',
    'phageboost' : 'PhageBoost',
    'phigaro' : 'Phigaro',
    'phispy' : 'PhiSpy',
    # 'phispy_trained' : 'PhiSpy (trained)',
    'phispy_pvog' : 'PhiSpy (+pVOGs)',
    'virsorter' : 'VirSorter',
    'virsorter2' : 'VirSorter2',
    'vibrant' : 'VIBRANT'
}


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


rule all:
    input:
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_benchmarks.tsv"),
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_tptn.tsv"),
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_dirSize.tsv")


rule out_dir_sizes:
    input:
        expand(os.path.join(testDir, '{tool}/{genome}_{tool}_tptn.tsv'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_dirSize.tsv")
    run:
        import subprocess
        out = open(output[0], 'w')
        for tool in TOOLS:
            dirSize = subprocess.check_output(['du', '-sb', os.path.join(testDir,tool)])
            l = dirSize.decode().split()
            out.write(f'{TOOLS[tool]}\t{l[0]}\n')
        out.close()


rule combine_tptn:
    input:
        expand(os.path.join(testDir, '{tool}/{genome}_{tool}_tptn.tsv'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_tptn.tsv")
    run:
        out = open(output[0],'w')
        out.write(f'Prophage_caller\tGenome\tTP\tTN\tFP\tFN\tAccuracy\tPrecision\tRecall\tSpecificity\tf1_score\n')
        for tool in TOOLS:
            for genome in GENOMES:
                s = {}
                tptn = open(os.path.join(testDir, f'{tool}/{genome}_{tool}_tptn.tsv'), 'r')
                for line in tptn:
                    if not line.startswith('#'):
                        l = line.split()
                        s[l[0]] = l[1]
                tptn.close()
                out.write('\t'.join( [ TOOLS[tool], genome, s["TP:"], s["TN:"], s["FP:"], s["FN:"], s["Accuracy:"],
                    s["Precision:"], s["Recall:"], s["Specificity:"], s["f1_score:"] + "\n" ] ))
        out.close()


rule combine_benchmarks:
    input:
        expand(os.path.join(testDir, '{tool}/benchmarks/{genome}_{tool}.txt'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_benchmarks.tsv")
    run:
        out = open(output[0], 'w')
        print_head = True
        for tool in TOOLS:
            for genome in GENOMES:
                bench = open(os.path.join(testDir, f'{tool}/benchmarks/{genome}_{tool}.txt'), 'r')
                h = bench.readline().split()
                b = bench.readline().split()
                bench.close()
                if print_head:
                    out.write('\t'.join(['Prophage_caller','Genome'] + h ))
                    out.write("\n")
                    print_head = False
                out.write('\t'.join([TOOLS[tool], genome] + b ))
                out.write("\n")
        out.close()
