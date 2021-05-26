
# CONFIG
outDirName = ''
TOOLS = {
    'dbscan-swa' : 'DBSCAN-SWA',
    'phage_finder' : 'Phage Finder',
    'phageboost' : 'PhageBoost',
    'phigaro' : 'Phigaro',
    'phispy' : 'PhiSpy',
    'phispy_trained' : 'PhiSpy (trained)',
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
            dirSize = subprocess.check_output(['du', '-sb', os.path.join(testDir, tool)], text=True)
            out.write(dirSize)
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
        bmk = {}
        benchmarks = ['User time (seconds)', 'Maximum resident set size (kbytes)', 'File system outputs']
        for tool in TOOLS:
            try:
                bmk[tool]
            except KeyError:
                bmk[tool] = {}
            for genome in GENOMES:
                try:
                    bmk[tool][genome]
                except KeyError:
                    bmk[tool][genome] = {}
                bench = open(os.path.join(testDir, f'{tool}/benchmarks/{genome}_{tool}.txt'), 'r')
                for line in bench:
                    l = line.strip().split(': ')
                    try:
                        benchmarks.index(l[0])
                        bmk[tool][genome][l[0]] = l[1]
                    except ValueError:
                        continue
                bench.close()
        out = open(output[0],'w')
        out.write('Prophage_caller\tGenome')
        for bench in benchmarks:
            out.write(f'\t{bench}')
        out.write('\n')
        for tool in TOOLS:
            for genome in GENOMES:
                out.write(f'{tool}\t{genome}')
                for bench in benchmarks:
                    out.write(f'\t{bmk[tool][genome][bench]}')
                out.write('\n')
        out.close()