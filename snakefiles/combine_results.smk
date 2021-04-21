
# CONFIG
outDirName = ''
TOOLS = ['dbscan-swa',
         'phage_finder',
         'phageboost',
         'phigaro',
         'phispy',
         'phispy_trained',
         'phispy_pvog',
         'virsorter',
         'virsorter2']


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")


rule all:
    input:
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_benchmarks.tsv"),
        os.path.join(workflow.basedir, "../jupyter_notebooks/all_tptn.tsv")


rule combine_tptn:
    input:
        expand(os.path.join(testDir, '{tool}/{genome}_{tool}_tptn.tsv'), tool=TOOLS, genome=GENOMES)
    output:
        os.path.join(workflow.basedir,"../jupyter_notebooks/all_tptn.tsv")
    run:
        out = open(output[0],'w')
        out.write(f'Prophage Caller\tGenome\tTP\tTN\tFP\tFN\tAccuracy\tPrecision\tRecall\tSpecificity\tf1 score\n')
        for tool in TOOLS:
            for genome in GENOMES:
                s = {}
                tptn = open(os.path.join(testDir, f'{tool}/{genome}_{tool}_tptn.tsv'), 'r')
                for line in tptn:
                    if not line.startswith('#'):
                        l = line.split()
                        s[l[0]] = l[1]
                tptn.close()
                out.write('\t'.join( [ tool, genome, s["TP:"], s["TN:"], s["FP:"], s["FN:"], s["Accuracy:"],
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
                out.write('\t'.join([tool, genome] + b ))
                out.write("\n")
        out.close()
