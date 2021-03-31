"""
Phispy (v2)

Manuscript: (v1) https://academic.oup.com/nar/article/40/16/e126/1027055
Software: (v2) https://github.com/linsalrob/PhiSpy

This snakefile runs Phispy with HMM searches of the pVOG database

"""

import os
import sys


# CONFIG
outDirName = "phispy_pvog"
use_pvogs = '--phmms ' + os.path.join(psBuild, 'pVOGs.hmm')
training_set = "data/testSet_genericAll.txt"


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")
include: os.path.join(workflow.basedir, "../rules/phispy_rules.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, f"{{genome}}_{outDirName}_tptn.tsv"), genome=GENOMES)
