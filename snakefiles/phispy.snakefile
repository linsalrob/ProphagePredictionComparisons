"""
Phispy (v2)

Manuscript: (v1) https://academic.oup.com/nar/article/40/16/e126/1027055
Software: (v2) https://github.com/linsalrob/PhiSpy

"""

import os
import sys


# CONFIG
outDirName = "phispy"
psBuild = os.path.join(workflow.basedir, "../build/phispy")
use_pvogs = ''
training_set = ''


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")
include: os.path.join(workflow.basedir, "../rules/phispy_rules.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, f"{{genome}}_{outDirName}_tptn.tsv"), genome=GENOMES)
