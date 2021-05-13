"""
Phispy (v2)

Manuscript: (v1) https://academic.oup.com/nar/article/40/16/e126/1027055
Software: (v2) https://github.com/linsalrob/PhiSpy

This snakefile runs Phispy using the training sets for each genbank file.

"""

import os
import sys


# CONFIG
outDirName = "phispy_trained"
psBuild = os.path.join(workflow.basedir, "../build/phispy")
use_pvogs = ''

def training_set(wildcards):
    ts = {
        "Achromobacter_denitrificans_strain_PR1": "-t data/trainSet_32002.17.txt",
        "Bacillus_halodurans_C-125": "-t data/trainSet_272558.23.txt",
        "Bacillus_subtilis_subsp._subtilis_str._168": "-t data/trainSet_224308.360.txt",
        "Bacteroides_uniformis_ATCC_8492_strain_81A2": "-t data/trainSet_411479.31.txt",
        "Bifidobacterium_longum_NCC2705": "-t data/trainSet_206672.37.txt",
        "Brucella_melitensis_16M": "-t data/trainSet_224914.79.txt",
        "Caulobacter_crescentus_CB15": "-t data/trainSet_190650.21.txt",
        "Clostridium_perfringens_str._13": "-t data/trainSet_195102.53.txt",
        "Clostridium_tetani_E88": "-t data/trainSet_212717.31.txt",
        "Deinococcus_radiodurans_R1": "-t data/trainSet_243230.96.txt",
        "Enterococcus_faecalis_strain_V583": "-t data/trainSet_1351.557.txt",
        "Enterococcus_faecalis_V583": "-t data/trainSet_226185.9.txt",
        "Escherichia_coli_CFT073": "-t data/trainSet_199310.168.txt",
        "Escherichia_coli_K12": "-t data/trainSet_83333.998.txt",
        "Escherichia_coli_O157-H7_EDL933": "-t data/trainSet_155864.289.txt",
        "Escherichia_coli_O157-H7": "-t data/trainSet_83334.295.txt",
        "Haemophilus_influenzae_Rd_KW20": "-t data/trainSet_71421.45.txt",
        "Lactococcus_lactis_subsp._lactis_Il1403": "-t data/trainSet_272623.42.txt",
        "Listeria_innocua_Clip11262": "-t data/trainSet_272626.22.txt",
        "Listeria_monocytogenes_EGD-e": "-t data/trainSet_169963.176.txt",
        "Mesorhizobium_loti_MAFF303099": "-t data/trainSet_266835.41.txt",
        "Mycobacterium_tuberculosis_CDC1551": "-t data/trainSet_83331.121.txt",
        "Mycobacterium_tuberculosis_H37Rv": "-t data/trainSet_83332.460.txt",
        "Neisseria_meningitidis_MC58": "-t data/trainSet_122586.26.txt",
        "Neisseria_meningitidis_Z2491": "-t data/trainSet_122587.18.txt",
        "Paracoccus_aminophilus_JCM_7686": "-t data/trainSet_1367847.3.txt",
        "Paracoccus_denitrificans_PD1222": "-t data/trainSet_318586.5.txt",
        "Paracoccus_sanguinis_5503": "-t data/trainSet_1525717.3.txt",
        "Paracoccus_sp._SCN_68-21": "-t data/trainSet_1660154.3.txt",
        "Paracoccus_yeei_TT13": "-t data/trainSet_147645.106.txt",
        "Pasteurella_multocida_subsp._multocida_str._Pm70": "-t data/trainSet_272843.53.txt",
        "Pseudomonas_aeruginosa_PAO1": "-t data/trainSet_208964.452.txt",
        "Pseudomonas_putida_KT2440": "-t data/trainSet_160488.79.txt",
        "Ralstonia_solanacearum_GMI1000": "-t data/trainSet_267608.42.txt",
        "Salmonella_enterica_subsp._enterica_serovar_Typhi_str._CT18": "-t data/trainSet_220341.87.txt",
        "Shewanella_oneidensis_MR-1": "-t data/trainSet_211586.69.txt",
        "Shigella_flexneri_2a_str._301": "-t data/trainSet_198214.txt",
        "Staphylococcus_aureus_strain_Sa_Newman_UoM": "-t data/trainSet_1280.10152.txt",
        "Staphylococcus_aureus_subsp._aureus_Mu50": "-t data/trainSet_158878.38.txt",
        "Staphylococcus_aureus_subsp._aureus_MW2": "-t data/trainSet_196620.15.txt",
        "Streptococcus_pyogenes_M1_GAS": "-t data/trainSet_160490.61.txt",
        "Streptococcus_pyogenes_MGAS315": "-t data/trainSet_198466.10.txt",
        "Streptococcus_pyogenes_MGAS8232": "-t data/trainSet_186103.26.txt",
        "Vibrio_cholerae_O1_biovar_eltor_str._N16961": "-t data/trainSet_243277.252.txt",
        "Xanthomonas_axonopodis_pv._citri_str._306": "-t data/trainSet_190486.46.txt",
        "Xylella_fastidiosa_9a5c": "-t data/trainSet_160492.65.txt",
        "Xylella_fastidiosa_Temecula1": "-t data/trainSet_183190.38.txt",
        "Yersinia_pestis_CO92": "-t data/trainSet_214092.200.txt",
        "Yersinia_pestis_KIM": "-t data/trainSet_187410.24.txt",
    }
    if wildcards.genome in ts:
        return ts[wildcards.genome]
    else:
        sys.stderr.write(f"FATAL: No training set for {wildcards.genome}\n")
        sys.exit(-2)


# GENERIC CONFIG/RECIPES
include: os.path.join(workflow.basedir, "../rules/preflight.smk")
include: os.path.join(workflow.basedir, "../rules/phispy_rules.smk")


# TARGETS
rule all:
    input:
        expand(os.path.join(outputdir, f"{{genome}}_{outDirName}_tptn.tsv"), genome=GENOMES)
