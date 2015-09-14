# 
import subprocess
'''
# Software requirements
bowtie2/2.2.3
'''

RAWREADS = ["Pmucro1_S1_L001", "Pmucro1_S1_L002", "Pmucro3_S1_L001", "Pmucro3_S1_L002"]
KMERS = 96
BASEDIR = "/work/MikheyevU/sasha/mucrosquamatus-genome/src/"

rule all:
	input: "../ref/GCF_000211495.1_ASM21149v1.1.bt2"

#build bowtie2 index
rule bowtie2_build:
     input: "../ref/GCF_000211495.1_ASM21149v1_genomic.fna"
     output: "../ref/GCF_000211495.1_ASM21149v1.1.bt2"
     shell: "module load bowtie2/2.2.3 ; bowtie2-build {input} ../ref/GCF_000211495.1_ASM21149v1"

