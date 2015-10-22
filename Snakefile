# Workflow for RNA-seq analysis of Paul's bacteria (Photobacterium leiognathi subsp. mandapamensis) 

'''
# Software requirements
pyfasta
'''

BASEDIR = "/work/MikheyevU/sasha/paul-rnaseq"
BACTERIA = ['Bacteria_29a_S23', 'Bacteria_29a_S23']

rule all:
	input: expand(BASEDIR + "/data/trinity/blast/FPKMfiltered_Trinity.{n}.xml", n = map(lambda x: str(x).zfill(4), range(0,1000)))

#split trinity into 1000 chunks
rule splitTrinity:
	input: BASEDIR + "/data/trinity/FPKMfiltered_Trinity.fasta"
	output: temp(expand(BASEDIR + "/data/trinity/FPKMfiltered_Trinity.{n}.fasta", n = map(lambda x: str(x).zfill(4), range(0,1000))))
	params: outdir = BASEDIR + "/data/trinity/blast"
	shell: "pyfasta split -n 1000 {input}; mkdir -p {params.outdir}"

rule blastTrinity:
	input: BASEDIR + "/data/trinity/FPKMfiltered_Trinity.{n}.fasta"
	output: BASEDIR + "/data/trinity/blast/FPKMfiltered_Trinity.{n}.xml"
	params: outfile = BASEDIR + "/data/trinity/blast/FPKMfiltered_Trinity.{n}.xml"
	shell: "module load ncbi-blast/2.2.30+; blastx -num_threads 8 -query {input} -db /work/MikheyevU/ncbi/nr -evalue 1e-4 -max_target_seqs 1 -outfmt 5 -out {params.outfile}"