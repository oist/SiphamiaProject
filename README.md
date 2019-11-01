# SiphamiaProject
RNA-seq Differential Expression project with _Siphamia tubifer_ and _Photobacterium mandapamensis_

# Workflow
## Bacteria
1. Add Spike-in sequences to _P. mandapamensis_ genome
2. Create reference genome index with RSEM 
      `rsem-prepare-reference --bowtie2 refwspike.fa bt2ref/out/resem_bt2spike_ref`
3. Trim barcodes from sample files with Flexbar 
      `flexbar -t trimmed/B30_trim_ -r DATA/B30_1.fastq -p DATA/B30_2.fastq -f fastq -a adapterfiles/B_adapters.fasta -ao 1`

4. Calculate differential expression with RSEM 
      `rsem-calculate-expression --bowtie2 --bowtie2-sensitivity-level very_sensitive --seed-length 19 --num-threads 12 --paired-end trimmed/B30_trim__1.fastq /trimmed/B30_trim__2.fastq  bt2ref/out/rsem_bt2spike_ref RSEM/B30`
5. Check ERCC results `bact.ERCC.R`
5. Analyse differential expression with edgeR (``bacteria.genes.R``)
  * Results are `combResults_tgw_ex1.csv`
6. GO enrichment
  * gene ID to GO term reference is `geneGOref.csv`, created with `BactXML.ipynb` (utilizing `geneidref.txt`, `geneids.txt`, and `genenume.txt`)
  * preparing `geneGOref.csv` for GOstats: `GoTermsR.R`
