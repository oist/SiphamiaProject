# SiphamiaProject
RNA-seq Differential Expression project with _Siphamia tubifer_ and _Photobacterium mandapamensis_

# Workflow
## Bacteria
1. Add Spike-in sequences to _P. mandapamensis_ genome
2. Create reference genome index with RSEM 
      `rsem-prepare-reference --bowtie2 $DATA/siphgenome/spike/refwspike.fa $DATA/siphgenome/bt2ref/out/resem_bt2spike_ref`
3. Trim barcodes from sample files with Flexbar 
      `flexbar -t $HOME/SiphComp/trimmed/B30_trim_ -r $DATA/B30_1.fastq -p $DATA/B30_2.fastq -f fastq -a $HOME/SiphComp/adapterfiles/B_adapters.fasta -ao 1`

4. Calculate differential expression with RSEM (``bactRSEM.slurm``, ``loRSEM.slurm``)
5. Check ERCC results `bact.ERCC.R`
5. Analyse differential expression with edgeR (``bacteria.genes.R``)
  * Results are `combResults_tgw_ex1.csv`
6. GO enrichment
  * gene ID to GO term reference is `geneGOref.csv`, created with `BactXML.ipynb` (utilizing `geneidref.txt`, `geneids.txt`, and `genenume.txt`)
  * preparing `geneGOref.csv` for GOstats: `GoTermsR.R`
