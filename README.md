# SiphamiaProject
RNA-seq Differential Expression project with _Siphamia tubifer_ and _Photobacterium mandapamensis_

# Workflow
## Bacteria
1. Add Spike-in sequences to _P. mandapamensis_ genome
2. Create reference genome index with RSEM (```rsem_ref.slurm```)
3. Trim barcodes from sample files with Flexbar (``trim_bact.slurm``, ``trim_lo.slurm``)
4. Calculate differential expression with RSEM (``bactRSEM.slurm``, ``loRSEM.slurm``)
5. Check ERCC results `bact.ERCC.R`
5. Analyse differential expression with edgeR (``bacteria.genes.R``)
  * Results are `combResults_tgw_ex1.csv`
6. GO enrichment
  * gene ID to GO term reference is `geneGOref.csv`, created with `BactXML.ipynb` (utilizing `geneidref.txt`, `geneids.txt`, and `genenume.txt`)
  * preparing `geneGOref.csv` for GOstats: `GoTermsR.R`

## Fish
2. Trim barcodes from fish muscle sequences with Flexbar (``flexM.slurm``)
1. Create bowtie2 reference index for _P. mandapamensis_ genome
2. Align mixed fish and bacteria (Light Organ) samples to bacteria genome with bowtie2 (``align.slurm``)
4. Create file of unaligned sequences using samtools - these are the fish LO sequences  
  * ``sambam.slurm`` - convert ``.sam`` file from ``align.slurm`` to a ``.bam`` file
  * ``bamsort.slurm`` - sort the ``.bam`` file
  * ``unmapped.slurm`` - remove unmapped reads 
5. Convert LO files to fastq files with bedtools (``sort_unmapped.slurm``, ``bamtofast.slurm``)
6. Construct a reference transcriptome with Trinity (``trinity.slurm``)
7. Add spike-in sequences to Trinity file
8. Creat reference index with RSEM (``trans_genemap.slurm``, ``reftrans.slurm``)
9. Calculate differential expression with RSEM (``floRSEM.slurm``, ``fmRSEM.slurm``)
10. Analyse differential expression with edgeR
11. Filter Trinity.fasta file to remove genes with less than 1 FPKM averaged across samples and keep only the longest isoform for each gene with Biopython (``lowFPKM_filter.ipynb``, ``isoform_filter.ipynb``)
12. GO enrichment
  * created `FishNotFish.csv` with `multiblastXML.ipynb`
