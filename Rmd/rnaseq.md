# Analysis of biomAP mRNA-Seq dataset

## Methods

## Results
* Intermediate files are all under MSI scratch space directory:
    `/scratch.global/zhoux379/rnaseq/me99c`, with the following subdirectories:
  * `10_fastq`, `15_trim`: raw and trimmed fastq files
  * `21_star`, `22_bam`: raw and coordinate-sorted BAM files
  * `31_mmquant`: raw read counts

## QC:
- [read trimming and mapping statitics](/data/41_qc/01.readmapping.pdf)
- [PCA plot](/data/41_qc/12.pca.pdf)
- [t-SNE plot](/data/41_qc/12.tsne.pdf): note the much improved separation of tissues compared to PCA plot, and the separation of leaf samples into two distinct clusters

## ASE analysis:
- [proportion conflicting reads for each sample](/data/41_qc/15.ase.pcft.pdf): 
  genes in most samples have very low proportion of conflicting reads
  (typically less than 2%); abnormally high rates of conflicting reads
  might indicate sample genotype mis-labelling
- [proportion reads w. paternal allele in each sample](/data/41_qc/15.ase.pref.pdf):
  numbers indicate sample size (i.e., number of genes). most inbred 
  samples have 0 paternal allele proportion while hybrid samples 
  have 0.5 parternal allele proportion, exceptions may indicate 
  sample genotype mis-labelling