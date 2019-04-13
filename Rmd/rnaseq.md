# Analysis of biomAP mRNA-Seq dataset

## Methods

## Results
- [Raw sample meta table](https://github.com/orionzhou/maize.expression/blob/master/data/05_read_list/me99c.tsv)
- [corrected sample meta table](https://github.com/orionzhou/maize.expression/blob/master/data/05_read_list/me99c.c.tsv):
  - sample `bm252` corrected from `Root` to `Leaf`
  - msi location: `/home/springer/zhoux379/projects/rnaseq/data/05_read_list/me99c.c.tsv`
- Intermediate files are all under MSI scratch space directory:
    `/scratch.global/zhoux379/rnaseq/me99c/`, with the following subdirectories:
  - `10_fastq`, `15_trim`: raw and trimmed fastq files
  - `21_star`, `22_bam`: raw and coordinate-sorted BAM files
  - `31_featurecounts`: read count table
  - `31_mmquant`: raw read counts
  - `33_ase`: intermediate files in the allele count pipeline

### QC:
- [read trimming and mapping statitics](/data/41_qc/01.readmapping.pdf)
- [PCA plot](/data/41_qc/12.pca.pdf)
- [t-SNE plot](/data/41_qc/12.tsne.pdf): note the much improved separation of tissues compared to PCA plot, and the separation of leaf samples into two distinct clusters

### ASE analysis:
- [proportion conflicting reads for each sample](/data/41_qc/15.ase.pcft.pdf): 
  genes in most samples have very low proportion of conflicting reads
  (typically less than 2%); abnormally high rates of conflicting reads
  might indicate sample genotype mis-labelling
- [proportion reads w. paternal allele in each sample](/data/41_qc/15.ase.pref.pdf):
  numbers indicate sample size (i.e., number of genes). most inbred 
  samples have 0 paternal allele proportion while hybrid samples 
  have 0.5 parternal allele proportion, exceptions may indicate 
  sample genotype mis-labelling

### [Mapping stats table](https://github.com/orionzhou/rnaseq/blob/master/data/raw_output/me99c/bamstats.tsv)
  contains trimming, mapping and counting statistics for each sample with the following columns:
  - msi path: `/home/springer/zhoux379/projects/rnaseq/data/raw_output/bamstats.tsv`
  - Sample meta data: `SampleID`, `Tissue`, `Genotype`, `Treatment`, `Replicate`:
  - Trimming statistics: `total`, `surviving`, `surviving_f`, `surviving_r`, `dropped`
  - Mapping statistics:
    - `pair`: read pairs
      - `pair_bad`, `pair_dup`: pairs that failed QC or duplicates
      - `pair_map`: mapped pairs (both ends)
      - `pair_orphan`: pairs with one end mapped
      - `pair_unmap`: unmapped pairs
    - `unpair`: singletons (single-end reads or pairs with one end failing QC)
      - `unpair_bad`, `unpair_dup`: singletons that failed QC or duplicates
      - `unpair_map`: mapped reads
      - `unpair_unmap`: unmapped reads
    - `pair_map_hq`, `pair_orphan_hq`, `unpair_map_hq`: pairs/reads mapped with high quality (i.e., unique)
    - `pair_map0`, `pair_orphan0`, `unpair_map0`: pairs/reads mapped with 0 mismatch
    - `pair_map_hq0`, `pair_orphan_hq0`, `unpair_map_hq0`: pairs/reads mapped with high quality (i.e., unique) and with 0 mismatch
  - Read counting statistics:
    - `Assigned`: reads assigned to exonic regions and thus counted
    - `Unassigned_MultiMapping`, `Unassigned_NoFeatures`, `Unassigned_Ambiguity`, `Unassigned_Unmapped`: reads not counted due to various reasons

### R data file (msi path: `/home/springer/zhoux379/projects/biomap/data/41_qc/10.rc.ase.rda`)
  containing raw read count tables, normalized expression values and allele-specific read counts:
* th - tibble for library (sample), with columns:
  * `SampleID`: bm001 - bm467
  * ` Tissue`: Leaf, Internode, Root, etc.
  * `Genotype`: B73, Mo17xPH207, etc
  * `Treatment`: replicate 1 or 2
  * `inbred`: whether this is the inbred parent (TRUE) or hybrid (FALSE)
  * `sizeFactor`, `libSize`: library size and normalization factor calculated using the median log ratio approach by DESeq2, accounts for library size
  * `normFactor`: library normalzation factor computed by edgeR using the TMM approach, does NOT account for library size
* tm - tibble for biomap expression data
  * `gid`: Gene ID (AGP_v4, Ensembl Plants v37, 46,117 in total)
  * `SampleID`: bm001 - bm467
  * `ReadCount`: raw read count
  * `nRC`: normalized read count (`nRC = ReadCount / sizeFactor`)
  * `rCPM`: raw CPM (adds up to 1,000,000 for each sample/library)
  * `rFPKM`: raw FPKM calculated using rCPM and gene exon length
  * `CPM`: CPM calculated by edgeR (`CPM = rCPM / normFactor`)
  * `FPKM`: FPKM calculated using CPM and gene exon length
* ta - tibble with allele-specific read counts
  * `SampleID`: bm001 - bm467
  * `gid`: AGP_v4 Gene ID
  * `n0`, `n1`: number of reads supporting the paternal and maternal allele
  * `ncft`: number of reads with 'conflicting' evidence supporting both the paternal and maternal alleles, this is rare and most often caused by mis-alignment around InDel regions.
