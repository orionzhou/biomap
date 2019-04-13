# Short variant calling of the biomAP resequencing dataset

## Methods
* Read were trimmed using [fastp], aligned using [BWA]() against the maize 
  [B73 AGP_v4]() reference.
* Alignment BAM files were duplicate marked and base recalibrated using [GATK4]
* Per-sample GVCF were first called using [GATK4] HaplotypeCaller, followed by 
  a joint variant calling among all samples by [GATK4] GenotypeGVCF to 
  generate the raw variant set
* [GATK4] VariantRecalibration was then performed to filter the raw set of 
  variants to generate the final variant calling set.

## Results
* [Per-sample SNP statistics](../data/variants/06.stat.snp.txt)
* [Per-sample Indel statistics](../data/variants/06.stat.indel.txt)
* Sample genotype table with variant effect annotation:

    /home/springer/zhoux379/projects/biomap/data/variants/10.tsv.gz

* Sample phylogeny (to be created)
