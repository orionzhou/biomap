# biomAP resequencing-based variants

- convert vcf to bcf
  - `bcftools view biomap62.vcf.gz -O b -o biomap62.bcf`
  - `bcftools view biomap62.snp.vcf.gz -O b -o biomap62.snp.bcf`
- filter raw resequencing VCF to 35 inbred lines used in the study
```
  bcftools view -S ../01_exp_design/11.inbred.txt --force-samples biomap62.snp.bcf -Ou |\
  bcftools view -m2 -M2 -v snps -i 'INFO/AC[0] > 1' -Oz -o biomap29.1.snp.vcf.gz
```
  - 29 out of the 35 inbred lines were resequenced and have genotypic information
  - 5 (Ny821, NKS8326, H99, DKFAPW, W64A) were not resequenced
  - 1 (Q381) was treated as identical to PH207
  
- genotype VCF (`biomap29.2.expanded.vcf.gz`) was created for 30 inbreds (`1|1`) and 96 hybrids (`0|1` or `1|0`) by running `biomap.py`.
- lift over to new genome coordinates:
```
  zcat biomap29.2.expanded.vcf.gz |\
  vcf-liftover.pl $genome/B73/08_seq_map/mapf.chain unmapped > biomap29.3.lifted.vcf
```
- create one vcf and bed for each genotype by running `biomap.py`

# biomAP RNA-Seq-based variants

- convert vcf to bcf
  - `bcftools view biomap35.snp.vcf.gz -O b -o biomap35.snp.bcf`
- genotype VCF (`biomap35.1.expanded.vcf.gz`) was created for 35 inbreds (`1|1`) and 96 hybrids (`0|1` or `1|0`) by running biomap.py.
- lift over to new genome coordinates:
```
  zcat biomap35.1.expanded.vcf.gz |\
  vcf-liftover.pl $genome/B73/08_seq_map/mapf.chain unmapped > biomap35.2.lifted.vcf
```
- create one vcf and bed for each genotype by running `biomap.py`