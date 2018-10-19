#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column
import vcf

from maize.formats.base import must_open
from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which
from maize.utils.location import locAry2Str, locStr2Ary

def merge_refs(dirw, genomes):
    pgenomes = ['Mo17a', 'Mo17b', 'Mo17c', 'Mo17d', 'Mo17e']
    agenomes = ['B73', 'Mo17a', 'Mo17b', 'Mo17c', 'Mo17d', 'Mo17e']
    dirw = '/home/springer/zhoux379/data/genome/pseudo'
    merge_refs(dirw, pgenomes)
    hgenomes = ["B73+%s" % pa for pa in pgenomes]
    #build_db(dirw, hgenomes)
    #build_db(dirw, agenomes)
    pa1 = "B73"
    for pa2 in genomes:
        pas = [pa1, pa2]
        genome = "+".join(pas)
        dirp = op.join(dirw, genome)
        print("mkdir -p %s" % dirp)
        print("cd %s" % dirp)
        
        cfgstr = "\\n".join(["%s,../%s/11_genome.fas" % (pa,pa) for pa in pas])
        print("echo -e \"%s\" > tmp.csv" % cfgstr)
        print("fasta merge tmp.csv > 11_genome.fas")
        
        cfgstr = "\\n".join(["%s,../%s/51.gff" % (pa,pa) for pa in pas])
        print("echo -e \"%s\" > tmp.csv" % cfgstr)
        print("gff merge tmp.csv > 51.gff")
        print("gff 2gtf 51.gff > 51.gtf")
        
        print("rm tmp.csv")

def build_db(dirw, genomes):
    for genome in genomes:
        dirp = op.join(dirw, genome)
        #print("genome bowtie %s" % dirp)
        #print("genome hisat %s --p 24 --overwrite" % dirp)
        print("genome star %s --p 24 --overwrite" % dirp)

def expand_vcf(dirw, fv, fs, diro):
    os.chdir(dirw)
    vcfr = vcf.Reader(filename = fv)
    sids = vcfr.samples

    t = Table.read(fs, format = 'ascii.tab')
    gdic = dict()
    gts = set(t['Genotype'])
    for gt in gts:
        pas = gt.split("x")
        if len(pas) > 2:
            gt0 = gt
            pas = pas[:2]
            logging.debug("%s -> %s" % (gt0, "x".join(pas)))
        
        flag = ''
        for pa in pas:
            if pa not in sids:
                flag = pa
        if flag != '':
            logging.debug("%s [%s] not in VCF samples" % (flag, gt))
            continue 
        
        fo = "%s/%s.vcf" % (diro, gt)
        fho = open(fo, "w")
        fho.write("##fileformat=VCFv4.2\n")
        fho.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fho.write("\t".join("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT DUMMY".split())+"\n")
        gdic[gt] = {'gts': pas, 'fho': fho}
    
    for rcd in vcfr:
        ary = [rcd.CHROM, rcd.POS, rcd.ID, rcd.REF, rcd.ALT[0], 999, '.', '.', 'GT']
        ary = [str(x) for x in ary]
        #print(len(rcd.samples))
        for gt, idic in gdic.items():
            gts = idic['gts']
            if len(gts) == 1:
                if rcd.genotype(gt)['GT'] == '1|1':
                    idic['fho'].write("\t".join(ary + ["0|1"])+"\n")
            else:
                assert len(gts) == 2, "error gts"
                gt1, gt2 = rcd.genotype(gts[0])['GT'], rcd.genotype(gts[1])['GT']
                if gt1 == "0|0" and gt2 == "1|1":
                    idic['fho'].write("\t".join(ary + ["0|1"])+"\n")
                elif gt1 == "1|1" and gt2 == "0|0":
                    idic['fho'].write("\t".join(ary + ["1|0"])+"\n")
                elif gts[0] == 'B73' and gt2 == "1|1":
                    idic['fho'].write("\t".join(ary + ["0|1"])+"\n")
                elif gts[1] == 'B73' and gt1 == "1|1":
                    idic['fho'].write("\t".join(ary + ["1|0"])+"\n")
    
    for gt, idic in gdic.items():
        idic['fho'].close()
        cmd = "vcf hybrid_bed %s/%s.vcf > %s/%s.bed" % (diro, gt, diro, gt)
        sh(cmd)

def expand_biomap_vcf(dirw, fs, fv, fo):
    os.chdir(dirw)
    vcfr = vcf.Reader(filename = fv)
    sids = vcfr.samples

    t = Table.read(fs, format = 'ascii.tab')
    gts = set(t['Genotype'])
    gdic = dict()
    for gt in gts:
        pas = gt.split("x")
        if len(pas) > 2:
            gt0 = gt
            pas = pas[:2]
            logging.debug("%s -> %s" % (gt0, "x".join(pas)))
        
        flag = ''
        for i in range(len(pas)):
            pa = pas[i]
            if pa == 'Q381':
                pa = 'PH207'
                pas[i] = pa
            if pa not in sids:
                flag = pa

        if flag != '':
            logging.debug("%s [%s] not in VCF samples" % (flag, gt))
            continue
        gdic[gt] = pas
        
    fho = open(fo, "w")
    fho.write("##fileformat=VCFv4.2\n")
    fho.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    fds = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
    fho.write("\t".join(fds + list(gdic.keys()))+"\n")
    
    for rcd in vcfr:
        ary = [rcd.CHROM, rcd.POS, rcd.ID, rcd.REF, rcd.ALT[0], 999, '.', '.', 'GT']
        alleles = []
        for gt, pas in gdic.items():
            allele = './.' 
            if len(pas) == 1:
                al = rcd.genotype(pas[0])['GT']
                if al == '1/1':
                    allele = "1|1"
                elif al == '0/0':
                    allele = "0|0"
            else:
                assert len(pas) == 2, "error gts"
                al1, al2 = rcd.genotype(pas[0])['GT'], rcd.genotype(pas[1])['GT']
                if pas[0] == 'B73':
                    al1 = "0/0"
                if pas[1] == 'B73':
                    al2 = "0/0"
                if al1 == '0/0' and al2 == '0/0':
                    allele = '0|0'
                elif al1 == '1/1' and al2 == '1/1':
                    allele = '1|1'
                elif al1 == "0/0" and al2 == "1/1":
                    allele = "0|1"
                elif al1 == "1/1" and al2 == "0/0":
                    allele = "1|0"
            alleles.append(allele)
        ary = [str(x) for x in ary + alleles]
        fho.write("\t".join(ary) + "\n")
    fho.close()
 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'biomap variant'
    )
    #parser.add_argument('fi', help = 'input gtb')
    #args = parser.parse_args()
    
    dirw = '/home/springer/zhoux379/projects/biomap/data/variants'
    fs = "../10.reads.tsv"
    
    fv = 'biomap35.snp.vcf.gz'
    fo = "biomap35.1.expanded.vcf"
    #expand_biomap_vcf(dirw, fs, fv, fo)
    ff = "biomap35.2.lifted.vcf"
    #expand_vcf(dirw, ff, fs, 'ase35')

    fv = "biomap29.1.snp.vcf"
    fo = "biomap29.2.expanded.vcf"
    #expand_biomap_vcf(dirw, fs, fv, fo)
    ff = "biomap29.3.lifted.vcf"
    expand_vcf(dirw, ff, fs, 'ase29')

