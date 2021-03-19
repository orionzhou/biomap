source("functions.R")
dirw = file.path(dird, 'variants/widiv942')

#widiv942

fi = file.path(dirw, '00.raw.tsv')
ti = read_tsv(fi)
fr = file.path(dirw, '01.ref.tsv')
tr = read_tsv(fr, col_names=c('chrom','start','pos','ref')) %>% select(-start)

ti2 = ti %>% select(-alleles,-strand,-assembly,-center,-protLSID,-assayLSID,
    -panel,-QCcode) %>%
    gather(gt, nt, -rs, -chrom, -pos)

ti3 = ti2 %>% group_by(rs, chrom, pos, gt) %>%
    summarise(n_allele = length(unique(nt))) %>% ungroup()
