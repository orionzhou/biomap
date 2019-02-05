source("functions.R")
dirw = '/home/springer/zhoux379/projects/biomap/data/variants'

fv = file.path(dirw, "db.snp")
tv = read_tsv(fv, col_names = c("vid", "type", "chrom", "pos", "note"))

#fh = file.path(dirw, "db.haplotype")
#th = read_tsv(fh, col_names = c("hid", "chrom", "start", "end", "vids"))

tv2 = tv %>% mutate(gpos = ceiling(chrom * 100000000 + pos/10))
tvs = tv2 %>% count(gpos) %>% arrange(-n)

nrow(tv)

gpos_rm = tvs$gpos[tvs$n >= 3]
to = tv2 %>% filter(! gpos %in% gpos_rm) %>% select(-gpos)
nrow(to)
fo = file.path(dirw, 'db.d3.snp')
write_tsv(to, fo, col_names = F)

gpos_rm = tvs$gpos[tvs$n >= 2]
to = tv2 %>% filter(! gpos %in% gpos_rm, type == 'single') %>% select(-gpos)
nrow(to)
fo = file.path(dirw, 'db.d2.snp')
write_tsv(to, fo, col_names = F)
