source("bm.fun.r")
dirw = file.path(dird, '03_raw_data')

if(1 == 2) {
#{{{ collect featurecounts data
fi = file.path(dirp, 'cache/31_featurecounts/01.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw
fo = file.path(dirw, "10.RawReadCount.RData")
save(t_rc, file = fo)
#}}}

#{{{ collect ASE
fh = file.path(dirw, '../01_exp_design/10.reads.tsv')
th = read_tsv(fh)
t_ase = tibble()
for (i in 1:nrow(th)) {
    fi = sprintf("%s/cache/33_ase/%s.tsv", dirp, th$SampleID[i])
    ti = read_tsv(fi)
    ti2 = ti %>% mutate(sid = sid) %>%
        select(sid, everything())
    t_ase = rbind(t_ase, ti2)
    cat(i, sid, "\n")
}
fo = file.path(dirw, "11.ase.RData")
save(t_ase, file = fo)
#}}}
}


