source("functions.R")
dirw = file.path(dird, '71_share')

fi = '~/projects/rnaseq/data/raw_output/me99c/cpm.rds'
res = readRDS(fi)
tm = res$tm_m

th = res$th_m %>% select(-Treatment) %>%
    mutate(Tissue = str_sub(Tissue, 1, 1))
th %>% count(Tissue)
th %>% count(Genotype)

to = tm %>% select(gid, SampleID, CPM) %>%
    inner_join(th, by = 'SampleID') %>%
    mutate(cond = sprintf("Zm_%s_%s_x_mRNA", Genotype, Tissue)) %>%
    select(gid, cond, CPM) %>%
    spread(cond, CPM)

fo = file.path(dirw, 'cpm_for_pete.rds')
saveRDS(to, file=fo)
