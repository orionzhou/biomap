source("functions.R")
dirw = file.path(dird, '41_qc')
fi = file.path(dirw, "10.rc.ase.rds")
res = readRDS(fi)
th = res$th; tm = res$tm; ta = res$ta

#{{{ remove samples, merge reps & create trios
sids_rm0 = 'bm466' # 'W64AxMo17xB73'
sids_rm1 = c('bm050','bm252','bm136','bm125','bm113')
sids_rm2 = c('bm446','bm440','bm447','bm442','bm448','bm441','bm445','bm429')
sids_rm3 = c('bm319','bm443','bm255','bm124')
sids_rm = c(sids_rm0, sids_rm1)
th = th %>% filter(!SampleID %in% sids_rm)
th1 = th %>% filter(inbred) %>% count(Tissue, Genotype)
th2 = th %>% filter(!inbred) %>% distinct(Tissue, Genotype) %>%
    separate(Genotype, c('pa1','pa2'), sep = 'x', remove=F) %>%
    left_join(th1, by = c('Tissue'='Tissue','pa1'='Genotype')) %>%
    rename(n1 = n) %>%
    left_join(th1, by = c('Tissue'='Tissue','pa2'='Genotype')) %>%
    rename(n2 = n) %>% filter(!is.na(n1), !is.na(n2)) %>%
    mutate(Tissue=factor(Tissue,levels=tissues5)) %>%
    mutate(pa2 = factor(pa2, levels = gts5)) %>% select(-n1,-n2)

#{{{ create study design table
to = th2 %>% count(Tissue, pa2) %>% spread(pa2,n)
colnames(to)[-1] = sprintf("[]x%s", colnames(to)[-1])
fo = file.path(dirw, '20.trio.tsv')
#write_tsv(to,fo, na = '')
#}}}

# create trios
tm = tm %>% filter(!SampleID %in% sids_rm) %>%
    inner_join(th[,c('SampleID','Tissue','Genotype')], by = 'SampleID') %>%
    group_by(Tissue, Genotype, gid) %>%
    summarise(ReadCount=sum(ReadCount),
              nRC = sum(nRC), CPM = mean(CPM), FPKM = mean(FPKM)) %>%
    ungroup()

ta = ta %>% filter(!SampleID %in% sids_rm) %>%
    inner_join(th[,c('SampleID','Tissue','Genotype','sizeFactor')], by = 'SampleID') %>%
    mutate(n0 = n0/sizeFactor, n1=n1/sizeFactor, ncft=ncft/sizeFactor) %>%
    group_by(Tissue, Genotype, gid) %>%
    summarise(n0=mean(n0),n1=mean(n1),ncft=mean(ncft)) %>%
    ungroup()

tm0 = tm %>% select(Tissue,gid,Genotype,nRC,CPM)
tt = tm0 %>% inner_join(th2, by = c("Tissue","Genotype")) %>%
    rename(nRC.h = nRC, CPM.h = CPM) %>%
    inner_join(tm0, by = c("Tissue"="Tissue","gid"="gid","pa1"="Genotype")) %>%
    rename(nRC.p1 = nRC, CPM.p1 = CPM) %>%
    inner_join(tm0, by = c("Tissue"="Tissue","gid"="gid","pa2"="Genotype")) %>%
    rename(nRC.p2 = nRC, CPM.p2 = CPM)

res = list(th1=th1,th2=th2,tt=tt,ta=ta)
fo = file.path(dirw, '21.trio.rds')
saveRDS(res file = fo)
#}}}


