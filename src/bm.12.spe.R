source("functions.R")
dirw = file.path(dird, '41_qc')
fi = file.path(dirw, '21.trio.rda')
x = load(fi)
cpm_min = .5

#{{{ biomap
td1 = tt %>%
    mutate(ep1 = (CPM.p1 >= cpm_min),
           ep2 = (CPM.p2 >= cpm_min),
           eh = (CPM.h >= cpm_min)) %>%
    select(Tissue,Genotype,gid,ep1,ep2,eh)
#}}}

#{{{ check Baldauf2018 data
fi = '~/projects/rnaseq/data/08_raw_output/me18b/cpm.rds'
x = readRDS(fi)
tm = x$tm
th = read_tsv('~/projects/rnaseq/data/05_read_list/me18b.c.tsv') %>%
    mutate(Tissue = sprintf("Root %s", Treatment))
tmm = tm %>% inner_join(th, by = 'SampleID') %>%
    group_by(Genotype,Tissue,gid) %>%
    summarise(CPM=mean(CPM))

th0 = th %>% distinct(SampleID, Tissue, Genotype) %>%
    separate(Genotype, c("pa1","pa2"), sep = 'x', remove = F) %>%
    mutate(inbred = is.na(pa2)) %>%
    distinct(Tissue,Genotype,pa1,pa2,inbred)
th1 = th0 %>% filter(inbred) %>% select(-inbred,-pa1, -pa2)
th2 = th0 %>% filter(!inbred) %>% select(-inbred)
tm1 = tmm %>% inner_join(th1, by = c("Tissue","Genotype"))
tm2 = tmm %>% inner_join(th2, by = c("Tissue","Genotype")) %>%
    mutate(cpm.h = CPM) %>% select(-CPM) %>%
    inner_join(tm1, by = c('Tissue'='Tissue','pa1'='Genotype','gid'='gid')) %>%
    rename(cpm.p1 = CPM) %>%
    inner_join(tm1, by = c('Tissue'='Tissue','pa2'='Genotype','gid'='gid')) %>%
    rename(cpm.p2 = CPM)

td2 = tm2 %>% mutate(ep1 = (cpm.p1 >= cpm_min),
                   ep2 = (cpm.p2 >= cpm_min),
                   eh = (cpm.h >= cpm_min)) %>%
    select(Tissue,Genotype,gid,ep1,ep2,eh)
#}}}

#{{{ check briggs
fi = '~/projects/briggs/data/41_qc/10.rda'
x = load(fi)
x

tm1 = tmm %>% select(Tissue,Genotype,gid, CPM) %>%
    spread(Genotype, CPM) %>%
    transmute(Tissue = Tissue, gid = gid, cpm.h = BxM, cpm.p1 = B73, cpm.p2 = Mo17)

td3 = tm1 %>% mutate(Genotype = 'BxM') %>%
    mutate(ep1 = (cpm.p1 >= cpm_min),
           ep2 = (cpm.p2 >= cpm_min),
           eh = (cpm.h >= cpm_min)) %>%
    select(Tissue,Genotype,gid,ep1,ep2,eh)
#}}}


#{{{ plot
datasets = c("BiomAP", "Baldauf2018", "Zhou2018")
td = mutate(td1, dataset = 'BiomAP') %>%
    bind_rows(mutate(td2, dataset = 'Baldauf2018')) %>%
    bind_rows(mutate(td3, dataset = 'Zhou2018')) %>%
    mutate(dataset = factor(dataset, levels = datasets))

tp = td %>% group_by(dataset, Tissue, Genotype) %>%
    summarise(unexp_diff = sum(!ep1&!ep2&eh) - sum(ep1&ep2&!eh),
              spe_diff = sum(ep1&!ep2&eh) + sum(!ep1&ep2&eh) -
                  sum(ep1&!ep2&!eh) - sum(!ep1&ep2&!eh)) %>%
    ungroup() %>% gather(diff_type, ng, -Tissue, - Genotype, -dataset)
tps = tp %>% group_by(dataset, Tissue, diff_type) %>%
    summarise(ng = max(ng)) %>% ungroup()
p = ggplot(tp, aes(x = Tissue, y = ng)) +
    geom_boxplot() +
    geom_text_repel(data = tps, aes(x=Tissue, y=ng, label=Tissue), size = 3) +
    geom_hline(yintercept = 0) +
    facet_grid(dataset~diff_type, scale = 'free') +
    otheme(ytext = T)
fo = file.path(dirw, 'test.pdf')
ggsave(p, file = fo, width = 10, height = 10)

tp = td %>% group_by(dataset, Tissue, Genotype) %>%
    summarise(ng = sum(eh) - pmax(sum(ep1), sum(ep2))) %>%
    ungroup()
tps = tp %>% group_by(dataset, Tissue) %>%
    summarise(ng = max(ng), txt = sprintf("N=%d", n())) %>% ungroup() %>%
    mutate(txt = ifelse(dataset=='Zhou2018', '', txt))
p = ggplot(tp, aes(x = Tissue, y = ng)) +
    geom_boxplot() +
    scale_y_continuous(name='Num. More Genes expressed in hybrid than inbred parents') +
    geom_text_repel(data = tps, aes(x=Tissue, y=ng, label=txt), size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0) +
    facet_wrap(~dataset, scale = 'free', ncol=1) +
    otheme(xtext=T, ytext=T, xtitle=T, xticks=T, yticks=T, ygrid=T)
fo = file.path(dirw, '81.expressed.genes.pdf')
ggsave(p, file = fo, width = 6, height = 10)



tp = td %>% group_by(dataset, Tissue, Genotype) %>%
    summarise(een = sum(ep1&ep2&!eh),
              nne = sum(!ep1&!ep2&eh),
              ene = sum(ep1&!ep2&eh),
              enn = sum(ep1&!ep2&!eh),
              nee = sum(!ep1&ep2&eh),
              nen = sum(!ep1&ep2&!eh)) %>%
    ungroup() %>% gather(diff_type, ng, -Tissue, -Genotype, -dataset)
tps = tp %>% group_by(dataset, diff_type) %>%
    summarise(ng = max(ng)) %>% ungroup()
p = ggplot(tp, aes(x = diff_type, y = ng)) +
    geom_violin() +
    geom_text_repel(data = tps, aes(x=diff_type, y=ng, label=diff_type), size = 3) +
    geom_hline(yintercept = 0) +
    facet_wrap(~dataset, nrow = 1, scale = 'free') +
    otheme(ytext = T)
fo = file.path(dirw, 'test.pdf')
ggsave(p, file = fo, width = 10, height = 6)
#}}}

