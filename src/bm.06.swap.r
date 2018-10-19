source("bm.fun.r")
dirw = file.path(dirp, 'analysis/06_sample_swap')
fi = file.path(dirw, "../03_raw_data/20.rc.norm.RData")
x = load(fi)
x
fh = file.path(dirw, "../01_exp_design/10.reads.tsv")
th = read_tsv(fh) %>%
    separate(Genotype, c("pa1", "pa2"), "x", remove = F) %>%
    replace_na(list(pa2=''))

#{{{ collect gtcheck results
tck = tibble()
for (i in 1:nrow(th)) {
    sid = th$SampleID[i]
    fck = sprintf("%s/cache/27_gtcheck/%s.txt", dirp, sid)
    if(file.exists(fck)) {
    tck1 = read_tsv(fck, col_names = F, comment = "#") %>%
        transmute(SampleID = sid,
                  dist.tot = X2,
                  dist.avg = X3,
                  nSite = X4,
                  tGenotype = X5,
                  sIdx = X6)
    tck = rbind(tck, tck1)
    }
}
tck %>% count(SampleID)
fo = file.path(dirw, "01.gtcheck.RData")
save(tck, file = fo)
#}}}

#{{{ read & process gtcheck results
fck = file.path(dirw, "01.gtcheck.RData")
x = load(fck)
gts_unk = c('Ny821', 'NKS8326', 'H99', 'DKFAPW', 'W64A')
#
tk1 = tck %>% group_by(SampleID) %>%
     top_n(-5, dist.avg) %>% select(SampleID, tGenotype, dist.avg)
th2 = th %>% filter(!pa1 %in% gts_unk & !pa2 %in% gts_unk) %>%
    inner_join(tk1, by = 'SampleID') %>%
    arrange(SampleID, dist.avg) %>%
    group_by(SampleID) %>%
    summarise(Tissue = Tissue[1], Genotype = Genotype[1],
              fc = dist.avg[2]/dist.avg[1],
              ml.gt = tGenotype[1],
              fc.conf = fc >= 1.5,
              swap = Genotype != ml.gt,
              swap2 = Genotype %in% tGenotype,
              scol = ifelse(fc.conf,
                            ifelse(swap, '3', '2'), 
                            ifelse(swap2, '1', '4')))
th2 %>% count(fc.conf)
#
th2 %>% filter(swap, fc.conf) %>% print(n=50)
tck %>% filter(SampleID == 'bm005') %>% top_n(-10, dist.avg)
#}}}

#{{{ PCA
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tissues = unique(th$Tissue)
tismap = LETTERS[1:length(tissues)]
names(tismap) = tissues
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(label = tismap[Tissue])
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Tissue), size = 1.5) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = ) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/07.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}

#{{{ hclust tree
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = th %>% mutate(sid = SampleID, Rep = as.character(Treatment),
    lab = sprintf("%s %s %s %s", SampleID, Tissue, Genotype, Treatment)) %>%
    left_join(th2[,c('SampleID','ml.gt','scol')], by = 'SampleID')
cols1 = c('gray80','black','red','seagreen3',pal_d3()(5))
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    scale_x_continuous(expand = c(0,0), limits=c(-1,65)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = sid), size = 2, nudge_x = .2, hjust = 0) +
    geom_text(aes(label = Tissue, color = Tissue), size = 2, nudge_x = 2.2, hjust = 0) +
    geom_text(aes(label = Genotype), size = 2, nudge_x = 5, hjust = 0) +
    geom_text(aes(label = Rep, color = Rep), size = 2, nudge_x = 9, hjust = 0) +
    geom_text(aes(label = ml.gt, color = scol), size = 2, nudge_x = 10, hjust = 0) +
    scale_color_manual(values = cols1)
fo = sprintf("%s/11.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 12, height = 40)
#}}}


#{{{ sample list -> analysis/01.sample.RData
fh = file.path(dirw, "../01_exp_design/10.reads.tsv")
th = read_tsv(fh) %>%
    select(-paired) %>%
    separate(Genotype, c("pa1", "pa2"), "x", remove = F) %>%
    replace_na(list(pa2=''))

thg = th %>% distinct(inbred, Genotype) %>%
    arrange(desc(inbred), Genotype)
gts = thg %>% pull(Genotype)
gts_i = thg %>% filter(inbred) %>% pull(Genotype)
gts_h = thg %>% filter(!inbred) %>% pull(Genotype)
tissues5 = c("Root", "Internode", "Leaf", "Seedling", "Endosperm")
tl = th %>% select(-pa1, -pa2) %>%
    mutate(Tissue = factor(Tissue, levels = tissues5),
           Genotype = factor(Genotype, levels = gts))
fo = file.path(dirw, '01.sample.RData')
save(tl, tissues5, gts, gts_i, gts_h, file = fo)
#}}}


