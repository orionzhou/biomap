source("bm.fun.r")
dirw = file.path(dirp, "41.qc")

#{{{ compute nRC, rCPM, CPM & FPKM and save to 10.norm.RData
fi = '/home/springer/zhoux379/data/genome/B73/v32/t5.gtb'
gids = read_tsv(fi) %>% distinct(par) %>% pull(par)
fi = file.path('/home/springer/zhoux379/data/genome/B73', "51.tsv")
t_gs = read_tsv(fi, col_types = 'ccccciic') %>% 
    filter(etype == 'exon') %>% 
    group_by(gid, tid) %>% 
    summarise(size = sum(end - beg + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))

fi = file.path(dirw, '../10.RawReadCount.RData')
x = load(fi)
x

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, tl, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '10.norm.RData')
save(tl, tm, file = fo)
#}}}

#{{{ hclust tree
fi = file.path(dirw, "10.norm.RData")
x = load(fi)
x

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

tp = tl %>% mutate(sid = SampleID, Rep = as.character(Treatment),
    lab = sprintf("%s %s %s %s", SampleID, Tissue, Genotype, Treatment)) %>%
    left_join(tl2[,c('SampleID','ml.gt','scol')], by = 'SampleID')
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
fo = sprintf("%s/11.cpm.hclust.pdf", dirw, cor_opt, hc_opt)
ggsave(p1, filename = fo, width = 12, height = 40)
#}}}

#{{{ ASE qc
fi = file.path(dirw, '../11.ase.RData')
x = load(fi)
x
#}}}

#{{{ PCA
fi = file.path(dirw, "10.norm.RData")
x = load(fi)
x

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
tissues = unique(tl$Tissue)
tismap = LETTERS[1:length(tissues)]
names(tismap) = tissues
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(tl, by = 'SampleID') %>%
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
fp = sprintf("%s/12.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}


