source("bm.fun.r")
dirw = file.path(dird, '41_qc')
fi = file.path(dirw, '10.rc.norm.ase.rda')
x = load(fi)

#{{{ prepare for hclust and pca 
ths = th %>% distinct(Tissue, Genotype)
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
cols5 = pal_d3()(5)
shapes2 = c(15,1)
#}}}

#{{{ hclust 
cor_opt = "pearson"
#cor_opt = "spearman"
hc_opt = "ward.D"
edist <- as.dist(1-cor(e, method = cor_opt))
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]
#
tp = th %>% mutate(taxa = SampleID, lab = SampleID) 
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
if(length(reps)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Replicate), lab)
t_hc = tp %>% select(taxa, everything())

thp = mutate(th, 
    lab = sprintf("%s %s %s %s", SampleID, Tissue, Genotype, Replicate))
p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    ggplot2::xlim(0, 15) + 
    theme_tree2()
p1 = p1 %<+% thp + 
    geom_tiplab(aes(label = lab), size = 4, offset = 0.04) + 
    #geom_text(aes(color = as.character(gt_ok), label = gt), size = 4, nudge_x = 6, hjust = 0) + 
    scale_color_manual(values = c("black", "royalblue", "tomato"))
fo = sprintf("%s/08.hclust.pdf", dirw, cor_opt, hc_opt)
ggsave(p1, filename = fo, width = 12, height = 30)
#}}}

#{{{
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tismap = LETTERS[1:length(tissues5)]
names(tismap) = tissues5
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(label = tismap[Tissue]) %>%
    mutate(Tissue = factor(Tissue, levels = tissues5))
p1 = ggplot(tp) +
    geom_point(aes(x=PC1, y=PC2, color=Tissue, shape=inbred), size = 1.5) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_manual(values = cols5) +
    scale_shape_manual(values = shapes2, labels = c("inbred", "hybrid")) +
    otheme(legend.pos = 'bottom.right', legend.dir = 'v',
           xtitle = T, ytitle = T, xtext = T, ytext = T,
           xgrid = T, ygrid = T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/12.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 7, height = 7)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims = 2, perplexity=30, verbose=T, 
              pca = T, max_iter = 500)

cols23 = c(pal_ucscgb()(23), pal_uchicago()(5))
tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID') %>%
    mutate(Tissue = factor(Tissue, levels = tissues5))
#tps = tp %>% mutate(txt = ifelse(Genotype == 'BxM' & Replicate == 1, Tissue, ''))
p_tsne = ggplot(tp) +
    #geom_text_repel(data=tps, aes(x=V1,y=V2,label=txt), size = 3, alpha = .8) +
    geom_point(aes(x = V1, y = V2, shape=inbred, color=Tissue), size = 2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = shapes2, labels = c("inbred", "hybrid")) +
    scale_color_manual(values = cols5) +
    otheme(legend.pos = 'bottom.right', legend.dir = 'v',
           xtitle = T, ytitle = T, xtext = T, ytext = T,
           xgrid = T, ygrid = T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
    #theme(panel.border = element_blank()) +
fp = sprintf("%s/12.tsne.pdf", dirw)
ggsave(p_tsne, filename = fp, width = 7, height = 7)
#}}}

#{{{ heatmap
tps = th %>% select(-Treatment, -paired) %>%
    mutate(Tissue = factor(Tissue, levels = tissues23),
           Genotype = factor(Genotype, levels = gts)) %>%
    arrange(Tissue, Genotype) %>%
    mutate(x = length(Tissue):1)
mat = 1 - as.matrix(edist)
tp = mat %>% as.data.frame() %>% 
    rownames_to_column(var='s1') %>% as_tibble() %>%
    gather(s2, sim, -s1) %>%
    inner_join(tps[,c('SampleID','x')], by = c('s1'='SampleID')) %>%
    inner_join(tps[,c('SampleID','x')], by = c('s2'='SampleID')) %>%
    rename(x = x.x, y = x.y)
tpx = tps %>% group_by(Tissue) %>%
    summarise(xmin = min(x), xmax = max(x), x = median(x),
              lab = Tissue[1]) %>% ungroup()
cols3 = pal_npg()(3)
names(cols3) = gts
tpg = tps %>% group_by(Tissue, Genotype) %>%
    summarise(xmin = min(x), xmax = max(x), x = median(x),
              lab = Genotype[1]) %>%
    ungroup() %>%
    mutate(col.gt = cols3[Genotype]) %>%
    mutate(lab = ifelse(lab=='BxM', 'F1', str_sub(lab,1,1)))
xt = -7; xg = -.5

p = ggplot(tp) +
    geom_tile(aes(x = x, y = y, fill = sim)) +
    geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(breaks = tpx$x, labels = tpx$lab, expand = c(0,0)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    #scale_fill_viridis() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           ytext = T) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 7))
fo = sprintf("%s/08.heatmap.pdf", dirw)
ggsave(p, file = fo, width = 9, height = 8.5)
#}}}

