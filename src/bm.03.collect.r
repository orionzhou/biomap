source("bm.fun.r")
dirw = file.path(dird, '41_qc')

sid = 'me99c'
diri = '~/projects/maize.expression/data/11_qc'

#{{{ read mapping stats & plot
fi = file.path(diri, sid, '10.mapping.stat.tsv') 
ti = read_tsv(fi)
types = c("FailedQC", "Unmap", "Map_LowQual", "Map_HighQual", 'Total')
tx = ti %>%
    mutate(Total = total,
           FailedQC = dropped,
           Unmap = pair_unmap + unpair_unmap,
           Map_LowQual = pair_map-pair_map_hq + pair_orphan-pair_orphan_hq + unpair_map-unpair_map_hq,
           Map_HighQual = pair_map_hq + pair_orphan_hq + unpair_map_hq) %>%
    select(Tissue, Genotype, inbred, FailedQC, Unmap, Map_LowQual, Map_HighQual, Total) %>%
    gather(type, ReadPairCount, -Tissue, -Genotype, -inbred) %>%
    group_by(Tissue, Genotype, inbred, type) %>%
    summarise(ReadPairCount = sum(ReadPairCount)) %>% ungroup() %>%
    mutate(Tissue = factor(Tissue, levels = tissues5)) %>%
    mutate(type = factor(type, levels = rev(types)))

tx1 = tx %>% mutate(mreads = ReadPairCount/1000000) %>% select(-ReadPairCount)
tps = tx1 %>% filter(type == 'Total') %>% select(-type) %>%
    rename(total = mreads)
tp = tx1 %>% 
    filter(type != 'Total') %>%
    inner_join(tps, by = c('Tissue','Genotype','inbred')) %>%
    mutate(prop = mreads/total) %>%
    mutate(lab = sprintf("%.02f", mreads/total)) %>%
    mutate(lab = str_replace(lab, "^0[\\.]", "."))
#
cols5 = pal_simpsons()(5)[c(1,2,4,5,3)]
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = mreads, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_manual(values = cols5) +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 11, strip.position = 'top', scales = 'free_x') +
    otheme(xtitle = T, xtext = T, ytext = T, xgrid = T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines"))
fp = file.path(dirw, "01.readmapping.pdf")
ggsave(p1, filename = fp, width = 10, height = 12)
#}}}

#{{{ save to 
fi = file.path(diri, sid, '20.rc.norm.rda') 
x = load(fi)
fi = file.path(diri, sid, '10.mapping.stat.tsv') 
ti = read_tsv(fi)
th = ti %>% select(SampleID, Tissue, Genotype, Replicate) %>%
    inner_join(tl, by = 'SampleID')
fa = file.path(diri, '../08_raw_output', sid, 'ase.tsv') 
ta = read_tsv(fa)

fo = file.path(dirw, '10.rc.dom.ase.rda')
save(th, tm, ta, file = fo)
#}}}

#{{{ ## save to 00.table.rda
tissues6 = c('coleoptile_tip', 'radicle_root', 'embryo_imbibedseed', 'seedlingleaf_11DAS', 'seedlingroot_11DAS', 'seedlingmeristem_11DAS')
tt = ti %>%
    mutate(Condition = ifelse(Tissue %in% tissues6, 'Growth chamber', 'Field')) %>%
    mutate(TotalReadPair = total,
           TrimmedReadPair = surviving,
           MappingRate = pair_map/surviving,
           UniqueMappingRate = pair_map_hq/surviving,
           AssignedRate = Assigned/surviving) %>%
    mutate(MappingRate = sprintf("%.1f%%", MappingRate*100), 
           UniqueMappingRate = sprintf("%.1f%%", UniqueMappingRate*100), 
           AssignedRate = sprintf("%.1f%%", AssignedRate*100)) %>%
    select(SampleID, Tissue, Genotype, Replicate, Condition,
           TotalReadPair, TrimmedReadPair,
           MappingRate, UniqueMappingRate, AssignedRate)
ft = file.path(dirw, "00.table.rda")
save(tt, file = ft)
#}}}
