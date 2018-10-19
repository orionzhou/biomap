source("bm.fun.r")
dirw = "~/projects/biomap/analysis"
diri = sprintf("%s/multiqc_data", dirw)


#{{{ trim - number read pairs
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% mutate(sid2 = Sample) %>%
    select(sid2, input_read_pairs, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) %>%
    separate(sid2, c("SampleID","suf"), sep="_") %>%
    filter(suf == 1) %>% select(-suf)
sum(ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>% pull(ndiff))
ti3 = ti2 %>% select(-input_read_pairs) %>%
    gather(type, nseq, -SampleID)
types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = tl %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(Tissue = factor(Tissue, levels = rev(tissues5)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 11, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/41.qc/01.reads.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 12)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 10, height = 10)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(sid2 = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    separate(sid2, c("SampleID", "suf"), sep = "_") %>% select(-suf) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(tl, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(Tissue = factor(Tissue, levels = rev(tissues5)),
           type = factor(type, levels = types))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 11, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/41.qc/02.mapping.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 10)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 10, height = 10)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)

types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(tl, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues5)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 11, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/41.qc/03.assigned.pdf", dirw)
#ggsave(p1, filename = fp, width = 10, height = 10)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 10, height = 10)
#}}}


