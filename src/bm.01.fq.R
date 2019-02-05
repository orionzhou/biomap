source("functions.R")
dirw = file.path(dird, '01_exp_design')

#{{{ rename samples
fi = file.path(dirw, "01_exp_design/04.typo.corrected.tsv")
ti = read_tsv(fi) %>%
    mutate(raw_idx = sid,
           sid = sprintf("bm%03d", idx)) %>%
    select(sid, raw_idx, sample, everything()) %>%
    select(-idx)
fo = file.path(dirw, "06.sample.renamed.tsv")
write_tsv(ti, fo)
#}}}

#{{{ create 11.inbred.txt
fi = file.path(dirw, '10.reads.tsv')
ti = read_tsv(fi)

t_in = ti %>% filter(inbred) %>% distinct(Genotype)
t_hy = ti %>% filter(!inbred) %>% distinct(Genotype) %>%
    separate(Genotype, c("pa1","pa2"), sep = "x", remove = F)
fo = file.path(dirw, "11.inbred.txt")
write_tsv(t_in, fo, col_names = F)
fo = file.path(dirw, "11.hybrid.txt")
write_tsv(t_hy, fo, col_names = F)
#}}}

# characterize pairwise SNPs
get_line_count <- function(fx) {
    #{{{ wc -l fx
    if(!file.exists(fx))
        0
    else
        as.integer(system2("wc", args = c("-l", fx, " | awk '{print $1}'"), stdout = T))
    #}}}
}
fi = file.path(dirw, "10.reads.tsv")
ti = read_tsv(fi)
tg = ti %>% distinct(type, genotype) %>%
    arrange(desc(type), genotype) %>%
    mutate(fvr = sprintf("%s/data/variants_rnaseq/%s.bed", dirw, genotype),
           fvg = sprintf("%s/data/variants_reseq/%s.bed", dirw, genotype),
           num_snp_rnaseq = sapply(fvr, get_line_count),
           num_snp_reseq = sapply(fvg, get_line_count))

to = tg %>% select(-fvr, -fvg)
fo = file.path(dirw, '../05.num_snps.tsv')
write_tsv(to, fo)


