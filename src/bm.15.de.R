source("functions.R")
dirw = file.path(dird, '41_qc')
fi = file.path(dirw, "10.rc.ase.rds")
res = readRDS(fi)
th = res$th; tm = res$tm; ta = res$ta

require(DESeq2)
require(edgeR)

tissue = tissues5[1]
sids = th %>% filter(inbred, Tissue==tissue) %>% pull(SampleID)
th1 = th %>% filter(SampleID %in% sids)
tm1 = tm %>% filter(SampleID %in% sids)

#{{{ prepare data
vh = th1 %>% mutate(Genotype = factor(Genotype)) %>% arrange(SampleID)
vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
    filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
vm = tm1 %>% filter(gid %in% gids) %>%
    select(SampleID, gid, ReadCount)
x = readcount_norm(vm)
mean.lib.size = mean(x$tl$libSize)
vm = x$tm
vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
stopifnot(identical(rownames(vh.d), colnames(vm.d)))
#}}}

#{{{ DESeq2
dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design = ~Genotype)
sizeFactors(dds) = vh$sizeFactor
dds = estimateDispersions(dds, fitType = 'parametric')
disp = dispersions(dds)
dds = nbinomLRT(dds, reduced = ~1)
resultsNames(dds)

fo = file.path(dirw, 'de.rds')
saveRDS(dds, file=fo)

fi = file.path(dirw, 'de.rds')
dds = readRDS(fi)
res = results(dds)

#dds = nbinomWaldTest(dds)
#res1 = results(dds, contrast=c(-1,0,1), pAdjustMethod="fdr")
#res2 = results(dds, contrast=c(-1,1,0), pAdjustMethod="fdr")
#res3 = results(dds, contrast=c(0,1,-1), pAdjustMethod="fdr")
#res4 = results(dds, contrast=c(-.5,1,-.5), pAdjustMethod="fdr")
#stopifnot(rownames(res1) == gids)
#stopifnot(rownames(res2) == gids)
#stopifnot(rownames(res3) == gids)
#stopifnot(rownames(res4) == gids)
## hvm
#dds = DESeqDataSetFromMatrix(countData=hm.d, colData=hh.d, design = ~ Genotype)
#sizeFactors(dds) = rep(1, nrow(hh))
#dds = DESeq(dds, fitType = 'parametric')
#resultsNames(dds)
#res5 = results(dds, contrast=c("Genotype","Hybrid","MidParent"), pAdjustMethod="fdr")
#stopifnot(rownames(res5) == gids)
##
#t_ds = tibble(gid = gids, disp = disp,
            #padj.mb = res1$padj, log2mb = res1$log2FoldChange,
            #padj.hb = res2$padj, log2hb = res2$log2FoldChange,
            #padj.hm = res3$padj, log2hm = res3$log2FoldChange,
            ##padj.fm = res4$padj, log2fm = res4$log2FoldChange
            #padj.fm = res5$padj, log2fm = res5$log2FoldChange
            #) %>%
    #replace_na(list(padj.mb = 1, padj.hb = 1, padj.hm = 1, padj.fm = 1))
##}}}

if(FALSE) {
#{{{ edgeR
y = DGEList(counts = vm.d, group = vh$Genotype)
y = calcNormFactors(y, method = 'TMM') #RLE
t_nf = y$samples %>% as_tibble() %>%
    mutate(SampleID = rownames(y$samples)) %>%
    select(SampleID = SampleID, libSize = lib.size, normFactor = norm.factors)
design = model.matrix(~0 + Genotype, data = vh)
colnames(design) = levels(vh$Genotype)
#y = estimateDisp(y, design)
y = estimateGLMCommonDisp(y, design, verbose = T)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)
t_cpm_merged = cpmByGroup(y) %>% as_tibble() %>%
    transmute(cpm.b = B73, cpm.m = Mo17, cpm.h = BxM)
t_cpm = cpm(y, normalized.lib.sizes = T) %>% as_tibble() %>%
    mutate(gid = gids) %>%
    gather(sid, cpm, -gid) %>% select(sid, gid, cpm)
# mb, hb, hm, fm 
lrt1 = glmLRT(fit, contrast = c(-1, 0, 1))
lrt2 = glmLRT(fit, contrast = c(-1, 1, 0))
lrt3 = glmLRT(fit, contrast = c(0, -1, 1))
lrt4 = glmLRT(fit, contrast = c(-.5, -.5, 1))
stopifnot(identical(gids, rownames(lrt1$table)))
stopifnot(identical(gids, rownames(lrt2$table)))
stopifnot(identical(gids, rownames(lrt3$table)))
stopifnot(identical(gids, rownames(lrt4$table)))
#tags = decideTestsDGE(lrt, adjust.method = "BH", p.value = .05, lfc = 1)
#stopifnot(identical(gids, rownames(tags)))
#{{{ fm
y = DGEList(counts = hm.d, lib.size = hh$libSize, norm.factors = hh$normFactor)
design = model.matrix(~0 + Genotype, data = hh)
colnames(design) = unique(hh$Genotype)
y = estimateGLMCommonDisp(y, design, verbose = T)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)
lrt5 = glmLRT(fit, contrast = c(-1, 1))
stopifnot(identical(gids, rownames(lrt5$table)))
#}}}
tr1 = lrt1$table %>% rownames_to_column("gid") %>% as_tibble() %>% 
    mutate(padj.mb = p.adjust(PValue, method = 'BH')) %>%
    select(gid, log2mb = logFC, padj.mb)
tr2 = lrt2$table %>% 
    mutate(padj.hb = p.adjust(PValue, method = 'BH')) %>%
    select(log2hb = logFC, padj.hb)
tr3 = lrt3$table %>%
    mutate(padj.hm = p.adjust(PValue, method = 'BH')) %>%
    select(log2hm = logFC, padj.hm)
tr4 = lrt4$table %>%
    mutate(padj.fm = p.adjust(PValue, method = 'BH')) %>%
    select(log2fm = logFC, padj.fm)
tr5 = lrt5$table %>%
    mutate(padj.fm = p.adjust(PValue, method = 'BH')) %>%
    select(log2fm = logFC, padj.fm)
t_eg = tr1 %>% bind_cols(tr2) %>% bind_cols(tr3) %>% bind_cols(tr5)
#}}}
}

