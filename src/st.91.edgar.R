source('functions.R')
dirw = file.path(dird, '91_edgar')
get_ds <- function(condC, condB, dds, gids) {
    #{{{
    res1 = results(dds, contrast = c("condC",condC,condB), pAdjustMethod='fdr')
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}
run_deseq2 <- function(th, tm, base.gt = 'B73') {
    #{{{
    th1 = th %>% mutate(condC = str_c(cond, Genotype, sep='.'))
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    ct = th1 %>% distinct(condC, cond, Genotype) %>% filter(Genotype != base.gt) %>%
        mutate(condB = str_c(cond, base.gt, sep="."))
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
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~condC)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = ct %>% mutate(ds = map2(condC, condB, get_ds, dds = dds, gids = gids))
    res %>% select(cond, condC, condB, ds)
    #}}}
}

yid = 'rn17b'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
th = th %>% mutate(lab = sprintf("%s_%s_%d", Genotype, Treatment, Replicate))
th %>% count(Tissue, Genotype, Treatment)

gts = c("B73",'Mo17'); treats=c('cold','control')
th_m1 = th_m %>% filter(Genotype %in% gts, Treatment %in% treats) %>%
    mutate(cond = str_c(Treatment, Genotype, sep='.'))
tm_m1 = tm_m %>% select(SampleID, gid, CPM) %>%
    filter(SampleID %in% th_m1$SampleID) %>%
    inner_join(th_m1[,c("SampleID",'cond')], by='SampleID') %>%
    select(-SampleID) %>% mutate(cond = str_c(cond, 'CPM', sep='.')) %>%
    spread(cond, CPM)

sids = th %>% filter(Genotype %in% c("B73",'Mo17'), Treatment %in% c('cold','control')) %>% pull(SampleID)
th1 = th %>% filter(SampleID %in% sids) %>% mutate(cond = Treatment)

x = run_deseq2(th1, tm)

x1 = x %>% select(cond, ds) %>% unnest(ds) %>%
    pivot_wider(names_from=cond, values_from=c(padj, log2fc), names_sep='.')
x2 = x1 %>% inner_join(tm_m1, by='gid') %>%
    select(gid, control.B73.CPM, control.Mo17.CPM,
           control.padj=padj.control, control.log2fc=log2fc.control,
           cold.B73.CPM, cold.Mo17.CPM,
           cold.padj=padj.cold, cold.log2fc=log2fc.cold)

fo = file.path(dirw, '01.DE.BvM.cold.control.tsv')
write_tsv(x2, fo)

