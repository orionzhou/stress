source('functions.R')
dirw = file.path(dird, '91_edgar')


yid = 'rn17b'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
th = th %>% mutate(lab = sprintf("%s_%s_%d", Genotype, Treatment, Replicate))
th %>% dplyr::count(Tissue, Genotype, Treatment)

gts = c("B73",'Mo17'); treats=c('cold','control')
th_m1 = th_m %>% filter(Genotype %in% gts, Treatment %in% treats) %>%
    mutate(cond = str_c(Genotype, Treatment, sep='.'))
tm_m1 = tm_m %>% select(SampleID, gid, CPM) %>%
    filter(SampleID %in% th_m1$SampleID) %>%
    inner_join(th_m1[,c("SampleID",'cond')], by='SampleID') %>%
    select(-SampleID) %>% mutate(cond = str_c(cond, 'CPM', sep='.')) %>%
    spread(cond, CPM)

sids = th %>% filter(Genotype %in% gts, Treatment %in% treats) %>% pull(SampleID)
th1 = th %>% filter(SampleID %in% sids) %>%
    dplyr::select(SampleID, cond1=Genotype, cond2=Treatment)

x = run_deseq2(th1, tm, base.cond='control')

x1 = x %>% select(cond1, ds) %>% unnest(ds) %>%
    pivot_wider(names_from=cond1, values_from=c(padj, log2fc), names_sep='.')
x2 = x1 %>% inner_join(tm_m1, by='gid') %>%
    select(gid, B73.control.CPM, B73.cold.CPM,
           B73.padj=padj.B73, B73.log2fc=log2fc.B73,
           Mo17.control.CPM, Mo17.cold.CPM,
           Mo17.padj=padj.Mo17, Mo17.log2fc=log2fc.Mo17)

fo = file.path(dirw, '01.DE.BvM.cold.control.tsv')
write_tsv(x2, fo)

