source('functions.R')
dirw = glue('{dird}/71_share')
setwd(dirw)
tga = read_loci() %>% select(gid,symbol,note)
txdb = load_txdb('Zmays_B73', primary=T)

#{{{ Ruben's gene
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm
#
th1 = th %>% filter(Experiment=='NM') %>%
    select(sid=SampleID, gt=Genotype, cond=Treatment, time=Timepoint) %>%
    mutate(cond = str_to_lower(cond))

gid = 'Zm00001d017584'
gid = 'Zm00001d039542'
rr = tm %>% filter(gid == !!gid) %>%
    select(gid, sid=SampleID, cpm=CPM, fpkm=FPKM) %>%
    inner_join(th1, by ='sid') %>%
    select(-sid)

to = rr %>% mutate(cond=factor(cond,levels=c("control",'cold'))) %>%
    arrange(gid, gt, cond, time) %>%
    mutate(cond=glue("{cond}_{time}h")) %>%
    select(gid, genotype=gt, condition=cond, cpm, fpkm)

fo = glue('{dirw}/01.cpm.tsv')
write_tsv(to, fo)
#}}}

#{{{ Bridget's CPM values
yid = 'zm.rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th1 = th %>% filter(Experiment=='HY') %>%
    group_by(Genotype,Treatment,Timepoint) %>%
    mutate(Replicate = 1:n()) %>% ungroup() %>%
    mutate(cond = glue("{Genotype}_{Timepoint}h_{Treatment}_Rep{Replicate}")) %>%
    select(SampleID, cond)

to = tm %>% filter(SampleID %in% th1$SampleID) %>%
    inner_join(th1, by='SampleID') %>%
    select(gid, cond, CPM) %>%
    spread(cond, CPM)

fo = glue("{dird}/71_share/02.hy.cpm.tsv.gz")
write_tsv(to, fo)
#}}}

