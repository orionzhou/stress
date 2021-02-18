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


