source('functions.R')
dirw = file.path(dird, '17_cluster')

yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m


th1 = th %>% filter(Experiment=='TC') %>% select(SampleID, Genotype, Treatment, Timepoint)
th1b = th1 %>% filter(Timepoint==0) %>% mutate(Treatment='Cold')
th1c = th1 %>% filter(Timepoint==0) %>% mutate(Treatment='Heat')
th2 = rbind(th1, th1b, th1c) %>% filter(Timepoint != 8)

rc = tm %>% select(gid, SampleID, CPM) %>%
    inner_join(th2, by ='SampleID') %>%
    mutate(has = sprintf("h%03d", Timepoint*10)) %>%
    select(-SampleID, -Timepoint) %>% spread(has, CPM) %>%
    mutate(h015 = if_else(is.na(h015), (h010+h020)/2, h015)) %>%
    mutate(h030 = if_else(is.na(h030), (h020+h040)/2, h015))
#    mutate(h080 = if_else(is.na(h080), h040*.8+h250*.2, h080))
sum(is.na(rc))
rc %>% dplyr::count(Genotype, Treatment)

ra = rc %>% gather(has, rc, -gid, -Genotype, -Treatment) %>%
    spread(Treatment, rc) %>%
    mutate(rac = log2(Cold/Control), rah = log2(Heat/Control)) %>%
    select(gid,Genotype,has,Cold=rac, Heat=rah) %>%
    gather(Treatment, ra, -gid,-Genotype,-has) %>%
    spread(has, ra) %>% arrange(gid, Genotype,Treatment)

rd = rc %>% gather(has, rc, -gid, -Genotype, -Treatment) %>%
    spread(Treatment, rc) %>%
    mutate(rdc = Cold-Control, rdh = Heat-Control) %>%
    select(gid,Genotype,has,Cold=rdc, Heat=rdh) %>%
    gather(Treatment, rd, -gid,-Genotype,-has) %>%
    spread(has, rd) %>% arrange(gid, Genotype,Treatment)

res = list(raw=rc, ratio=ra, diff=rd)
fo = file.path(dirw, '01.tc.rds')
saveRDS(res, fo)

