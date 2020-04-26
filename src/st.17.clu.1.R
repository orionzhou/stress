source('functions.R')
require(WGCNA)
dirw = file.path(dird, '17_cluster')
enableWGCNAThreads()

fi = file.path(dirw, '../15_de/10.rds')
deg = readRDS(fi)

fi = file.path(dirw, '01.tc.rds')
tc = readRDS(fi)

f_cfg = file.path(dirw, 'clu_options.xlsx')
cfg = read_xlsx(f_cfg)

#cfg %>% mutate(r1 = pmap_lgl(list(cid,cond,opt_deg,opt_clu), run_softPower, deg=!!deg, dirw=!!dirw))
cfg = read_xlsx(f_cfg)
res = cfg %>%# filter(cid >= 'c21') %>%
    mutate(x=pmap(list(cid,cond,drc,opt_deg,opt_clu,optQ,
                       softPower,deepSplit,MEDissThres,minGap),
                  run_wgcna_pipe, tc=!!tc, deg=!!deg, dirw=!!dirw))

fo = sprintf("%s/12.wgcna.rds", dirw)
saveRDS(res, fo)

