source('functions.R')
require(WGCNA)
dirw = file.path(dird, '17_cluster')
enableWGCNAThreads()

#{{{ read DE, TC, config
fi = file.path(dirw, '../15_de/05.rds')
deg = readRDS(fi)$deg48 %>%
    select(Genotype,Treatment,Timepoint,cond2,up,down) %>%
    gather(drc, gids, -Treatment,-Genotype,-Timepoint,-cond2) %>%
    spread(cond2, gids) %>%
    dplyr::rename(gids0 = time0, gids1 = timeM) %>%
    mutate(gids = map2(gids0, gids1, intersect)) %>%
    select(Genotype,Treatment,Timepoint,drc, gids)
#
fi = file.path(dirw, '01.tc.rds')
tc = readRDS(fi)
f_cfg = file.path(dirw, 'config.xlsx')
cfg = read_xlsx(f_cfg)
#}}}

#cfg %>% mutate(r1 = pmap_lgl(list(cid,cond,opt_deg,opt_clu), run_softPower, deg=!!deg, dirw=!!dirw))
cfg = read_xlsx(f_cfg)
res = cfg %>%# filter(cid >= 'c21') %>%
    mutate(x=pmap(list(cid,cond,drc,opt_deg,opt_clu,optQ,
                       softPower,deepSplit,MEDissThres,minGap),
                  run_wgcna_pipe, gt_map=!!gt_map, tc=!!tc, deg=!!deg, dirw=!!dirw))

fo = sprintf("%s/12.wgcna.rds", dirw)
saveRDS(res, fo)

