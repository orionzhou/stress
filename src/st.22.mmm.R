source('functions.R')
dirw = file.path(dird, '22_mmm')

#{{{ get known cold/heat TF motifs
fi = file.path(dirw, '../06_tf_list/10.stress.tf.tsv')
ti = read_tsv(fi) %>% select(stress,gid=At_gid,name_At)
#
fm = '/home/springer/zhoux379/projects/cre/data/01_tfbs/01.motifs.rds'
tm0 = readRDS(fm)$gmtf %>% select(motif, gid, name)
#
tm = tm0 %>% left_join(ti, by='gid') %>%
    mutate(sname=ifelse(!is.na(stress) & !is.na(name_At), str_c(stress, name_At, sep='_'), NA)) %>%
    mutate(name = ifelse(is.na(sname), name, sprintf("%s(%s)", name, sname))) %>%
    select(-stress,-name_At,-sname) %>%
    group_by(motif) %>%
    summarise(name = str_c(name, collapse=' ')) %>% ungroup()
#}}}

#{{{ extract motif enrichment scores [large-mem]
fg = file.path(dirw, '../21_seq/02.cre.tsv')
tg = read_tsv(fg) %>% select(epi, bin, gid, sid)
tgs = tg %>% distinct(epi, bin, gid) %>% count(epi, bin) %>%
    rename(ng_total = n)

fl = file.path(dirw, '../21_seq/10.rds')
tl = readRDS(fl) #%>%
#    mutate(db_size=map_dbl(tg, get_db_size <- function(x) sum(x$size))) %>%
#    select(-tg)

dirr = '/home/springer/zhoux379/projects/stress/nf/raw'
fi = sprintf("%s/15.fimo.tsv", dirr)
ti = read_tsv(fi, col_names=c('mid','sid','hit','bp','score','eval'))
min_eval = 1e-5
ti2 = ti %>% filter(eval <= min_eval)
tis = ti2 %>%
    inner_join(tg, by='sid') %>% distinct(epi, bin, mid, gid) %>%
    count(epi, bin, mid) %>% rename(ng_hit = n) %>%
    inner_join(tgs, by=c('epi','bin')) %>%
    select(epi, bin, mid, ng_hit, ng_total)

tx = tl %>% #filter(epi=='umr',bin=='-1k') %>%
    unnest(tg) %>% select(-bat,-mid,-size) %>%
    inner_join(ti2, by='sid') %>%
    distinct(lid, ng0, pick, note, bat_mid, bin_epi,epi,bin, ng, mid, gid) %>%
    count(lid, pick, note, bat_mid, bin_epi,epi,bin, ng, mid) %>%
    inner_join(tis, by=c('epi','bin','mid')) %>%
    rename(hitInSample=n, sampleSize=ng, hitInPop=ng_hit, popSize=ng_total) %>%
    mutate(pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail=F)) %>%
    mutate(ratioInSample = sprintf("%d/%d", hitInSample, sampleSize)) %>%
    mutate(ratioInPop = sprintf("%d/%d", hitInPop, popSize)) %>%
    select(-hitInSample,-sampleSize,-hitInPop,-popSize)
tx2 = tx %>% rename(motif=mid) %>% left_join(tm, by='motif') %>%
    select(note,pick,bat_mid,bin_epi,pval=pval.raw,ratioS=ratioInSample,
           ratioP=ratioInPop, motif,name) %>%
    arrange(pval) %>% print(n=30,width=Inf)
tx3 = tx2 %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid) %>% slice(1) %>% ungroup() %>%
    print(n=15, width=Inf)

#{{{ save tx2
to = tx2 %>% select(bat_mid, bin_epi, pick, note, ratioS, ratioP, pval, motif, name)
fo = file.path(dirw, '01.motif.rds')
saveRDS(to, fo)
#}}}
#}}}

fi = file.path(dirw, '01.motif.rds')
mtf_enrich = readRDS(fi)
mtf_enrich1 = mtf_enrich %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid) %>% slice(1) %>% ungroup() %>%
    print(n=15, width=Inf)

#{{{ heatmap all
tp = mtf_enrich1 %>% select(bat_mid, motif) %>%
    inner_join(tx2, by=c('bat_mid','motif')) %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(mtf = str_sub(name, 1, 50)) %>%
    mutate(bat_mid = factor(bat_mid, levels=rev(levels(tl$bat_mid)))) %>%
    mutate(ycol = ifelse(pick, "brown", 'black'),
           ytag = ifelse(pick, "b", 'p'), note=ifelse(pick, note, ''),
           ylab = glue("<{ytag} style='color:{ycol}'>{bat_mid} {note}</{ytag}>"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid)))
tpy = tp %>% distinct(ylab, mtf)
# plot heatmap
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=ylab)) +
    geom_tile(aes(fill=score)) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_text(data=tpy, aes(x=24.5,y=ylab,label=mtf), hjust=0, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,.4)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7.5)) +
    theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/03.pdf", dirw), width = 10, height = 12)
#}}}

#{{{ picked heatmap
tp = tp %>% filter(pick)
tpy = tp %>% distinct(ylab, mtf)
# plot heatmap
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=ylab)) +
    geom_tile(aes(fill=score)) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_text(data=tpy, aes(x=24.5,y=ylab,label=mtf), hjust=0, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,.4)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7.5)) +
    theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/03.picked.pdf", dirw), width = 10, height = 6)
#}}}

#{{{ all enriched motif for +/-2k:raw
tp = mtf_enrich %>% filter(bin_epi=='+/-2k:raw', pick, pval<1e-5) %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(mtf = str_sub(name, 1, 50)) %>%
    mutate(ylab = glue("{bat_mid} {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(ylab, mtf) %>%
    summarise(score = max(score), lab = str_c(lab, collapse=', ')) %>% ungroup()
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=ylab,y=mtf)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/04.picked.all_motif.pdf", dirw), width = 6, height = 12)
#}}}



#{{{ # old
fr = file.path(dird, '../nf/raw/fimo.tsv')
tr = read_tsv(fr, col_names=c("lid",'motif','bp','score')) %>%
    inner_join(tl, by='lid') %>%
    separate(bat, c('stress','drc'), sep="_", remove=F) %>%
    mutate(pscore=score/db_size) %>%
    mutate(pbp=bp/db_size) %>%
    inner_join(tm1, by=c('motif'))
tr %>% distinct(name) %>% arrange(name) %>% pull(name)

names = c('DREB1A', 'DREB1F', 'DREB2A',
          'CRF3', 'ERF5','CRF2','CRF3','ERF011','ERF054',
          'HSFA6A','HSFA1A,HSFA1B','HSFA2,HSFA6B,HSFA7B',
          'HSFA6B','HSFA8','HSFB1','HSFB4,HSFA3', 'HSFC1',
          'HY5')
#{{{
for (name in names) {
tp = tr %>% mutate(score = pscore) %>%
    filter(name == !!name) %>%
    mutate(motif = sprintf("%s %s %s", name, src_type, src_id)) %>%
    mutate(bat_mid = factor(bat_mid, levels=rev(levels(tl$bat_mid))))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=bat_mid)) +
    geom_tile(aes(fill=score)) +
    #geom_text(aes(label=score, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='number genes in cluster',colors=cols100) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(.~motif, nrow=2, dir='v', scale='free_y') +
    otheme(legend.pos='right', legend.dir='v', legend.title=F,
           margin = c(.3,.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7.5)) +
    theme(axis.text.y = element_text(size=7.5)) +
    guides(color=F)
#
fo = sprintf("%s/03.%s.pdf", dirw, name)
p %>% ggexport(filename = fo, width = 12, height = 9)
}
#}}}

name='HSFA1A,HSFA1B'
name='DREB1F'
bat='heat_up'; mid='m18'
tm1 %>% filter(name == !!name)
tl %>% filter(bat==!!bat, mid==!!mid)

mid = 'M06948_2.00'
mid = 'M06564_2.00'
dirr = '/home/springer/zhoux379/projects/stress/nf/raw'
fo = sprintf("%s/12_fimo_sum/%s.tsv", dirr, mid)
to = read_tsv(fo, col_names=c('sid','hit','bp','score','eval'))

to1 = to %>% filter(eval <= 1e-5)
to2 = to1 %>% inner_join(tg, by='sid') %>% distinct(gid)
bg_prop = nrow(to2) / length(unique(tg$gid))
bg_prop

tx = tl %>% filter(epi=='raw',bin=='-500') %>% unnest(tg) %>% select(-size) %>%
    inner_join(to1, by='sid') %>%
    distinct(lid, bat,mid, ng0, note, bat_mid, bin_epi, ng, gid) %>%
    count(lid, note, bat_mid, bin_epi, ng) %>%
    mutate(prop = n / ng, fc=prop/bg_prop) %>%
    arrange(-prop) %>% print(n=40)
#}}}




