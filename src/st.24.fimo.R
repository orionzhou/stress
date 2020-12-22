source('functions.R')
require(universalmotif)
dirw = file.path(dird, '24_fimo')
fm = "~/projects/cre/data/01_tfbs/10.fam.rds"
mtf_fam = readRDS(fm)
tm = mtf_fam$fam %>% select(fid,conseq,fname)

#{{{ get known cold/heat TF motifs
fi = file.path(dirw, '../06_tf_list/10.stress.tf.tsv')
ti = read_tsv(fi) %>% select(stress,gid=At_gid,name_At)

tm2 = mtf_fam$gmtf %>% select(mid, gid, name)
    left_join(ti, by='gid') %>%
    mutate(sname=ifelse(!is.na(stress) & !is.na(name_At), str_c(stress, name_At, sep='_'), NA)) %>%
    mutate(name = ifelse(is.na(sname), name, sprintf("%s(%s)", name, sname))) %>%
    select(-stress,-name_At,-sname) %>%
    group_by(mid) %>%
    summarise(name = str_c(name, collapse=' ')) %>% ungroup()
#}}}

#{{{ extract motif enrichment scores [large-mem]
fg = file.path(dirw, '../21_seq/05.cre.loc.tsv')
tg = read_tsv(fg, col_types='ccccii') %>% select(epi, bin, gid, sid)
tgs = tg %>% distinct(epi, bin, gid) %>% count(epi, bin) %>%
    rename(ng_total = n)
#
fl = file.path(dirw, '../21_seq/15.rds')
tl = readRDS(fl)

fi = sprintf("%s/fimo.rds", dirr)
fimo = readRDS(fi)

min_eval = 1e-4
ti2 = fimo %>% filter(!is.na(cnt)) %>% select(cnt) %>% unnest(cnt) %>%
    filter(pval <= min_eval) %>% rename(fid=mid)
tis = ti2 %>%
    inner_join(tg, by='sid') %>% distinct(epi, bin, fid, gid) %>%
    count(epi, bin, fid) %>% rename(ng_hit = n) %>%
    inner_join(tgs, by=c('epi','bin')) %>%
    select(epi, bin, fid, ng_hit, ng_total)

tx = tl %>% #filter(epi=='umr',bin=='-1k') %>%
    unnest(tg) %>% select(-bat,-mid,-size,-cond,-clid) %>%
    inner_join(ti2, by='sid') %>%
    distinct(lid, ng0, pick, note, bat_mid, bin_epi,epi,bin, ng, fid, gid) %>%
    count(lid, pick, note, bat_mid, bin_epi,epi,bin, ng, fid) %>%
    inner_join(tis, by=c('epi','bin','fid')) %>%
    rename(hitInSample=n, sampleSize=ng, hitInPop=ng_hit, popSize=ng_total) %>%
    mutate(pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail=F)) %>%
    mutate(ratioInSample = sprintf("%d/%d", hitInSample, sampleSize)) %>%
    mutate(ratioInPop = sprintf("%d/%d", hitInPop, popSize)) %>%
    select(-hitInSample,-sampleSize,-hitInPop,-popSize)
tx2 = tx %>% left_join(tm, by='fid') %>%
    select(lid,note,pick,bat_mid,bin_epi,pval=pval.raw,ratioS=ratioInSample,
           ratioP=ratioInPop, fid,fname) %>%
    arrange(pval) %>% print(n=30,width=Inf)
tx3 = tx2 %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid) %>% slice(1) %>% ungroup() %>%
    print(n=15, width=Inf)

to = tx2 %>% select(lid,bat_mid, bin_epi, pick, note,ratioS,ratioP,pval,fid,fname)
fo = file.path(dirw, '01.motif.enrich.rds')
saveRDS(to, fo)
#}}}

fi = file.path(dirw, '01.motif.enrich.rds')
mtf_enrich = readRDS(fi)

#{{{ ## best1 motif - picked heatmap
mtf_enrich1 = mtf_enrich %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid) %>% slice(1) %>% ungroup() %>%
    print(n=15, width=Inf)
tp = mtf_enrich1 %>% filter(pick, pval<1e-4) %>%
    mutate(score = -log10(pval)) %>%
    mutate(ylab = glue("{bat_mid} {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(bat_mid,bin_epi, ylab, fid, fname) %>%
    summarise(score = max(score)) %>% ungroup() %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(name = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpy = tp %>% distinct(ylab, fname)
# plot heatmap
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=bat_mid)) +
    geom_tile(aes(fill=score)) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_text(data=tpy, aes(x=28.5,y=ylab,label=fname), hjust=0, size=2) +
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


#{{{ all enriched motif: +/-2k:raw
tp = mtf_enrich %>% filter(bin_epi=='+/-2k:raw', pick, pval<1e-4) %>%
    mutate(score = -log10(pval)) %>%
    mutate(ylab = glue("{bat_mid} {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(ylab, fid, fname) %>%
    summarise(score = max(score)) %>% ungroup() %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tps = tp %>% distinct(fid, name_s) %>% arrange(desc(fid))
tp = tp %>% mutate(fid=factor(fid, levels=tps$fid))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/05.all.motif.grps.pdf", dirw)
p %>% ggexport(filename=fo, width=5, height=7)
#}}}

#{{{ case study of known motif enrichment
bat_mids = c("cold_up:m16","cold_up:m18","cold_up:m21","heat_up:m08","heat_up:m18")
tp = mtf_enrich %>% #filter(pval<1e-5) %>%
    #inner_join(mtfs, by='fid') %>%
    mutate(score = -log10(pval)) %>%
    filter((bat_mid %in% bat_mids[c(1,3)] & fid %in% c('f032','f003','f002','f111','f013','f035','f068','f093')) |
           (bat_mid %in% bat_mids[4:5] & fid == 'f137')) %>%
    mutate(ylab = glue("{bat_mid} | {name}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(ylab, bin_epi) %>%
    summarise(score = max(score)) %>% ungroup() %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(ylab_s = str_sub(ylab, 1, 60))
tps = tp %>% distinct(ylab, ylab_s)
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=ylab)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$ylab, labels=tps$ylab_s, expand=c(0,0)) +
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
p %>% ggexport(filename =sprintf("%s/08.case.motif.pdf", dirw), width = 8, height = 7)

bin_epis = c("+/-500:raw","+/-500:umr","+/-1k:raw","+/-1k:umr","+/-2k:raw","+/-2k:umr")
tp1 = tp %>% filter(bin_epi %in% bin_epis)
p = ggplot(tp1, aes(x=ylab_s,y=score)) +
    geom_point(aes(color=bin_epi, shape=bin_epi), size=2, na.rm = F) +
    scale_x_discrete(name='-log10(pvallue)', expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.03,.03))) +
    scale_color_aaas() +
    scale_shape_manual(values=0:6) +
    coord_flip() +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(size=7.5)) +
    theme(axis.text.y = element_text(size=7))
p %>% ggexport(filename =sprintf("%s/09.case.motif.pdf", dirw), width = 8, height = 6)
#}}}

#{{{ meta plot of selected TFBS motifs
bat_mids = c("cold_up:m16","cold_up:m18","cold_up:m21","heat_up:m08","heat_up:m18")
tm0 = mtf_enrich %>% filter(bin_epi=='+/-2k:raw', pick, pval<1e-4) %>%
    filter((bat_mid %in% bat_mids[c(1,3)] & fid %in% c('f032','f003','f002','f111','f013','f035','f068','f093')) |
           (bat_mid %in% bat_mids[4:5] & fid == 'f137'))
tl0 = tm0 %>% select(fid, lid) %>%
    inner_join(tl, by='lid') %>%
    select(fid, lid, bat_mid, bin_epi, note, ng, ng_c, tg, tg_c)
tl0a = tl0 %>% select(fid, lid, ng, tg) %>% unnest(tg) %>% select(-size)
tl0b = tl0 %>% select(fid, lid, ng=ng_c, tg=tg_c) %>% unnest(tg) %>% select(-size)

#{{{ selected groups + background
tp0 = fimo %>% filter(lid %in% tm0$fid) %>% select(raw) %>% unnest(raw) %>%
    rename(fid=mid) %>%
    inner_join(tl0b, by=c('fid','sid')) %>%
    mutate(pos = (start+end)/2) %>%
    arrange(lid, fid, gid, desc(score)) %>%
    group_by(lid, fid, gid) %>% slice(1) %>% ungroup() %>%
    mutate(bin = cut(pos, breaks=seq(0,4000,by=200))) %>%
    mutate(bin = as.numeric(bin)) %>%
    count(lid, fid, bin, srd, ng) %>%
    complete(nesting(lid, fid, srd, ng), bin, fill=list(n=0)) %>%
    mutate(p = n/ng) %>% select(-n, -ng) %>% rename(p2=p)
tp1 = fimo %>% filter(lid %in% tm0$fid) %>% select(raw) %>% unnest(raw) %>%
    rename(fid=mid) %>%
    inner_join(tl0a, by=c('fid','sid')) %>%
    #inner_join(tm0[,c('bat_mid','fid')], by=c('bat_mid','fid')) %>%
    mutate(pos = (start+end)/2) %>%
    arrange(lid, fid, gid, desc(score)) %>%
    group_by(lid, fid, gid) %>% slice(1) %>% ungroup() %>%
    mutate(bin = cut(pos, breaks=seq(0,4000,by=200))) %>%
    mutate(bin = as.numeric(bin)) %>%
    count(lid, fid, bin, srd, ng) %>%
    complete(nesting(lid, fid, srd, ng), bin, fill=list(n=0)) %>%
    mutate(p = n/ng) %>% select(-n, -ng) %>% rename(p1=p)

types = c("cold/heat-responsive",'control')
tp = tp1 %>%
    left_join(tp0, by=c('lid','fid','bin','srd')) %>%
    replace_na(list(p2=0)) %>%
    gather(type, prop, -lid, -fid, -bin, -srd) %>%
    inner_join(tl0[,c('lid','bat_mid')], by='lid') %>%
    inner_join(tm, by='fid') %>%
    mutate(ylab = glue("{bat_mid} | {fname}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    mutate(ylab_s = str_sub(ylab, 1, 30))
tpx = tibble(x=c(.5,10.5,20.5),lab=c('-2kb','TSS','+2kb'))
p = ggplot(tp, aes(x=bin,y=prop,col=srd)) +
    geom_line(aes(linetype=type), size=.5, na.rm = F) +
    geom_point(aes(shape=type), size=1, na.rm = F) +
    scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
    scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    scale_shape(labels=types) +
    scale_linetype(labels=types) +
    facet_wrap(~ylab_s, scale='free', nrow=4) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           strip.style='white',margin = c(.3,.3,.3,.3),
           xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
    guides(fill=F)
p %>% ggexport(filename =sprintf("%s/15.tfbs.metaplot.pdf", dirw), width = 10, height = 10)
#}}}
#}}}


#{{{ # obsolete
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


