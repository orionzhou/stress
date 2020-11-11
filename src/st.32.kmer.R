source('functions.R')
dirw = glue('{dird}/32_kmer')
#{{{ functions
filter_kmer_loc <- function(gids, gids_c, loc) {
    #{{{
    t1 = loc %>% filter(gid %in% gids) %>% mutate(type = 'cold/heat-responsive')
    t2 = loc %>% filter(gid %in% gids_c) %>% mutate(type = 'control')
    t1 %>% bind_rows(t2)
    #}}}
}
loc_freq <- function(tl, ng, ng_c) {
    #{{{
    to1 = tibble(type=c('cold/heat-responsive','control'), nc=c(ng,ng_c))
    to = tl %>% mutate(pos = round((start + 1 + end)/2)) %>%
        mutate(bin = cut(pos, breaks=seq(0,4000,by=200))) %>%
        mutate(bin = as.numeric(bin)) %>%
        #distinct(gid,type,srd,bin,umr)
        distinct(gid,type,bin,umr)
    tc1 = to %>% count(type, bin) %>% mutate(epi = 'raw')
    tc2 = to %>% filter(umr > 0) %>% count(type, bin) %>% mutate(epi = 'umr')
    #tc3 = to %>% filter(acrL > 0) %>% count(type, srd, bin) %>% mutate(epi = 'acrL')
    #tc4 = to %>% filter(acrE > 0) %>% count(type, srd, bin) %>% mutate(epi = 'acrE')
    rbind(tc1,tc2) %>% complete(type, bin, epi, fill=list(n=0)) %>%
        inner_join(to1, by='type') %>%
        mutate(p=n/nc)
    #}}}
}
loc_ftype <- function(tl) {
    #{{{
    tl %>% mutate(ftype = 'intergenic') %>%
        mutate(ftype = ifelse(intron > 0, 'intron', ftype)) %>%
        mutate(ftype = ifelse(utr3 > 0, 'utr3', ftype)) %>%
        mutate(ftype = ifelse(utr5 > 0, 'utr5', ftype)) %>%
        mutate(ftype = ifelse(cds > 0, 'cds', ftype)) %>%
        filter(umr > 0) %>%
        count(type, ftype) %>%
        group_by(type) %>% mutate(p = n / sum(n)) %>% ungroup()
    #}}}
}
assess_signif <- function(tl, ng, ng_c) {
    #{{{
    nh = tl %>% filter(type!='control') %>% distinct(gid) %>% pull(gid) %>% length()
    nh_c = tl %>% filter(type=='control') %>% distinct(gid) %>% length()
    nhu = tl %>% filter(type!='control',umr>0) %>% distinct(gid) %>% length()
    nhu_c = tl %>% filter(type=='control',umr>0) %>% distinct(gid) %>% length()
    fc = (nh/ng) / (nh_c/ng_c)
    fcu = (nhu/ng) / (nhu_c/ng_c)
    pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)
    pval.u = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F)
    tibble(nh=nh,nh_c=nh_c,nhu=nhu,nhu_c=nhu_c,fc=fc,fcu=fcu,pval=pval,pval.u=pval.u)
    #}}}
}
#}}}

# run kmer finding job for B/M/W

gt = 'Zmays_B73'
gt = 'Zmays_Mo17'
gt = 'Zmays_W22'

#{{{ read in kmer/module info
fr = file.path(dird, '23_mmd/01.kmer.grp.rds')
x = readRDS(fr)
km = x$kmer; mtf=x$mtf; grp=x$grp
#
fm = file.path(dird, "31_promoter", gt, "15.module.rds")
md = readRDS(fm)
mds = md %>% distinct(bat_mid, ng0, note, ctag)
md2 = md %>%
    select(bat_mid, ctag, ng0, ng, tg, ng_c, tg_c) %>%
    mutate(gids = map(tg, 'gid'), gids_c = map(tg_c,'gid')) %>%
    select(-tg, -tg_c, -ng0)
# estabilish kmer-module for M/W
tk = md %>%# filter(bat %in% c("cold_up",'heat_up'), pick) %>%
    select(note, bat, mid, ctag) %>%
    mutate(bat_mid = str_c(bat, mid, sep=":")) %>%
    select(-mid) %>%
    inner_join(km, by=c('bat','note')) %>%
    inner_join(md2, by=c('bat_mid','ctag')) %>%
    select(mid, kmers, bat_mid, note, ctag, ng, ng_c, gids, gids_c)
tk %>% count(bat_mid, note, ctag) %>% spread(ctag, n) %>% print(n=21)
#}}}

#{{{ [run only when gt='B73'] calc motif significance & order in B73
#{{{ read in kmer location
fi = file.path(dirw, gt, 'x.bed')
ti = read_tsv(fi, col_names=c('gid','start','end','kmer','idx','srd','kmer2','idx2','umr'))
identical(ti$kmer,ti$kmer2)
identical(ti$idx,ti$idx2)
kl = ti %>% select(-idx2,-kmer2) %>% group_by(kmer) %>% nest() %>% rename(loc=data)
#}}}

#{{{ assess significance in c2 control
x3 = tk %>% filter(ctag=='c2') %>% select(-ctag) %>%
    unnest(kmers) %>% rename(kmer=kmers) %>%
    inner_join(kl, by='kmer') %>%
    mutate(tl = pmap(list(gids, gids_c, loc), filter_kmer_loc)) %>%
    select(mid, kmer, bat_mid, note, ng, ng_c, tl) %>% unnest(tl) %>%
    group_by(mid, bat_mid, note, ng, ng_c) %>%
    nest() %>% rename(tl=data) %>% ungroup() %>%
    mutate(tp = pmap(list(tl, ng, ng_c), loc_freq))
x4a = x3 %>% select(mid,tl) %>% unnest(tl) %>%
    distinct(mid,gid,type,umr) %>%
    group_by(mid) %>%
    summarise(nh = length(unique(gid[type!='control'])),
              nh_c = length(unique(gid[type=='control'])),
              nhu = length(unique(gid[type!='control' & umr > 0])),
              nhu_c = length(unique(gid[type=='control' & umr >0]))) %>%
    ungroup()
x4 = x3 %>% inner_join(x4a, by='mid') %>%
    mutate(fc = (nh/ng)/(nh_c/ng_c), fcu = (nhu/ng)/(nhu_c/ng_c)) %>%
    mutate(pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)) %>%
    mutate(pval.umr = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F))
#
mer = x4 %>% inner_join(mtf[,c('mid','fid')], by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid')
#
tmb = mer %>%
  arrange(bat_mid, note, fid, pval.umr) %>%
  group_by(bat_mid, note, fid) %>%
  slice(1) %>%
  arrange(bat_mid, note, pval.umr) %>%
  group_by(bat_mid, note) %>% mutate(i = 1:n()) %>% ungroup() %>%
  separate(bat_mid, c('bat','mid0'), sep=":") %>%
  inner_join(km[,c('mid','kmers')], by='mid') %>%
  mutate(bat = factor(bat, levels=bats)) %>%
  arrange(bat, note, i) %>%
  select(bat, note, mid, kmers, fid, fname, i, ng, ng_c, pval, pval.umr)
#}}}

fo = file.path(dirw, '03.best.motifs.rds')
saveRDS(tmb, fo)

tmb2 = tmb %>% mutate(kmers=map_chr(kmers, str_c, collapse=',')) %>%
    select(bat,note,i,mid,fid,fname,ng,ng_c,pval,kmers)
fo = file.path(dirw, '03.best.motifs.tsv')
write_tsv(tmb2, fo)
#}}}

#{{{ prepare for ML input
#{{{ read in kmer location
fi = file.path(dirw, gt, 'x.bed')
#ti = read_tsv(fi, col_names=c('sid','start','end','kmer','idx','srd','kmer2','idx2','umr','acrL','acrE','cds','utr5','utr3','intron'))
ti = read_tsv(fi, col_names=c('gid','start','end','kmer','idx','srd','kmer2','idx2','umr'))
identical(ti$kmer,ti$kmer2)
identical(ti$idx,ti$idx2)
kl = ti %>% select(-idx2,-kmer2) %>% group_by(kmer) %>% nest() %>% rename(loc=data)
#}}}
fmb = file.path(dirw, '03.best.motifs.rds')
tmb = readRDS(fmb)

#{{{ filter kmer loc to 100 best motifs
cnt_mtf <- function(x) length(unique(x$i))
mloc = tmb %>% select(bat, note, mid, fid, i, kmers) %>%
    filter(i <= 100) %>%
    unnest(kmers) %>% rename(kmer=kmers) %>%
    inner_join(kl, by='kmer') %>%
    unnest(loc) %>%
    mutate(pos = round((start+end)/2)) %>%
    select(bat, note, mid, fid, i, gid, kmer, pos, srd, umr) %>%
    group_by(bat, note) %>%
    nest() %>% rename(tl=data) %>% ungroup() %>%
    mutate(mtf = map_int(tl, cnt_mtf))

fo = file.path(dirw, gt, '05.mtf.loc.rds')
saveRDS(mloc, fo)

#{{{ heatmap
tp = x5 %>% filter(nh_c >= 25) %>%
    mutate(score = -log10(pval.umr)) %>%# filter(!is.infinite(score)) %>%
    filter(score >= 10) %>%
    mutate(score2 = -log10(pval)) %>%# filter(!is.infinite(score2)) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} ({ng0}) {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    arrange(ylab, fid, fname, desc(score)) %>%
    group_by(ylab, fid, fname) %>% slice(1) %>% ungroup() %>%
    mutate(lab = sprintf("%.0f %.0f %.1f %.1f", score2, score, fc, fcu)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
#fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
fids = rev(tpm$gid)
tp = tp %>% mutate(fid=factor(fid, levels=fids))
tps = tp %>% distinct(fid, name_s) %>% arrange(fid)
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = file.path(dirw, gt, "05.mtf.grp.enrich.pdf")
p %>% ggexport(filename=fo, width=10, height=80)
#}}}
#}}}



#{{{ second pass to use module-wide best motif
x6 = mer %>% inner_join(mtf[,c('mid','fid')], by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid') %>%
    filter(nh_c >= 25, pval.umr < 5) %>%
    arrange(fid, fname, pval.umr) %>%
    group_by(fid, fname) %>% slice(1) %>% ungroup() %>%
    select(fid, mid)
x = tk %>% filter(!is.na(re)) %>% select(lid, mid) %>%
    inner_join(mtf[,c('mid','fid')], by='mid') %>% rename(mid0=mid) %>%
    inner_join(x6, by='fid') %>%
    inner_join(tk0, by='mid') %>%
    inner_join(tl0[,c('lid','bat_mid')], by='lid') %>%
    distinct(bat_mid, fid, mid, kmers)
#
x2 = x %>%
    inner_join(tln2, by='bat_mid') %>%
    select(fid, mid, kmers, bat_mid, ng, ng_c, sids, sids_c)

x3 = x2 %>% unnest(kmers) %>% rename(kmer=kmers) %>%
    inner_join(kl, by='kmer') %>%
    mutate(tl = pmap(list(sids, sids_c, loc), filter_kmer_loc)) %>%
    select(fid, mid, kmer, bat_mid, ng, ng_c, tl) %>% unnest(tl) %>%
    group_by(fid, mid, bat_mid, ng, ng_c) %>%
    nest() %>% rename(tl=data) %>% ungroup() %>%
    mutate(tp = pmap(list(tl, ng, ng_c), loc_freq)) %>%
    mutate(tf = map(tl, loc_ftype))

x4a = x3 %>% select(bat_mid,mid,tl) %>% unnest(tl) %>%
    distinct(bat_mid,mid,sid,type,umr) %>%
    group_by(bat_mid,mid) %>%
    summarise(nh = length(unique(sid[type!='control'])),
              nh_c = length(unique(sid[type=='control'])),
              nhu = length(unique(sid[type!='control' & umr > 0])),
              nhu_c = length(unique(sid[type=='control' & umr >0]))) %>%
    ungroup()
x4 = x3 %>% inner_join(x4a, by=c('bat_mid','mid')) %>%
    mutate(fc = (nh/ng)/(nh_c/ng_c), fcu = (nhu/ng)/(nhu_c/ng_c)) %>%
    mutate(pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)) %>%
    mutate(pval.umr = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F)) %>%
    inner_join(grp[,c('fid','fname')], by='fid')
mer = x4

fo = file.path(dirw, gt, '07.mtf.grp.enrich.rds')
saveRDS(mer, fo)

#{{{ heatmap
tp = mer %>% filter(nh_c >= 25) %>%
    mutate(score = -log10(pval.umr)) %>%# filter(!is.infinite(score)) %>%
    filter(score >= 10) %>%
    mutate(score2 = -log10(pval)) %>%# filter(!is.infinite(score2)) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} ({ng0}) {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    arrange(ylab, fid, fname, desc(score)) %>%
    group_by(ylab, fid, fname) %>% slice(1) %>% ungroup() %>%
    mutate(lab = sprintf("%.0f %.0f %.1f %.1f", score2, score, fc, fcu)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
#fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
fids = rev(tpm$gid)
tp = tp %>% mutate(fid=factor(fid, levels=fids))
tps = tp %>% distinct(fid, name_s) %>% arrange(fid)
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/%s/07.mtf.grp.enrich.pdf", dirw, gt)
p %>% ggexport(filename=fo, width=10, height=80)
#}}}
#}}}
#}}}



#{{{ #single bat_mid motif plot
i = 1
tg0 = tl0 %>% filter(bat_mid == bat_mids[i]) %>% pull(tg)
tg0 = tg0[[1]] %>% select(-size)
tp0 = r$mtf %>% inner_join(ts, by=c('bat_mid','bin_epi')) %>%
    filter(bat_mid == bat_mids[i]) %>%
    select(grp,grpname,conseq,cnt)
tp1 = tp0 %>% mutate(nhits = map_int(cnt, nrow)) %>%
    arrange(grp,desc(nhits)) %>% group_by(grp) %>%
    summarise(conseq = conseq[1]) %>% ungroup()
tp = tp0 %>% select(-conseq) %>% inner_join(tp1, by='grp') %>%
    unnest(cnt) %>%
    group_by(grp, grpname, conseq, sid) %>%
    summarise(nseq=max(nseq)) %>% ungroup() %>%
    inner_join(tg0, by='sid') %>%
    mutate(score=nseq) %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(grpname = ifelse(grpname=='', str_c(grp,conseq,sep=' '), grpname)) %>%
    mutate(grpname_s = str_sub(grpname, 1, 15))
#
tps = tp %>% distinct(grpname, grpname_s)
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=grpname,y=gid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=1.5) +
    scale_x_discrete(breaks=tps$grpname, labels=tps$grpname_s, expand=c(0,0), position='top') +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,3.3,.3,.3),
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=25, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/09.%d.pdf", dirw, i), width = 8, height = 12)
#}}}

#{{{ motif group meta plot
fi = file.path(dirw, '07.mtf.grp.enrich.rds')
mer = readRDS(fi)

fids = c(
'g0437',
'g0098'
'g0054',
'g0645',
'g0024',
'g0130',
'g0014',
'g0050',
'g0009',
'g0290',
'g0722',
'g0191',
'g0154'
)
mtf_meta_plot <- function(fid, fname, mer, dirw=dirw) {
# {{{ plot
fname=grp %>% filter(grp == grp1) %>% pull(fname)
tp = mer %>% filter(fid==grp1) %>%
    mutate(pan = sprintf("%s n=%d fc=%.1f fcu=%.1f pval=%.2g pal.umr=%.2g", bat_mid, ng, fc, fcu, pval, pval.umr)) %>%
    select(pan, tp) %>% unnest(tp) #%>%
    #mutate(p = ifelse(srd == '-', -p, p))
#tp1 = tp %>% filter(srd=='+')
#tp2 = tp %>% filter(srd=='-')
types = c("cold/heat-responsive",'control')
tpx = tibble(x=c(.5,10.5,20.5),lab=c('-2kb','TSS','+2kb'))
p = ggplot(tp, aes(x=bin,y=p,col=epi)) +
    geom_line(aes(linetype=type), size=.5, na.rm = F) +
    geom_point(aes(shape=type), size=1, na.rm = F) +
    scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
    scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    #scale_shape(labels=types) +
    #scale_linetype(labels=types) +
    facet_wrap(~pan, scale='free', ncol=2) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           strip.style='white',margin = c(.3,.3,.3,.3),
           xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
    guides(fill=F)
fp = sprintf("%s/%s.%s.pdf", dirw, grp1, fname)
p %>% ggexport(filename=fp, width = 10, height = 10)
#}}}
}
tibble(fid=fids) %>% inner_join(grp, by='fid') %>%
    mutate(p=pmap(list(fid,fname,mtf_meta_plot, mer=mer, dirw=dirw)))
#}}}



#{{{ ##assess significance of all found motifs [obsolete]
#{{{ run st.23.mmd.1.R
fi = file.path(dirr, '23_kmer_count.rds')
kmer_cnt = readRDS(fi)
kmer_hits <- function(kmers, kmer_cnt)
    kmer_cnt %>% filter(kmer %in% kmers) %>% pull(sids) %>% unlist() %>% unique()

fi = file.path(dirw, '01.motif.grp.rds')
res = readRDS(fi)

rename_mid <- function(lid, mid) sprintf("%s_%03d", lid, as.numeric(str_replace(mid,'DREME-','')))
fi = sprintf("%s/dreme_kmer.tsv", dirr)
tk = read_tsv(fi) %>% mutate(mid = map2_chr(lid, mid, rename_mid))
tk0 = tk %>% filter(is.na(re)) %>% select(-re,-seq_rc) %>%
    group_by(mid) %>% summarise(kmers = list(seq)) %>% ungroup()

x = tk %>% filter(!is.na(re)) %>% select(lid, mid) %>%
    inner_join(tk0, by='mid') %>%
    mutate(hits = map(kmers,kmer_hits, kmer_cnt=kmer_cnt))
x1 = tl0 %>% filter(pick) %>% select(bat_mid, bin_epi, lid)

x2 = x %>% inner_join(x1, by='lid') %>% rename(lid0=lid,bin_epi0=bin_epi) %>%
    inner_join(x1, by='bat_mid') %>%
    select(bat_mid, mid, lid0, bin_epi0, lid, bin_epi, hits) %>%
    inner_join(tl0[,c('lid','ng0','ng','tg','ng_c','tg_c')], by='lid') %>%
    mutate(sids = map(tg, 'sid'), sids_c = map(tg_c, 'sid')) %>%
    select(-tg, -tg_c)

x3 = x2 %>%# slice(1:300) %>%
    mutate(sids_p = map2(sids, hits, intersect)) %>%
    mutate(sids_n = map2(sids_c, hits, intersect)) %>%
    mutate(pos = map_int(sids_p, length)) %>%
    mutate(neg = map_int(sids_n, length)) %>%
    mutate(pval = phyper(pos-1, neg, ng_c-neg, ng, lower.tail=F)) %>%
    select(-sids,-sids_c,-hits)
#}}}

fi = file.path(dirw, '11.phyper.rds')
t_hy = readRDS(fi)

fi = file.path(dirw, '01.motif.grp.rds')
tg = readRDS(fi)
mtf = tg$mtf; grp=tg$grp
#
te = t_hy %>% inner_join(mtf[,c('mid','fid')],by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid')

#{{{ all enriched motif: +/-2k:umr
tls = tl0 %>% distinct(bat_mid, note)
tp0 = te %>% filter(bin_epi=='+/-2k:umr') %>%
    mutate(score = -log10(pval)) %>% filter(!is.infinite(score))
tp = tp0 %>%
    filter(score >= 10) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(ylab, fid, fname) %>%
    summarise(score = max(score)) %>% ungroup() %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
tps = tp %>% distinct(fid, name_s) %>% arrange(desc(fid))
tp = tp %>% mutate(fid=factor(fid, levels=fids))
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/14.all.motif.grps.pdf", dirw)
p %>% ggexport(filename=fo, width=7, height=80)
#}}}

grp %>% filter(grp=='g0119')
grp %>% filter(str_detect(mid,'l0790_001')) %>% pull(mids)
grp %>% filter(str_detect(mid,'f137'))
mtf0 = rk$fam %>% filter(fid=='f137') %>% pull(pwm)
mtf1 = mtf %>% filter(mid=='l0790_001') %>% pull(pwm)
mtf0a = rk$mtf %>% filter(mid=='M06942_2.00') %>% pull(pwm)
mtf2 = mtf %>% filter(mid=='l1086_001') %>% pull(pwm)
mtf3 = mtf %>% filter(mid=='l0585_012') %>% pull(pwm)
compare_motifs(c(mtf0,mtf1,mtf2), method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)
compare_motifs(c(mtf3,mtfs_known), compare.to=1, method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)
compare_motifs(c(mtf2,mtfs_all), compare.to=1, method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)

fo=file.path(dirw,'tmp.pdf')
pdf(file=fo, width=6,height=16)
midsd = grp %>% filter(str_detect(mid,'l0790_001')) %>% pull(mids)
mtfsd = mtf %>% filter(str_detect(midsd, mid)) %>% pull(pwm)
view_motifs(c(mtf0,mtf0a,mtfsd), method="PCC",min.mean.ic=0,min.overlap=6)
dev.off()
#}}}



