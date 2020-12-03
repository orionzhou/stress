source('functions.R')
require(universalmotif)
dirw = glue('{dird}/25_dreme')

fi = glue("{dird}/21_seq/regions.xlsx")
bins = read_xlsx(fi)$bin
epis = c("raw",'umr','acrE','acrL')
#{{{ read module lists  & write seq lists
tag = 'degA'
tag = 'degB'
tag = 'dmodA'
tag = 'dmodB'
tag = 'var2'
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)
tls = crossing(bin=bins, epi=epis) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>% arrange(bin, epi)

tcp = md %>% select(cid, gids)
tcn = md %>% select(gids=gids_c) %>% distinct(gids) %>% mutate(cid=glue("cc{1:n()}"))
tcl = tcp %>% bind_rows(tcn)
tc = md
#
tlp = tcp %>% crossing(tls) %>%
    arrange(cid, bin, epi) %>%
    mutate(lid = sprintf("l%03d", 1:n())) %>%
    select(lid,cid,bin,epi,gids)
tln = tcn %>% crossing(tls) %>%
    arrange(cid, bin, epi) %>%
    mutate(lid = sprintf("cl%03d", 1:n())) %>%
    select(clid=lid,ccid=cid,bin,epi,gids_c=gids)
tl = tc %>% inner_join(tlp, by=c('cid','gids')) %>%
    inner_join(tln, by=c('bin','epi','gids_c')) %>%
    select(lid,cid,cond,note,bin,epi,clid,ccid,ng,ng_c)
#
#{{{ write seq lists and fas
diro = glue('{dirw}/00_nf/{tag}')
if(!dir.exists(diro)) dir.create(diro)
fo = glue("{diro}/01.tc.tsv")
write_tsv(tc %>% select(cid,cond,note,ng,ng_c), fo)
#
diro1 = glue('{dirw}/00_nf/{tag}/02_gene_lists')
if(!dir.exists(diro1)) dir.create(diro1)
tcl %>% mutate(fo = glue("{diro1}/{cid}.txt")) %>% mutate(j=map2(gids, fo, write))
diro1 = glue('{dirw}/00_nf/{tag}/03_gene_status')
if(!dir.exists(diro1)) dir.create(diro1)
tc %>% mutate(fo = glue("{diro1}/{cid}.tsv")) %>% mutate(j=map2(ts, fo, write_tsv))
#
fo = glue("{diro}/05.tl.tsv")
write_tsv(tl, fo)
r = list(tc=tc, tl=tl)
fo = glue("{diro}/08.tc.tl.rds")
saveRDS(r, fo)
#}}}

#{{{ [old] WGCNA-based seq list pairs
tag = 'wgcna'
fi = glue("{dird}/17_cluster/65.modules.ctrl.rds")
i_tibble <- function(gids) tibble(gid=gids)
mdp = md %>% select(cid, gid=gids)
mdn = md %>% select(cond,gid=gids_c) %>% distinct(cond,gid)
#tls = tls0 %>% filter(str_detect(bin, "\\+/\\-"),!epi %in% c('acrL','acrE'))
tls = tls0 %>% filter(!epi %in% c('acrL','acrE'))
#}}}
#}}}


#{{{ read DREME motifs/kmers and save
tag = 'degA'
tag = 'degB'
tag = 'dmodA'
tag = 'dmodB'
tag = 'var2'
diri = glue("{dirw}/00_nf")
fk = glue("{diri}/{tag}/23.kmer.tsv")
fm = glue("{diri}/{tag}/23.kmer.motif.tsv")
tkk = read_tsv(fk) %>% select(-seq_rc) %>% rename(kmer=seq)
tkm = read_tsv(fm) %>% select(-seq_rc) %>% rename(kmer=seq)
make_kmer_motif <- function(kmer) create_motif(kmer, name=kmer, type='PPM', alphabet='DNA')
kmers = tkk %>% distinct(kmer) %>% mutate(mtf = map(kmer, make_kmer_motif))
tk0 = tkk %>% group_by(mid) %>% summarise(kmers = list(kmer)) %>% ungroup()
tk = tkm %>% inner_join(tk0, by='mid')
#
fl = glue('{diri}/{tag}/08.tc.tl.rds')
r8 = readRDS(fl)
tc=r8$tc; tl=r8$tl
#
fi = glue("{diri}/{tag}/23.dreme.rds")
rd = readRDS(fi)
#
km = tl %>% select(lid, cid, ng,ng_c=ng_c) %>%
    inner_join(tk, by='lid') %>%
    inner_join(rd, by=c('lid','mid')) %>%
    select(lid,cid,mid,kmer,kmers,pos,pos_c,neg,neg_c,ng,ng_c,pval,mtf)

fo = glue("{diri}/{tag}/25.rds")
saveRDS(km, fo)
#}}}

#{{{ merge clusters, lists, modules
#{{{ read in
diri = glue("{dirw}/00_nf")
tag = 'degB'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc1 = r$tc %>% mutate(tag = tag)
tl1 = r$tl %>% mutate(tag = tag)
fi = glue("{diri}/{tag}/25.rds")
tk1 = readRDS(fi) %>% mutate(tag = tag)
#
tag = 'dmodB'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc2 = r$tc %>% mutate(tag = tag)
tl2 = r$tl %>% mutate(tag = tag)
fi = glue("{diri}/{tag}/25.rds")
tk2 = readRDS(fi) %>% mutate(tag = tag)
#}}}

f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg) %>% fill(tag, .direction='down') %>% select(tag,ocid,cid)

tc = tc1 %>% bind_rows(tc2) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>% arrange(cid) %>%
    select(cid,cond,note, ng,ng_c,gids,gids_c,ts)
tl = tl1 %>% bind_rows(tl2) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>% arrange(cid, lid) %>%
    select(-ocid,-tag) %>% select(cid, everything())
tk = tk1 %>% bind_rows(tk2) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>% arrange(cid, lid) %>%
    select(-ocid,-tag) %>% mutate(mid = glue("{cid}_{mid}")) %>%
    mutate(mtf = map2(mtf, mid, rename_mtf)) %>%
    select(lid, cid, everything())

fo = glue("{dirw}/02.rds")
r = list(tc=tc, tl=tl, tk=tk)
saveRDS(r, fo)

to = tc %>% select(cid, cond, note, gids) %>% unnest(gids) %>%
    filter(str_detect(gids, '^B73_')) %>%
    mutate(gids=str_replace(gids, 'B73_','')) %>%
    group_by(cid, cond, note) %>%
    summarise(ng = n(), gids=str_c(gids, collapse=',')) %>% ungroup()
fo = glue("{dird}/17_cluster/25.modules.B73.tsv")
write_tsv(to, fo)
#}}}

#{{{ collapse DREME motifs w. known cisbp motifs
diri = '~/projects/cre/data/01_tfbs'
fi = glue('{diri}/10.fam.rds')
fi = glue('{diri}/05.motifs.rds')
rk = readRDS(fi)
tm1 = rk %>% select(ctag, mid, name, mtf=pwm)
#
fi = glue("{dirw}/02.rds")
r = readRDS(fi)
tm2 = r$tk %>% mutate(ctag='dreme',name=map_chr(mtf,'consensus')) %>%
    select(ctag,mid,name,mtf)
tm = tm1 %>% bind_rows(tm2) %>% mutate(conseq=map_chr(mtf,'consensus'))

x = cluster_motifs(tm, cutHeight1=.1, cutHeight2=.1)

ctags = c('maize','cisbp','dreme')
mtf = tm %>% inner_join(x, by='mid') %>%
    mutate(ctag = factor(ctag, levels=ctags)) %>%
    select(ctag,mid,name,conseq,icscore,grp,mtf)
grp = mtf %>%
    arrange(grp, ctag, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mid = mid[1], name=name[1], n_mid = n(),
              n_known = sum(ctag %in% c('cisbp','maize')),
              n_dreme = sum(ctag=='dreme'),
              mids = list(mid),
              conseq = conseq[1], mtf=mtf[1]) %>%
    ungroup() %>%
    arrange(desc(n_known), grp) %>%
    mutate(fid = sprintf("f%04d", 1:n())) %>%
    mutate(mtf = map2(mtf, fid, rename_mtf)) %>%
    mutate(known = n_known > 0) %>%
    select(fid, fname=name, mid, n_mid, n_known, n_dreme, known, conseq, mtf, grp)
mtf = mtf %>% inner_join(grp %>% select(grp,fid), by='grp') %>% select(-grp)
grp = grp %>% select(-grp)
#
grp %>% filter(n_known>0)
grp %>% filter(n_known==0) %>% count(n_dreme)

tk = r$tk %>% inner_join(mtf %>% select(mid,fid), by='mid') %>%
    inner_join(grp %>% select(fid,fname,known), by='fid')
tl = r$tl; tc = r$tc
fo = glue("{dirw}/03.mtf.grp.rds")
r3 = list(tc=tc, tl=tl, tk=tk)
saveRDS(r3, fo)
#}}}


#{{{ summerize kmers found in each module/searching parameter
fi = glue("{dirw}/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk

#{{{ num. motifs found in each module
tps = tk %>% group_by(cid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup() %>%
    inner_join(tc %>% select(cid,cond,note,ng), by='cid')
tk1 = tk %>% group_by(cid,lid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup()
tp0 = tl %>% select(lid,cid, bin, epi, ng) %>%
    inner_join(tk1, by=c('cid','lid')) %>% rename(score = n_grp) %>%
    mutate(lab = score)
tpy = tps %>%
    mutate(ylab = glue("{cond}: {note} ({ng})")) %>%
    arrange(desc(cid)) %>% mutate(i=1:n())
tpy_l = tibble(o=cumsum(c(0,11))+.5)
tpy_l2 = tibble(o=c(2,5,8)+.5)
ymax = max(tpy$i)

#{{{ all motifs
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=72,by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
pa = ggplot(tp, aes(x=x,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    #scale_x_continuous(breaks=tpx$x, labels=tpx$epi, expand=expansion(mult=c(0,0)), position='top') +
    scale_x_continuous(expand=expansion(mult=c(.001,.005))) +
    scale_y_continuous(breaks=tpy$i, labels=tpy$ylab,
                       expand=c(0,0), limits=c(.4,ymax+3.5),
            sec.axis=sec_axis(~., breaks=as.numeric(tpy$i), labels=tpy$n_grp)) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,.6,.3,.3), panel.border=F,
           xtick=F, ytick=F, xtitle=F, xtext=F, ytext=T) +
    annotate(geom='text',x=tpx$x,y=ymax+.7,label=tpx$epi,size=2.5,vjust=0,angle=30,hjust=0) +
    annotate(geom='segment',x=tpx1$xmin,xend=tpx1$xmax,y=ymax+1.5,yend=ymax+1.5,size=.5) +
    annotate(geom='text',x=tpx1$x,y=ymax+2,label=tpx1$bin2,size=3,vjust=0) +
    annotate(geom='segment',x=tpx2$xmin,xend=tpx2$xmax,y=ymax+2.8,yend=ymax+2.8,size=.5) +
    annotate(geom='text',x=tpx2$x,y=ymax+3,label=tpx2$opt,size=3,vjust=0) +
    annotate(geom='segment',x=.5,xend=max(tp$x)+.5,y=tpy_l$o,yend=tpy_l$o, color='blue', size=.5) +
    annotate(geom='segment',x=.5,xend=max(tp$x)+.5,y=tpy_l2$o,yend=tpy_l2$o, color='black', size=.3) +
    annotate(geom='segment',x=tpx_l$o, xend=tpx_l$o,y=.5,yend=ymax+.5, color='blue', size=.5)
#}}}
ggsave(pa, file=glue("{dirw}/10.n.mtf.pdf"), width=14, height=4.5)
ggsave(pa, file=glue("{dirf}/sf07.pdf"), width=14, height=4.5)

#{{{ selected motifs
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    filter(bin2=='-/+2k') %>% mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=8,by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
pb = ggplot(tp, aes(x=x,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    #scale_x_continuous(breaks=tpx$x, labels=tpx$epi, expand=expansion(mult=c(0,0)), position='top') +
    scale_x_continuous(expand=expansion(mult=c(.001,.005))) +
    scale_y_continuous(breaks=tpy$i, labels=tpy$ylab,
                       expand=c(0,0), limits=c(.4,ymax+3.5),
            sec.axis=sec_axis(~., breaks=as.numeric(tpy$i), labels=tpy$n_grp)) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,.6,.3,.3), panel.border=F,
           xtick=F, ytick=F, xtitle=F, xtext=F, ytext=T) +
    annotate(geom='text',x=tpx$x,y=ymax+.7,label=tpx$epi,size=2.5,vjust=0,angle=30,hjust=0) +
    annotate(geom='segment',x=tpx1$xmin,xend=tpx1$xmax,y=ymax+1.5,yend=ymax+1.5,size=.5) +
    annotate(geom='text',x=tpx1$x,y=ymax+2,label=tpx1$bin2,size=3,vjust=0) +
    annotate(geom='segment',x=tpx2$xmin,xend=tpx2$xmax,y=ymax+2.8,yend=ymax+2.8,size=.5) +
    annotate(geom='text',x=tpx2$x,y=ymax+3,label=tpx2$opt,size=3,vjust=0) +
    annotate(geom='segment',x=.5,xend=max(tp$x)+.5,y=tpy_l$o,yend=tpy_l$o, color='blue', size=.5) +
    annotate(geom='segment',x=.5,xend=max(tp$x)+.5,y=tpy_l2$o,yend=tpy_l2$o, color='black', size=.3) +
    annotate(geom='segment',x=tpx_l$o, xend=tpx_l$o,y=.5,yend=ymax+.5, color='blue', size=.5)
#}}}
ggsave(pb, file=glue("{dirw}/10.n.mtf.2.pdf"), width=4.5, height=4.5)
saveRDS(pb, glue("{dirf}/f3b.rds"))
#}}}

#{{{ top 40 motifs found in each module - fig S8
r1 = tk %>% select(mid,fid,fname,known, lid,pval) %>%
    inner_join(tl %>% select(lid,cid,bin,epi,ng), by='lid') %>%
    inner_join(tc %>% select(cid, cond, note), by='cid')
tp = r1 %>% arrange(cid, pval) %>%
    separate(bin,c('opt','bin'),sep=":") %>%
    group_by(cid, cond, note, ng, opt, fid, fname, known) %>%
    summarise(pval = min(pval)) %>% ungroup() %>%
    arrange(cid, opt, pval) %>%
    group_by(cid, cond, note, opt) %>%
    slice(1:40) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = ifelse(known, fname, '')) %>%
    mutate(i = factor(i, levels=50:1))
tps = tp %>% group_by(cid,fid) %>%
    summarise(opts = str_c(sort(opt),collapse=',')) %>% ungroup() %>%
    mutate(tts_only = opts=='TTS') %>% select(-opts)
tp = tp %>% inner_join(tps, by=c('cid','fid')) %>%
    mutate(lab = ifelse(tts_only, glue("{lab}*"), lab))
tpx = tp %>% distinct(cid, cond, ng, note, opt) %>%
    mutate(xlab = glue("{cond}: {note} ({ng})"))
tp1 = tibble(o=cumsum(c(11))+.5)

#{{{ plot
swit = (min(tp$score) + max(tp$score)) / 2
p4 = ggplot(tp, aes(x=cid,y=i)) +
    geom_tile(aes(fill=score), na.rm = F, size=.5, color='white') +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_vline(xintercept=tp1$o, color='blue') +
    scale_x_discrete(breaks=tpx$cid, labels=tpx$xlab, expand=c(0,0), position='top') +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='-log10(Pval)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(opt ~ ., nrow=2) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           margin = c(.3,4.3,.3,.3), legend.vjust=-1.7,
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(legend.position = c(1,0), legend.justification = c(0,0)) +
    theme(axis.text.x = element_text(angle=25, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
#}}}
p4 %>% ggexport(filename=glue("{dirw}/11.top.mtf.pdf"), width=10, height=10)
p4 %>% ggexport(filename=glue("{dirf}/sf08.pdf"), width=10, height=7)
#}}}
#}}}

#{{{ pick best motifs for each gene list
fi = glue("{dirw}/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk

t1 = tk %>% select(mid,fid,fname,known, cid,lid,pval) %>%
    inner_join(tl %>% select(lid, cid, bin, epi), by=c('cid','lid'))
t2 = t1 %>% arrange(cid, pval) %>%
    separate(bin,c('opt','bin'),sep=":", remove=F) %>%
    arrange(cid, pval) %>%
    group_by(cid, fid, fname) %>%
    slice(1) %>% ungroup() %>%
    arrange(cid, pval) %>%
    group_by(cid) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    inner_join(tk %>% select(mid,kmer,kmers), by='mid') %>%
    mutate(kmers = map_chr(kmers, str_c, collapse=',')) %>%
    select(cid, i, opt,bin,epi, pval,fid,fname,kmers)
t3 = t2 %>% group_by(cid) %>% nest() %>% rename(kmer=data) %>% ungroup() %>%
    mutate(n_mtf = map_int(kmer, nrow)) %>%
    select(cid, n_mtf, kmer)

r5 = list(tc=tc, tl=tl, tk=t3)
fo = glue("{dirw}/05.best.mtfs.rds")
saveRDS(r5, fo)
#}}}


### check enrichment of kmers in variable gene lists ###
#{{{
tag = 'var2'
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)
#
fi = glue("{dird}/41_ml/06.tk.tc.rds")
r6 = readRDS(fi)
km = r6$tk %>% select(kmer) %>% unnest(kmer) %>% distinct(fid,fname)

read_mtf_status <- function(fi) read_tsv(fi) %>% gather(fid,mcp,-gid,-status)
ts = md %>% select(cid,cond,note) %>%
    mutate(fi = glue("{dirw}/31_mtf_status/{cid}.tsv")) %>%
    mutate(ts = map(fi, read_mtf_status)) %>%
    unnest(ts) %>% select(-fi)

ts2 = ts %>% mutate(status = ifelse(status==1, 'p', 'n')) %>%
    mutate(m_status = ifelse(mcp>0, 1, 0)) %>%
    mutate(st = str_c(status, m_status, sep='')) %>%
    count(cid,cond,note,fid,st) %>%
    spread(st, n)
ts3 = ts2 %>%
    separate(fid, c('fid','suf'), sep="_") %>% select(-suf) %>%
    left_join(km, by='fid') %>%
    mutate(rate.p = p1/(p0+p1), rate.n = n1/(n0+n1)) %>%
    mutate(fc = rate.p / rate.n) %>%
    mutate(pval = phyper(p1-1, n1+p1, n0+p0, p0+p1, lower.tail=F)) %>%
    select(cid,cond,note,fname,n0,n1,p0,p1,rate.p,rate.n,pval) %>%
    arrange(cid,pval) %>%
    group_by(cid) %>% slice(1:10) %>% ungroup() %>% print(n=50)
#}}}



