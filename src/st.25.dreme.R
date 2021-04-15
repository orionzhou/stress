source('functions.R')
require(universalmotif)
dirw = glue('{dird}/25_dreme')

fi = glue("{dird}/21_seq/regions.xlsx")
bins = read_xlsx(fi)$bin
epis = c("raw",'umr','acrE','acrL')
#{{{ read module lists  & write seq lists
tag = 'degB'
tag = 'degA'
tag = 'dmodB'
tag = 'dmodA'
#tag = 'var2'
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)
tls = crossing(bin=bins, epi=epis) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>% arrange(bin, epi)
#
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
#tag = 'var2'
tag = 'degB'
tag = 'degA'
tag = 'dmodA'
tag = 'dmodB'
diri = glue("{dirw}/00_nf")
#
fi = glue("{diri}/{tag}/20.nseq.txt")
tn = read_delim(fi, delim=':', col_names=c('fn','n')) %>%
    mutate(lid = str_replace_all(fn,'\\..*','')) %>% select(lid,n)
#
fi = glue("{diri}/{tag}/23.dreme.rds")
rd = readRDS(fi)
#
fl = glue('{diri}/{tag}/08.tc.tl.rds')
r8 = readRDS(fl)
tc=r8$tc
tl=r8$tl %>% rename(ng0=ng,ng0_c=ng_c) %>%
    inner_join(tn, by='lid') %>% rename(ng=n) %>%
    inner_join(tn, by=c('clid'='lid')) %>% rename(ng_c=n)
#
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)
#
fi = glue("{diri}/{tag}/23.fimo.rds")
fm = readRDS(fi) %>% unnest(loc) %>% distinct(lid,mid,gid)
#
md2 = md %>% select(cid,ts) %>% unnest(ts)
fm2 = fm %>%
    inner_join(tl %>% select(lid,cid), by='lid') %>%
    inner_join(md2, by=c('cid','gid')) %>%
    count(lid,mid,status) %>%
    spread(status,n) %>% rename(pos=`1`,neg=`0`) %>%
    replace_na(list(pos=0,neg=0))
#
tm = fm2 %>% inner_join(tl, by='lid') %>%
    arrange(lid) %>%
    select(cid,lid,mid,pos,ng,neg,ng_c,ng0,ng0_c) %>%
    mutate(pval = phyper(pos-1, pos+neg, ng+ng_c-pos-neg, ng, lower.tail=F)) %>%
    inner_join(rd, by=c('lid','mid')) %>%
    print(n=20)
#
fo = glue("{diri}/{tag}/25.rds")
saveRDS(tm, fo)
#}}}

#{{{ merge clusters, lists, modules
#{{{ read in B
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
#
scope = 'B'
B_tc = tc1 %>% bind_rows(tc2) %>% mutate(scope=scope)
B_tl = tl1 %>% bind_rows(tl2) %>% mutate(scope=scope)
B_tk = tk1 %>% bind_rows(tk2) %>% mutate(scope=scope)
#}}}
#{{{ read in BMW
diri = glue("{dirw}/00_nf")
tag = 'degA'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc1 = r$tc %>% mutate(tag = tag)
tl1 = r$tl %>% mutate(tag = tag)
fi = glue("{diri}/{tag}/25.rds")
tk1 = readRDS(fi) %>% mutate(tag = tag)
#
tag = 'dmodA'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc2 = r$tc %>% mutate(tag = tag)
tl2 = r$tl %>% mutate(tag = tag)
fi = glue("{diri}/{tag}/25.rds")
tk2 = readRDS(fi) %>% mutate(tag = tag)
#
scope = 'BMW'
BMW_tc = tc1 %>% bind_rows(tc2) %>% mutate(scope=scope) %>% select(-me)
BMW_tl = tl1 %>% bind_rows(tl2) %>% mutate(scope=scope)
BMW_tk = tk1 %>% bind_rows(tk2) %>% mutate(scope=scope)
#}}}

#{{{ merge B and BMW
f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg) %>% fill(tag, .direction='down') %>% select(tag,ocid,cid)
#
scopes = c('B','BMW')
tc = rbind(B_tc, BMW_tc) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>%
    mutate(scope = factor(scope, levels=scopes)) %>%
    arrange(scope, cid) %>%
    mutate(bid = str_c('b', str_pad(1:n(), width=2,pad='0'))) %>%
    select(bid, everything()) %>%
    select(bid,scope,cid,cond,note, ng,ng_c,gids,gids_c,ts)
tl = rbind(B_tl, BMW_tl) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>%
    inner_join(tc %>% select(bid,scope,cid), by=c('scope','cid')) %>%
    arrange(bid, lid) %>%
    select(-ocid,-tag,-scope,-cid) %>% select(bid, everything())
tk = rbind(B_tk, BMW_tk) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>%
    inner_join(tc %>% select(bid,scope,cid), by=c('scope','cid')) %>%
    arrange(bid, lid, mid) %>%
    select(-ocid,-tag,-scope,-cid) %>%
    mutate(mid = glue("{bid}_{mid}")) %>%
    mutate(mtf = map2(mtf, mid, rename_mtf)) %>%
    select(bid, everything())
#}}}

fo = glue("{dirw}/02.rds")
r = list(tc=tc, tl=tl, tk=tk)
saveRDS(r, fo)

to = tc %>% filter(scope=='BMW') %>%
    select(cid, cond, note, gids) %>% unnest(gids) %>%
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
tl = r$tl %>% filter(str_detect(bin, "^TSS"))
tm2 = r$tk %>% inner_join(tl %>% select(bid,lid), by=c('bid','lid')) %>%
    filter(pval < .05) %>%
    mutate(ctag='dreme',name=map_chr(mtf,'consensus')) %>%
    select(ctag,mid,name,mtf)
tm = tm1 %>% bind_rows(tm2) %>% mutate(conseq=map_chr(mtf,'consensus'))

x = cluster_motifs(tm, cutHeight1=.1, cutHeight2=.1)

#{{{ tmp hack
x1 = x %>% filter(str_detect(mid,'^B')) %>%
    separate(mid,c('scope','cid','lid','idx'),sep="_", remove=F) %>%
    inner_join(tc %>% select(bid,scope,cid), by=c('scope','cid')) %>%
    mutate(mid1 = str_c(bid,lid,idx, sep='_')) %>%
    select(mid,mid1)
x2 = x %>% left_join(x1, by='mid') %>%
    mutate(mid = ifelse(is.na(mid1), mid, mid1)) %>%
    select(-mid1)
#}}}

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
tc = r$tc
fo = glue("{dirw}/03.mtf.grp.rds")
r3 = list(tc=tc, tl=tl, tk=tk)
saveRDS(r3, fo)
#}}}


#{{{ summerize motifs found in each module/searching parameter - f4b-c
fi = glue("{dirw}/03.mtf.grp.rds")
r3 = readRDS(fi)
#{{{ table stats
tc0 = r3$tc %>% select(bid,scope)
r3$tk %>% inner_join(tc0, by='bid') %>% count(scope)
r3$tk %>% inner_join(tc0, by='bid') %>% distinct(scope,fid) %>% count(scope)
r3$tk %>% inner_join(tc0, by='bid') %>% distinct(scope,fid,known) %>% count(scope,known)
#}}}

#{{{ num. motifs found in each module - f4b
#{{{ f4b
tc = r3$tc %>% filter(scope == 'B')
tc1 = tc %>% select(bid,cid,scope)
tl = r3$tl %>% filter(epi=='raw', !str_detect(bin,'1k')) %>% inner_join(tc1,by='bid')
tk = r3$tk %>% filter(lid %in% tl$lid) %>% inner_join(tc1,by='bid')
#{{{ select plot - prepare plot data
tps = tk %>% group_by(bid,scope,cid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup() %>%
    inner_join(tc %>% select(bid,cond,note,ng), by='bid')
tk1 = tk %>% group_by(bid,lid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup()
tp0 = tl %>% select(scope,lid,bid,cid, bin, epi, ng) %>%
    inner_join(tk1, by=c('bid','lid')) %>% rename(score = n_grp) %>%
    mutate(lab = score)
tpy = tps %>%
    mutate(cond_note = glue("{cond}: {note}")) %>%
    mutate(ylab = glue("{cond}: {note} ({ng})")) %>%
    mutate(ylab1 = glue("({ng})")) %>%
    mutate(ylab2 = glue("({n_grp})")) %>%
    arrange(desc(bid)) %>% mutate(i=1:n())
tpp = tpy %>% group_by(scope) %>% summarise(ymin=min(i), ymax=max(i)) %>% ungroup()
tpy_l = tibble(o=c(0,2,5,8,11)) %>% crossing(tpp) %>% mutate(o=o+ymin)
#
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n()) %>% crossing(tpp)
xmax = max(tpx$x)
tpx1 = tpx %>% group_by(scope,ymax,opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(scope,ymax,opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('scope','bin','epi')) %>%
    inner_join(tpy %>% select(bid,i), by='bid')
tpx_l = tibble(o=seq(0,to=6, by=3)+.5) %>%
    crossing(tpp)
swit = (min(tp$score) + max(tp$score)) / 2
tp = tp %>% mutate(lab.col=score>swit)
#}}}

#{{{ select plot - ps
ps = ggplot(tp, aes(x=x,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=lab.col), hjust=.5, size=2.5) +
    geom_segment(data=tpx_l,aes(x=o,xend=o,y=ymin-.5,yend=ymax+.5),color='blue',size=.5)+
    geom_segment(data=tpy_l,aes(x=.5,xend=xmax+.5,y=o-.5,yend=o-.5),color='blue',size=.5)+
    geom_text(data=tpx,aes(x=x,y=ymin-.6,label=bin2),size=2.5,vjust=1,angle=30,hjust=.8) +
    geom_text(data=tpy,aes(x=.3,y=i,label=ylab),size=2.5,hjust=1) +
    geom_text(data=tpy,aes(x=xmax+.7,y=i,label=ylab2),size=2.5,hjust=0) +
    #geom_text(data=tpy,aes(x=xmax+5,y=i,label=cond_note),size=2.8,hjust=.5) +
    scale_x_continuous(expand=expansion(add=c(5,1.2))) +
    scale_y_continuous(expand=expansion(add=c(.5,.05))) +
    scale_fill_gradientn(name='# motifs',colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    #facet_wrap(scope~., scale='free_y', ncol=2) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           margin = c(1.3,.3,.3,.3), panel.border=F,
           xtick=F, ytick=F, xtitle=F, xtext=F, ytext=F) +
    #theme(legend.position = c(0,0), legend.justification = c(1,0)) +
    #theme(legend.title = element_text(size=10, margin=margin(b=.8,unit="lines"))) +
    theme(strip.text.x=element_text(hjust=.3,size=10,face='bold')) +
    guides(color=F)
#}}}
fo = glue("{dirw}/11.n.mtf.pdf")
ggexport(ps, filename=fo, width=4, height=4)
saveRDS(ps, glue("{dirf}/f4b.rds"))
#}}}

#{{{ ## f4b (old w. B and BMW models)
tc = r3$tc %>% mutate(scope = glue("{scope} model"))
tc1 = tc %>% select(bid,cid,scope)
tl = r3$tl %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
tk = r3$tk %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
#{{{ full plot - prepare plot data
tps = tk %>% group_by(bid,scope,cid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup() %>%
    inner_join(tc %>% select(bid,cond,note,ng), by='bid')
tk1 = tk %>% group_by(bid,lid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup()
tp0 = tl %>% select(scope,lid,bid,cid, bin, epi, ng) %>%
    inner_join(tk1, by=c('bid','lid')) %>% rename(score = n_grp) %>%
    mutate(lab = score)
tpy = tps %>%
    mutate(cond_note = glue("{cond}: {note}")) %>%
    mutate(ylab = glue("{cond}: {note} ({ng})")) %>%
    mutate(ylab1 = glue("({ng})")) %>%
    mutate(ylab2 = glue("({n_grp})")) %>%
    arrange(desc(bid)) %>% mutate(i=1:n())
tpp = tpy %>% group_by(scope) %>% summarise(ymin=min(i), ymax=max(i)) %>% ungroup()
tpy_l = tibble(o=c(0,2,5,8,11)) %>% crossing(tpp) %>% mutate(o=o+ymin)
#
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n()) %>% crossing(tpp)
xmax = max(tpx$x)
tpx1 = tpx %>% group_by(scope,ymax,opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(scope,ymax,opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('scope','bin','epi')) %>%
    inner_join(tpy %>% select(bid,i), by='bid')
tpx_l = tibble(o=seq(0,to=4*length(unique(tp$bin)),by=4)+.5) %>%
    crossing(tpp)
swit = (min(tp$score) + max(tp$score)) / 2
#}}}
#{{{ full plot - pf
pf = ggplot(tp, aes(x=x,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_segment(data=tpx_l,aes(x=o, xend=o,y=ymin-.5,yend=ymax+.5),color='blue',size=.5)+
    geom_segment(data=tpy_l,aes(x=.5,xend=xmax+.5,y=o-.5,yend=o-.5),color='blue',size=.5)+
    geom_text(data=tpx,aes(x=x,y=ymax+.7,label=epi),size=2.5,vjust=0,angle=30,hjust=0) +
    geom_segment(data=tpx1,aes(x=xmin,xend=xmax,y=ymax+1.5,yend=ymax+1.5),size=.5) +
    geom_text(data=tpx1,aes(x=x,y=ymax+2,label=bin2),size=3,vjust=0) +
    #geom_segment(data=tpx2,aes(x=xmin,xend=xmax,y=ymax+2.8,yend=ymax+2.8),size=.5) +
    #geom_text(data=tpx2,aes(x=x,y=ymax+3,label=opt),size=3,vjust=0) +
    geom_text(data=tpy,aes(x=.3,y=i,label=ylab),size=2.5,hjust=1) +
    geom_text(data=tpy,aes(x=xmax+.7,y=i,label=ylab2),size=2.5,hjust=0) +
    scale_x_continuous(expand=expansion(mult=c(.001,.001),add=c(10,2))) +
    scale_y_continuous(expand=expansion(mult=c(.001,.001),add=c(0,.5))) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(scope~., scale='free_y', ncol=1) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,.3,.3,.3), panel.border=F,
           xtick=F, ytick=F, xtitle=F, xtext=F, ytext=F)
#}}}

fo = glue("{dirw}/10.n.mtf.pdf")
ggexport(pf, filename=fo, width=7, height=7)

tc = r3$tc %>% mutate(scope = glue("{scope} model"))
tc1 = tc %>% select(bid,cid,scope)
tl = r3$tl %>% filter(epi=='raw', !str_detect(bin,'1k')) %>% inner_join(tc1,by='bid')
tk = r3$tk %>% filter(lid %in% tl$lid) %>% inner_join(tc1,by='bid')
#{{{ select plot - prepare plot data
tps = tk %>% group_by(bid,scope,cid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup() %>%
    inner_join(tc %>% select(bid,cond,note,ng), by='bid')
tk1 = tk %>% group_by(bid,lid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup()
tp0 = tl %>% select(scope,lid,bid,cid, bin, epi, ng) %>%
    inner_join(tk1, by=c('bid','lid')) %>% rename(score = n_grp) %>%
    mutate(lab = score)
tpy = tps %>%
    mutate(cond_note = glue("{cond}: {note}")) %>%
    mutate(ylab = glue("{cond}: {note} ({ng})")) %>%
    mutate(ylab1 = glue("({ng})")) %>%
    mutate(ylab2 = glue("({n_grp})")) %>%
    arrange(desc(bid)) %>% mutate(i=1:n())
tpp = tpy %>% group_by(scope) %>% summarise(ymin=min(i), ymax=max(i)) %>% ungroup()
tpy_l = tibble(o=c(0,2,5,8,11)) %>% crossing(tpp) %>% mutate(o=o+ymin)
#
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n()) %>% crossing(tpp)
xmax = max(tpx$x)
tpx1 = tpx %>% group_by(scope,ymax,opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(scope,ymax,opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('scope','bin','epi')) %>%
    inner_join(tpy %>% select(bid,i), by='bid')
tpx_l = tibble(o=seq(0,to=6, by=3)+.5) %>%
    crossing(tpp)
swit = (min(tp$score) + max(tp$score)) / 2
#}}}
#{{{ select plot - ps
ps = ggplot(tp, aes(x=x,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_segment(data=tpx_l,aes(x=o,xend=o,y=ymin-.5,yend=ymax+.5),color='blue',size=.5)+
    geom_segment(data=tpy_l,aes(x=.5,xend=xmax+.5,y=o-.5,yend=o-.5),color='blue',size=.5)+
    geom_text(data=tpx,aes(x=x,y=ymax+.6,label=bin2),size=2.5,vjust=0,angle=30,hjust=.2) +
    #geom_segment(data=tpx1,aes(x=xmin,xend=xmax,y=ymax+1.5,yend=ymax+1.5),size=.5) +
    #geom_text(data=tpx1,aes(x=x,y=ymax+2,label=bin2),size=3,vjust=0) +
    geom_text(data=tpy,aes(x=.3,y=i,label=ylab1),size=2.5,hjust=1) +
    geom_text(data=tpy,aes(x=xmax+.7,y=i,label=ylab2),size=2.5,hjust=0) +
    geom_text(data=tpy,aes(x=xmax+5,y=i,label=cond_note),size=2.8,hjust=.5) +
    scale_x_continuous(expand=expansion(mult=c(.001,.001),add=c(1.5,2.5))) +
    scale_y_continuous(expand=expansion(mult=c(.001,.001),add=c(.1,.5))) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(scope~., scale='free_y', ncol=2) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,-6,.3,.3), panel.border=F,
           xtick=F, ytick=F, xtitle=F, xtext=F, ytext=F) +
    theme(strip.text.x=element_text(hjust=.3,size=10,face='bold'))
#}}}

fo = glue("{dirw}/11.n.mtf.s.pdf")
ggexport(ps, filename=fo, width=6, height=5)
#saveRDS(ps, glue("{dirf}/f4b.rds"))
#}}}

#{{{ ## below are old
scopes = c('B','BMW')
scope=scopes[1]
#{{{ scope B: pf1 * ps1
tc = r3$tc %>% filter(scope==!!scope)
tc1 = tc %>% select(bid,cid)
tl = r3$tl %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
tk = r3$tk %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
#{{{ full plot: pf1
#{{{
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
#
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=4*length(unique(tp$bin)),by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
#}}}
#{{{ plot
pf1 = ggplot(tp, aes(x=x,y=i)) +
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
#}}}
#{{{ selected motifs - ps1
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    filter(bin2=='-/+2k') %>% mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=4,by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
ps1 = ggplot(tp, aes(x=x,y=i)) +
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
#}}}
scope=scopes[2]
#{{{ scope BMW: pf2 & ps2
tc = r3$tc %>% filter(scope==!!scope)
tc1 = tc %>% select(bid,cid)
tl = r3$tl %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
tk = r3$tk %>% filter(bid %in% tc1$bid) %>% inner_join(tc1,by='bid')
#{{{ full plot: pf2
#{{{
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
#
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=4*length(unique(tp$bin)),by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
#}}}
#{{{ plot
pf2 = ggplot(tp, aes(x=x,y=i)) +
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
#}}}
#{{{ selected motifs - ps2
tpx = tl %>% distinct(bin,epi) %>% arrange(bin,epi) %>%
    separate(bin, c('opt','bin2'), sep=':', remove=F) %>%
    filter(bin2=='-/+2k') %>% mutate(x = 1:n())
tpx1 = tpx %>% group_by(opt, bin2) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tpx2 = tpx %>% group_by(opt) %>% summarise(xmin=min(x), xmax=max(x), x=(xmin+xmax)/2) %>% ungroup()
tp = tp0 %>% inner_join(tpx, by=c('bin','epi')) %>%
    inner_join(tpy %>% select(cid,i), by='cid')
tpx_l = tibble(o=seq(0,to=4,by=4)+.5)
swit = (min(tp$score) + max(tp$score)) / 2
ps2 = ggplot(tp, aes(x=x,y=i)) +
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
#}}}
#}}}
#}}}

#{{{ top 40 motifs found in each module - prepare data 
tc=r3$tc; tl=r3$tl; tk=r3$tk
r1 = tk %>% select(mid,fid,fname,known, bid,lid,pval) %>%
    inner_join(tl %>% select(lid,bid,bin,epi,ng), by=c('bid','lid')) %>%
    inner_join(tc %>% select(bid,scope,cid,cond,note), by='bid')
tp = r1 %>%
    separate(bin,c('opt','bin'),sep=":") %>%
    arrange(bid,cid, opt, fid, pval) %>%
    group_by(bid,scope,cid, cond, note, ng, opt, fid, fname, known) %>%
    summarise(pval = pval[1]) %>% ungroup() %>%
    arrange(bid,cid, opt, pval) %>%
    group_by(bid,scope,cid, cond, note, opt) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = ifelse(known, fname, ''))
tps = tp %>% group_by(cid,fid) %>%
    summarise(scopes = str_c(sort(scope),collapse=',')) %>% ungroup() %>%
    mutate(bmw_only = scopes=='BMW') %>% select(-scopes)
tp = tp %>% inner_join(tps, by=c('cid','fid')) %>%
    mutate(lab = ifelse(bmw_only, glue("{lab}*"), lab)) %>%
    filter(i <= 40) %>%
    mutate(i = factor(i, levels=50:1))
tpx = tp %>% distinct(bid,scope,cid, cond, ng, note, opt) %>%
    mutate(xlab = glue("{cond}: {note} ({ng})"))
tp1 = tibble(o=cumsum(c(11))+.5)
tp = tp %>% filter(scope=='B')
tpx = tpx %>% filter(scope=='B')
swit = (min(tp$score) + max(tp$score)) / 2
tp = tp %>% mutate(lab.col=score>swit)
#}}}

#{{{ plot f4c
p4 = ggplot(tp, aes(x=cid,y=i)) +
    geom_tile(aes(fill=score), na.rm = F, size=.5, color='white') +
    geom_text(aes(label=lab, color=lab.col), hjust=.5, size=2.5) +
    geom_vline(xintercept=tp1$o, color='blue') +
    scale_x_discrete(breaks=tpx$cid, labels=tpx$xlab, expand=expansion(add=c(3,0))) +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='-log10(Pval)',colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    annotate('text', x=0, y=20, label='Enrichment Rank', size=3, angle=90) +
    annotate('segment', x=0, xend=0, y=17, yend=15, arrow=arrow(length=unit(.1,'inches'), type='closed')) +
    annotate('segment', x=0, xend=0, y=23, yend=25, arrow=arrow(length=unit(.1,'inches'), type='closed')) +
    annotate('text', x=0, y=26, label='Higher', size=2.5) +
    annotate('text', x=0, y=14, label='Lower', size=2.5) +
    otheme(legend.pos='bottom.left', legend.dir='v', legend.title=T,
           margin = c(.3,.3,.3,.3), legend.vjust=-1.7, panel.border=F,
           strip.size=12, strip.compact=F, panel.spacing = .5,
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(legend.position = c(0,0), legend.justification = c(0,0)) +
    theme(legend.title = element_text(size=10, margin=margin(b=.8,unit="lines"))) +
    theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
#}}}
#
fo = glue("{dirw}/12.top.mtf.pdf")
p4 %>% ggexport(filename=fo, width=6, height=8)
saveRDS(p4, glue("{dirf}/f4c.rds"))
#}}}

#{{{ check fold enrichment in all versus umr
#{{{ prepare
t1 = tk %>% select(mid,fid,fname,known, bid,lid,pval,pos,ng,neg,ng_c) %>%
    inner_join(tl %>% select(lid, bid, bin, epi), by=c('bid','lid')) %>%
    arrange(bid,bin,epi,fid, pval) %>%
    group_by(bid,bin,epi,fid,fname,known) %>%
    slice(1) %>% ungroup()
t2 = t1 %>% arrange(bid, pval) %>%
    separate(bin,c('opt','bin'),sep=":", remove=F) %>%
    arrange(bid, pval) %>%
    group_by(bid, fid, fname,known) %>%
    slice(1) %>% ungroup() %>%
    arrange(bid, pval) %>%
    group_by(bid) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    select(bid,i,fid)
t3 = t2 %>% inner_join(t1,by=c('bid','fid')) %>%
    filter(bin == 'TSS:-/+2k', epi %in% c("raw",'umr')) %>%
    select(bid,i,fid,fname,known,epi,pval)
#}}}

tc0  = tc %>% filter(scope=='B') %>%
    mutate(pnl=glue("{cond}:{note}")) %>% select(bid,pnl)
tp = t3 %>%
    inner_join(tc0, by='bid') %>%
    spread(epi,pval) %>%
    filter(!is.na(raw), !is.na(umr)) %>%
    arrange(bid,i) %>% group_by(bid) %>% mutate(i = 1:n()) %>% ungroup() %>%
    gather(epi,pval, -bid,-pnl,-i,-fid,-fname,-known) %>%
    mutate(score = -log10(pval)) %>%
    mutate(txt = number(score, accuracy=1)) %>%
    mutate(lab = ifelse(known, fname, '')) %>%
    select(bid,pnl,i,txt,lab,known,epi,score)
#
tpy = tp %>% distinct(pnl,i,lab)
#{{{ plot
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=epi,y=i)) +
    geom_tile(aes(fill=score), na.rm=F, size=0, color='white') +
    geom_text(aes(label=txt, color=score>swit), hjust=.5, size=2.5) +
    geom_text(data=tpy, aes(x=.4,y=i,label=lab), hjust=1, size=2) +
    #geom_vline(xintercept=tp1$o, color='blue') +
    scale_x_discrete(expand=expansion(add=c(2,.5))) +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='-log10(Pval)',colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(pnl ~ ., ncol=3) +
    otheme(legend.pos='none', legend.dir='v', legend.title=T,
           margin = c(.3,.3,.3,.3), legend.vjust=-1.7,
           strip.compact=F, panel.spacing = .5,
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    guides(color=F)
#}}}
fo = glue("{dirw}/z.pdf")
ggsave(p, file=fo, width=6, height=6)
#}}}

#{{{ pick best motifs for each gene list
fi = glue("{dirw}/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk

t1 = tk %>% select(mid,fid,fname,known, bid,lid,pval,pos,ng,neg,ng_c) %>%
    inner_join(tl %>% select(lid, bid, bin, epi), by=c('bid','lid'))
t2 = t1 %>% arrange(bid, pval) %>%
    separate(bin,c('opt','bin'),sep=":", remove=F) %>%
    arrange(bid, pval) %>%
    group_by(bid, fid, fname,known) %>%
    slice(1) %>% ungroup() %>%
    arrange(bid, pval) %>%
    group_by(bid) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    inner_join(tk %>% select(mid,mtf), by='mid') %>%
    mutate(conseq = map_chr(mtf,'consensus')) %>%
    select(bid, i, opt,bin,epi, pval,mid,fid,fname,known,conseq,mtf,pos,ng,neg,ng_c)
t3 = t2 %>% group_by(bid) %>% nest() %>% rename(mtfs=data) %>% ungroup() %>%
    mutate(n_mtf = map_int(mtfs, nrow)) %>%
    select(bid, n_mtf, mtfs) %>% print(n=30)

r5 = list(tc=tc, tl=tl, tk=t3)
fo = glue("{dirw}/05.best.mtfs.rds")
saveRDS(r5, fo)
tk2 = t3 %>% unnest(mtfs) %>% select(-mtf)
fo = glue("{dirw}/05.best.mtfs.tsv")
write_tsv(tk2, fo)
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


#{{{ read DREME motifs/kmers and save [old]
tag = 'degA'
tag = 'dmodA'
tag = 'dmodB'
tag = 'degB'
#tag = 'var2'
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

