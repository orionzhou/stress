source('functions.R')
dirw = file.path(dird, '17_cluster')
bats = c('cold_up', 'heat_up', 'cold_down', 'heat_down')

#{{{ prepare TC matrix
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
#}}}
#run st.17.clu.1.R

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

#{{{ test run wgcna
require(WGCNA)
enableWGCNAThreads()

cfg = read_xlsx(f_cfg)
i=9
cid=cfg$cid[i];cond=cfg$cond[i];drc=cfg$drc[i];
opt_deg=cfg$opt_deg[i];opt_clu=cfg$opt_clu[i];optQ=cfg$optQ[i];
softPower=cfg$softPower[i];deepSplit=cfg$deepSplit[i];
MEDissThres=cfg$MEDissThres[i]; minGap=cfg$minGap[i]
r = run_wgcna_pipe(cid,cond,drc,opt_deg,opt_clu,optQ,
    softPower,deepSplit,MEDissThres,minGap, gt_map, tc, deg, dirw)

r$clu %>% filter(str_detect(gid,gid0))
r$me %>% filter(ctag=='raw',clu==14) %>% pull(me)

cfg = read_xlsx(f_cfg)
res = cfg %>% dplyr::filter(cid == 'c09') %>%
    mutate(x=pmap(list(cid,cond,drc,opt_deg,opt_clu,optQ,
                       softPower,deepSplit,MEDissThres,minGap),
                  run_wgcna_pipe, tc=!!tc, deg=!!deg, dirw=!!dirw))
#}}}

#{{{ process wgcna results and save as module
fi = sprintf("%s/12.wgcna.rds", dirw)
res = readRDS(fi) %>% rename(stress=cond) %>%
    mutate(ng=map_int(x, 'ng'), np=map_int(x, 'np')) %>%
    mutate(nm1=map_int(x, 'nm1'), nm2=map_int(x, 'nm2')) %>%
    mutate(me=map(x, 'me'), clu=map(x,'clu')) %>% select(-x)

#{{{ testing cpm_diff w. cpm_ratio / cpm_raw
tx = res %>% filter(cid >= 'c01', cid <= 'c12') %>%
    select(-clu) %>% unnest(me) %>%
    group_by(stress,drc,opt_deg,opt_clu,softPower,deepSplit,MEDissThres,minGap, ng, np) %>%
    nest() %>% rename(me = data)

plot_me <- function(stress,drc,opt_deg,opt_clu,
                    softPower,deepSplit,MEDissThres,minGap, ng, np, me, dirw) {
    #{{{
    tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
        number(ng), opt_deg, stress, drc, number(np), opt_clu)
    labx = ifelse(minGap==0, 'auto', minGap)
    lab = sprintf("deepSplit=%g minGap=%s", deepSplit, labx)
    #
    me0 = me %>% filter(ctag == 'merged') %>% select(-ctag) %>%
        arrange(cid, optQ, nm1, nm2, n, clu) %>%
        group_by(cid,optQ,nm1,nm2,n,clu) %>%
        mutate(val = scale_by_first(val)) %>%
        ungroup() %>%
        mutate(val = ifelse(val > 1, 1, val)) %>%
        mutate(val = ifelse(val < -1, -1, val))
    tpy0 = me0 %>% distinct(optQ, clu, n)
    tpy = me0 %>% select(optQ,clu,cond,val) %>%
        spread(cond, val) %>% rename(id = clu) %>%
        group_by(optQ) %>% nest() %>%
        mutate(clu = map(data, order_id_by_hclust)) %>%
        select(-data) %>% unnest(clu) %>%
        mutate(y = str_c(optQ, clu, sep="_")) %>%
        inner_join(tpy0, by=c("optQ",'clu')) %>%
        mutate(lab = sprintf("%s (%d)", clu, n))
    optQs = c('raw','diff','ratio')
    tp = me0 %>%
        mutate(y = str_c(optQ, clu, sep="_")) %>%
        mutate(y = factor(y, levels=rev(tpy$y))) %>%
        mutate(optQ = factor(optQ, levels=optQs)) %>%
        mutate(pan=sprintf("cpm [%s]: %s raw -> %d merged", optQ, nm1, nm2))
    pans = tp %>% distinct(optQ, pan) %>% arrange(optQ) %>% pull(pan)
    tp = tp %>% mutate(pan=factor(pan, levels=pans))
    p = ggplot(tp, aes(x=cond,y=y)) +
        geom_tile(aes(fill=val)) +
        #geom_text(aes(label=lab, color=fc>swit), hjust=.5, size=2.5) +
        scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
        scale_y_discrete(breaks=tpy$y, labels=tpy$lab, expand=c(0,0)) +
        scale_fill_gradientn(name='normalized eigengene value',colors=cols100,limits=c(-1,1)) +
        #scale_fill_viridis(name='normalized eigengene value') +
        #scale_color_manual(values=c('black','white')) +
        facet_wrap(~pan, nrow=1, scale='free_y') +
        otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
               margin = c(.2,.1,.2,.1),
               ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
        guides(color = F)
    p = annotate_figure(p, top=text_grob(tit, size=9.5, hjust=.5, color='maroon'),
        fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
    fo = sprintf("%s/17.%s.%s.pdf", dirw, stress, drc)
    p %>% ggexport(filename = fo, width = 10, height = 7)
    T
    #}}}
}

tx %>% mutate(r=pmap_lgl(list(stress,drc,opt_deg,opt_clu,
                              softPower,deepSplit,MEDissThres,minGap, ng,np,me),
    plot_me, dirw=!!dirw))
#}}}

process_wgcna <- function(opt_deg, opt_clu, ti, gt_map, tc, dirw) {
    #{{{
    gts_deg = gt_map[[opt_deg]]; gts_clu = gt_map[[opt_clu]]
    diro = sprintf("%s/20_%s_%s", dirw, opt_deg, opt_clu)
    if(!dir.exists(diro)) dir.create(diro)
    ti1 = ti %>% group_by(bat, stress, drc) %>%
        nest() %>%
        mutate(r = pmap(list(bat, stress, drc, data), process_wgcna_bat,
                        diro, opt_deg, opt_clu)) %>%
        mutate(toc=map(r, 'toc'), tom=map(r, 'tom')) %>%
        select(bat,stress,drc,toc,tom)
    ti1
    #}}}
}
process_wgcna_bat <- function(bat, stress, drc, ti, diro, opt_deg, opt_clu) {
#{{{ different module merging parameters
#{{{ read, filter & preprocess
me = ti %>% select(-clu) %>% unnest(me) %>%
    filter(ctag != 'raw' | (ctag=='raw' & MEDissThres==.1)) %>%
    select(cid,optQ,deepSplit,MEDissThres,minGap,
           np,ng,ctag,n,gids,clu,me)
mex = me %>% distinct(ctag, MEDissThres, clu) %>%
    count(ctag, MEDissThres) %>%
    mutate(tag = sprintf("%s_%s (%s)", ctag, MEDissThres, n)) %>%
    select(-n) %>% arrange(desc(ctag), MEDissThres) %>% mutate(x = 1:n())
#}}}
#
#{{{ process tc
clu = ti %>% select(MEDissThres, clu) %>% mutate(ctag = 'merged') %>%
    unnest(clu)
clu0 = clu %>% filter(MEDissThres==0.1) %>%
    mutate(col2=col1,clu2=clu1,ctag='raw')
tc = clu %>% bind_rows(clu0) %>% rename(clu = clu2, raw_clu = clu1) %>%
    count(ctag, MEDissThres, clu, raw_clu) %>%
    arrange(ctag, MEDissThres, clu, raw_clu) %>%
    group_by(ctag, MEDissThres, clu) %>%
    summarise(n_clu = n(), raw_clus = str_c(raw_clu,collapse='+'), n = sum(n), raw_clu=raw_clu[1]) %>%
    ungroup() %>%
    arrange(desc(ctag), MEDissThres) %>%
    group_by(raw_clus,n_clu,raw_clu,n) %>% nest() %>% ungroup() %>%
    mutate(ctag1 = map_chr(data, mf<-function(df) df$ctag[[1]])) %>%
    mutate(MEDissThres1 = map_dbl(data, mf<-function(df) df$MEDissThres[[1]])) %>%
    mutate(clu = map_dbl(data, mf<-function(df) df$clu[[1]])) %>%
    arrange(n_clu, raw_clu)  %>%
    mutate(mid = sprintf("m%02d", 1:n())) %>%
    select(mid, raw_clus, n_clu, n, ctag1, MEDissThres1, raw_clu, clu, data)
    #%>% print(n=40)
#}}}
#
#{{{ plot cluster merging
tp = tc %>%
    #rename(ctag=ctag1, MEDissThres=MEDissThres1) %>%
    select(mid, raw_clus, n_clu, n, data) %>%
    unnest(data) %>%# filter(! (ctag=='merged' & n_clu == 1)) %>%
    mutate(raw_clus = str_split(raw_clus, '[\\+]')) %>%
    unnest(raw_clus) %>% rename(raw_clu=raw_clus) %>%
    mutate(raw_clu = as.integer(raw_clu)) %>%
    mutate(fcol = ifelse(n_clu==1, 'm01', mid)) %>%
    mutate(lab = sprintf("%s (%d)", mid, n)) %>%
    inner_join(mex, by=c('ctag','MEDissThres'))
tpy = tp %>% filter(ctag=='raw') %>% distinct(raw_clu) %>%
    arrange(desc(raw_clu)) %>% mutate(y=1:n())
tp = tp %>% inner_join(tpy, by='raw_clu')
#
pc = ggplot(tp, aes(x, y)) +
    geom_tile(aes(fill=fcol), color='white', size=.5) +
    geom_text(aes(label=lab, color=fcol), size=2.5) +
    scale_x_continuous(breaks=mex$x, labels=mex$tag, expand=expansion(mult=c(0,0))) +
    scale_y_continuous(breaks=tpy$y, expand=c(0,0)) +
    scale_fill_manual(values=cols36) +
    scale_color_manual(values=brights36) +
    otheme(ygrid=F, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F,
           margin = c(.3,.5,.3,.5)) +
    guides(fill=F,color=F) +
    theme(axis.text.x = element_text(angle=15, hjust=.5, vjust=.8, size=8))
#fo = file.path(dirw, '19.1.pdf')
#ggsave(fo, p, width=6, height=6)
#}}}
#
#{{{ save module IDs,  gene IDs and MEs
toc = tp %>% distinct(ctag, MEDissThres, mid)
tom = tc %>% select(-n) %>% rename(ctag=ctag1,MEDissThres=MEDissThres1) %>%
    inner_join(me, by=c('ctag','MEDissThres','clu')) %>%
    select(mid, n, gids, me)
#}}}
#
#{{{ plot title and label
y = me[1,]
tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
    number(y$ng), opt_deg, stress, drc, number(y$np), opt_clu)
labx = ifelse(y$minGap==0, 'auto', y$minGap)
lab = sprintf("deepSplit=%g minGap=%s", y$deepSplit, labx)
#
tp = tc %>% select(-n) %>% rename(ctag=ctag1,MEDissThres=MEDissThres1) %>%
    mutate(fcol = ifelse(n_clu==1, 'm01', mid)) %>%
    inner_join(me, by=c('ctag','MEDissThres','clu')) %>%
    select(mid, n, clu, fcol, me) %>%
    unnest(me) %>% gather(cond, val, -mid,-n, -clu,-fcol) %>%
    group_by(mid,clu) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    mutate(val = ifelse(val > 1, 1, val)) %>%
    mutate(val = ifelse(val < -1, -1, val))
tpy = tp %>% distinct(mid,n) %>% arrange(desc(mid)) %>%
    mutate(lab=sprintf("%s (%d)", mid, n)) %>% mutate(y=1:n())
tpx = tp %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n())
tp = tp %>% inner_join(tpy[,c('mid','lab','y')], by='mid') %>%
    inner_join(tpx, by='cond')
#{{{ heatmap ph
p = ggplot(tp, aes(x=cond,y=y)) +
    geom_tile(aes(fill=val)) +
    #geom_text(aes(label=lab, color=fc>swit), hjust=.5, size=2.5) +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_continuous(breaks=tpy$y, labels=tpy$lab, expand=c(0,0)) +
    scale_fill_gradientn(name='normalized eigengene value',colors=cols100,limits=c(-1,1)) +
    #scale_fill_viridis(name='normalized eigengene value') +
    #scale_color_manual(values=c('black','white')) +
    #facet_wrap(~, nrow=2, dir='v', scale='free_y') +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           margin = c(1,.3,.3,.3),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=1, size=8)) +
    guides(color = F)
ph = annotate_figure(p, top=text_grob(tit, size=9, hjust=.5, color='maroon'),
    fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
#fo = sprintf("%s/18.1.pdf", dirw)
#p %>% ggexport(filename = fo, width = 5, height = 6)
#}}}
#{{{ lineplot pl
p = ggplot(tp, aes(x=x,y=val)) +
    geom_line(aes(color=fcol), size=.5) +
    geom_point(color='black', size=.5) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$cond, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(values=cols36) +
    facet_wrap(~lab, ncol=3, dir='v', scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    guides(color = F)
pl = annotate_figure(p, top=text_grob(tit, size=9, hjust=.5, color='maroon'),
    fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
#}}}
#}}}
#
fo = sprintf("%s/%s.pdf", diro, bat)
ggarrange(pc, pl,
    widths=c(1,1.7), heights=c(1), labels=LETTERS[1:6], nrow=1, ncol=2) %>%
ggexport(filename = fo, width = 10,  height = 8)
list(toc=toc, tom=tom)
#}}}
}
tr = res %>% rename(bat=batch) %>% filter(bat %in% bats) %>%
    group_by(opt_deg, opt_clu) %>% nest() %>% ungroup() %>%
    mutate(r = pmap(list(opt_deg, opt_clu, data), process_wgcna,
                    gt_map=!!gt_map, tc=!!tc, dirw=!!dirw))

to = tr %>% select(-data) %>% unnest(r)
fo = file.path(dirw, "25.modules.rds")
saveRDS(to, fo)
#}}}

#{{{ plot MEs for selected merging parameters
#{{{ read modules and config
fi = file.path(dirw, "25.modules.rds")
x = readRDS(fi)
#
ff = file.path(dirw, 'config.xlsx')
tf = read_xlsx(ff, sheet='picked') %>%
    fill(opt_deg, .direction='down') %>%
    fill(opt_clu, .direction='down') %>%
    fill(bat, .direction='down') %>%
    group_by(opt_deg, opt_clu, bat) %>%
    mutate(i = 1:n()) %>% ungroup() %>% mutate(bat=factor(bat,levels=bats))
tfs = tf %>% count(opt_deg, opt_clu, bat) %>%
    mutate(y0 = c(0,cumsum(n)[1:(n()-1)]+(1:(n()-1))))
tf = tf %>% inner_join(tfs, by=c('opt_deg','opt_clu','bat')) %>%
    mutate(y=y0+i) %>% select(-n,-y0)
#
ti1 = tf %>% group_by(opt_deg, opt_clu,bat) %>% nest() %>% rename(picked=data)
ti = x %>% mutate(bat=factor(bat,levels=bats)) %>%
    inner_join(ti1, by=c("opt_deg",'opt_clu','bat')) %>%
    group_by(opt_deg,opt_clu) %>% nest()
#}}}
plot_me <- function(opt_deg, opt_clu, x, tc, dirw) {
    #{{{
    diro = sprintf("%s/20_%s_%s", dirw, opt_deg, opt_clu)
    if(!dir.exists(diro)) dir.create(diro)
    xc = x %>% select(bat, toc) %>% unnest(toc)
    xm = x %>% select(bat, tom) %>% unnest(tom)
    xp = x %>% select(bat, picked) %>% unnest(picked)
#
    #{{{ prop. of genes picked
    xm1 = xm %>% select(bat,gids) %>% unnest(gids) %>%
        distinct(bat, gids) %>% count(bat) %>% rename(ng_tot = n)
    xp1 = xp %>% select(bat, mid) %>% inner_join(xm, by=c('bat','mid')) %>%
        select(bat, gids) %>% unnest(gids) %>%
        distinct(bat, gids) %>% count(bat) %>% rename(ng_pick = n)
    bat_lab = xm1 %>% inner_join(xp1, by=c('bat')) %>% mutate(pg_pick=ng_pick/ng_tot) %>%
        mutate(lab = glue("{bat}\n{ng_pick}/{ng_tot}: {percent(pg_pick, accuracy=.1)}")) %>%
        select(bat, lab)
    #}}}
#{{{ heatmap of all modules
tp = xm %>% select(bat, mid, n, me) %>%
    unnest(me) %>% gather(cond, val, -bat,-mid,-n) %>%
    group_by(bat,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    mutate(bat = factor(bat, levels=bats))
tpx = tp %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n()) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tpy = tp %>% distinct(bat, mid, n) %>%
    mutate(lab=sprintf("%s (%d)", mid, n)) %>%
    mutate(y=str_c(bat, mid, sep=' ')) %>% arrange(bat, desc(mid)) %>% select(-n)
tp = tp %>% inner_join(tpx, by='cond') %>%
    inner_join(tpy, by=c('bat','mid')) %>%
    mutate(y = factor(y, levels=tpy$y))
p = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.2) +
    scale_x_continuous(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.01,.01))) +
    scale_y_discrete(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01))) +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    facet_wrap(~bat, ncol=4, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T) +
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '01.heatmap.all.pdf')
ggsave(p, file=fo, width=10, height=8)
#}}}
#
#{{{ heatmap for picked modules
tp0 = xp %>% inner_join(xm, by=c('bat','mid')) %>%
    select(y,bat,i,mid,n,me,note) %>%
    unnest(me) %>% gather(cond, val, -bat,-i,-y,-mid,-n, -note) %>%
    group_by(bat,i,y,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    #mutate(val = ifelse(val > 1, 1, val)) %>%
    #mutate(val = ifelse(val < -1, -1, val)) %>%
    mutate(lab=sprintf("%s (%d): %s", mid, n, note)) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(i = factor(i))
tpx = tp0 %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n()) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tp0 = tp0 %>% inner_join(tpx, by='cond')
tp = tp0
tpy = tp %>% distinct(y, bat, lab) %>% arrange(y)
tpys = tpy %>% group_by(bat) %>%
    summarise(ymin=min(y), ymax=max(y), ymid=(ymin+ymax)/2) %>% ungroup() %>%
    inner_join(bat_lab, by='bat')
p = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.3) +
    geom_segment(data=tpys, aes(x=.3,xend=.3,y=ymin,yend=ymax),color='dodgerblue',size=1) +
    geom_text(data=tpys, aes(x=-1,y=ymid,label=lab),size=3,angle=0,hjust=.5,vjust=.5) +
    scale_x_continuous(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.15,.01))) +
    scale_y_reverse(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01)), position='right') +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T)
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    #theme(strip.text.y = element_blank()) +
    #theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '03.heatmap.picked.pdf')
ggsave(p, file=fo, width=6, height=6)
#}}}
#
#{{{ violin & line plot
# scale asinh expression
scale2 <- function(ti) {r=max(ti$val)-min(ti$val); ti %>% mutate(val=val/r)}
te = tc$diff %>% filter(Genotype %in% gt_map[[opt_clu]]) %>%
    mutate(gid = str_c(gid, Genotype, sep="_")) %>%
    mutate(stress=str_to_lower(Treatment)) %>%
    select(-Genotype, -Treatment) %>%
    gather(has, val, -stress, -gid) %>% mutate(val=asinh(val)) %>%
    group_by(stress, gid) %>% nest() %>% ungroup() %>%
    mutate(ndata = map(data, scale2)) %>%
    select(stress, gid, r=ndata) %>% unnest(r)
#
tp = tp0 %>% mutate(x=factor(x))
tpx = tp %>% distinct(cond, x) %>% arrange(cond, x) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tpl = tp %>% distinct(bat,i,mid,lab)
te1 = tp0 %>% distinct(bat,mid) %>% inner_join(xm, by=c('bat','mid')) %>%
    select(bat, mid, gids) %>% unnest(gids) %>%
    rename(gid=gids) %>%
    #mutate(gid = str_replace(gid, "_.*$", '')) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    separate(bat, c('stress','drc'), remove=F) %>% select(-drc) %>%
    inner_join(te, by=c('stress','gid')) %>%
    rename(cond=has) %>%
    inner_join(tpx, by='cond') %>%
    inner_join(tpl, by=c('bat','mid'))
#
p = ggviolin(te1, x='x', y='val', color='grey', fill='gray') +
    #geom_boxplot(te1, mapping=aes(x=x,y=val,group=x), width=.6, color='grey') +
    geom_line(tp, mapping=aes(x=x,y=val,group=y,color=bat), size=.5) +
    geom_point(tp, mapping=aes(x=x,y=val), color='gray35', size=1) +
    geom_text(tpl, mapping=aes(x=1,y=1.25,label=lab), color='royalblue', hjust=0, vjust=0, size=2.5) +
    scale_x_discrete(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.01,.08))) +
    scale_color_aaas() +
    facet_grid(i ~ bat, scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F) +
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '05.violin.picked.pdf')
ggsave(p, file=fo, width=8, height=6)
#}}}
#
#{{{ multi-line plot
# scale asinh expression
scale2 <- function(ti) {r=max(ti$val)-min(ti$val); ti %>% mutate(val=val/r)}
tp = tp0 %>% mutate(x=factor(x))
tpx = tp %>% distinct(cond, x) %>% arrange(cond, x) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tpl = tp %>% filter(cond=='h250') %>% select(bat, mid,i, lab,val)
#
p = ggplot(tp) +
    geom_line(aes(x=x,y=val,group=y,color=as.character(i)), size=.5) +
    geom_point(tp, mapping=aes(x=x,y=val), color='gray35', size=1) +
    geom_text_repel(tpl, mapping=aes(x=8.1,y=val,label=lab,color=as.character(i)), hjust=0, vjust=0, size=2.5, direction='y',nudge_x=.5, nudge_y=0, segment.size=.2) +
    scale_x_discrete(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.02,.5))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.01,.08))) +
    scale_color_aaas() +
    facet_wrap(~bat, ncol=1, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3), strip.style='white',
           ygrid=T, xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F) +
    theme(strip.text = element_text(hjust=0, size=11))
fo = file.path(diro, '05.violin.picked.combine.pdf')
ggsave(p, file=fo, width=5, height=7)
#}}}
#
#{{{ write xm for sharing
to = xm %>% mutate(bat=factor(bat, levels=bats)) %>% arrange(bat,mid) %>%
    mutate(gids = map_chr(gids, str_c, collapse=',')) %>%
    mutate(me = map_chr(me, str_c, collapse=','))
fo = file.path(diro, '09.modules.tsv')
write_tsv(to, fo)
#}}}
    #}}}
}

ti %>% mutate(r = pmap(list(opt_deg, opt_clu, data), plot_me,
                       tc=!!tc, dirw=!!dirw))
#}}}


#{{{ [obsolete] test plot members within individual module
deg1 = deg$deg48 %>% filter(Genotype=='B73',Treatment=='Cold',cond2=='timeM') %>%
    select(Treatment,Timepoint,ds) %>% unnest(ds)
xm1 = xm %>% select(-me) %>% unnest(gids) %>% rename(gid=gids) %>%
    mutate(gid = str_replace(gid, "_.*$", ''))
gid0 = 'Zm00001d002065'
xm1 %>% filter(gid==gid0)
deg1 %>% filter(gid==gid0)

bat = 'cold_up'; mid = 'm20'
gids = xm %>% filter(bat == !!bat, mid == !!mid) %>% unnest(gids) %>%
    mutate(gid = str_replace(gids, "_.*$", '')) %>% pull(gid)

p = ggplot(tp, aes(x=x,y=val)) +
    #geom_line(aes(group=gid),size=.3, color='grey') +
    geom_boxplot(aes(group=x)) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$has, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized expression value',expand=expansion(mult=c(.05,.05))) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(dirw, 't.pdf')
ggsave(p, file=fo, width=8, height=6)
#
#}}}


