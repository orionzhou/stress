source('functions.R')
require(WGCNA)
dirw = file.path(dird, '17_cluster')
enableWGCNAThreads()

#{{{ read
fi = file.path(dirw, '../15_de/10.rds')
deg = readRDS(fi)
#
fi = file.path(dirw, '01.tc.rds')
tc = readRDS(fi)
#
f_cfg = file.path(dirw, 'clu_options.xlsx')
cfg = read_xlsx(f_cfg)
#}}}

#{{{ ad hoc pipe run
cfg = read_xlsx(f_cfg)
res = cfg %>% filter(cid == 'c21') %>%
    mutate(x=pmap(list(cid,cond,drc,opt_deg,opt_clu,optQ,
                       softPower,deepSplit,MEDissThres,minGap),
                  run_wgcna_pipe, tc=!!tc, deg=!!deg, dirw=!!dirw))
#}}}

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

bats = c('cold_up', 'heat_up', 'cold_down', 'heat_down')
bat = bats[1]
#{{{ different module merging parameters
#{{{ read, filter & preprocess
res0 = res %>% filter(batch==bat)
me = res0 %>% select(-clu) %>% unnest(me) %>%
    filter(ctag != 'raw' | (ctag=='raw' & cid %in% c('c21','c26','c31','c36'))) %>%
    select(cid,stress,drc,opt_deg,opt_clu,optQ,deepSplit,MEDissThres,minGap,
           np,ng,ctag,n,gids,clu,me)
mex = me %>% distinct(ctag, MEDissThres, clu) %>%
    count(ctag, MEDissThres) %>%
    mutate(tag = sprintf("%s_%s (%s)", ctag, MEDissThres, n)) %>%
    select(-n) %>% arrange(desc(ctag), MEDissThres) %>% mutate(x = 1:n())
#}}}
#
#{{{ process tc
clu = res0 %>% select(MEDissThres, clu) %>% mutate(ctag = 'merged') %>%
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
    select(mid, raw_clus, n_clu, n, ctag1, MEDissThres1, raw_clu, clu, data) %>% print(n=40)
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
toc1 = tp %>% distinct(ctag, MEDissThres, mid)
tom1 = tc %>% select(-n) %>% rename(ctag=ctag1,MEDissThres=MEDissThres1) %>%
    inner_join(me, by=c('ctag','MEDissThres','clu')) %>%
    select(mid, n, gids, me)
if (!exists('toc')) toc = list()
if (!exists('tom')) tom = list()
toc[[bat]] = toc1
tom[[bat]] = tom1
#}}}
#
#{{{ plot title and label
y = me[1,]
tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
    number(y$ng), y$opt_deg, y$stress, y$drc, number(y$np), y$opt_clu)
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
fo = sprintf("%s/20.%s.pdf", dirw, bat)
ggarrange(pc, pl,
    widths=c(1,1.7), heights=c(1), labels=LETTERS[1:6], nrow=1, ncol=2) %>%
ggexport(filename = fo, width = 10,  height = 8)

to1 = tibble(bat=names(toc)) %>% mutate(clu=map(bat, myf <- function(x) toc[[x]]))
to2 = tibble(bat=names(tom)) %>% mutate(me=map(bat, myf <- function(x) tom[[x]]))
to = to1 %>% inner_join(to2, by='bat')
fo = file.path(dirw, "25.modules.rds")
saveRDS(to, fo)
#}}}

#{{{ plot MEs for selected merging parameters
fi = file.path(dirw, "25.modules.rds")
x = readRDS(fi)
xc = x %>% select(bat, clu) %>% unnest(clu)
xm = x %>% select(bat, me) %>% unnest(me)

ff = file.path(dirw, 'config.xlsx')
tf = read_xlsx(ff, sheet='picked')
ti = tf %>% inner_join(xm, by=c('bat','mid'))
#
tp0 = ti %>% select(y,bat,i,mid,n,me,note) %>%
    unnest(me) %>% gather(cond, val, -bat,-i,-y,-mid,-n, -note) %>%
    group_by(bat,i,y,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    #mutate(val = ifelse(val > 1, 1, val)) %>%
    #mutate(val = ifelse(val < -1, -1, val)) %>%
    mutate(lab=sprintf("%s (%d): %s", mid, n, note)) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(i = factor(i))
tpx = tp0 %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n())
tp0 = tp0 %>% inner_join(tpx, by='cond')

#{{{ heatmap of all modules
tp = xm %>% select(bat, mid, n, me) %>%
    unnest(me) %>% gather(cond, val, -bat,-mid,-n) %>%
    group_by(bat,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    mutate(bat = factor(bat, levels=bats))
tpx = tp %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n())
tpy = tp %>% distinct(bat, mid, n) %>%
    mutate(lab=sprintf("%s (%d)", mid, n)) %>%
    mutate(y=str_c(bat, mid, sep=' ')) %>% arrange(bat, desc(mid)) %>% select(-n)
tp = tp %>% inner_join(tpx, by='cond') %>%
    inner_join(tpy, by=c('bat','mid')) %>%
    mutate(y = factor(y, levels=tpy$y))
p = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.2) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$cond, expand=expansion(mult=c(.01,.01))) +
    scale_y_discrete(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01))) +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    facet_wrap(~bat, ncol=4, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(dirw, '27.all.pdf')
ggsave(p, file=fo, width=10, height=8)
#}}}

#{{{ heatmap for picked modules
tp = tp0
tpy = tp %>% distinct(y, bat, lab) %>% arrange(y)
tpys = tpy %>% group_by(bat) %>%
    summarise(ymin=min(y), ymax=max(y), ymid=(ymin+ymax)/2) %>% ungroup()
p = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.3) +
    geom_segment(data=tpys, aes(x=.3,xend=.3,y=ymin,yend=ymax),color='dodgerblue',size=1) +
    geom_text(data=tpys, aes(x=.1,y=ymid,label=bat),size=3,angle=0,hjust=1,vjust=.5) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$cond, expand=expansion(mult=c(.2,.01))) +
    scale_y_reverse(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01)), position='right') +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    #theme(strip.text.y = element_blank()) +
    #theme(strip.text.x = element_text(size=10))
fo = file.path(dirw, '27.picked.pdf')
ggsave(p, file=fo, width=6, height=5)
#}}}

#{{{ lineplot
tp = tp0
tpx = tp %>% distinct(cond, x) %>% arrange(cond, x)
p = ggplot(tp, aes(x=x,y=val)) +
    geom_line(size=.5) +
    geom_point(color='black', size=.5) +
    geom_text(aes(x=1,y=1,label=lab), color='royalblue', hjust=0, vjust=1, size=3) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$cond, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.05,.05))) +
    facet_grid(y ~ bat) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(dirw, '28.picked.pdf')
ggsave(p, file=fo, width=8, height=6)
#}}}

#{{{ write xm for sharing
to = xm %>% mutate(bat=factor(bat, levels=bats)) %>% arrange(bat,mid) %>%
    mutate(gids = map_chr(gids, str_c, collapse=',')) %>%
    mutate(me = map_chr(me, str_c, collapse=','))
fo = file.path(dirw, '25.modules.tsv')
write_tsv(to, fo)
#}}}
#}}}


i=7
cid=cfg$cid[i];cond=cfg$cond[i];drc=cfg$drc[i];
opt_deg=cfg$opt_deg[i];opt_clu=cfg$opt_clu[i];optQ=cfg$optQ[i];
softPower=cfg$softPower[i];deepSplit=cfg$deepSplit[i];
MEDissThres=cfg$MEDissThres[i]; minGap=cfg$minGap[i]
