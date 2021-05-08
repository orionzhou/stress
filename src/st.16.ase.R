source('functions.R')
dirw = glue('{dird}/16_ase')
regs = c('conserved','unexpected','cis','trans','cis+trans')
conds = c("Control0",'Control1','Control25','Cold1','Cold25','Heat1','Heat25')
crosses = c("B73xMo17",'W22xB73','W22xMo17')

#{{{ demo cis/trans modes
fi = file.path(dirw, 'modes_demo.xlsx')
ti = read_xlsx(fi) %>% fill(mode)

tp = ti %>% gather(opt, lfc, -mode, -i) %>%
    mutate(mode=factor(mode,levels=unique(ti$mode))) %>%
    separate(opt, c('gen','cond','gt'), sep=c(1,2)) %>%
    mutate(gt = ifelse(gt == '2', 'B', 'M')) %>%
    mutate(x.gen = ifelse(gen == 'h', 3, 0)) %>%
    mutate(x.cond = ifelse(cond == 's', 2, 0)) %>%
    mutate(x = x.gen + x.cond) %>% select(-x.gen, -x.cond) %>%
    mutate(gen = ifelse(gen == 'h', 'hybrid', 'parent')) %>%
    mutate(gen = factor(gen, levels=c('parent','hybrid')))
tp2 = tp %>% pivot_wider(names_from=cond, values_from=c(x,lfc), names_sep='.')
tp3 = tp %>% filter(i==1, cond=='c')
tp4 = tp2 %>% distinct(mode,i, gen, x.c, x.s)
tp4s = tp4 %>% filter(i==1)
#
xnbreaks = tp %>% distinct(x) %>% arrange(x) %>% pull(x)
p = ggplot(tp) +
    geom_rect(data=tp4, aes(xmin=x.c,xmax=x.s, ymin=-Inf,ymax=Inf, fill=gen), alpha=.3) +
    geom_text(data=tp4s, aes(x=(x.c+x.s)/2, y=5, label=gen), size=3) +
    geom_segment(data=tp2, aes(x=x.c,xend=x.s,y=lfc.c,yend=lfc.s,color=gen), size=1) +
    geom_point(aes(x=x, y=lfc, shape=gt), size = 2) +
    geom_text(data=tp3, aes(x=x-.05,y=lfc+.2,label=gt), hjust=1,vjust=0, size=3) +
    scale_x_continuous(breaks=xbreaks, labels=rep(c('control','stress'),2), expand=expansion(mult=c(.1,.1))) +
    scale_y_continuous(name='log2(ReadCount+1)', expand=expansion(mult=c(.15,.15))) +
    scale_linetype_manual(values=c('solid','dashed')) +
    scale_shape_manual(name='Genotype', values=c(15,16)) +
    scale_color_manual(name='Generation', values=pal_npg()(3)) +
    scale_fill_manual(values=pal_simpsons()(3)) +
    facet_grid(i ~ mode) +
    otheme(xtick=T,xtext=T,ytick=T,ytext=T,ytitle=T,xgrid=F,ygrid=T,
           legend.pos='none', legend.dir='v', legend.title=T,
           legend.box='h') +
    theme(axis.text.x = element_text(angle=15))
fo = file.path(dirw, 'modes_demo.pdf')
ggsave(p, file=fo, width=8, height=5)
#}}}

#{{{ organize data for tests
yid = 'rn20a'; res = rnaseq_cpm(yid)

#{{{ get sizeFactor and dispersions using DESeq2
require(DESeq2)
th = res$th %>% filter(Experiment=='HY') %>%
    mutate(cond=str_c(Treatment,Timepoint,sep='')) %>%
    select(SampleID, gt=Genotype, cond)
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    select(gid, SampleID, ReadCount)
tw = tm %>%
    spread(SampleID, ReadCount) %>%
    replace(., is.na(.), 0)
gids = tw$gid
twd = data.frame(tw[,-1])
rownames(twd) = tw$gid
th2 = th %>% mutate(sid = SampleID, gt=factor(gt), cond=factor(cond))
thd = column_to_rownames(as.data.frame(th2), var = 'sid')
stopifnot(identical(rownames(thd), colnames(twd)))
dds = DESeqDataSetFromMatrix(countData=twd, colData=thd, design = ~ gt + cond)
dds = estimateSizeFactors(dds)
sf = sizeFactors(dds)
t_sf = tibble(SampleID = names(sf), sizeFactor = as.numeric(sf))
dds = estimateDispersions(dds, fitType = 'parametric')
stopifnot(identical(gids, names(dds)))
t_disp = tibble(gid=gids, disp = dispersions(dds))
tl = th %>% inner_join(t_sf, by = 'SampleID')
#}}}

th1 = th %>% filter(str_detect(gt, 'x')) %>% distinct(cond, gt) %>%
    mutate(cross = gt) %>%
    dplyr::rename(h=gt) %>% separate(h, c("p1",'p2'), sep='x', remove=F) %>%
    gather(gen, gt, -cond, -cross) %>%
    inner_join(th, by=c("cond",'gt')) %>% arrange(cond,cross,gen,gt,SampleID)
th2a = th1 %>% distinct(cond) %>% filter(!str_detect(cond, '^Control')) %>%
    mutate(condB = str_replace(cond, "Cold", "Control")) %>%
    mutate(condB = str_replace(condB, "Heat", "Control"))
th2b = th1 %>% distinct(cond) %>% filter(!str_detect(cond, '^Control')) %>%
    mutate(condB = 'Control0')
th2 = th2a %>% bind_rows(th2b)

tm_p = th1 %>% inner_join(t_sf, by='SampleID') %>%
    inner_join(tm, by=c('SampleID')) %>%
    mutate(nRC = ReadCount / sizeFactor) %>%
    group_by(cond,cross,gen,gid) %>%
    summarise(rc=list(nRC), trc=sum(nRC)) %>% ungroup() %>%
    pivot_wider(names_from=gen, values_from=c(rc,trc), names_sep='.')
#
tm_ase = th1 %>% filter(gen=='h') %>%
    inner_join(t_sf, by='SampleID') %>%
    inner_join(res$ase_gene, by=c("SampleID")) %>%
    mutate(nRC1 = allele1 / sizeFactor, nRC2 = allele2 / sizeFactor) %>%
    group_by(cond,cross,gid) %>%
    summarise(rc.h1=list(nRC1), rc.h2=list(nRC2),
              trc.h1=sum(nRC1), trc.h2=sum(nRC2)) %>% ungroup()
#
tmt = tm_p %>% inner_join(t_disp, by='gid') %>%
    inner_join(tm_ase, by=c("cond",'cross','gid')) %>%
    select(cond,cross,gid,disp,trc.p1,trc.p2,trc.h,trc.h1,trc.h2,
           rc.p1,rc.p2,rc.h,rc.h1,rc.h2)

ra = list(tl=tl, disp=t_disp, th=th1, th.cmp = th2, tmt = tmt)
fo = file.path(dirw, '01.raw.rds')
saveRDS(ra, fo)
#}}}

# run st.16.ase.1.R
# run st.16.ase.2.R

#{{{ simple cis/trans
fi = file.path(dirw, '05.modes.rds')
ti = readRDS(fi)
ti %>% count(cond,cross,reg) %>% spread(reg,n) %>% print(n=40)
ct_basic = ti

nconds = crossing(x=conds, y=crosses) %>% mutate(cond=str_c(x,y,sep='_')) %>%
    mutate(x=factor(x,levels=conds), y=factor(y,levels=crosses)) %>%
    arrange(x,y) %>% pull(cond)
tpx = tibble(reg=regs) %>% mutate(x=1:n())
tp = ti %>% mutate(cond = str_c(cond, cross, sep='_')) %>%
    inner_join(tpx, by='reg') %>%
    mutate(cond = factor(cond, levels=nconds)) %>%
    mutate(reg=factor(reg, levels=regs))
linecol = 'azure3'; lty = 'solid'
cols9 = pal_npg()(10)
#{{{ scatter plot
p = ggplot(tp, aes(x=prop.p, y=prop.h)) +
    #geom_vline(xintercept=0, linetype=lty, color=linecol) +
    #geom_hline(yintercept=0, linetype=lty, color=linecol) +
    geom_point(aes(color=reg, shape=reg), size=1) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    scale_x_continuous(name='Allele 1 proportion in Parent', limits=c(0,1),expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Allele 1 proportion in F1', limits=c(0,1),expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='Mode:', values=cols9) +
    scale_shape_manual(name='Mode:', values=c(1:5)) +
    facet_wrap(~cond, ncol=6) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=1))
fp = sprintf('%s/10.modes.pdf', dirw)
ggsave(p, file = fp, width = 12, height = 9)
#}}}

#{{{ bar plot showing counts
tp1 = tp %>% count(cond, x, reg)
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
p = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=reg), width=.7, stat='identity') +
    geom_text(data=tp1s,aes(x=5.6,y=6200,label=lab), size=2.5, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+100, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$reg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.1))) +
    scale_fill_manual(name='Modes:', values=cols9) +
    facet_wrap(~cond, ncol=6) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           legend.vjust = -.4, legend.box='v',
           xtick=F, xtext=F, xtitle=F, ytitle=T, ytick=T, ytext=T)
fp = sprintf('%s/10.modes.cnt.pdf', dirw)
ggsave(p, file = fp, width = 8, height = 6)
#}}}
#}}}

#{{{ complex cis/trans
get_reg <- function(ti) ti$reg[1]
fmax = .5
fi = file.path(dirw, '05.modes.x.rds')
tr = readRDS(fi) %>%
    mutate(prop.p.s = mrc.p1/(mrc.p1+mrc.p2)) %>%
    mutate(prop.p.c = mrc0.p1/(mrc0.p1+mrc0.p2)) %>%
    mutate(prop.h.s = mrc.h1/(mrc.h1+mrc.h2)) %>%
    mutate(prop.h.c = mrc0.h1/(mrc0.h1+mrc0.h2)) %>%
    mutate(prop.p = prop.p.s - prop.p.c) %>%
    mutate(prop.h = prop.h.s - prop.h.c) %>%
    mutate(prop.p = map_dbl(prop.p, set_bound, minV=-fmax, maxV=fmax)) %>%
    mutate(prop.h = map_dbl(prop.h, set_bound, minV=-fmax, maxV=fmax))
tr2 = tr %>% mutate(reg = map_chr(reg, get_reg))
tr2 %>% count(cond,condB,cross,reg) %>% spread(reg,n) %>% print(n=40)
ct_stress = tr2

nconds = crossing(x=conds, y=crosses) %>% mutate(cond=str_c(x,y,sep='_')) %>%
    mutate(x=factor(x,levels=conds), y=factor(y,levels=crosses)) %>%
    arrange(x,y) %>% pull(cond)
tpx = tibble(reg=regs) %>% mutate(x=1:n())
tp = tr2 %>% filter(condB!='Control0') %>%
    mutate(cond = str_c(cond, cross, sep='_')) %>%
    inner_join(tpx, by='reg') %>%
    mutate(cond = factor(cond, levels=nconds)) %>%
    mutate(reg=factor(reg, levels=regs))
t_ase = tp
linecol = 'azure3'; lty = 'solid'
cols9 = pal_npg()(10)

#{{{ scatter plot
limits = c(-fmax, fmax)
p = ggplot(tp, aes(x=prop.p, y=prop.h)) +
    geom_point(aes(color=reg, shape=reg), size=1) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    scale_x_continuous(name= 'Proportion of Allele 1 Change in Parent', limits=limits, ,expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Proportion of Allele 1 Change in F1', limits=limits, expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='Mode:', values=cols9) +
    scale_shape_manual(name='Mode:', values=c(1:5)) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=1))
fp = sprintf('%s/11.modex.pdf', dirw)
ggsave(p, file = fp, width = 10, height = 14)
#}}}

#{{{ bar plot showing counts
tp1 = tp %>% count(cond, x, reg)
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
p = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=reg), width=.7, stat='identity') +
    geom_text(data=tp1s,aes(x=5.6,y=6200,label=lab), size=2.5, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+100, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$reg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.1))) +
    scale_fill_manual(name='Modes:', values=cols9) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=F, xtext=F, xtitle=F, ytitle=T, ytick=T, ytext=T)
fp = sprintf('%s/11.modex.cnt.pdf', dirw)
ggsave(p, file = fp, width = 4.5, height = 6)
#}}}
#}}}

#{{{ cis/trans for dDEGs - f2d-e, sf04
#{{{ read
fi = file.path(dird, '15_de/05.rds')
r = readRDS(fi)
td = r$td
tp1 = td %>% filter(ddrc != 0, deg != 'A=B=') %>%
    mutate(cond = str_c(stress,Timepoint, sep='')) %>%
    select(cond,qry,tgt,gid,deg,x)
#
get_reg <- function(ti) ti$reg[1]
fr = file.path(dirw, '05.modes.x.rds')
tr = readRDS(fr) %>% mutate(reg = map_chr(reg, get_reg)) #%>%
    #mutate(prop.p.s = mrc.p1/(mrc.p1+mrc.p2)) %>%
    #mutate(prop.p.c = mrc0.p1/(mrc0.p1+mrc0.p2)) %>%
    #mutate(prop.h.s = mrc.h1/(mrc.h1+mrc.h2)) %>%
    #mutate(prop.h.c = mrc0.h1/(mrc0.h1+mrc0.h2)) %>%
    #mutate(prop.p = prop.p.s - prop.p.c) %>%
    #mutate(prop.h = prop.h.s - prop.h.c) %>%
    #mutate(prop.p = map_dbl(prop.p, set_bound, minV=-fmax, maxV=fmax)) %>%
    #mutate(prop.h = map_dbl(prop.h, set_bound, minV=-fmax, maxV=fmax))
tr %>% count(cond,condB,cross,reg) %>% spread(reg,n) %>% print(n=40)
tp2 = tr %>%
    mutate(cross = ifelse(cross=='B73xMo17', 'Mo17xB73', cross)) %>%
    separate(cross, c('qry','tgt'), sep='x')
#
ddeg = tp1 %>% inner_join(tp2, by=c("cond","qry",'tgt','gid'))
ddeg %>% count(cond,qry,tgt,deg, reg) %>%
    spread(reg,n) %>% print(n=40)
#
limits = c(0,1)
nconds = crossing(x=conds, y=crosses) %>%
    mutate(cond=str_c(x,y,sep='_')) %>%
    mutate(x=factor(x,levels=conds), y=factor(y,levels=crosses)) %>%
    arrange(x,y) %>% pull(cond)
tpx = tibble(reg=regs) %>% mutate(x=1:n())
tp = ddeg %>% filter(condB!='Control0', deg %in% c("A+B=",'A=B+')) %>% select(-x) %>%
    mutate(cross = str_c(qry, tgt, sep='x')) %>%
    mutate(cross = ifelse(cross=='Mo17xB73', 'B73xMo17', cross)) %>%
    mutate(cond = str_c(cond, cross, sep='_')) %>%
    inner_join(tpx, by='reg') %>%
    mutate(cond = factor(cond, levels=nconds)) %>%
    mutate(reg=factor(reg, levels=regs))
linecol = 'azure3'; lty = 'solid'
cols9 = pal_simpsons()(10)
#}}}

#{{{ scatter plot - sf04a
p = ggplot(tp, aes(x=prop.p, y=prop.h)) +
    geom_point(aes(color=reg, shape=reg), size=1) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    scale_x_continuous(name= 'Proportion of Allele 1 Change in Parent', limits=limits,expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Proportion of Allele 1 Change in F1', limits=limits,expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='Mode:', values=cols9) +
    scale_shape_manual(name='Mode:', values=c(1:5)) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -1.2, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=1))
fo = glue('{dirf}/sf04a.rds')
saveRDS(p, fo)
#}}}
fp = glue('{dirw}/13.modex.pdf')
ggsave(p, file = fp, width = 5, height = 7)
#{{{ bar plot showing counts - sf04b
tp1 = tp %>% count(cond, x, reg)
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
ymax = 100
p = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=reg), width=.7, stat='identity') +
    geom_text(data=tp1s,aes(x=5.6,y=ymax,label=lab), size=2.5, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+10, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$reg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.1))) +
    scale_fill_manual(name='Mode:', values=cols9) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -1.2, legend.box='h',
           xtick=F, xtext=F, xtitle=F, ytitle=T, ytick=T, ytext=T)
fo = glue('{dirf}/sf04b.rds')
saveRDS(p, fo)
#}}}
fp = glue('{dirw}/13.modex.cnt.pdf')
ggsave(p, file = fp, width = 4.5, height = 6)

#{{{ scatter plot f2d
tp1 = tp %>% filter(cond=='Cold25_B73xMo17')
p = ggplot(tp1, aes(x=prop.p, y=prop.h)) +
    geom_point(aes(color=reg, shape=reg), size=1) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    scale_x_continuous(name= 'Proportion of Allele 1 Change in Parent', limits=limits,expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Proportion of Allele 1 Change in F1', limits=limits,expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='Mode:', values=cols9) +
    scale_shape_manual(name='Mode:', values=c(1:5)) +
    facet_wrap(~cond, ncol=1) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=F,
           legend.spacing.x=.05, legend.vjust = -.6, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=2))
fo = glue('{dirf}/f2d.rds')
saveRDS(p, fo)
#}}}
fp = glue("{dirw}/14.modex.1.pdf")
ggsave(p, file = fp, width = 4, height = 4)
#{{{ bar plot f2e
tp1 = tp %>% count(cond, x, reg) %>%
    filter(cond=='Cold25_B73xMo17')
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
p = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=reg), width=.7, stat='identity') +
    geom_text(data=tp1s,aes(x=5.6,y=120,label=lab), size=2.5, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+10, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$reg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.1))) +
    scale_fill_manual(name='Modes:', values=cols9) +
    facet_wrap(~cond, ncol=1, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           legend.vjust = -.6, legend.box='h',
           xtick=T, xtext=T, xtitle=F, ytitle=T, ytick=T, ytext=T) +
    theme(axis.text.x = element_text(angle=20, size=7.5, vjust=1, hjust=1))
fo = glue('{dirf}/f2e.rds')
saveRDS(p, fo)
#}}}
#}}}

res = list(ct_basic=ct_basic, ct_stress=ct_stress, ddeg=ddeg)
fo = glue('{dirw}/20.rds')
saveRDS(res, fo)

#{{{ showcase dDEG and cis/trans
#{{{ read in
fi = glue('{dirw}/20.rds')
res = readRDS(fi)
ddeg = res$ddeg

#{{{ prepare & filtering
fi = file.path(dirw, '01.raw.rds')
ra = readRDS(fi)
tmt = ra$tmt
#
min_rc = 10
t_rc = ra$th.cmp %>% inner_join(tmt, by=c('condB'='cond')) %>% select(-trc.h) %>%
    rename(trc0.p1=trc.p1,trc0.p2=trc.p2,trc0.h1=trc.h1,trc0.h2=trc.h2,
           rc0.p1=rc.p1,rc0.p2=rc.p2,rc0.h1=rc.h1,rc0.h2=rc.h2) %>%
    inner_join(tmt, by=c("cond",'cross','gid','disp')) %>%
    filter(trc.p1 + trc.p2 >= 2*min_rc, trc.h1+trc.h2 >= min_rc) %>%
    filter(trc0.p1 + trc0.p2 >= 2*min_rc, trc0.h1+trc0.h2 >= min_rc) %>%
    mutate(n.p1 = map_int(rc.p1, length), n.p2=map_int(rc.p2, length),
           n.h1 = map_int(rc.h1, length), n.h2=map_int(rc.h2, length),
           n0.p1 = map_int(rc0.p1, length), n0.p2=map_int(rc0.p2, length),
           n0.h1 = map_int(rc0.h1, length), n0.h2=map_int(rc0.h2, length)) %>%
    mutate(mrc.p1 = trc.p1/n.p1, mrc.p2 = trc.p2/n.p2,
           mrc.h1 = trc.h1/n.h1, mrc.h2 = trc.h2/n.h2,
           mrc0.p1 = trc0.p1/n0.p1, mrc0.p2 = trc0.p2/n0.p2,
           mrc0.h1 = trc0.h1/n0.h1, mrc0.h2 = trc0.h2/n0.h2) %>%
    mutate(dmrc.p1 = mrc.p1 - mrc0.p1, dmrc.p2 = mrc.p2 - mrc0.p2,
           dmrc.h1 = mrc.h1 - mrc0.h1, dmrc.h2 = mrc.h2 - mrc0.h2) %>%
    filter((dmrc.p1>=0 & dmrc.p2>=0 & dmrc.p1+dmrc.p2>0 &
            dmrc.h1>=0 & dmrc.h2>=0 & dmrc.h1+dmrc.h2>0) |
           (dmrc.p1<=0 & dmrc.p2<=0 & dmrc.p1+dmrc.p2<0 &
            dmrc.h1<=0 & dmrc.h2<=0 & dmrc.h1+dmrc.h2<0)) %>%
    mutate(prop.p=dmrc.p1/(dmrc.p1+dmrc.p2), prop.h=dmrc.h1/(dmrc.h1+dmrc.h2)) %>%
    select(cond,condB,cross,gid,mrc.p1,mrc.p2,mrc.h1,mrc.h2,
           mrc0.p1,mrc0.p2,mrc0.h1,mrc0.h2,
           dmrc.p1,dmrc.p2,dmrc.h1,dmrc.h2,prop.p,prop.h,
           rc.p1,rc.p2,rc.h1,rc.h2,
           rc0.p1,rc0.p2,rc0.h1,rc0.h2,disp) %>%
    mutate(i= 1:n())
#}}}
t_rc = t_rc %>%
    mutate(cross = ifelse(cross=='B73xMo17', 'Mo17xB73', cross)) %>%
    separate(cross, c('qry','tgt'), sep='x')

ti = ddeg %>% filter(condB != 'Control0') %>%
    select(cond,condB,qry,tgt,gid,deg,reg) %>%
    inner_join(t_rc, by=c("cond",'condB','qry','tgt','gid'))
#}}}

#{{{ plot
tis = ti %>% filter(cond=='Cold25', str_detect(deg, '[+]')) %>%
    arrange(reg) %>% group_by(reg) %>% slice(1:3) %>% ungroup()
tp = ti %>% filter(i %in% tis$i) %>%
    mutate(pan = sprintf("%s: %s\n%s %s_%s %s", reg, str_replace(gid,'0000',''), cond, qry, tgt, deg)) %>%
    select(pan, reg,
        pc1=mrc0.p1, ps1=mrc.p1, pc2=mrc0.p2, ps2=mrc.p2,
        hc1=mrc0.h1, hs1=mrc.h1, hc2=mrc0.h2, hs2=mrc.h2) %>%
    #mutate(pc1=asinh(pc1), ps1=asinh(ps1), pc2=asinh(pc2), ps2=asinh(ps2),
           #hc1=asinh(hc1), hs1=asinh(hs1), hc2=asinh(hc2), hs2=asinh(hs2)) %>%
    mutate(reg = factor(reg, levels = regs)) %>%
    gather(opt, lfc, -pan, -reg) %>%
    separate(opt, c('gen','cond','gt'), sep=c(1,2)) %>%
    mutate(gt = ifelse(gt == '2', 'B', 'A')) %>%
    mutate(x.gen = ifelse(gen == 'h', 3, 0)) %>%
    mutate(x.cond = ifelse(cond == 's', 2, 0)) %>%
    mutate(x = x.gen + x.cond) %>% select(-x.gen, -x.cond) %>%
    mutate(gen = ifelse(gen == 'h', 'hybrid', 'parent')) %>%
    mutate(gen = factor(gen, levels=c('parent','hybrid')))
pans =  tp %>% distinct(reg, pan) %>% arrange(reg) %>% pull(pan)
tp = tp %>% mutate(pan = factor(pan, levels=pans))
tp2 = tp %>% pivot_wider(names_from=cond, values_from=c(x,lfc), names_sep='.')
tp3 = tp %>% filter(cond=='c')
tp4 = tp2 %>% distinct(pan, gen, x.c, x.s)
tp4s = tp4
#
xbreaks = tp %>% distinct(x) %>% arrange(x) %>% pull(x)
p = ggplot(tp) +
    geom_rect(data=tp4, aes(xmin=x.c,xmax=x.s, ymin=-Inf,ymax=Inf, fill=gen), alpha=.3) +
    geom_text(data=tp4s, aes(x=(x.c+x.s)/2, y=5, label=gen), size=3) +
    geom_segment(data=tp2, aes(x=x.c,xend=x.s,y=lfc.c,yend=lfc.s,color=gen), size=1) +
    geom_point(aes(x=x, y=lfc, shape=gt), size = 2) +
    geom_text(data=tp3, aes(x=x-.05,y=lfc+.2,label=gt), hjust=1,vjust=0, size=3) +
    scale_x_continuous(breaks=xbreaks, labels=rep(c('control','stress'),2), expand=expansion(mult=c(.1,.1))) +
    scale_y_continuous(name='Normalized Read Count', expand=expansion(mult=c(.15,.15))) +
    scale_linetype_manual(values=c('solid','dashed')) +
    scale_shape_manual(name='Genotype', values=c(15,16)) +
    scale_color_manual(name='Generation', values=pal_npg()(3)) +
    scale_fill_manual(values=pal_simpsons()(3)) +
    facet_wrap(~pan, scale='free_y', nrow=3, dir='v') +
    otheme(xtick=T,xtext=T,ytick=T,ytext=T,ytitle=T,xgrid=F,ygrid=T,
           legend.pos='none', legend.dir='v', legend.title=T,
           legend.box='h') +
    theme(axis.text.x = element_text(angle=15))
fo = file.path(dirw, '21.cases.pdf')
ggsave(p, file=fo, width=9, height=5)
#}}}

#}}}


