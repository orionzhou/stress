source('functions.R')
#require(UpSetR)
dirw = file.path(dird, '15_de')
#{{{ functions
stresses = c("Control","Cold","Heat")
drcs = c('up', 'down')
run_matplotlib_venn3 <- function(s1, s2, s3, labs=LETTERS[1:3]) {
    #{{{
    require(venneuler)
    ti = tibble(a=c(str_c(s1,collapse=','), str_c(s2,collapse=','), str_c(s3,collapse=',')))
    fi = file.path('tmpi.txt')
    write_tsv(ti, fi, col_names=F)
    system("venn3.py coord tmpi.txt tmpo1.tsv tmpo2.tsv")
    tc = read_tsv('tmpo1.tsv', col_names=c('x','y','r')) %>% mutate(lab=labs)
    tl = read_tsv('tmpo2.tsv', col_names=c('x','y','cnt'))
    system("rm tmpi.txt tmpo1.tsv tmpo2.tsv")
    list(tc = tc, tl = tl)
    #}}}
}
get_venncounts <- function(v1, v2, v3) {
    #{{{
    require(limma)
    x1 = tibble(id = v1, v1 = T)
    x2 = tibble(id = v2, v2 = T)
    x3 = tibble(id = v3, v3 = T)
    x = x1 %>% full_join(x2, by='id') %>% full_join(x3, by='id') %>%
        replace_na(list(v1=F,v2=F,v3=F))
    x0 = as.data.frame(x[,-1])
    rownames(x0) = x$id
    vdc <- vennCounts(x0)
    class(vdc) = 'matrix'
    as_tibble(vdc[-1,]) %>%
        mutate(x=c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
            y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))
    #}}}
}
run_venn3 <- function(s1, s2, s3) {
    #{{{
    s1 = unique(s1); s2 = unique(s2); s3 = unique(s3)
    n1 = length(s1); n2 = length(s2); n3 = length(s3)
    n12 = sum(s1 %in% s2); n13=sum(s1 %in% s3); n23=sum(s2 %in% s3)
    n123 = sum((s1 %in% s2) && (s1 %in% s3))
    x = c(A=n1, B=n2, C=n3, "A&B"=n12, "A&C"=n13, "B&C"=n23 ,"A&B&C"=n123)
    l <- venneuler(x)
    tc = as_tibble(l$centers) %>%
        mutate(r = as.double(l$diameters)/2) %>%
        mutate(lab = l$labels)
    x1=tc$x[1];y1=tc$y[1];x2=tc$x[2];y2=tc$y[2];x3=tc$x[3];y3=tc$y[3]
    r1=tc$r[1];r2=tc$r[2];r3=tc$r[3]
    x123 = (x1+x2+x3)/3; y123 = (y1+y2+y3)/3
    tl = tibble()
    #}}}
}
#}}}

#{{{ read in
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
#
th = res$th
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ call DEGs
th1 = th %>% filter(Experiment=='HY') %>%
    distinct(SampleID, Genotype, Treatment, Timepoint) %>%
    mutate(cond = str_c(Genotype, Treatment, Timepoint, sep="_"))

tc1 = th1 %>% distinct(Genotype, Treatment, Timepoint, cond) %>%
    filter(Treatment != 'Control') %>%
    dplyr::rename(cond1 = cond) %>%
    mutate(TreatmentB = 'Control', TimepointB = Timepoint) %>%
    mutate(cond2 = str_c(Genotype, TreatmentB, TimepointB, sep="_")) %>%
    select(Genotype, Treatment, Timepoint, TreatmentB, TimepointB, cond1, cond2)
tc2 = th1 %>% distinct(Genotype, Treatment, Timepoint, cond) %>%
    filter(Treatment != 'Control') %>%
    dplyr::rename(cond1 = cond) %>%
    mutate(TreatmentB = 'Control', TimepointB = 0) %>%
    mutate(cond2 = str_c(Genotype, TreatmentB, TimepointB, sep="_")) %>%
    select(Genotype, Treatment, Timepoint, TreatmentB, TimepointB, cond1, cond2)
tc = rbind(tc1, tc2)

x = call_deg(th1, tm, tc)
fo = file.path(dirw, '01.de.rds')
saveRDS(x, fo)
#}}}

#{{{ call DEGs w. interaction effect
#{{{ prepare data, filter gids
conds = th %>% filter(Experiment=='HY') %>%
    distinct(Treatment, Timepoint) %>%
    mutate(Treatment = factor(Treatment, levels=stresses)) %>%
    mutate(Timepoint = factor(Timepoint, levels=c(0,1,25))) %>%
    mutate(cond = fct_cross(Treatment, Timepoint, sep='')) %>%
    mutate(cond = fct_reorder2(cond, Timepoint, Treatment, .desc = F)) %>%
    arrange(cond) %>% pull(cond)
th0 = th %>% filter(Experiment=='HY') %>%
    distinct(SampleID, Genotype, Treatment, Timepoint) %>%
    mutate(cond = str_c(Treatment, Timepoint, sep="")) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    dplyr::rename(gt = Genotype)
th0a = th0 %>% filter(Timepoint %in% c(0,1))
th0b = th0 %>% filter(Timepoint %in% c(0,25))
#
tm0 = tm %>% filter(SampleID %in% th0$SampleID)
gids = tm0 %>%
    group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
    filter(n.sam > .2 * nrow(th0)) %>% pull(gid)
tm0 = tm0 %>% filter(gid %in% gids)
#}}}
run_deseq_gt_cond <- function(th1, tm1) {
    #{{{ process data & run DESeq2
    require(DESeq2)
    vh = th1 %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    vm = tm1 %>% filter(SampleID %in% th1$SampleID) %>%
        select(SampleID, gid, ReadCount)
    vm.w = vm %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d,
                                 design = ~ gt + cond + gt:cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    dds
    #}}}
}
call_deg_gt3 <- function(dds, gts, stress, tp) {
    #{{{
    gt1 = gts[1]; gt2 = gts[2]; gt3 = gts[3]
    c1 = sprintf("gt%s.cond%s%s", gt2, stress, tp)
    c2 = sprintf("gt%s.cond%s%s", gt3, stress, tp)
    res1 = results(dds, name=c1)
    res2 = results(dds, name=c2)
    res3 = results(dds, contrast = list(c2, c1))
    stopifnot(rownames(res1)==gids && rownames(res2)==gids && rownames(res3)==gids)
    #
    r1 = tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1)) %>% mutate(qry=gt2, tgt=gt1)
    r2 = tibble(gid = gids, padj = res2$padj, log2fc = res2$log2FoldChange) %>%
        replace_na(list(padj = 1)) %>% mutate(qry=gt3, tgt=gt1)
    r3 = tibble(gid = gids, padj = res3$padj, log2fc = res3$log2FoldChange) %>%
        replace_na(list(padj = 1)) %>% mutate(qry=gt3, tgt=gt2)
    rbind(r1,r2,r3) %>% mutate(stress=!!stress,Timepoint=tp) %>%
        select(stress,Timepoint, qry, tgt, everything())
    #}}}
}

#{{{ call cold/heat 1/25h gt-specific DEGs vs control_0
th1 = th0a; tm1 = tm0
dds = run_deseq_gt_cond(th1, tm1)
res11 = call_deg_gt3(dds, gts3, 'Cold', 1)
res21 = call_deg_gt3(dds, gts3, 'Heat', 1)
#
th1 = th0b; tm1 = tm0
dds = run_deseq_gt_cond(th1, tm1)
res12 = call_deg_gt3(dds, gts3, 'Cold', 25)
res22 = call_deg_gt3(dds, gts3, 'Heat', 25)
#
res1 = rbind(res11,res12,res21,res22) %>% mutate(cond2 = 'control0')
#}}}
#{{{ call cold/heat 1/25h gt-specific DEGs vs control_1/25
th1 = th0a %>% mutate(cond=fct_relevel(cond, 'Control1')); tm1 = tm0
dds = run_deseq_gt_cond(th1, tm1)
res11 = call_deg_gt3(dds, gts3, 'Cold', 1)
res21 = call_deg_gt3(dds, gts3, 'Heat', 1)
#
th1 = th0b %>% mutate(cond=fct_relevel(cond, 'Control25')); tm1 = tm0
dds = run_deseq_gt_cond(th1, tm1)
res12 = call_deg_gt3(dds, gts3, 'Cold', 25)
res22 = call_deg_gt3(dds, gts3, 'Heat', 25)
#
res2 = rbind(res11,res12,res21,res22) %>% mutate(cond2 = 'controlM')
#}}}
#
res = rbind(res1, res2) %>%
    mutate(cond2 = ifelse(cond2 == 'control0', 'c0', 'cm')) %>%
    pivot_wider(names_from=cond2, values_from=c(padj,log2fc), names_sep=".") %>%
    group_by(stress,Timepoint,qry,tgt) %>% nest() %>% ungroup()
res %>% unnest(data) %>%
    dplyr::rename(padj=padj.c0, log2fc=log2fc.c0) %>%
    #dplyr::rename(padj=padj.cm, log2fc=log2fc.cm) %>%
    mutate(d = ifelse(padj < 0.05, ifelse(log2fc <= -1, "-",
               ifelse(log2fc >= 1, "+", "ND")), "ND")) %>%
    dplyr::count(stress, Timepoint, qry,tgt, d) %>%
    spread(d, n)

fo = file.path(dirw, '02.de.gt.rds')
saveRDS(res, fo)
#}}}

#{{{ merge called DEGs
#{{{ deg48
fi = file.path(dirw, '01.de.rds')
x = readRDS(fi)
get_de_gids <- function(ti, padj=.05, log2fc=1, drc='all') {
    #{{{
    to = ti %>% filter(padj < !!padj, abs(log2fc) >= !!log2fc)
    if( drc == 'up' ) to = to %>% filter(log2fc > 0)
    else if( drc == 'down' ) to = to %>% filter(log2fc < 0)
    to %>% pull(gid)
    #}}}
}
#
deg48 = x %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(Treatment = factor(Treatment, levels=stresses)) %>%
    mutate(cond2 = ifelse(TimepointB==0, 'time0', 'timeM')) %>%
    select(-cond1,-TreatmentB,-TimepointB) %>%
    mutate(up = map(ds, get_de_gids, padj=.05, log2fc=1, drc='up')) %>%
    mutate(down = map(ds, get_de_gids, padj=.05, log2fc=1, drc='down')) %>%
    arrange(Treatment,Genotype,Timepoint,cond2)
#}}}

#{{{ deg12
merge_ds <- function(ds.q, ds.t, ds) {
    #{{{
    ds.q %>% rename(padj.q=padj, log2fc.q=log2fc) %>%
    inner_join(ds.t, by='gid') %>% rename(padj.t=padj, log2fc.t=log2fc) %>%
    inner_join(ds, by='gid')
    #}}}
}
x1 = deg48 %>% filter(cond2=='timeM') %>%
    dplyr::rename(stress = Treatment, gt=Genotype) %>%
    select(stress, gt, Timepoint, ds)
fi = file.path(dirw, '02.de.gt.rds')
x2 = readRDS(fi) %>%
    mutate(qry = factor(qry, levels=gts3)) %>%
    mutate(tgt = factor(tgt, levels=gts3)) %>%
    mutate(stress = factor(stress, levels=stresses))
#
deg12 = x2 %>%
    rename(gt=qry) %>% inner_join(x1, by=c("gt","stress","Timepoint")) %>%
    rename(ds.q = ds) %>% rename(qry = gt) %>%
    rename(gt=tgt) %>% inner_join(x1, by=c("gt",'stress','Timepoint')) %>%
    rename(ds.t = ds) %>% rename(tgt = gt) %>%
    mutate(ds = pmap(list(ds.q, ds.t, data), merge_ds)) %>%
    select(qry, tgt, stress,Timepoint,ds) %>%
    mutate(qry = factor(qry, levels=gts3)) %>%
    mutate(tgt = factor(tgt, levels=gts3))
#}}}

res = list(deg48=deg48, deg12=deg12)
fo = file.path(dirw, '05.rds')
saveRDS(res, fo)
#}}}

____plot DEGs___

fi = file.path(dirw, '05.rds')
x = readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12

#{{{ bar plot
td1 = deg48 %>%
    select(Genotype,Treatment,Timepoint,cond2,up,down) %>%
    gather(drc, gids, -Treatment,-Genotype,-Timepoint,-cond2) %>%
    spread(cond2, gids) %>%
    dplyr::rename(gids0 = time0, gids1 = timeM) %>%
    mutate(gids = map2(gids0, gids1, intersect)) %>%
    mutate(n0 = map_int(gids0, length)) %>%
    mutate(n1 = map_int(gids1, length)) %>%
    mutate(n = map_int(gids, length))
#
tp = td1 %>% filter(Genotype %in% gts3) %>% mutate(n0 = n0-n, n1=n1-n) %>%
    select(Genotype,Treatment,Timepoint,drc,ctrl0_uniq=n0,ctrlx_uniq=n1,
           ctrl0_share=n, ctrlx_share=n) %>%
    gather(opt, n, -Genotype,-Treatment,-Timepoint,-drc) %>%
    separate(opt, c('opt1','opt2'), sep='_') %>%
    mutate(x = ifelse(opt1=='ctrlx', 2, 1)) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    mutate(y = ifelse(drc=='down', -sqrt(n), sqrt(n))) %>%
    mutate(opt2 = factor(opt2, levels=c("uniq","share")))
off = 6
tpt1 = tp %>% filter(opt2=='share') %>%
    mutate(yt = y/2) %>%
    mutate(lab=number(abs(n), big.mark=',', accuracy=1))
tpt2 = tp %>% group_by(Genotype,cond,drc,x) %>%
    summarise(n=sum(n), y=sum(y)) %>% ungroup() %>%
    mutate(lab=number(abs(n), big.mark=',', accuracy=1)) %>%
    mutate(yt = ifelse(drc=='down', y-off, y+off))
tpa = tibble(x=.35, y=c(20,-20), color=pal_d3()(2)[c(2,1)], lab=c('up','down')) %>%
    mutate(yt = ifelse(lab=='down', y-off, y+off)) %>%
    mutate(cond = 'Cold_1h') %>% crossing(Genotype=gts3)
labs = c('time0_control', 'timeM_control')

pv = ggplot(tp) +
    geom_bar(aes(x=x, y=y, fill=opt2), stat='identity', position='stack', alpha=1, size=0, color='white', width=.8) +
    geom_hline(yintercept = 0, size=.1) +
    geom_text(data=tpt1, aes(x=x, y=yt, label=lab), color='white', size=2.5) +
    geom_text(data=tpt2, aes(x=x, y=yt, label=lab), color='black', size=2.5) +
    geom_segment(data=tpa, aes(x=x, xend=x, y=0, yend=y, color=lab), arrow=arrow(length=unit(.1,'inches'), type='closed')) +
    geom_text(data=tpa, aes(x=x, y=yt, label=lab, color=lab), size=2.5, hjust=0) +
    scale_x_continuous(breaks=c(1,2),labels=labs, expand=expansion(mult=c(.1,.1))) +
    scale_fill_manual(values=pal_jco()(4)[c(3,4)]) +
    scale_color_manual(values=pal_d3()(4)[c(1,2)]) +
    facet_grid(Genotype ~ cond, switch='y', scales='free_x', space='free_x') +
    #otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4,
    otheme(legend.pos='none', strip.style='white', panel.spacing='.1',
           xtext=T, xtick=T) +
    theme(axis.text.x = element_text(angle = 20, hjust=.8, vjust=1.1)) +
    guides(color=F,fill=F)
fo = file.path(dirw, '10.deg.pdf')
ggsave(pv, filename=fo, width=5, height=8)
saveRDS(pv, file.path(dirf, 'f.1b.rds'))
#}}}

#{{{ ## time0 v timeM: venn graph
tv = tibble(x=c(-1,1), y=c(0,0), lab=labs) %>% mutate(lab=factor(lab,levels=labs))
tp = td1 %>% mutate(cond = sprintf("%s%d_%s", Treatment,Timepoint,drc)) %>%
    select(Genotype,cond,n0,n1,n) %>% crossing(tv)
pv = ggplot(tp) +
    geom_circle(aes(x0=x, y0=y, r=1.5, fill=lab), alpha=.3, size=0, color='grey') +
    geom_text(tp, mapping=aes(x=-1,y=0,label=n0-n), size=2.5) +
    geom_text(tp, mapping=aes(x=1,y=0,label=n1-n), size=2.5) +
    geom_text(tp, mapping=aes(x=0,y=0,label=n), size=2.5) +
    facet_grid(Genotype ~ cond, switch = 'y') +
    scale_fill_manual(values=pal_npg()(2)) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4)
fo = file.path(dirw, '10.deg.venn.pdf')
ggsave(pv, filename=fo, width=12, height=8)
#}}}

#{{{ simple gt venn graph
drcs = c('up','down')
td2 = td1 %>% mutate(cond=sprintf("%s_%dh", Treatment, Timepoint)) %>%
    select(Genotype, cond, drc,  gids) %>%
    filter(Genotype %in% gts3) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    #inner_join(td0, by=c('Genotype','cond','gid')) %>%
    group_by(Genotype, cond) %>%
    summarise(gidsU=list(gid[drc=='up']), gidsD=list(gid[drc=='down'])) %>%
    ungroup() %>%
    dplyr::rename(up=gidsU, down=gidsD) %>%
    gather(drc, gids, -Genotype, -cond) %>%
    mutate(drc = factor(drc, levels=c(drcs)))
td3 = td2 %>% select(Genotype,cond,drc, gids) %>%
    spread(Genotype, gids) %>%
    dplyr::rename(gidsB=B73, gidsM=Mo17, gidsW=W22) %>%
    mutate(vdc = pmap(list(gidsB, gidsM, gidsW), run_matplotlib_venn3, labs=gts3)) %>%
    mutate(tc = map(vdc, 'tc'), tl = map(vdc, 'tl'))

cols3 = pal_tron()(5)[c(1,2,3)]
#{{{
tp = td3 %>% filter(drc=='up') %>%
    mutate(pan = str_c(cond,str_to_upper(drc), sep='_')) %>%
    arrange(drc,cond)
tp = tp %>% mutate(pan = factor(pan, levels=tp$pan))
tpc = tp %>% select(pan, tc) %>% unnest(tc)
tpl = tp %>% select(pan, tl) %>% unnest(tl)
p3a = ggplot(tpc) +
    geom_circle(aes(x0=x, y0=y, r=r, fill=lab, color=lab), alpha=.2, size=1) +
    geom_text(tpl, mapping=aes(x=x,y=y,label=cnt), size=2.5) +
    facet_wrap(~ pan, scale='free', ncol=2) +
    scale_color_manual(values=cols3) +
    scale_fill_manual(values=cols3) +
    #scale_fill_manual(values=pal_tron()(5)) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4)
fo = sprintf('%s/11.gt.venn.pdf', dirw)
ggsave(p3a, filename=fo, width=4, height=4)
#}}}
#{{{
tp = td3 %>% mutate(pan = str_c(cond,str_to_upper(drc), sep='_')) %>%
    arrange(drc,cond)
tp = tp %>% mutate(pan = factor(pan, levels=tp$pan))
tpc = tp %>% select(pan, tc) %>% unnest(tc)
tpl = tp %>% select(pan, tl) %>% unnest(tl)
pv = ggplot(tpc) +
    geom_circle(aes(x0=x, y0=y, r=r, fill=lab, color=lab), alpha=.2, size=1) +
    geom_text(tpl, mapping=aes(x=x,y=y,label=cnt), size=2.5) +
    facet_wrap(~ pan, scale='free', ncol=4) +
    scale_color_manual(values=cols3) +
    scale_fill_manual(values=cols3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4)
fo = sprintf('%s/11.gt.venn.all.pdf', dirw)
ggsave(pv, filename=fo, width=8, height=4.5)
#}}}
#}}}

#{{{ log2fc heatmap of DEGs in 1/2 genotypes
drcs = c('up','down')
ty = td1 %>% mutate(cond=sprintf("%s_%dh", Treatment, Timepoint)) %>%
    select(Genotype, cond, drc,  gids) %>%
    filter(Genotype %in% gts3) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    group_by(Genotype, cond) %>%
    summarise(gidsU=list(gid[drc=='up']), gidsD=list(gid[drc=='down'])) %>%
    ungroup() %>%
    dplyr::rename(up=gidsU, down=gidsD) %>%
    gather(drc, gids, -Genotype, -cond) %>%
    mutate(drc = factor(drc, levels=c(drcs))) %>%
    unnest(gids) %>% rename(gid=gids) %>%
    mutate(DE=T) %>%
    spread(Genotype, DE) %>% replace_na(list(B73=F,Mo17=F,W22=F)) %>%
    mutate(nDE = B73+Mo17+W22) %>%
    mutate(tag = str_c(as.integer(B73),as.integer(Mo17),as.integer(W22))) %>%
    select(-B73,-Mo17,-W22)

tf = deg48 %>% mutate(cond=sprintf("%s_%dh", Treatment, Timepoint)) %>%
    filter(cond2 == 'time0', Genotype %in% gts3) %>%
    select(Genotype, cond, ds) %>% unnest(ds) %>%
    select(-padj)
#{{{ order by hc
tyw = ty %>% inner_join(tf, by=c("cond",'gid')) %>%
    select(-nDE) %>% spread(Genotype, log2fc) %>%
    group_by(cond,drc,tag) %>%
    nest() %>% ungroup() %>%
    mutate(gid = map(data, hc_order_row, cor.opt='euclidean')) %>%
    select(cond, drc, tag, gid) %>% unnest(gid) %>%
    mutate(y = n():1)
ty2 = ty %>%
    inner_join(tyw, by=c('cond','drc','tag','gid')) %>%
    arrange(cond, drc, desc(nDE), desc(tag), y) %>%
    group_by(cond,drc) %>% mutate(y = 1:n()) %>% ungroup()
ty2 %>% count(cond, drc, nDE) %>% print(n=4)
#}}}

fmax=5
drc = 'up'
tp = ty2 %>% filter(drc==!!drc) %>%
    inner_join(tf, by=c("cond",'gid')) %>%
    mutate(log2fc = ifelse(log2fc < -fmax, -fmax, log2fc)) %>%
    mutate(log2fc = ifelse(log2fc > fmax, fmax, log2fc))
tpy = ty2 %>% group_by(cond, drc, tag) %>%
    summarise(ymin=min(y), ymax=max(y), ymid=(ymin+ymax)/2, n=n()) %>%
    ungroup() %>%
    mutate(tag = sprintf("%s (%d)", tag, n)) %>% filter(drc == !!drc)
#
p1 = ggplot() +
    geom_tile(tp, mapping=aes(x=Genotype, y=y, fill=log2fc)) +
    geom_segment(tpy, mapping=aes(x=.4,xend=.4,y=ymin,yend=ymax), lineend='round',size=.4, color='royalblue') +
    geom_segment(tpy, mapping=aes(x=.32,xend=.4,y=ymin,yend=ymin), lineend='round', size=.3, color='royalblue') +
    geom_segment(tpy, mapping=aes(x=.32,xend=.4,y=ymax,yend=ymax), lineend='round', size=.3, color='royalblue') +
    geom_text(tpy, mapping=aes(x=.3,y=ymid,label=tag), size=2.5, hjust=1, color='gray30') +
    scale_x_discrete(expand=expansion(mult=c(.9,.05))) +
    scale_y_reverse(expand=expansion(mult=c(.005,.005))) +
    #scale_fill_viridis(name = 'log2fc', direction=-1) +
    scale_fill_gradientn(name='log2(FoldChange)', colors=cols100) +
    facet_wrap(~cond, scale='free',nrow=1) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           panel.border=F, xtick=T, xtext=T, legend.vjust=-.4,
           margin=c(1,.1,.1,0), strip.style = 'dark') +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7))
fp = sprintf('%s/12.heatmap.pdf', dirw)
ggsave(p1, file = fp, width = 8, height = 8)

#{{{ old
tp0 = td3 %>% filter(!is.na(BMW))
dmap = c('-1'='-', '0'='.', '1'='+')
tpw = tp0 %>% select(cond,gid,gt,deg,log2fc, BMW0, BMW) %>%
    mutate(deg = dmap[as.character(deg)]) %>%
    pivot_wider(names_from=gt, values_from=c(deg, log2fc), names_sep='.') %>%
    mutate(grp=str_c(deg.B73, deg.Mo17, deg.W22)) %>%
    select(-deg.B73, -deg.Mo17, -deg.W22)
tpw %>% count(cond, grp) %>%
    complete(cond, nesting(grp), fill=list(n=0)) %>% spread(cond, n) %>%
    arrange(desc(grp)) %>% print(n=30)

tpw %>% count(cond, grp, BMW0) %>%
    filter(cond == 'Heat25') %>% select(-cond) %>%
    complete(grp, nesting(BMW0), fill=list(n=0)) %>%
    spread(BMW0, n) %>%
    arrange(desc(grp)) %>% print(n=30)

e = tpw %>% select(-cond)
mat = data.frame(e)
rownames(mat) = tpw$cond
colnames(mat) = colnames(e)
opt.dist = 'euclidean'; opt.hc = 'complete'#ward.D'
edist = idist(mat, method = opt.dist)
ehc = hclust(edist, method = opt.hc)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]

lnames = tpw$gid
cols100b=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)
tp = tp0 %>% mutate(gid = str_c(cond, grp, gid, sep='_')) %>%
    mutate(gid = factor(gid, levels=rev(lnames)))
p1 = ggplot(tp, aes(x=gt, y=gid, fill=log2fc)) +
    geom_tile() +
    scale_x_discrete(expand=expansion(mult=c(.05,.05))) +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    #scale_fill_viridis(name = 'log2fc', direction=-1) +
    scale_fill_gradientn(name='log2fc', colors=cols100b) +
    facet_wrap(~cond, scale='free') +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           xtick=T, xtext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7))
fp = sprintf('%s/21.heatmap.pdf', dirw)
ggsave(p1, file = fp, width = 8, height = 10)

    tp = th %>% mutate(taxa = SampleID) %>%
        select(taxa, everything())
    p1 = ggtree(tree, layout = 'rectangular') +
        scale_x_continuous(expand = expansion(mult=c(1e-2,expand.x))) +
        scale_y_discrete(expand = expansion(mult=c(.01,.01)))
    p1 = p1 %<+%
        tp + geom_tiplab(aes(label=get(var.lab), color=get(var.col)), size=lab.size) +
        get(str_c('scale_color', pal.col, sep="_"))() +
        guides(color=F)
#}}}
#}}}

#{{{ genotypic difference in cold/heat response
#{{{ characterize pairwise contrasts
fmax=10
drcs = c("\u2191",'=',"\u2193")
drcs = c("+",'=',"-")
tx = tibble(deg = c("A+B+",'A+B=','A=B+', 'A-B-','A-B=','A=B-', 'A+B-','A-B+','A=B='),
             x = c(1:3, 5:7, 9:11))
td = deg12 %>%
    unnest(ds) %>%
    mutate(deg.q = padj.q < .05 & abs(log2fc.q) >= 1) %>%
    mutate(deg.t = padj.t < .05 & abs(log2fc.t) >= 1) %>%
    mutate(ddeg = padj.c0<.05 & abs(log2fc.c0)>=1 & padj.cm<.05 & abs(log2fc.cm)>=1) %>%
    mutate(drc.q = ifelse(log2fc.q < 0, "-", "+")) %>%
    mutate(drc.t = ifelse(log2fc.t < 0, "-", "+")) %>%
    mutate(deg.qt = ifelse(deg.q, ifelse(deg.t, 'q+t', 'q'),
                    ifelse(deg.t, 't', 'none'))) %>%
    mutate(ddrc = ifelse(ddeg, ifelse(log2fc.cm<0, -1, 1), 0)) %>%
    mutate(deg.qt = factor(deg.qt, levels=c("q",'t','q+t','none'))) %>%
    mutate(deg = NA) %>%
    mutate(deg = ifelse(deg.qt=='none', 'A=B=', deg)) %>%
    mutate(deg = ifelse(deg.qt=='q+t', ifelse(drc.q=='-',
                        ifelse(drc.t=='-', 'A-B-', 'A-B+'),
                        ifelse(drc.t=='-', 'A+B-', 'A+B+')), deg)) %>%
    mutate(deg = ifelse(deg.qt=='q', ifelse(drc.q=='-','A-B=','A+B='), deg)) %>%
    mutate(deg = ifelse(deg.qt=='t', ifelse(drc.t=='-','A=B-','A=B+'), deg)) %>%
    inner_join(tx, by='deg') %>%
    mutate(deg = factor(deg, levels=tx$deg)) %>%
    mutate(log2fc.q = map_dbl(log2fc.q, set_bound, minV=-fmax, maxV=fmax)) %>%
    mutate(log2fc.t = map_dbl(log2fc.t, set_bound, minV=-fmax, maxV=fmax)) %>%
    mutate(cond = str_c(stress,Timepoint,qry,tgt, sep="_")) #%>%
    #filter(! (deg=='0' & ddeg == 'no-diff-deg'))
td %>% count(qry,tgt,stress,Timepoint,ddrc,drc.q,drc.t) %>%
    print(n=50)
#}}}
fo = file.path(dirw, '08.de.rds')
#saveRDS(td, fo)

#{{{ multi-panel scatter plot
tp = td %>% filter(ddrc != 0, deg != 'A=B=')
linecol = 'azure3'; lty = 'solid'
cols9 = brewer.pal(10, "Paired")[c(2,4,6,8,9,7,5,3,1)]
cols9 = pal_npg()(10)[c(1:7,10)]
#{{{ scatter plot
p = ggplot(tp, aes(x=log2fc.q, y=log2fc.t)) +
    geom_vline(xintercept=0, linetype=lty, color=linecol) +
    geom_hline(yintercept=0, linetype=lty, color=linecol) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    geom_point(aes(color=deg, shape='k'), size=1.5) +
    scale_x_continuous(name='Genotype (A) log2fc', limits=c(-fmax,fmax),expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Genotype (B) log2fc', limits=c(-fmax,fmax),expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='DEG status:', values=cols9) +
    scale_shape_manual(name='Genotype Effect:', labels=c("negative","positive"), values=c(1,4)) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=1), shape=F)
fp = sprintf('%s/15.ddeg.pdf', dirw)
ggsave(p, file = fp, width = 8, height = 9.5)
#}}}
#{{{ bar plot showing counts
tp1 = tp %>% count(cond, x, deg)
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
p3b = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=deg), width=.8, stat='identity') +
    geom_text(data=tp1s,aes(x=11,y=500,label=lab), size=3, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+10, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tx$x, labels=tx$deg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.1))) +
    scale_fill_manual(name='DEG status:', values=cols9) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=F, ytitle=T, ytick=T, ytext=T) +
    theme(axis.text.x = element_text(angle=40, size=7, vjust=1.2, hjust=1))
fp = sprintf('%s/15.ddeg.cnt.pdf', dirw)
ggsave(p3b, file = fp, width = 6, height = 6)
#}}}

#{{{ plot mean CPM for each category of genes
te = th %>% filter(Experiment=='HY', Genotype %in% gts3) %>%
    select(SampleID, Genotype, Treatment, Timepoint) %>%
    inner_join(tm[,c('gid','SampleID','CPM')], by='SampleID') %>%
    mutate(cond = str_c(Treatment,Timepoint,sep='')) %>%
    group_by(Genotype,cond, gid) %>%
    summarise(CPM=mean(CPM)) %>% ungroup() %>% rename(gt=Genotype)

tp1 = td %>% filter(ddrc != 0) %>%
    mutate(cond=str_c(stress,Timepoint,sep='')) %>%
    select(qry,tgt,cond,gid,deg) %>%
    inner_join(te, by=c('cond'='cond','qry'='gt','gid'='gid')) %>%
    rename(cpm.q1 = CPM) %>%
    inner_join(te, by=c('cond'='cond','tgt'='gt','gid'='gid')) %>%
    rename(cpm.t1 = CPM)
tp2 = td %>% filter(ddrc != 0) %>%
    mutate(cond=str_c('Control',Timepoint,sep='')) %>%
    select(qry,tgt,cond,gid,deg) %>%
    inner_join(te, by=c('cond'='cond','qry'='gt','gid'='gid')) %>%
    rename(cpm.q0 = CPM) %>%
    inner_join(te, by=c('cond'='cond','tgt'='gt','gid'='gid')) %>%
    rename(cpm.t0 = CPM) %>% select(-cond)

tp = tp1 %>% inner_join(tp2, by=c("qry",'tgt','gid','deg')) %>%
    mutate(cpm = pmin(cpm.q1, cpm.t1, cpm.q0, cpm.t0)) %>%
    mutate(cpm = asinh(cpm)) %>%
    mutate(cond = str_c(cond, qry, tgt, sep="_"))
p = ggplot(tp, aes(x=deg, y=cpm)) +
    geom_violin(aes(fill=deg)) +
    scale_x_discrete(expand=expansion(mult=c(.10,.10))) +
    scale_y_continuous(name='asinh(CPM)', expand=expansion(mult=c(.05,.05))) +
    scale_fill_manual(name='DEG status:', values=cols9) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=F, ytitle=T, ytick=T, ytext=T) +
    theme(axis.text.x = element_text(angle=40, size=7, vjust=1.2, hjust=1))
fp = sprintf('%s/15.ddeg.cpm.pdf', dirw)
ggsave(p, file = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ # infer 3-way gt comparison
fc = file.path(dirw, '18.ddeg.cat.tsv')
tc = read_tsv(fc)
td2 = td %>% mutate(cond=str_c(stress,Timepoint)) %>%
    mutate(comp = str_c(qry,tgt,sep="_")) %>%
    select(cond,comp,gid,ddrc) %>%
    spread(comp, ddrc) %>%
    rename(MB=Mo17_B73, WB=W22_B73, WM=W22_Mo17) %>%
    inner_join(tc, by=c("MB","WB",'WM')) %>%
    select(cond,gid,MB,WB,WM,BMW0,BMW)
td2s = td2 %>% filter(!is.na(BMW0), BMW0 != '000') %>%
    select(cond, gid, BMW0, BMW)
td3 = deg48 %>% filter(cond2=='timeM', Genotype %in% gts3) %>%
    mutate(cond = str_c(Treatment, Timepoint, sep='')) %>%
    select(gt=Genotype,cond, ds) %>%
    unnest(ds) %>%
    mutate(deg = ifelse(padj < .05, ifelse(log2fc <= -1, -1,
                 ifelse(log2fc >= 1, 1, 0)), 0)) %>%
    inner_join(td2s, by=c("cond","gid"))

ddrcs = c(1,0,-1)
ac = crossing(cond=unique(td2$cond), MB=ddrcs, WB=ddrcs, WM=ddrcs)
y = td2 %>% count(cond, MB, WB, WM) %>%
    right_join(ac, by=c('cond','MB','WB','WM')) %>%
    replace_na(list(n=0)) %>%
    spread(cond,n) %>%
    arrange(desc(MB), desc(WB), desc(WM)) %>%
    print(n=40)
fo = file.path(dirw, '18.ddeg.cat.fill.tsv')
write_tsv(y, fo)

bmws = c("001",'010','100','011','101','110')
z = td2 %>% filter(!is.na(BMW)) %>%
    count(cond, BMW) %>%
    mutate(BMW = factor(BMW, levels=bmws)) %>%
    complete(cond, nesting(BMW), fill=list(n=0)) %>%
    spread(cond,n) %>% print(n=10)
#}}}
#}}}

#{{{ create gene status table for ML evaluation
fi = file.path(dirw, '05.rds')
x = readRDS(fi)
deg48 = x$deg48
#
td1 = deg48 %>%
    select(gt=Genotype,cond=Treatment,time=Timepoint,cond2,ds) %>%
    unnest(ds) %>%
    mutate(deg = ifelse(padj<.05 & abs(log2fc)>=1, ifelse(log2fc < 0, 'D','U'), 'N')) %>%
    select(-padj, -log2fc) %>%
    spread(cond2, deg) %>%
    mutate(st = 'N') %>%
    mutate(st = ifelse(time0=='U' & timeM=='U', 'U', st)) %>%
    mutate(st = ifelse(time0=='D' & timeM=='D', 'D', st)) %>%
    mutate(cond = str_to_lower(cond)) %>%
    select(cond,time,gt,gid,st)
ton = gcfg$gene %>% filter(! gid %in% td1$gid) %>%
    select(gid) %>% mutate(st = 'zero') %>%
    crossing(td1 %>% distinct(cond, time, gt))
td1 = td1 %>% bind_rows(ton) %>%
    mutate(st = factor(st, levels=c("U",'D','N', 'zero')))
td1 %>% count(cond,time,gt,st) %>% spread(st,n)

fi = file.path(dirw, '08.de.rds')
ti = readRDS(fi)
ti2 = ti %>% select(-cond) %>% rename(cond=stress, time=Timepoint) %>%
    mutate(cond=str_to_lower(cond))
ddegs=c("A+B+",'A+B=','A=B+', 'A-B-','A-B=','A=B-', 'A+B-','A-B+','A=B=')
ddegs1 = c(str_c('d',ddegs), str_c('n',ddegs))
td2 = ti2 %>% mutate(st1 = ifelse(ddeg, 'd', 'n')) %>%
    mutate(st = str_c(st1, deg)) %>%
    select(cond,time, qry,tgt, gid,st) %>%
    mutate(st = factor(st, levels = ddegs1))
td2 %>% count(cond,time, qry,tgt,st) %>% spread(st,n)

r = list(td1=td1, td2=td2)
fo = file.path(dirw, '09.gene.status.rds')
saveRDS(r, fo)
#}}}

#{{{ share DEG list
to = td2 %>% mutate(n_deg = map_int(gids, length)) %>%
    mutate(gid = map_chr(gids, str_c, collapse=',')) %>%
    select(-gids)
fo = file.path(dirw, '05.deg.tsv')
write_tsv(to, fo)
#}}}

________BELOW ARE DEPRECATED_________
#{{{ alluvial
require(ggalluvial)

gt = 'Mo17'
degs = c("up-DE",'down-DE','non-DE')
tp = t_ds %>% filter(gt == !!gt) %>%
    mutate(deg = ifelse(padj>.01, 'non-DE', ifelse(lfc<0, 'up-DE', 'down-DE'))) %>%
    select(gid,hac,deg) %>% spread(hac, deg) %>%
        count(hac005,hac029,hac048,hac053,hac077,hac149) %>%
        mutate(subject=1:length(n)) %>%
        gather(hac, deg, -subject, -n) %>%
        mutate(deg = factor(deg, levels=degs))
#
cols3 = c(pal_aaas()(2), 'gray')
p = ggplot(tp, aes(x=hac, stratum=deg, alluvium=subject, y=n, fill=deg, label=deg)) +
    scale_fill_manual(values=cols3) +
    #geom_flow(stat='alluvium', lode.guidance='rightleft', color='darkgray') +
        geom_flow() +
        geom_text(stat = 'stratum', size=2) +
    geom_stratum(alpha = .6) +
        otheme(xtext=T,ytext=T,xtick=T,ytick=T)
fo = sprintf('%s/43.%s.pdf', dirw, gt)
ggsave(p, file=fo, width=8, height=6)
#}}}

#{{{ [not working] call DEGs 2 [cold 1h vs. (control 1h + control 0h)]
th0 = th %>% filter(Experiment=='HY') %>%
    distinct(SampleID, Genotype, Treatment, Timepoint) %>%
    mutate(cond = str_c(Treatment, Timepoint, sep='_')) %>%
    mutate(cond_1h = ifelse(cond=='Control_0', 'Control_1', cond)) %>%
    mutate(cond_25h = ifelse(cond=='Control_0', 'Control_25', cond)) %>%
    mutate(cond_1h = str_c(Genotype, cond_1h, sep="_")) %>%
    mutate(cond_25h = str_c(Genotype, cond_25h, sep="_")) %>%
    select(-cond)
th1 = th0 %>% rename(cond=cond_1h)
th2 = th0 %>% rename(cond=cond_25h)

tc1 = th0 %>% filter(Timepoint==1) %>% mutate(cond1 = cond_1h) %>%
    distinct(Genotype, Treatment, Timepoint, cond1) %>%
    filter(Treatment != 'Control') %>%
    mutate(cond2 = str_c(Genotype, 'Control', Timepoint, sep="_")) %>%
    select(Genotype, Treatment, Timepoint, cond1, cond2)
tc2 = th0 %>% filter(Timepoint==25) %>% mutate(cond1 = cond_25h) %>%
    distinct(Genotype, Treatment, Timepoint, cond1) %>%
    filter(Treatment != 'Control') %>%
    mutate(cond2 = str_c(Genotype, 'Control', Timepoint, sep="_")) %>%
    select(Genotype, Treatment, Timepoint, cond1, cond2)

x1 = call_deg(th1, tm, tc1)
x2 = call_deg(th2, tm, tc2)


td2 = x2 %>% unnest(ds) %>% filter(padj < .05, abs(log2fc) >= 1) %>%
    mutate(drc = ifelse(log2fc < 0, 'down', 'up')) %>%
    group_by(Genotype,Treatment,Timepoint,drc) %>%
    summarise(gids = list(gid)) %>% ungroup() %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(Treatment = factor(Treatment, levels=stresses))

tx = td1 %>% select(Genotype,Treatment,Timepoint,drc,gids) %>%
    dplyr::rename(gids0=gids) %>%
    inner_join(td2, by=c("Genotype","Treatment","Timepoint",'drc')) %>%
    dplyr::rename(gids1=gids) %>%
    mutate(gids = map2(gids0, gids1, intersect)) %>%
    mutate(n0 = map_int(gids0, length)) %>%
    mutate(n1 = map_int(gids1, length)) %>%
    mutate(n = map_int(gids, length))

fo = file.path(dirw, '01.de2.rds')
saveRDS(x, fo)
#}}}



