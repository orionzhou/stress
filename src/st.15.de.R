source('functions.R')
require(UpSetR)
require(limma)
dirw = file.path(dird, '15_de')
#{{{ functions
stresses = c("Control","Cold","Heat")
drcs = c('up', 'down')
set_bound <- function(x, minV, maxV) min(max(x,minV),maxV)
get_venncounts <- function(v1, v2, v3) {
    #{{{
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
#}}}

#{{{ read in
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

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
conds = crossing(stress=c("Control",'Cold','Heat'), tp=c(1,25)) %>%
    mutate(stress = factor(stress, levels=stresses)) %>% arrange(stress) %>%
    mutate(cond=str_c(stress,tp)) %>% pull(cond)
th0 = th %>% filter(Experiment=='HY') %>%
    distinct(SampleID, Genotype, Treatment, Timepoint) %>%
    filter(Timepoint != 0) %>%
    mutate(cond = str_c(Treatment, Timepoint, sep="")) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    dplyr::rename(gt = Genotype)
th0a = th0 %>% filter(Timepoint == 1)
th0b = th0 %>% filter(Timepoint == 25)
#
#{{{ filter gids
tm0 = tm %>% filter(SampleID %in% th0$SampleID)
gids = tm0 %>%
    group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
    filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
tm0 = tm0 %>% filter(gid %in% gids)
#}}}

#{{{ call cold/heat 1h DEGs
th1 = th0a; tm1 = tm0
#{{{ process data & run DESeq2
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
resultsNames(dds)
#}}}
#
#{{{ cold 1
res1 = results(dds, name="gtMo17.condCold1")
res2 = results(dds, name="gtW22.condCold1")
res3 = results(dds, contrast = list("gtW22.condCold1", "gtMo17.condCold1"))
stopifnot(rownames(res1)==gids && rownames(res2)==gids && rownames(res3)==gids)
#
r1 = tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'M_B')
r2 = tibble(gid = gids, padj = res2$padj, log2fc = res2$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_B')
r3 = tibble(gid = gids, padj = res3$padj, log2fc = res3$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_M')
res11 = rbind(r1,r2,r3) %>% mutate(stress='Cold',Timepoint=1) %>%
    select(stress,Timepoint, comp, everything())
#}}}
#{{{ heat 1
res1 = results(dds, name="gtMo17.condHeat1")
res2 = results(dds, name="gtW22.condHeat1")
res3 = results(dds, contrast = list("gtW22.condHeat1", "gtMo17.condHeat1"))
stopifnot(rownames(res1)==gids && rownames(res2)==gids && rownames(res3)==gids)
#
r1 = tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'M_B')
r2 = tibble(gid = gids, padj = res2$padj, log2fc = res2$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_B')
r3 = tibble(gid = gids, padj = res3$padj, log2fc = res3$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_M')
res21 = rbind(r1,r2,r3) %>% mutate(stress='Heat',Timepoint=1) %>%
    select(stress,Timepoint, comp, everything())
#}}}
#}}}

#{{{ call cold/heat 25h DEGs
th1 = th0b; tm1 = tm0
#{{{ process data & run DESeq2
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
resultsNames(dds)
#}}}
#
#{{{ cold 25
res1 = results(dds, name="gtMo17.condCold25")
res2 = results(dds, name="gtW22.condCold25")
res3 = results(dds, contrast = list("gtW22.condCold25", "gtMo17.condCold25"))
stopifnot(rownames(res1)==gids && rownames(res2)==gids && rownames(res3)==gids)
#
r1 = tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'M_B')
r2 = tibble(gid = gids, padj = res2$padj, log2fc = res2$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_B')
r3 = tibble(gid = gids, padj = res3$padj, log2fc = res3$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_M')
res12 = rbind(r1,r2,r3) %>% mutate(stress='Cold',Timepoint=25) %>%
    select(stress,Timepoint, comp, everything())
#}}}
#{{{ heat 25
res1 = results(dds, name="gtMo17.condHeat25")
res2 = results(dds, name="gtW22.condHeat25")
res3 = results(dds, contrast = list("gtW22.condHeat25", "gtMo17.condHeat25"))
stopifnot(rownames(res1)==gids && rownames(res2)==gids && rownames(res3)==gids)
#
r1 = tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'M_B')
r2 = tibble(gid = gids, padj = res2$padj, log2fc = res2$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_B')
r3 = tibble(gid = gids, padj = res3$padj, log2fc = res3$log2FoldChange) %>%
    replace_na(list(padj = 1)) %>% mutate(comp = 'W_M')
res22 = rbind(r1,r2,r3) %>% mutate(stress='Heat',Timepoint=25) %>%
    select(stress,Timepoint, comp, everything())
#}}}
#}}}

res = rbind(res11,res12,res21,res22)
res %>% dplyr::count(stress,Timepoint,comp)

res %>% mutate(d = ifelse(padj < 0.05, ifelse(log2fc <= -1, "-",
                   ifelse(log2fc >= 1, "+", "ND")), "ND")) %>%
    dplyr::count(stress, Timepoint, comp, d) %>%
    spread(d, n)

fo = file.path(dirw, '02.de.gt.rds')
saveRDS(res, fo)
#}}}

#{{{ merge called DEGs to 10.rds and obtain gt-specific DEG sets
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
merge_ds <- function(ds1, ds2, ds) {
    #{{{
    ds1 %>% rename(padj.1=padj, log2fc.1=log2fc) %>%
    inner_join(ds2, by='gid') %>% rename(padj.2=padj, log2fc.2=log2fc) %>%
    inner_join(ds, by='gid')
    #}}}
}
fi = file.path(dirw, '02.de.gt.rds')
gmap = c("B"="B73","M"="Mo17","W"="W22")
x2 = readRDS(fi) %>% group_by(stress,Timepoint,comp) %>% nest() %>% ungroup() %>%
    separate(comp, c("qry", "tgt"), sep="_") %>%
    mutate(qry = gmap[qry], tgt= gmap[tgt]) %>%
    mutate(qry = factor(qry, levels=gts3)) %>%
    mutate(tgt = factor(tgt, levels=gts3)) %>%
    mutate(stress = factor(stress, levels=stresses))
x1 = deg48 %>% filter(cond2=='timeM') %>%
    rename(stress = Treatment, gt=Genotype) %>%
    select(stress, gt, Timepoint, ds)
#
deg12 = x2 %>%
    rename(gt=qry) %>% inner_join(x1, by=c("gt","stress","Timepoint")) %>%
    rename(ds.q = ds) %>% rename(qry = gt) %>%
    rename(gt=tgt) %>% inner_join(x1, by=c("gt",'stress','Timepoint')) %>%
    rename(ds.t = ds) %>% rename(tgt = gt) %>%
    mutate(ds = pmap(list(ds.q, ds.t, data), merge_ds)) %>%
    select(qry, tgt, ,stress,Timepoint,ds) %>%
    mutate(qry = factor(qry, levels=gts3)) %>%
    mutate(tgt = factor(tgt, levels=gts3))
#}}}

res = list(deg48=deg48, deg12=deg12)
fo = file.path(dirw, '05.rds')
saveRDS(res, fo)
#}}}

#{{{ plot DEGs
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

labs = c('ctrl_0', 'ctrl_M')
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
    otheme(legend.pos='none', xtext=T, xtick=T) +
    guides(color=F,fill=F)
fo = file.path(dirw, '10.deg.pdf')
ggsave(pv, filename=fo, width=5, height=8)
#}}}

#{{{ [messy] venn graph
tv = tibble(x=c(-1,1), y=c(0,0), lab=labs) %>% mutate(lab=factor(lab,levels=labs))
tp = td1 %>% select(Genotype,cond1,n0,n1,n) %>% crossing(tv)
pv = ggplot(tp) +
    geom_circle(aes(x0=x, y0=y, r=1.5, fill=lab), alpha=.3, size=0, color='grey') +
    geom_text(tp, mapping=aes(x=-1,y=0,label=n0-n), size=2.5) +
    geom_text(tp, mapping=aes(x=1,y=0,label=n1-n), size=2.5) +
    geom_text(tp, mapping=aes(x=0,y=0,label=n), size=2.5) +
    facet_grid(Genotype ~ cond1, switch = 'y') +
    scale_fill_manual(values=pal_npg()(2)) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4)
fo = file.path(dirw, '10.deg.venn.pdf')
ggsave(pv, filename=fo, width=12, height=8)
#}}}

#{{{ genotypic difference in cold/heat response
fmax=10
tp = deg12 %>% unnest(ds) %>%
    mutate(deg1 = padj.1 < .05 & abs(log2fc.1) >= 1) %>%
    mutate(deg2 = padj.2 < .05 & abs(log2fc.2) >= 1) %>%
    mutate(deg12 = padj < .05 & abs(log2fc) >= 1) %>%
    mutate(deg = ifelse(deg1, ifelse(deg2, '1+2', '1'), ifelse(deg2, '2', '0')),
        ddeg = ifelse(deg12, 'diff-deg', 'no-diff-deg')) %>%
    mutate(log2fc.1 = map_dbl(log2fc.1, set_bound, minV=-fmax, maxV=fmax)) %>%
    mutate(log2fc.2 = map_dbl(log2fc.2, set_bound, minV=-fmax, maxV=fmax)) %>%
    mutate(cond = str_c(stress,Timepoint,qry,tgt, sep="_")) #%>%
    #filter(! (deg=='0' & ddeg == 'no-diff-deg'))
tp %>% count(qry,tgt,stress,Timepoint,deg,ddeg)

p = ggplot(tp, aes(x=log2fc.1, y=log2fc.2)) +
    geom_point(aes(color=ddeg, shape=deg), size=1.5) +
    scale_x_continuous(limits=c(-fmax,fmax),expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(limits=c(-fmax,fmax),expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='deg', values=pal_aaas()(10)[c(2,9)]) +
    scale_shape_manual(name='ddeg', values=c(1:4)) +
    facet_wrap(~cond, ncol=3) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           xtick=T, xtext=T, ytick=T, ytext=T)
fp = sprintf('%s/15.pdf', dirw)
ggsave(p, file = fp, width = 12, height = 12)

y = tp %>% mutate(cond=str_c(stress,Timepoint)) %>%
    mutate(comp = str_c(qry,tgt,sep="_")) %>%
    mutate(drc = ifelse(deg12, ifelse(log2fc<0, -1, 1), 0)) %>%
    select(cond,comp,gid,drc) %>%
    spread(comp, drc) %>%
    select(cond, MB=Mo17_B73, WB=W22_B73, WM=W22_Mo17, gid)
y %>% count(cond, MB, WB, WM) %>% arrange(cond, desc(MB), desc(WB), desc(WM)) %>% print(n=40) 
#}}}
#}}}

#{{{ [not work] genotypic difference for cold/heat - up in B,M,W
stress='Heat'; opt = 'ctrl0'
stress='Heat'; opt = 'ctrl'
stress='Heat'; opt = 'ctrlx'
stress='Cold'; opt = 'ctrl'
stress='Cold'; opt = 'ctrl0'
stress='Cold'; opt = 'ctrlx'
gids = td1 %>% filter(Genotype %in% gts3, Treatment==stress, drc=='up') %>%
    select(Genotype,Treatment,Timepoint,drc, gids) %>%
    unnest(gids) %>% rename(gid = gids) %>% distinct(gid) %>% pull(gid)
#
tp0 = x %>% filter(Genotype %in% gts3, Treatment==stress) %>%
    mutate(cond=sprintf("%dh_%s", Timepoint, Genotype)) %>%
    #mutate(cond=sprintf("%s_%dh", Genotype, Timepoint)) %>%
    mutate(condB = ifelse(TimepointB==0, 'ctrl0', 'ctrlx')) %>%
    select(cond, condB, ds) %>% unnest(ds) %>%
    filter(gid %in% gids) %>% select(-padj) %>%
    spread(condB, log2fc) %>%
    mutate(ctrl = ifelse(abs(ctrl0) < abs(ctrlx), ctrlx, ctrl0))
tp = tp0 %>% select(cond, gid, log2fc=!!opt) %>%
    mutate(log2fc = map_dbl(log2fc, set_bound, minV=-3, maxV=3))
tpl = tp %>% spread(cond, log2fc)
tpw = tp %>% spread(gid, log2fc)

e = tpl %>% select(-gid)
mat = data.frame(e)
rownames(mat) = tpl$gid
colnames(mat) = colnames(e)
fp = sprintf('%s/10.lfc.%s.%s.pdf', dirw, stress, opt)
pdf(fp, 6,12)
pheatmap(mat, cluster_cols=F, show_rownames=F)
dev.off()

e = tpw %>% select(-cond)
mat = data.frame(e)
rownames(mat) = tpw$cond
colnames(mat) = colnames(e)
opt.dist = 'euclidean'; opt.hc = 'complete'#ward.D'
edist = idist(mat, method = opt.dist)
ehc = hclust(edist, method = opt.hc)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]
#
cols100b=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)
tp = tp %>% mutate(gid = factor(gid, levels=rev(lnames)))
p1 = ggplot(tp, aes(x=cond, y=gid, fill=log2fc)) +
    geom_tile() +
    scale_x_discrete(expand=expansion(mult=c(.05,.05))) +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    #scale_fill_viridis(name = 'log2fc', direction=-1) +
    scale_fill_gradientn(name='log2fc', colors=cols100b) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           xtick=T, xtext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7))
fp = sprintf('%s/10.lfc.%s.%s.pdf', dirw, stress, opt)
ggsave(p1, file = fp, width = 6, height = 12)


    tp = th %>% mutate(taxa = SampleID) %>%
        select(taxa, everything())
    p1 = ggtree(tree, layout = 'rectangular') +
        scale_x_continuous(expand = expansion(mult=c(1e-2,expand.x))) +
        scale_y_discrete(expand = expansion(mult=c(.01,.01)))
    p1 = p1 %<+%
        tp + geom_tiplab(aes(label=get(var.lab), color=get(var.col)), size=lab.size) +
        get(str_c('scale_color', pal.col, sep="_"))() +
        guides(color=F)

#{{{ [obsolete] venn graph
td0 = x %>% unnest(ds) %>% filter(TimepointB != 0, padj < 0.01, abs(log2fc)>=1) %>%
    mutate(drc = ifelse(log2fc < 0, 'down', 'up')) %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(Treatment = factor(Treatment, levels=stresses)) %>%
    mutate(cond1=str_c(Treatment,Timepoint, sep="_")) %>%
    select(Genotype, cond1, gid, drc)
#
td2 = td1 %>% mutate(cond1=str_c(Treatment,Timepoint,sep="_")) %>%
    select(Genotype,cond1, gids) %>%
    filter(Genotype %in% gts3) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    inner_join(td0, by=c('Genotype','cond1','gid')) %>%
    group_by(Genotype, cond1) %>%
    summarise(gids=list(gid), gidsU=list(gid[drc=='up']), gidsD=list(gid[drc=='down'])) %>%
    ungroup() %>%
    dplyr::rename(all=gids, up=gidsU, down=gidsD) %>%
    gather(drc, gids, -Genotype, -cond1) %>%
    mutate(drc = factor(drc, levels=c('all', drcs)))

tp = td2 %>% select(Genotype,cond1,drc, gids) %>%
    spread(Genotype, gids) %>%
    dplyr::rename(gidsB=B73, gidsM=Mo17, gidsW=W22) %>%
    mutate(vdc = pmap(list(gidsB, gidsM, gidsW), get_venncounts)) %>%
    select(cond1, drc, vdc) %>% unnest(vdc)
tv = tibble(x = c(0, 0.866, -0.866), y = c(1, -0.5, -0.5), lab=rev(gts3)) %>%
    mutate(lab=factor(lab, levels=gts3))
tp1 = tp %>% distinct(cond1) %>% crossing(tv)
pv = ggplot(tp1) +
    geom_circle(aes(x0=x, y0=y, r=1.5, fill=lab), alpha=.3, size=0, color='grey') +
    geom_text(tp, mapping=aes(x=x,y=y,label=Counts), size=2.5) +
    facet_grid(drc ~ cond1, switch='y') +
    scale_fill_manual(values=pal_tron()(5)) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-.4)
fo = sprintf('%s/05.gt.venn.pdf', dirw)
ggsave(pv, filename=fo, width=9, height=7)
#}}}
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



