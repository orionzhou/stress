#{{{ load & read
require(devtools)
load_all('~/git/rmaize')
require(progress)
require(ape)
require(ggtree)
require(ggforce)
require(Rtsne)
require(ggpubr)
require(lubridate)
options(dplyr.summarise.inform = F)
dirg = '~/data/genome'
dirp = "~/projects/stress"
dird = file.path(dirp, 'data')
dirr = '~/projects/stress/nf/raw'
dirf = file.path(dird, '95_figures', 'plots')
gcfg = read_genome_conf()
fh = file.path(dird, 'samples.xlsx')
th = read_xlsx(fh, sheet='merged') %>%
    mutate(Time = sprintf("%02d:%02d", hour(Time), minute(Time))) %>%
    mutate(Tissue = 'leaf') %>% rename(Replicate=Rep) %>%
    select(SampleID,Tissue,Genotype,Treatment,Timepoint,Experiment,Replicate)
th2 = th %>%
    mutate(txt=str_c(Experiment, Treatment, Timepoint, sep='_')) %>%
    select(SampleID, Tissue, Genotype, Treatment, Experiment, txt, Replicate)
gts3 = c("B73",'Mo17','W22')
gts6 = c("B73",'Mo17','W22','B73xMo17','W22xB73','W22xMo17')
gts25 = c("B73", "B97", "CML322", "CML333", "CML52", "CML69", "DK105",
    "EP1", "F7", "Il14H", "Ki11", "Ki3", "M162W", "M37W",
    "Mo17", "Mo18W", "MS71", "NC350", "NC358", "Oh43", "Oh7B",
    "P39", "PH207", "Tx303", "W22")
gt_map = list('B'='B73','M'='Mo17','W'='W22','BMW'=gts3)
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
cols100v = viridis_pal(direction=-1,option='magma')(100)
bats = c('cold_up', 'heat_up', 'cold_down', 'heat_down')
#}}}
colbright <- function(col) {x = col2rgb(col); as.integer(.2126*x[1] + .7152*x[2] + .0722*x[3])}
cols36 = c(pal_ucscgb()(18)[8], pal_igv()(18), pal_ucscgb()(18)[c(1:7,9:18)])
brights36 = tibble(col=cols36) %>% mutate(bright=map_int(col,colbright)) %>%
    mutate(b = ifelse(bright<128, 'white','black')) %>% pull(b)

get_ds <- function(cond, condB, dds, gids) {
    #{{{
    res1 = results(dds, contrast = c("cond",cond,condB), pAdjustMethod='fdr')
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}
call_deg_1 <- function(th, tm, base.cond = 'wt') {
    #{{{ th with 'cond1' and 'cond2' cols
    th1 = th %>% mutate(cond = str_c(cond1, cond2, sep='.'))
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    ct = th1 %>% distinct(cond, cond1, cond2) %>% filter(cond2 != base.cond) %>%
        mutate(condB = str_c(cond1, base.cond, sep="."))
    #{{{ prepare data
    vh = th1 %>% mutate(cond2 = factor(cond2)) %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = ct %>% mutate(ds = map2(cond, condB, get_ds, dds = dds, gids = gids))
    res %>% select(cond1, cond, condB, ds)
    #}}}
}
call_deg <- function(th, tm, comps) {
    #{{{ th with 'cond' column and comps with 'cond1' and 'cond2' columns
    th1 = th
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    #{{{ prepare data
    vh = th1 %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = comps %>% mutate(ds = map2(cond1, cond2, get_ds, dds=dds, gids=gids))
    res
    #}}}
}
#### WGCNA
order_id_by_hclust <- function(ti, cor.opt='pearson', hc.opt = 'ward.D') {
    #{{{
    ids = ti$id
    e = ti %>% select(-id) %>% t()
    edist = as.dist(1-cor(e, method = cor.opt))
    ehc = hclust(edist, method = hc.opt)
    ids[ehc$order]
    #}}}
}
scale_by_first <- function(x) {
    #{{{
    x1 = x - x[1]
    rge0 = max(x) - min(x)
    rge1 = max(abs(x))
    x1# / rge1
    #}}}
}
pick_soft_power <- function(ti) {
    #{{{
    datExpr = t(as.matrix(ti[,-1]))
    colnames(datExpr) = ti$gid
    #
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    powers = 1:20
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    #}}}
}
plot_soft_power <- function(sft, fo, wd=9, ht=5) {
    #{{{
    pdf(fo, width = 9, height = 5)
    #sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
        main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    #}}}
}
run_wgcna <- function(ti, softPower=8, type='signed hybrid', corFnc='cor',
                      TOM=T, TOMType='signed', hclust.opt='ward.D') {
    #{{{
    require(WGCNA)
    #enableWGCNAThreads()
    #
    datExpr = t(as.matrix(ti[,-1]))
    colnames(datExpr) = ti$gid
    datExpr[is.na(datExpr)] = 0
    #
    adj = adjacency(datExpr, power=softPower, type=type, corFnc=corFnc)
    dim(adj)
    if (TOM) {
        TOM = TOMsimilarity(adj, TOMType =TOMType)
        diss = 1-TOM
    } else {
        diss = 1-adj
    }
    tree = hclust(as.dist(diss), method = hclust.opt)
    list(datExpr=datExpr, diss=diss, tree=tree)
    #}}}
}
make_raw_modules <- function(datExpr, dissTOM, tree,
                             minModuleSize=20, deepSplit=2, minGap=0, pre='rc') {
    #{{{
    # Module identification using dynamic tree cut:
    if(minGap == 0) minGap = NULL
    dynamicMods = cutreeDynamic(dendro = tree, distM = dissTOM,
                  deepSplit = deepSplit, pamRespectsDendro = F,
                  minClusterSize = minModuleSize, minGap = minGap)
    dynamicColors = labels2colors(dynamicMods)
    clus1 = dynamicMods; cols1 = dynamicColors
    # Calculate eigengenes
    MEList = moduleEigengenes(datExpr, colors=dynamicColors, excludeGrey=T, subHubs=T, returnValidOnly = T)
    MEs = MEList$eigengenes
    #
    # Construct numerical labels corresponding to the colors
    tu = tibble(gid=colnames(datExpr), clu=clus1, col=cols1) %>%
        group_by(clu, col) %>%
        summarise(n=n(), gids=list(gid)) %>% ungroup()
    me = MEs %>% as_tibble() %>% mutate(cond=rownames(MEs)) %>%
        gather(col, val, -cond) %>% spread(cond, val) %>%
        mutate(col=str_replace(col, "^ME", "")) %>%
        group_by(col) %>% nest() %>% dplyr::rename(me = data) %>% ungroup() %>%
        inner_join(tu, by=c('col')) %>%
        arrange(desc(n)) %>%
        mutate(cid = str_c(pre, str_pad(1:n(),2,pad='0'))) %>%
        #select(ctag, col, clu, me, n, gids)
        select(cid, n, me, gids)
    # make clu
    clu = me %>% select(cid, gids) %>% unnest(gids) %>% dplyr::rename(gid=gids)
    clu = tibble(gid=colnames(datExpr)) %>% inner_join(clu,by='gid')
    #
    list(me=me, clu=clu)
    #}}}
}
merge_modules <- function(rc, datExpr, pre='mc', cutHeight=.2, mms=c()) {
    #{{{
    tMEs = rc$me %>% select(cid, me) %>% unnest(me) %>%
        gather(cond, val, -cid) %>% spread(cid, val)
    MEs = tMEs %>% select(-cond) %>% as.data.frame()
    rownames(MEs) = tMEs$cond
    # Cluster module eigengenes
    #MEDiss = 1-cor(MEs) #, use="pairwise.complete.obs")
    #METree = hclust(as.dist(MEDiss), method = 'ward.D')
    #
    mg = mergeCloseModules(datExpr, rc$clu$cid,# MEs=MEs,
                           cutHeight=cutHeight, verbose=0)
    MEs2 = mg$newMEs; mergedColors = mg$colors;
    #
    # Construct numerical labels corresponding to the colors
    tu = tibble(gid=colnames(datExpr), col=mergedColors) %>%
        group_by(col) %>% summarise(n=n(), gids=list(gid)) %>% ungroup() %>%
        arrange(desc(n)) %>%
        mutate(cid = str_c(pre, str_pad(1:n(),2,pad='0'))) %>%
        select(cid, n, gids)
    tu = tu %>% mutate(cid = ifelse(cid %in% names(mms), mms[cid], cid)) %>%
        unnest(gids) %>%
        group_by(cid) %>%
        summarise(n = n(), gids = list(gids)) %>% ungroup()
    cids = tu %>% select(-n) %>% unnest(gids) %>% rename(gid=gids) %>%
        arrange(match(gid, colnames(datExpr))) %>% pull(cid)
    r = moduleEigengenes(datExpr, cids)
    me = r$eigengenes %>% as_tibble() %>%
        mutate(cond=rownames(r$eigengenes)) %>%
        gather(cid, val, -cond) %>% spread(cond, val) %>%
        mutate(cid=str_replace(cid, "^ME", "")) %>%
        group_by(cid) %>% nest() %>% dplyr::rename(me = data) %>% ungroup()
    ae = r$averageExpr %>% as_tibble() %>%
        mutate(cond=rownames(r$averageExpr)) %>%
        gather(cid, val, -cond) %>% spread(cond, val) %>%
        mutate(cid=str_replace(cid, "^AE", "")) %>%
        group_by(cid) %>% nest() %>% dplyr::rename(ae = data) %>% ungroup()
    tu = tu %>% inner_join(me, by='cid') %>%
        inner_join(ae, by=c('cid')) %>%
        arrange(cid)
    tul = tu %>% select(cid, gids) %>% unnest(gids) %>%
        dplyr::rename(gid=gids) %>%
        arrange(match(gid, colnames(datExpr))) %>% select(gid,cid)
    #
    list(tu=tu, tul=tul)
    #}}}
}
cut_merge_modules <- function(datExpr, dissTOM, tree,
                              minModuleSize=20, deepSplit=2,
                              MEDissThres=.25, minGap=0) {
    #{{{
    #{{{ Module identification using dynamic tree cut:
    if(minGap == 0) minGap = NULL
    dynamicMods = cutreeDynamic(dendro = tree, distM = dissTOM,
                  deepSplit = deepSplit, pamRespectsDendro = FALSE,
                  minClusterSize = minModuleSize, minGap = minGap)
    dynamicColors = labels2colors(dynamicMods)
    clus1 = dynamicMods; cols1 = dynamicColors
    # Calculate eigengenes
    MEList = moduleEigengenes(datExpr, colors=dynamicColors, excludeGrey=T, subHubs=T, returnValidOnly = T)
    MEs = MEList$eigengenes
    MEs1 = MEs
    #}}}
    #{{{ Cluster module eigengenes
    MEDiss = 1-cor(MEs) #, use="pairwise.complete.obs")
    METree = hclust(as.dist(MEDiss), method = 'ward.D')
    #
    mg = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 0)
    mergedColors = mg$colors;
    mergedMEs = mg$newMEs;
    MEs2 = mergedMEs
    #}}}
    #{{{ Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(mergedColors, colorOrder)-1
    clus2 = moduleLabels; cols2 = mergedColors
    #
    clu = tibble(gid=colnames(datExpr), clu1=clus1, col1=cols1, clu2=clus2, col2=cols2)
    tu1 = clu %>% group_by(clu1, col1) %>%
        summarise(n=n(), gids=list(gid)) %>%
        ungroup() %>% mutate(ctag='raw') %>% rename(clu=clu1,col=col1)
    tu2 = clu %>% group_by(clu2, col2) %>%
        summarise(n=n(), gids=list(gid)) %>%
        ungroup() %>% mutate(ctag='merged') %>% rename(clu=clu2,col=col2)
    tu = tu1 %>% bind_rows(tu2)
    #
    tp1 = MEs1 %>% as_tibble() %>% mutate(cond=rownames(MEs1)) %>%
        gather(col, val, -cond) %>% spread(cond, val) %>% mutate(ctag = 'raw')
    tp2 = MEs2 %>% as_tibble() %>% mutate(cond=rownames(MEs2)) %>%
        gather(col, val, -cond) %>% spread(cond, val) %>% mutate(ctag = 'merged')
    me0 = tp1 %>% bind_rows(tp2) %>%
        mutate(col=str_replace(col, "^ME", "")) %>%
        group_by(ctag, col) %>% nest() %>% rename(me = data) %>% ungroup() %>%
        inner_join(tu, by=c('ctag','col')) %>%
        select(ctag, col, clu, me, n, gids)
    me = me0 %>% filter(ctag=='merged') %>%
        arrange(desc(n)) %>%
        mutate(cid = str_c('c', str_pad(1:n(),2,pad='0'))) %>%
        select(cid, n, me, gids)
    #}}}
    #clu = tu %>% select(-clu) %>% rename(clu=nclu) %>% select(ctag,clu,n,gids)
    list(me=me, me0=me0, clu=clu)
    #}}}
}
plot_me <- function(me, ncol=3, tit='', strip.compact=T) {
    #{{{
    tp = me %>%
        mutate(pnl = glue("{cid} [{number(n,accuracy=1)}]")) %>%
        select(pnl, n, me) %>% unnest(me) %>%
        gather(x, val, -pnl, -n) %>%
        mutate(hr=as.double(str_replace(x, 'h', ''))/10)
    tpx = tp %>% distinct(x, hr) %>% arrange(hr)
    tpp = tp %>% distinct(pnl, n) %>% arrange(desc(n))
    tp = tp %>% mutate(pnl=factor(pnl,levels=tpp$pnl))
    times = c(0,2,4,8,25)
    ggplot(tp, aes(x=x, y=val, group=1)) +
        geom_point(size=1.5) +
        geom_line() +
        #scale_x_continuous(name="Hours", breaks=times, expand=expansion(mult=c(.05,.05))) +
        scale_x_discrete(name="Hours", breaks=tpx$x, labels=tpx$hr, expand=expansion(mult=c(.05,.05))) +
        scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
        facet_wrap(pnl~., ncol=ncol, scale='free_y') +
        ggtitle(tit) +
        otheme(panel.spacing=.1, strip.compact=strip.compact,
               xtitle=T, xtext=T, xtick=T, margin=c(.2,.2,.2,.2)) +
        theme(plot.title=element_text(hjust=.5, size=10))
    #}}}
}
extract_avg_expr <- function(ti, tcr) {
    #{{{
    ti2 = ti %>% inner_join(tcr, by=c('cond','gid')) %>%
        gather(sid, val, -cond, -cid, -gid) %>%
        group_by(gid, cid, cond) %>%
        mutate(val = val/mean(val)) %>% ungroup() %>%
        select(cond, cid, gid, sid, val) %>%
        mutate(val = ifelse(is.na(val), 0, val))
    #
    ti3 = ti2 %>% group_by(cond, cid, sid) %>%
        summarise(avg=mean(val), std=sd(val),
                  y5=quantile(val,.05), y25=quantile(val,.25),
                  y50=quantile(val,.5), y75=quantile(val,.75),
                  y95=quantile(val,.95)) %>% ungroup()
    tis = ti3 %>% filter(sid=='h000') %>% select(cond,cid,off=y50)
    ti4 = ti3 %>% inner_join(tis, by=c('cond','cid')) %>%
        mutate(avg=avg-off,y5=y5-off,y25=y25-off,y50=y50-off,
               y75=y75-off,y95=y95-off) %>%
        select(-off) %>% rename(x= sid) %>%
        mutate(xn=as.double(str_replace(x, 'h', ''))/10) %>%
        select(cid,cond,x,xn,everything())
    ti4
    #}}}
}
plot_avg_expr <- function(tp, ncol=3, tit='', strip.compact=T) {
    #{{{
    tpx = tp %>% distinct(x, xn) %>% arrange(xn)
    #tpp = tp %>% distinct(pnl, n) %>% arrange(desc(n))
    #tp = tp %>% mutate(pnl=factor(pnl,levels=tpp$pnl))
    times = c(0,2,4,8,25)
    cols3 = pal_npg()(6)[c(2,1,4)]
    tpl = tp %>% arrange(cid, desc(y75)) %>%
        group_by(cid) %>% slice(1) %>% ungroup() %>%
        mutate(x='h250') %>% select(cid,pnl,txt,x,y=y75)
    p = ggplot(tp, aes(x=x)) +
        geom_ribbon(aes(ymin=y25,ymax=y75,fill=cond,group=cond), alpha=.2) +
        geom_line(aes(y=y50, color=cond, group=cond)) +
        geom_text(data=tpl, aes(x=x,y=y, label=txt), size=2, hjust=1, vjust=1) +
        #scale_x_continuous(name="Hours", breaks=times, expand=expansion(mult=c(.05,.05))) +
        scale_x_discrete(name="Hours", breaks=tpx$x, labels=tpx$xn, expand=expansion(mult=c(.05,.05))) +
        scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
        facet_wrap(pnl~., ncol=ncol, scale='free_y') +
        scale_color_manual(values=cols3) +
        scale_fill_manual(values=cols3) +
        otheme(legend.pos='bottom.right', legend.dir='v',
               panel.spacing=.1, strip.compact=strip.compact,
               xtitle=T, xtext=T, xtick=T, margin=c(.2,.2,.2,.2))
        #guides(color=F, fill=F)
    if (tit != '') p = p + ggtitle(tit) +
        theme(plot.title=element_text(hjust=.5, size=10))
    p
    #}}}
}
plot_cluster <- function(clu, tcr, ncol=3, tit='', strip.compact=T) {
    #{{{
    ti = tcr %>%# filter(cond %in% c(!!cond, 'control')) %>%
        inner_join(clu, by='gid') %>%
        gather(sid, val, -cond, -cid, -gid) %>%
        group_by(gid, cid, cond) %>%
        mutate(val = val/mean(val)) %>% ungroup() %>%
        select(cond, cid, gid, sid, val) %>%
        mutate(val = ifelse(is.na(val), 0, val))
    #
    tp = ti %>% group_by(cond, cid, sid) %>%
        summarise(n=n(), avg=mean(val), std=sd(val), y25=quantile(val,.25),
                  y50=quantile(val,.5), y75=quantile(val,.75)) %>% ungroup()
    tps = tp %>% filter(sid=='h000') %>% select(cond,cid,off=y50)
    tp = tp %>% inner_join(tps, by=c('cond','cid')) %>%
        mutate(avg=avg-off,y25=y25-off,y50=y50-off,y75=y75-off) %>%
        dplyr::rename(x=sid) %>%
        mutate(pnl = glue("{cid} [{number(n,accuracy=1)}]")) %>%
        mutate(hr=as.double(str_replace(x, 'h', ''))/10)
    tpx = tp %>% distinct(x, hr) %>% arrange(hr)
    tpp = tp %>% distinct(pnl, n) %>% arrange(desc(n))
    tp = tp %>% mutate(pnl=factor(pnl,levels=tpp$pnl))
    times = c(0,2,4,8,25)
    ggplot(tp, aes(x=x)) +
        geom_ribbon(aes(ymin=y25,ymax=y75,fill=cond,group=cond), alpha=.5) +
        geom_line(aes(y=y50, color=cond, group=cond)) +
        #scale_x_continuous(name="Hours", breaks=times, expand=expansion(mult=c(.05,.05))) +
        scale_x_discrete(name="Hours", breaks=tpx$x, labels=tpx$hr, expand=expansion(mult=c(.05,.05))) +
        scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
        facet_wrap(pnl~., ncol=ncol, scale='free_y') +
        scale_color_npg() +
        ggtitle(tit) +
        otheme(panel.spacing=.1, strip.compact=strip.compact,
               xtitle=T, xtext=T, xtick=T, margin=c(.2,.2,.2,.2)) +
        theme(plot.title=element_text(hjust=.5, size=10)) +
        guides(color=F, fill=F)
    #}}}
}
plot_me0 <- function(me, cond) {
    #{{{
    tp = me %>%
        mutate(time=as.double(str_replace(cond, !!cond, ''))/10) %>%
        mutate(pan = sprintf("%s [%s]", clu, number(n)))
    #times = th %>% distinct(Timepoint) %>% arrange(Timepoint) %>% pull(Timepoint)
    times = c(0,2,4,8,25)
    tpa = tp %>% filter(ctag=='raw')
    tpb = tp %>% filter(ctag=='merged')
    xtit = sprintf("Hour after %s Stress", cond)
    pe1 = ggplot(tpa, aes(x=time, y=val)) +
        geom_point(size=1.5) +
        geom_line() +
        scale_x_continuous(name=xtit, breaks=times, expand=expand_scale(mult=c(.05,.05))) +
        scale_y_continuous(expand=expand_scale(mult=c(.05,.05))) +
        facet_wrap(~pan, ncol=1) +
        otheme(xtitle=T, xtext=T, xtick=T, margin=c(.6,.2,.1,.1))
    pe2 = ggplot(tpb, aes(x=time, y=val)) +
        geom_point(size=1.5) +
        geom_line() +
        scale_x_continuous(name=xtit, breaks=times, expand=expand_scale(mult=c(.05,.05))) +
        scale_y_continuous(expand=expand_scale(mult=c(.05,.05))) +
        facet_wrap(~pan, ncol=1) +
        otheme(xtitle=T, xtext=T, xtick=T, margin=c(.6,.2,.1,.1))
    ggarrange(pe1, pe2,
        nrow = 1, ncol = 2, labels = c("Raw",'Merged'), widths = c(4,4),
        hjust=-.1, vjust=1.5)
    #}}}
}
run_softPower <- function(cid, cond, opt_deg, opt_clu, deg, dirw) {
    #{{{
    cat(cid, "\n")
    gts1 = 'B73'; gts3 = c("B73",'Mo17','W22')
    gts_deg = if(opt_deg == 'BMW') gts3 else gts1
    gts_clu = if(opt_clu == 'BMW') gts3 else gts1
    #
    gids = deg$td2 %>% filter(Genotype %in% gts_deg, str_detect(cond1, cond)) %>%
        select(gids) %>% unnest(gids) %>% distinct(gids) %>% pull(gids)
    #
    tt = tm %>% select(gid, SampleID, CPM) %>%
        inner_join(th[,c('SampleID','ExpID','Genotype','Treatment','Timepoint')], by='SampleID') %>%
        filter(ExpID=='TC', Treatment==!!cond | Timepoint==0, gid %in% gids) %>%
        mutate(cond = sprintf("%s%03d", cond, Timepoint*10)) %>%
        select(-ExpID, -SampleID, -Treatment, -Timepoint) %>%
        mutate(CPM = asinh(CPM)) %>%
        spread(cond, CPM)
    tiw = tt %>% filter(Genotype %in% gts_clu) %>% mutate(gid=str_c(gid, Genotype,sep="_")) %>% select(-Genotype)
    #
    fo = sprintf("%s/11.%s.pdf", dirw, cid)
    plot_soft_power(tiw, fo)
    TRUE
    #}}}
}
run_wgcna_pipe <- function(cid, cond, drc, opt_deg, opt_clu, optQ,
                           softPower, deepSplit, MEDissThres, minGap,
                           gt_map, tc, deg, dirw) {
    #{{{
    cat(cid, "\n")
    gts_deg = gt_map[[opt_deg]]; gts_clu = gt_map[[opt_clu]]
    gids = deg %>%
        filter(Genotype %in% gts_deg, Treatment==cond, drc == !!drc) %>%
        select(gids) %>% unnest(gids) %>% distinct(gids) %>% pull(gids)
    #
    ti = tc[[optQ]] %>%
        filter(Genotype %in% gts_clu, Treatment == cond, gid %in% gids) %>%
        mutate(gid=str_c(gid, Genotype,sep="_")) %>%
        select(-Genotype, -Treatment)
    #
    ng = length(gids); np = nrow(ti)
    tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
        number(ng), opt_deg, cond, drc, number(np), opt_clu)
    labx = ifelse(minGap==0, 'auto', minGap)
    lab = sprintf("deepSplit=%g minGap=%s", deepSplit, labx)
    cat(tit, '\n')
    cat(lab, '\n')
    x = run_wgcna(ti, optQ=optQ, softPower=softPower, deepSplit=deepSplit,
                  MEDissThres=MEDissThres, minGap=minGap)
    nm1 = nrow(x$clu %>% distinct(clu1))
    nm2 = nrow(x$clu %>% distinct(clu2))
    msg = sprintf("%d raw modules; %d merged modules", nm1, nm2)
    cat(msg, '\n')
    stats = list(ng=ng,np=np,nm1=nm1,nm2=nm2)
    c(x, stats)
    #}}}
}
get_tss_tts <- function(genome) {
    #{{{
    fi = file.path(dirg, genome, '50_annotation', '10.tsv')
    ti = read_tsv(fi) %>%
        filter(ttype == 'mRNA', etype == 'exon') %>%
        mutate(start=start-1) %>%
        group_by(gid, tid) %>%
        summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
        ungroup() %>% mutate(size = end-start) %>%
        arrange(gid, desc(size)) %>%
        group_by(gid) %>% slice(1) %>% ungroup() %>%
        mutate(tss = ifelse(srd=='-', end, start)) %>%
        mutate(tts = ifelse(srd=='-', start, end)) %>%
        filter(!str_detect(gid, "^(Zeam|Zema)"))
    ti %>% select(gid, chrom, tss, tts, srd)
    #}}}
}
get_tss <- function(genome) {
    #{{{
    fi = file.path(dirg, genome, '50_annotation', '10.tsv')
    ti = read_tsv(fi) %>% filter(etype == 'exon') %>% mutate(start=start-1) %>%
        filter(ttype=='mRNA') %>%
        group_by(gid, tid) %>%
        summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
        ungroup() %>%
        arrange(chrom, start, end) %>%
        mutate(tss = ifelse(srd=='-', end, start)) %>%
        mutate(tss0 = ifelse(srd == '-', -tss, tss)) %>%
        arrange(gid, desc(tss0)) %>%
        group_by(gid) %>% slice(1) %>% ungroup() %>% select(-tss0, -tid) %>%
        filter(!str_detect(gid, "^(Zeam|Zema)"))
    ti %>% select(gid, chrom, pos=tss, srd)
    #}}}
}
eval_pred <- function(ti, seed=1, downsample=F) { # status, pred, prob
  #{{{
  if(downsample) ti = downsample(ti, seed=seed, colname='status')
  metrics6 <- metric_set(sens,spec,precision,accuracy, f_meas, roc_auc, pr_auc)
  metrics6(ti, truth=status, estimate=pred, prob) %>%
    select(metric=.metric, estimate=.estimate)
  #}}}
}


