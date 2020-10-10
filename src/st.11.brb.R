source('functions.R')
dirw = file.path(dird, '11_brb')

#{{{ [obsolete] read1 UMI
#{{{ read umi_cnt and save
ti = crossing(x=1:12, y=1:8) %>%
    mutate(yl = LETTERS[y]) %>%
    mutate(fn = str_c(yl,x,sep='')) %>%
    rbind(tibble(x=13, y=9, yl='I',fn='undetermined')) %>%
    mutate(fi = sprintf("%s/04_umi_cnt/%s/%s.tsv", diri, sid, fn)) %>%
    mutate(res = map(fi, read_tsv, col_names=c('umi','cnt'))) %>%
    select(x, y, yl, fn, res) %>%
    unnest()

fo = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
saveRDS(ti, file=fo)
#}}}

fi = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
ti = readRDS(fi)
n_read = sum(ti$cnt); n_umi = nrow(ti); pct_dedup = percent(n_umi/n_read)
tit = sprintf("%s UMIs / %s total reads = %s", number(n_umi,big.mark=','), number(n_read,big.mark=','), pct_dedup)

#{{{ plot reads + umi
tp = ti %>% group_by(x,y,fn) %>%
    summarise(cnt = sum(cnt), umi=n()) %>% ungroup() %>%
    mutate(lab = str_c(number(cnt), number(umi), sep="\n"))
mid = (min(tp$cnt) + max(tp$cnt)) / 2
tp = tp %>% mutate(color = ifelse(cnt < mid, 'white','black'))
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=cnt), color='black') +
    geom_text(aes(x,y, label=lab, color=cnt<mid), hjust=1, size=2.5, nudge_x=.4) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/01.nseq.%s.pdf", dirw, sid)
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ fastqc
f_qc = sprintf("%s/02_trim/%s_2_fastqc.zip", diri, sid)
qc = qc_read(f_qc, modules=c("Basic Statistics", "Sequence Duplication Levels"))
n_read2 = as.integer(qc$basic_statistics$Value[[4]])
pct_dedup2 = qc$total_deduplicated_percentage
cat(str_c(number(n_read2), pct_dedup2, '\n', sep=', '))
#}}}
#}}}

#{{{ read DGE and store
read_dge <- function(ft)
    read_tsv(ft) %>% rename(gid=1) %>% gather(coord, cnt, -gid)

sids = str_c("batch", c('1','2a','2b','3'), sep="")
td = tibble(sid=sids) %>%
    mutate(fd1=sprintf("%s/cache/21_dge/%s/output.dge.reads.txt", dird, sid)) %>%
    mutate(fd2=sprintf("%s/cache/21_dge/%s/output.dge.umis.txt", dird, sid)) %>%
    mutate(td1 = map(fd1, read_dge)) %>%
    mutate(td2 = map(fd2, read_dge))
td1 = td %>% select(sid, td1) %>% unnest(td1) %>% rename(nr=cnt)
td2 = td %>% select(sid, td2) %>% unnest(td2) %>% rename(nu=cnt)
cnts = td1 %>% inner_join(td2, by=c('sid','coord','gid')) %>%
    select(sid, coord, gid, nr, nu)
tc = cnts %>% unite("SampleID", sid:coord, remove=T, sep="_")

fh = file.path(dird, '01_exp_design/05.brb.local.tsv')
thb = read_tsv(fh) %>%
    select(sid=Tissue,coord=file_prefix,SampleID=Treatment) %>%
    left_join(th, by='SampleID') %>%
    rename(ExpID = SampleID) %>%
    mutate(SampleID = str_c(sid, coord, sep="_")) %>%
    rename(batch = sid) %>%
    select(SampleID,Tissue,Genotype,Treatment,Timepoint,Replicate,everything())
thb %>% count(Experiment)

res = list(th = thb, tc = tc)
fo = file.path(dirw, '01.cnts.rds')
saveRDS(res, fo)

#{{{ ngene umi
l2n = 1:8; names(l2n) = LETTERS[l2n]
min_nr = 2; min_nu = 2 #min_nr = 10; min_nu = 10
tp = res$tc %>% separate(SampleID, c('sid','coord'), sep='_') %>%
    dplyr::rename(SampleID = coord) %>%
    #filter(sid == !!sid) %>%
    rename(n_read = nr, n_umi = nu) %>%
    group_by(sid, SampleID) %>%
    summarise(nr = sum(n_read >= min_nr), nu = sum(n_umi >= min_nu)) %>%
    ungroup() %>%
    mutate(yl = str_sub(SampleID, 1, 1), x = str_sub(SampleID, 2)) %>%
    mutate(y = as.integer(l2n[yl]), x = as.integer(x))

mid = (min(tp$nu) + max(tp$nu)) / 2
tit=sprintf('# genes with >= %d UMIs', min_nu)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nu), color='black') +
    geom_text(aes(x,y, label=number(nu), color=nu<mid), hjust=1, size=2.5, nudge_x=.3) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    facet_wrap(~sid, nrow=2) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/05.ngene.umi.pdf", dirw)
ggsave(fo, p, width=10, height=7)
#
mid = (min(tp$nr) + max(tp$nr)) / 2
tit=sprintf('# genes with >= %d reads', min_nr)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nr), color='black') +
    geom_text(aes(x,y, label=number(nr), color=nr<mid), hjust=1, size=2.5, nudge_x=.3) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    facet_wrap(~sid, nrow=2) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/05.ngene.read.pdf", dirw)
ggsave(fo, p, width=10, height=7)
#}}}
#}}}

#{{{ normalize and store
fi = file.path(dirw, '01.cnts.rds')
cnt = readRDS(fi)
thb = cnt$th

require(DESeq2); require(edgeR)
t_rc = cnt$tc %>% dplyr::select(SampleID, gid, ReadCount = nu) %>%
    filter(SampleID %in% thb$SampleID)
r = readcount_norm(t_rc)
tm_u = r$tm

t_rc = cnt$tc %>% dplyr::select(SampleID, gid, ReadCount = nr) %>%
    filter(SampleID %in% thb$SampleID)
r = readcount_norm(t_rc)
tm_r = r$tm

res = list(th=thb, tm_u=tm_u, tm_r = tm_r)
fo = file.path(dirw, '03.cpm.rds')
saveRDS(res, fo)
#}}}

fi = file.path(dirw, '03.cpm.rds')
cpm = readRDS(fi)

bid = 'b1'; batches = 'batch1'; wd=6; ht=8
bid = 'b2'; batches = c("batch2a","batch2b"); wd=6; ht=10
bid = 'b3'; batches = 'batch3'; wd=6; ht=10
th = cnt$th %>% filter(batch %in% batches) %>%
    mutate(lab = str_c(Genotype, Treatment, Timepoint, sep='_')) %>%
    mutate(grp = str_c(Genotype, Treatment, sep='_'))
tm = cnt$tm_r %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ QC
p1 = plot_hclust(tm,th,pct.exp=.4,cor.opt='pearson',var.col='Genotype',
    expand.x=.2)
fo = sprintf("%s/11.hclust.%s.pdf", dirw, bid)
ggsave(p1, filename = fo, width=wd, height=ht)

p2 = plot_tsne(tm,th,pct.exp=.4,perp=8,iter=1500,
    var.shape='Genotype', var.col='Treatment', var.lab='Timepoint',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
fp = sprintf("%s/12.tsne.%s.pdf", dirw, bid)
ggsave(p2, filename = fp, width=8, height=8)

bid = 'b23'; batches = c("batch2a","batch2b", 'batch3')
th = res$th %>% filter(batch %in% batches) %>%
    mutate(batch = str_replace(batch, '[ab]$', '')) %>%
    mutate(lab = str_c(batch, Genotype, Treatment, Timepoint, sep='_')) %>%
    mutate(Treatment = str_c(batch, Treatment, sep='_'))
tm = res$tm_r %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm,th,pct.exp=.4,cor.opt='pearson',var.col='Genotype',
    expand.x=.2)
fo = sprintf("%s/11.hclust.%s.pdf", dirw, bid)
ggsave(p1, filename = fo, width=8, height=20)

p2 = plot_tsne(tm,th,pct.exp=.4,perp=8,iter=1500,
    var.shape='Genotype', var.col='Treatment', var.lab='Timepoint',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
fp = sprintf("%s/12.tsne.%s.pdf", dirw, bid)
ggsave(p2, filename = fp, width=8, height=8)
#}}}

#{{{ share raw read count table w. Zach
fi = file.path(dirw, '01.cnts.rds')
cnt = readRDS(fi)

to1 = cnt$tc %>% select(-nu) %>% spread(SampleID, nr)
to2 = cnt$tc %>% select(-nr) %>% spread(SampleID, nu)

fo = file.path(dirw, '01.meta.tsv')
write_tsv(cnt$th, fo)
fo1 = file.path(dirw, '91.read.count.tsv.gz')
write_tsv(to1, fo1)
fo2 = file.path(dirw, '91.barcode.count.tsv.gz')
write_tsv(to2, fo2)
#}}}

#{{{ compare w. traditional RNA-Seq
fi = file.path(dirw, '03.cpm.rds')
cpm = readRDS(fi)
#
rn = rnaseq_cpm('rn20a')

th1 = cpm$th %>% select(sid1=SampleID, sid2=ExpID, batch)
th2 = rn$th %>% select(sid2 = SampleID) %>% distinct(sid2)
th = th1 %>% inner_join(th2, by='sid2')

tm1 = cpm$tm_r %>% select(sid1=SampleID, gid, cpm1=CPM)
tm2 = rn$tm %>% select(sid2=SampleID, gid, cpm2=CPM)
tm = th %>% inner_join(tm1, by='sid1') %>%
    inner_join(tm2, by=c('sid2','gid'))

get_cor <- function(ti, min_cpm=1) {
    #{{{
    ng1 = sum(ti$cpm1 >= min_cpm)
    ng2 = sum(ti$cpm2 >= min_cpm)
    tif = ti %>% filter(cpm1 >= min_cpm, cpm2 >= min_cpm)
    ng = nrow(tif)
    pcc = cor(tif$cpm1, tif$cpm2, method='pearson')
    spc = cor(tif$cpm1, tif$cpm2, method='spearman')
    #kdc = cor(tif$cpm1, tif$cpm2, method='kendall')
    tibble(ng1=ng1,ng2=ng2,ng=ng,pcc=pcc,spc=spc)
    #}}}
}

to = tm %>% group_by(sid1,sid2,batch) %>%
    nest() %>% rename(cpm = data) %>%
    mutate(x = map(cpm, get_cor)) %>%
    select(batch, sid1, sid2, x) %>% unnest(x)

to %>% group_by(batch) %>% skim(ng1,ng2,ng,pcc,spc)

fo = file.path(dirw, '61.cor.tsv')
write_tsv(to, fo)

p = ggboxplot(to, x = "batch", y = "ng1",
    color = "batch", palette = pal_npg()(5),
    add = "jitter", shape = "batch", ylab='number genes CPM > 1') +
    otheme(legend.pos='none',ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T)
fo = file.path(dirw, '62.ng1.pdf')
ggsave(p, file=fo, width=5, height=5)

p = ggboxplot(to, x = "batch", y = "ng2",
    color = "batch", palette = pal_npg()(5),
    add = "jitter", shape = "batch", ylab='number genes CPM > 1') +
    otheme(legend.pos='none',ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T)
fo = file.path(dirw, '62.ng2.pdf')
ggsave(p, file=fo, width=5, height=5)

p = ggboxplot(to, x = "batch", y = "ng",
    color = "batch", palette = pal_npg()(5),
    add = "jitter", shape = "batch", ylab='number genes CPM > 1') +
    otheme(legend.pos='none',ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T)
fo = file.path(dirw, '62.ng.pdf')
ggsave(p, file=fo, width=5, height=5)

p = ggboxplot(to, x = "batch", y = "pcc",
    color = "batch", palette = pal_npg()(5),
    add = "jitter", shape = "batch", ylab='number genes CPM > 1') +
    otheme(legend.pos='none',ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T)
fo = file.path(dirw, '62.pcc.pdf')
ggsave(p, file=fo, width=5, height=5)

p = ggboxplot(to, x = "batch", y = "spc",
    color = "batch", palette = pal_npg()(5),
    add = "jitter", shape = "batch", ylab='number genes CPM > 1') +
    otheme(legend.pos='none',ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T)
fo = file.path(dirw, '62.spc.pdf')
ggsave(p, file=fo, width=5, height=5)
#}}}
