source('functions.R')
dirw = glue('{dird}/10_qc')
require(kableExtra)

fi = "~/projects/s3/zhoup-nfo/archive/zm.rn20a/01.rds"
res = readRDS(fi) 
th = res$th %>%
    mutate(Genotype = ifelse(Genotype=='MS71', 'Ms71', Genotype))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
stresses = c("Control",'Cold','Heat')

#{{{ TC - sf07
ex = 'TC'
th1 = th %>% filter(Experiment==ex, Timepoint != 8) %>%
    mutate(Genotype = factor(Genotype, levels=gts3)) %>%
    mutate(Treatment = factor(Treatment, levels=c("Control",'Cold','Heat'))) %>%
    mutate(has = sprintf("h%03d", Timepoint*10)) %>%
    mutate(grp = sprintf("%s_%s_%s", Treatment, has, Genotype)) %>%
    mutate(cond = sprintf("%s_%s", Treatment, has)) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(SampleID, grp, sep=' '))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1a = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1a, width=6, height=8)
fo = glue('{dirf}/sf07a.rds')
saveRDS(p1a, fo)

p1b = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1b, width=6, height=8)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=2,iter=1000, seed=12,
    var.shape='Treatment',var.col='Genotype',var.lab='has',#var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.pdf")
ggsave(fo, p2, width=5, height=5)
fo = file.path(dirf, 'sf07b.rds')
saveRDS(p2, fo)

th2 = th1 %>% filter(Genotype=='B73')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2a = plot_tsne(tm2,th2,pct.exp=.7,perp=3,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.B73.pdf", dirw, ex), p2a, width=6, height=6)

th2 = th1 %>% filter(Genotype=='Mo17')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2b = plot_tsne(tm2,th2,pct.exp=.7,perp=3,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.Mo17.pdf", dirw, ex), p2b, width=6, height=6)

th2 = th1 %>% filter(Genotype=='W22')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2c = plot_tsne(tm2,th2,pct.exp=.8,perp=4,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.W22.pdf", dirw, ex), p2c, width=6, height=6)
#}}}

#{{{ HY - f1c, sf01a-f
#{{{ init
ex = 'HY'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = ifelse(Genotype %in% names(gt_map2), gt_map2[Genotype], Genotype)) %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(grp = sprintf("%s_%dh_%s", Treatment, Timepoint, Genotype)) %>%
    mutate(cond = sprintf("%s %d h", Treatment, Timepoint)) %>%
    mutate(Treatment = factor(Treatment, levels=stresses))
conds = th1 %>% distinct(Treatment, Timepoint, cond) %>%
    arrange(Treatment, Timepoint) %>% pull(cond)
th1 = th1 %>% mutate(cond = factor(cond, levels = conds)) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(Replicate = 1:n()) %>% ungroup() %>%
    mutate(lab = str_c(SampleID, grp, sep=' ')) %>%
    mutate(clab = ifelse(Replicate==1, cond, ''))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))
cols7 = c('gray', brewer.pal(n = 6, name = "Paired")[c(3,4,1,2,5,6)])
cols7 = pal_simpsons()(16)[c(3,4,14,10,7,1,5)]
#}}}

#{{{ hclust
p1a = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1a, width=8, height=10)

p1b = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1b, width=8, height=10)
#}}}

#{{{ f1c, sf01a-f
#{{{ tsne
p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=6,iter=800, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(0:2,3,4,8), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.pdf")
ggsave(fo, p2, width=6, height=6)
fo = glue('{dirf}/f1a.rds')
saveRDS(p2, fo)

th2 = th1 %>% filter(Genotype %in% gts3)
p2 = plot_tsne(tm1,th2,pct.exp=.7,perp=4,iter=1000, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(0:2,3,4,8), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.inbreds.pdf")
ggsave(fo, p2, width=5, height=5)
#}}}

#{{{ tSNE - indi gt 
gt = 'B73'
gt = 'Mo17'
gt = 'W22'
#{{{ f1c - indi
th2 = th1 %>% filter(Genotype %in% gt)
p2 = plot_tsne(tm1,th2,pct.exp=.7,perp=2.5,iter=1000, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='clab',var.ellipse='grp',
    legend.pos='top.right', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(0:2), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.{gt}.pdf")
ggsave(fo, p2, width=4, height=4)
fo = glue("{dirf}/f1a.{gt}.rds")
saveRDS(p2, fo)
#}}}

a = readRDS('f1a.B73.rds') + o_margin(.2,0,0,.2)
b = readRDS('f1a.Mo17.rds') + o_margin(.2,.2,0,0)
c = readRDS('f1a.W22.rds') + o_margin(0,0,.2,.2)
fo = glue("{dirw}/22.tsne.{ex}.BMW.pdf")
ggarrange(a, b, c, nrow=2, ncol=2,
          labels = c("B73","Mo17", "W22"),
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename=fo, width=6, height=6)
#}}}

#{{{ pca - f1c, sf01
#{{{ filter th and tm
th1b = th1 %>% mutate(grp = Genotype) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(grp) %>% mutate(i = 1:n()) %>% ungroup() %>%
    mutate(clab = ifelse(i==1, as.character(Genotype), ''))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ 6 genotypes
p2 = plot_pca(tm1,th1b, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='clab',var.ellipse='grp',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7)
fo = glue("{dirw}/22.pca.{ex}.pdf")
ggsave(fo, p2, width=6, height=6)
#}}}

#{{{ 3 genotypes - sf01a
th2 = th1b %>% filter(Genotype %in% gts3)
p2 = plot_pca(tm1,th2, min.value=1, pct.exp=.6, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='clab',var.ellipse='grp',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    theme(legend.position=c(0,.4))
fo = glue("{dirw}/22.pca.{ex}.inbreds.pdf")
ggsave(fo, p2, width=5, height=5)
#}}}
fo = glue("{dirf}/sf01a.rds")
saveRDS(p2, fo)

#{{{ B73 - f1c
gt = 'B73'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/f1c.rds')
saveRDS(p3, fo)

#{{{ sf01b
gt = 'Mo17'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    facet_wrap(~Genotype)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/sf01b.rds')
saveRDS(p3, fo)

#{{{ sf01c
gt = 'W22'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    facet_wrap(~Genotype)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/sf01c.rds')
saveRDS(p3, fo)

#{{{ sf01d
gt = 'B73xMo17'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.left', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    facet_wrap(~Genotype)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/sf01d.rds')
saveRDS(p3, fo)

#{{{ sf01e
gt = 'W22xB73'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.right', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    facet_wrap(~Genotype)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/sf01e.rds')
saveRDS(p3, fo)

#{{{ sf01f
gt = 'W22xMo17'
th2 = th1 %>% filter(Genotype == gt)
p3 = plot_pca(tm1,th2, min.value=1, pct.exp=.5, pca.center=T, pca.scale=T,
    var.col='cond',var.shape='cond',var.lab='',var.ellipse='',
    legend.pos='bottom.right', legend.dir='v', legend.box='v',legend.title=F,
    shapes=c(1,1,1,15,15,17,17), cols=cols7) +
    facet_wrap(~Genotype)
fo = glue("{dirw}/22.pca.{ex}.{gt}.pdf")
ggsave(fo, p3, width=4, height=4)
#}}}
fo = glue('{dirf}/sf01f.rds')
saveRDS(p3, fo)
#}}}
#}}}
#}}}

#{{{ NM
ex = 'NM'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts25)) %>%
    mutate(Treatment = factor(Treatment, levels=c("Control",'Cold'))) %>%
    mutate(grp = sprintf("%s_%dh_%s", Treatment, Timepoint, Genotype)) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    mutate(cond = factor(cond, levels=c("Control_1h",'Control_25h','Cold_1h','Cold_25h'))) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(SampleID, cond, Genotype, sep=' '))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.25)
fo = glue("{dirw}/21.hclust.{ex}.p.pdf")
ggsave(fo, p1, width=6, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.25)
fo = glue("{dirw}/21.hclust.{ex}.s.pdf")
ggsave(fo, p1, width=6, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=4,iter=1200, seed=12,
    var.shape='cond',var.col='cond',var.lab='Genotype',#var.ellipse='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:1,15:16), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.pdf")
ggsave(fo, p2, width=6, height=6)
#}}}

#{{{ sd1
exps = c('TC'="2 (time course)",'HY'='1 (inbreds + hybrids)','NM'='3 (genotype panel)')
conds = c("Control",'Cold','Heat')
gts = unique(c(gts6,gts25))
tp = res$bamstat %>%
    inner_join(th, by=c('sid'='SampleID')) %>%
    mutate(Experiment=factor(exps[Experiment])) %>%
    mutate(Treatment=factor(Treatment, levels=conds)) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    arrange(Experiment, Treatment, Timepoint, Genotype) %>%
    group_by(Experiment) %>% mutate(i = 1:n()) %>%
    ungroup() %>%
    mutate(Experiment = as.character(Experiment)) %>%
    mutate(Experiment=ifelse(i==1, Experiment, '')) %>%
    mutate(Timepoint = as.character(Timepoint)) %>%
    mutate(rate=percent(pair_map_hq/pair_map, accuracy=.1)) %>%
    mutate(pair=number(pair,big.mark=',')) %>%
    mutate(pair_map=number(pair_map,big.mark=',')) %>%
    mutate(pair_map_hq=number(pair_map_hq,big.mark=',')) %>%
    mutate(rate = glue("{pair_map_hq} ({rate})")) %>%
    select(Experiment,Treatment,Timepoint,Genotype,
        `Total Reads`=pair, `Mapped Reads`=pair_map,
        `Uniquely Mapped`=rate) %>%
    print(width=Inf)

fo = glue("{dird}/71_share/sd1.tsv")
write_tsv(tp, fo)

x = tp %>%
    kbl(format='latex', escape=T, longtable=T, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    #collapse_rows(1, latex_hline='major', valign='middle', longtable_clean_cut=F) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
#fo = file.path(dirf, 'st1.rds')
#saveRDS(x, file=fo)
#}}}

#{{{ zach's st6
fi = glue("{dirw}/st6.xlsx")
ti = read_xlsx(fi)

x = ti %>% select(-cid) %>%
    rename(GeneID=gid, `Response pattern`=`Heat Responsive`) %>%
    kbl(format='latex', escape=T, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    #collapse_rows(1, latex_hline='major', valign='middle', longtable_clean_cut=F) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirf, 'st6.rds')
saveRDS(x, file=fo)
#}}}

#{{{ share read count and CPM table
diro = glue("{dird}/71_share")
exps = c('TC'="2 (time course)",'HY'='1 (inbreds + hybrids)','NM'='3 (genotype panel)')
conds = c("Control",'Cold','Heat')
gts = unique(c(gts6,gts25))

to1 = th %>% 
    mutate(Experiment=factor(exps[Experiment])) %>%
    mutate(Treatment=factor(Treatment, levels=conds)) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    arrange(Experiment, Treatment, Timepoint, Genotype, SampleID) %>%
    #group_by(Experiment) %>% mutate(i = 1:n()) %>%
    #ungroup() %>%
    #mutate(Experiment = as.character(Experiment)) %>%
    #mutate(Experiment=ifelse(i==1, Experiment, '')) %>%
    mutate(Timepoint = as.character(Timepoint)) %>%
    group_by(Experiment,Treatment,Timepoint,Genotype) %>%
    mutate(Replicate=1:n()) %>% ungroup() %>%
    select(Experiment,Treatment,Timepoint,Genotype, Replicate, SampleID) %>%
    print(n=20)

tor = res$tm %>% filter(SampleID %in% to1$SampleID) %>%
    select(SampleID, gid, ReadCount) %>%
    spread(SampleID, ReadCount)
toc = res$tm %>%# filter(SampleID %in% to1$SampleID) %>%
    select(SampleID, gid, CPM) %>%
    spread(SampleID, CPM)
tof = res$tm %>%# filter(SampleID %in% to1$SampleID) %>%
    select(SampleID, gid, FPKM) %>%
    spread(SampleID, FPKM)
tot = res$tm %>%# filter(SampleID %in% to1$SampleID) %>%
    select(SampleID, gid, TPM) %>%
    spread(SampleID, TPM)

fo = glue("{diro}/01.meta.tsv")
write_tsv(to1, fo)
fo = glue("{diro}/05.readcount.tsv.gz")
write_tsv(tor, fo)
fo = glue("{diro}/05.cpm.tsv.gz")
write_tsv(toc, fo)
fo = glue("{diro}/05.fpkm.tsv.gz")
write_tsv(tof, fo)
fo = glue("{diro}/05.tpm.tsv.gz")
write_tsv(tot, fo)
#}}}

