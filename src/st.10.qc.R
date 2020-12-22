source('functions.R')
dirw = glue('{dird}/10_qc')

yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ TC
ex = 'TC'
th1 = th %>% filter(Experiment==ex) %>%
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
fo = glue('{dirf}/sf02a.rds')
saveRDS(p1a, fo)

p1b = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1b, width=6, height=8)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=3,iter=1000, seed=12,
    var.shape='Treatment',var.col='Genotype',var.lab='has',#var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
fo = glue("{dirw}/22.tsne.{ex}.pdf")
ggsave(fo, p2, width=5, height=5)
fo = file.path(dirf, 'sf02b.rds')
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

#{{{ HY
ex = 'HY'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(grp = sprintf("%s_%dh_%s", Treatment, Timepoint, Genotype)) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(Replicate = 1:n()) %>% ungroup() %>%
    mutate(lab = str_c(SampleID, grp, sep=' ')) %>%
    mutate(clab = ifelse(Replicate==1, cond, ''))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ f1a
p1a = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1a, width=8, height=10)

p1b = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1b, width=8, height=10)

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

gt = 'B73'
gt = 'Mo17'
gt = 'W22'
#{{{ f1a - indi
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


