source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
a = readRDS('f1a.rds') + o_margin(.5,0,0,0)
b = readRDS('f1b.rds') + o_margin(.1,.2,0,.4) + theme(axis.title.x=element_text(size=7))
c = readRDS('f1c.rds') + o_margin(.1,.2,.5,.5)
d = readRDS('f1d.rds')
pab = ggarrange(a, b, nrow=1, ncol=2, labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2))
pcd = ggarrange(c,d, nrow=1, ncol=2, labels = LETTERS[3:4],
          widths=c(2,2), heights=c(2,2))
pabc = ggarrange(a, b, c, nrow=3, ncol=1, labels = LETTERS[1:3],
          widths=c(2,2), heights=c(2.3,1.7,4))
ggarrange(pabc, d, nrow=1, ncol=2, labels=c('','D'),
          widths=c(2,2), heights=c(1,2.5)) %>%
    ggexport(filename="f1.pdf", width=9, height=8)
#}}}

#{{{ Fig 2
a = readRDS('f2a.rds')
b = readRDS('f2b.rds')
c = readRDS('f2c.rds')
d = readRDS('f2d.rds')
e = readRDS('f2e.rds')
ggarrange(a, b, c, d, e, nrow=2, ncol=3, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="f2.pdf", width=8, height=6)
#}}}

#{{{ Fig S1
    #theme(plot.margin = margin(1.2,0,.2,.2, "lines")) +
a = readRDS('sf01a.rds') + theme(legend.position=c(0,.4))
b = readRDS('sf01b.rds') + theme(legend.position='none')
c = readRDS('sf01c.rds') + theme(legend.position='none')
d = readRDS('sf01d.rds') + theme(legend.position='none')
e = readRDS('sf01e.rds') + theme(legend.position='none')
f = readRDS('sf01f.rds') + theme(legend.position='none')
pbc = ggarrange(b,c, nrow=2, ncol=1, labels=LETTERS[2:3])
pabc = ggarrange(a, pbc, nrow=1, ncol=2, labels=c('A',''), widths=c(2,1))
pdef = ggarrange(d,e,f, nrow=1, ncol=3, labels=LETTERS[4:6])
ggarrange(pabc, pdef, nrow=2, ncol=1, labels=c("A",''),
          widths=c(2,2), heights=c(2,1)) %>%
ggexport(filename="sf01.pdf", width=8, height=8)
#}}}

# Fig S2 and S3 are done

#{{{ Fig S4
a = readRDS('sf02a.rds') +
    theme(plot.margin = margin(1.2,0,.2,.2, "lines")) +
    theme(legend.justification=c(.5,-.5))
b = readRDS('sf02b.rds') +
    theme(plot.margin = margin(1.2,0,.2,0, "lines")) +
    theme(legend.justification=c(.5,-.5))
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(1.1,2), heights=c(2,2)) %>%
ggexport(filename="sf02.pdf", width=10.5, height=8)
#}}}

#{{{ Fig S5
a = readRDS('sf05a.rds')
b = readRDS('sf05b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf05.pdf", width=11, height=8)
#}}}

#{{{ Fig S6
a = readRDS('sf06a.rds')
b = readRDS('sf06b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf06.pdf", width=10, height=8)
#}}}

#{{{ Fig S7
a = readRDS('sf05a.rds') +
    scale_x_continuous(expand = expansion(mult=c(1e-2,.4))) +
    theme(plot.margin = margin(.2,0,.2,.2, "lines"))
b = readRDS('sf05b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="sf05.pdf", width=9, height=6)
#}}}

#{{{ Fig 3
a = image_read_pdf(glue("{dirf}/regions.pdf"))
b = readRDS(glue("{dird}/25_dreme/11.n.mtf.rds"))
c = readRDS(glue("{dird}/25_dreme/12.top.mtf.rds")) +
    theme(legend.position = c(.05,0), legend.justification = c(0,0))
ab = ggarrange(as.grob(a), b, nrow=2, ncol=1, heights=c(2,2), labels=LETTERS[1:2])
ggarrange(ab, c, nrow=1, ncol=2,
          labels = c('', 'C'),
          widths=c(2,2)) %>%
ggexport(filename="f3.pdf", width=8, height=7)
#}}}

#{{{ Fig 4
a = image_read(glue("{dirf}/scheme.png"))
b = readRDS(glue("{dird}/41_ml/31.acc.b.rds"))
c = readRDS(glue("{dird}/41_ml/31.acc.c.rds"))
d = readRDS(glue("{dird}/41_ml/41.fea.imp.rds"))
fo = glue("{dirf}/f4.pdf")
ggarrange(as.grob(a),b,c,d, nrow=4, ncol=1,
          heights=c(1.5,1,1,2), labels=LETTERS[1:4]) %>%
    ggexport(filename=fo, width=8, height=12)
#}}}

#{{{ Fig 5
a = readRDS(glue("{dird}/17_cluster/22.all.expr.rds"))
b = readRDS(glue("{dird}/41_ml/28.coexp.acc.rds"))
fo = glue("{dirf}/f5.pdf")
ggarrange(a,b, nrow=2, ncol=1,
          heights=c(2,1), labels=LETTERS[1:2]) %>%
    ggexport(filename=fo, width=8, height=11)
#}}}

#{{{ Fig 6
a = readRDS(glue("{dird}/41_ml/32.auroc.rds")) +
    scale_y_continuous(name='AUROC', expand=expansion(mult=c(0,.05)))
b = readRDS(glue("{dird}/41_ml/36.dde.bar.rds")) +
    theme(legend.title=element_text(size=8),legend.text=element_text(size=7)) +
    theme(axis.text.x=element_text(size=7,angle=10,vjust=1,hjust=.7))
c = readRDS(glue("{dird}/41_ml/38.dde.cis.rds")) +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) +
    theme(axis.text.x=element_text(size=7,angle=10,vjust=1,hjust=.7))
ggarrange(a, b, c, nrow=1, ncol=3, labels=LETTERS[1:3],
          widths=c(3,2,2), heights=c(2,2)) %>%
ggexport(filename="f6.pdf", width=8, height=7)
#}}}

# f7 is done

#{{{ Fig S8
a = readRDS('sf08a.rds')
b = readRDS('sf08b.rds') + o_margin(.2,.0,.0,.3)
c = readRDS('sf08c.rds') + o_margin(.0,.0,.3,.6)
bc = ggarrange(b, c, nrow=2, ncol=1, labels=LETTERS[2:3], heights=c(1,1.7))
ggarrange(a, bc, nrow=1, ncol=2, labels=LETTERS[1],
          widths=c(2,2.5), heights=c(2,2)) %>%
ggexport(filename="sf08.pdf", width=7, height=7)
#}}}

# Fig S9 is done
