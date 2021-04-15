source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
b = readRDS('f1a.rds') + o_margin(.1,.2,0,.4)
a = readRDS('f1b.rds') + o_margin(.5,.2,0,.4)
c = readRDS('f1c.rds') + o_margin(.1,.2,.5,.5)
d = readRDS('f1d.rds')
pab = ggarrange(a, b, nrow=1, ncol=2, labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2))
pcd = ggarrange(c,d, nrow=1, ncol=2, labels = LETTERS[3:4],
          widths=c(2,2), heights=c(2,2))
pabc = ggarrange(a, b, c, nrow=3, ncol=1, labels = LETTERS[1:3],
          widths=c(2,2), heights=c(2,1.8,4))
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
a = readRDS('sf01a.rds')
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

#{{{ Fig S2
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

#{{{ Fig S3
a = readRDS('sf03a.rds')
b = readRDS('sf03b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf03.pdf", width=10, height=8)
#}}}

#{{{ Fig S4
a = readRDS('sf04a.rds')
b = readRDS('sf04b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf04.pdf", width=10, height=8)
#}}}

#{{{ Fig S5
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
a = readRDS('f3a.rds') +
    theme(plot.margin = margin(.2,0,.2,.2, "lines"))
b = image_read_pdf("Fig3B.pdf")
ggarrange(a, as.grob(b), ncol=2,
          labels = LETTERS[1:2],
          widths=c(1,1.4), heights=c(1,1)) %>%
ggexport(filename="f3.pdf", width=8, height=5)
#}}}

# Fig S6 is already done
# Fig S7 is Zach's TF enrichment figure
# Fig S8-9 is Erika's HSF/CBF expression pattern figure

#{{{ Fig 4
a = image_read_pdf("regions.pdf")
b = readRDS('f4b.rds')
c = readRDS('f4c.rds') +
    theme(legend.position = c(.05,0), legend.justification = c(0,0))
ab = ggarrange(as.grob(a), b, nrow=2, ncol=1,
               heights=c(2,3), labels=LETTERS[1:2])
ggarrange(ab, c, nrow=1, ncol=2,
          labels = c('', 'C'),
          widths=c(1.2,2)) %>%
ggexport(filename="f4.pdf", width=10, height=8)
#}}}

# Fig 5 is done
# Fig S10 is done
# Fig S11 is done

#{{{ Fig 6
a = readRDS('f6a.rds') +
    scale_y_continuous(name='AUROC', expand=expansion(mult=c(0,.05)))
b = readRDS('f6b.rds') +
    theme(legend.title=element_text(size=8),legend.text=element_text(size=7)) +
    theme(axis.text.x=element_text(size=7,angle=10,vjust=1,hjust=.7))
c = readRDS('f6c.rds') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) +
    theme(axis.text.x=element_text(size=7,angle=10,vjust=1,hjust=.7))
ggarrange(a, b, c, nrow=1, ncol=3, labels=LETTERS[1:3],
          widths=c(3,2,2), heights=c(2,2)) %>%
ggexport(filename="f6.pdf", width=8, height=7)
#}}}

#{{{ Fig 7
a = readRDS('f7a.rds')
b = readRDS('f7b.rds')
c = readRDS('f7c.rds')
fo = glue("{dirf}/f7.pdf")
ab = ggarrange(a, b, nrow=1, widths=c(1,2), labels=LETTERS[1:2])
ggarrange(ab, c, nrow=2, heights=c(1,1), labels=c('','C')) %>%
    ggexport(filename=fo, width=7, height=7)
#}}}

#{{{ Fig S12
a = readRDS('sf12a.rds')
b = readRDS('sf12b.rds') + o_margin(.3,.3,.3,.3)
b = ggarrange(b, NULL, nrow=2, ncol=1, heights=c(1,1.8))
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf12.pdf", width=6, height=7)
#}}}

#{{{ Fig S13
a = readRDS('sf13a.rds')
b = readRDS('sf13b.rds')
ggarrange(a, b, nrow=2, ncol=1, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf13.pdf", width=6, height=6)
#}}}

