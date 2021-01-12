source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
a = readRDS('f1a.rds') + o_margin(.2,0,.2,.4) +
    theme(legend.position = c(1,1), legend.justification = c(1,1))
b = readRDS('f1b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="f1r.pdf", width=9, height=6)
#}}}

#{{{ Fig S2
a = readRDS('sf02a.rds') +
    scale_x_continuous(expand = expansion(mult=c(1e-2,.4))) +
    theme(plot.margin = margin(.2,0,.2,.2, "lines"))
b = readRDS('sf02b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="sf02.pdf", width=9, height=6)
#}}}

#{{{ Fig 2
a = readRDS('f2a.rds') +
    theme(plot.margin = margin(.2,0,.2,.2, "lines"))
b = image_read_pdf("Fig2B.pdf")
ggarrange(a, as.grob(b), ncol=2,
          labels = LETTERS[1:2],
          widths=c(1,1.4), heights=c(1,1)) %>%
ggexport(filename="f2.pdf", width=8, height=5)
#}}}

# Fig S3 is already done
# Fig S4 is Zach's TF enrichment figure
# Fig S5 is Erika's HSF/CBF expression pattern figure

#{{{ Fig 3
a = image_read_pdf("regions.pdf")
b = readRDS('f3b.rds')
ggarrange(as.grob(a), b, nrow=2, ncol=1,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,3)) %>%
ggexport(filename="f3.pdf", width=5.5, height=7)
#}}}

# Fig S6 is done
# Fig S7 is done
# Fig 4 is done

#{{{ Fig 5
a = readRDS('f5a.rds')
b = readRDS('f5b.rds')
c = readRDS('f5c.rds')
d = readRDS('f5d.rds')
e = readRDS('f5e.rds')
ggarrange(a, b, c, d, e, nrow=2, ncol=3, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="f5.pdf", width=8, height=6)
#}}}

# Fig S11 is done

#{{{ Fig S12
a = readRDS('sf12a.rds') +
    theme(plot.margin = margin(1.2,0,.2,.2, "lines")) +
    theme(legend.justification=c(.5,-.5))
b = readRDS('sf12b.rds') +
    theme(plot.margin = margin(1.2,0,.2,0, "lines")) +
    theme(legend.justification=c(.5,-.5))
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(1.1,2), heights=c(2,2)) %>%
ggexport(filename="sf12.pdf", width=10.5, height=8)
#}}}

#{{{ Fig S13
a = readRDS('sf13a.rds')
b = readRDS('sf13b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf13.pdf", width=10, height=8)
#}}}

#{{{ Fig S14
a = readRDS('sf14a.rds')
b = readRDS('sf14b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf14.pdf", width=10, height=8)
#}}}

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

#{{{ Fig S15
a = readRDS('sf15a.rds')
b = readRDS('sf15b.rds') + o_margin(.3,.3,.3,.3)
b = ggarrange(b, NULL, nrow=2, ncol=1, heights=c(1,1.8))
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf15.pdf", width=6, height=7)
#}}}

#{{{ Fig S16
a = readRDS('sf16a.rds')
b = readRDS('sf16b.rds')
ggarrange(a, b, nrow=2, ncol=1, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf16.pdf", width=6, height=6)
#}}}

