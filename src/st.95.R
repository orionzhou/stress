source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
a = readRDS('f1a.rds') +
    theme(plot.margin = margin(.2,0,.2,.4, "lines")) +
    theme(legend.position = c(1,1), legend.justification = c(1,1))
b = readRDS('f1b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="f1.pdf", width=9, height=6)
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
          widths=c(1.2,2), heights=c(1,1)) %>%
ggexport(filename="f2.pdf", width=9, height=5)
#}}}

# Fig S3 is already done
# Fig S4 is Zach's TF enrichment figure
# Fig S5 is Erika's HSF/CBF expression pattern figure

#{{{ Fig 3
a = image_read_pdf("regions.pdf")
b = readRDS('f3b.rds') + o_margin(.2,.2,0,.2)
ggarrange(as.grob(a), b, nrow=2, ncol=1,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,3)) %>%
ggexport(filename="f3.pdf", width=5, height=7)
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

#{{{ Fig S11
a = readRDS('sf11a.rds') +
    theme(plot.margin = margin(.2,0,.2,.2, "lines"))
b = readRDS('sf11b.rds') +
    theme(plot.margin = margin(.2,0,.2,0, "lines"))
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(1.1,2), heights=c(2,2)) %>%
ggexport(filename="sf11.pdf", width=10.5, height=8)
#}}}

#{{{ Fig S12
a = readRDS('sf12a.rds')
b = readRDS('sf12b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf12.pdf", width=10, height=8)
#}}}

#{{{ Fig S13
a = readRDS('sf13a.rds')
b = readRDS('sf13b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:5],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="sf13.pdf", width=10, height=8)
#}}}

#{{{ Fig 6
a = readRDS('f6a.rds')
b = readRDS('f6b.rds')
ggarrange(a, b, nrow=1, ncol=2, labels=LETTERS[1:2],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename="f6.pdf", width=8, height=7)
#}}}




