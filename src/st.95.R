source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
a = readRDS('f.1a.rds') +
    theme(plot.margin = unit(c(.2,0,.2,.4), "lines")) +
    theme(legend.position = c(1,1), legend.justification = c(1,1))
b = readRDS('f.1b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="f1.pdf", width=9, height=6)
#}}}

#{{{ Fig S2
a = readRDS('sf.2a.rds') +
    scale_x_continuous(expand = expansion(mult=c(1e-2,.4))) +
    theme(plot.margin = unit(c(.2,0,.2,.2), "lines"))
b = readRDS('sf.2b.rds')
ggarrange(a, b, nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="sf2.pdf", width=9, height=6)
#}}}

#{{{ Fig 2
a = readRDS('f.2a.rds') +
    theme(plot.margin = unit(c(.2,0,.2,.2), "lines"))
b = image_read_pdf("Fig2B.pdf")
ggarrange(a, as.grob(b), nrow=1, ncol=2,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,2)) %>%
ggexport(filename="f2.pdf", width=9, height=6)
#}}}

# Fig S3 is already done
# Fig S4 is Zach's TF enrichment figure
# Fig S5 is Erika's HSF/CBF expression pattern figure

#{{{ Fig 3
a = image_read_pdf("regions.pdf")
b = readRDS('f.3b.rds') +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines"))
ggarrange(as.grob(a), b, nrow=2, ncol=1,
          labels = LETTERS[1:2],
          widths=c(2,2), heights=c(2,3)) %>%
ggexport(filename="f3.pdf", width=5, height=9)
#}}}

# Fig S6 is done
# Fig S7 is done

#{{{ Fig 4
a = readRDS('f.4a.rds')
b = readRDS('f.4b.rds')
c = readRDS('f.4c.rds')
d = readRDS('f.4d.rds')
ggarrange(a, b, c, d, nrow=2, ncol=2,
          labels = c("A - F1", "B - AUROC", "C - F1", "D - AUROC"),
          widths=c(2,2), heights=c(2,1.5)) %>%
ggexport(filename="f4.pdf", width=6, height=6)
#}}}


