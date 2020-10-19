source('functions.R')
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



