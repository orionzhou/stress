source('functions.R')
dirw = file.path(dird, '01_exp_design')

#{{{ create snk/nf sample list
fo = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fo)
#}}}

#{{{ BRB-seq design
fi = file.path(dirw, 'brbseq_layout.xlsx')
ti1 = read_xlsx(fi, 'batch1')
ti2 = read_xlsx(fi, 'batch2')
ti3 = read_xlsx(fi, 'batch3')

to = rbind(ti1, ti2, ti3)
to %>% count(batch)
to %>% count(SampleID) %>% count(n)
to2 = to %>% mutate(Genotype='',Treatment=SampleID, Replicate=1) %>%
    mutate(directory=file.path('/scratch.global/zhoux379/stress/data',batch,'demul')) %>%
    mutate(SampleID=str_c(batch,coord,sep='-')) %>%
    select(SampleID, Tissue=batch, Genotype, Treatment, Replicate, directory,
        file_prefix=coord)

fo = file.path(dirw, '05.brb.local.tsv')
write_tsv(to2, fo)
#}}}

#{{{ plot design
load_fonts()
dirw = dird
fas = c('\uf4d8', '\uf00d')
fa1 = fas[1]; fa2 = fas[2]

plot_design <- function(tp, tix, tiy, box.wd=.5, box.ht=.5, leg.pos.x=0, mar.left=3, x.off=.5, x3=F, leg.nrow=3, xlab=F) {
    #{{{
    ff='fas'; fa1='\uf4d8'; fa2='\uf00d'
    ff='icm'; fa1='\ue900'; fa2='\uf05e'
    tfa = tibble(Treatment=c('Control','Cold','Heat'),
      faw = c('\uf6c4','\uf2dc','\uf185'),
      fat = c('\uf2c9','\uf2cb','\uf2c7')) %>%
        mutate(fawt = str_c(faw,fat,sep='')) %>% inner_join(tiy, by='Treatment')
    tpb = tiy %>% mutate(ymin = y-.5, ymax=y+.5)
    tpb2 = tp %>% distinct(x0, y0) %>%
        mutate(xmin=x0-box.wd, xmax=x0+box.wd, ymin=y0-box.ht, ymax=y0+box.ht)
    tpl = tp %>% filter(lost)
    tpt = tp %>% filter(Timepoint==min(tp$Timepoint), Treatment=='Control')
    xb = min(tpb2$xmin); xe = max(tpb2$xmax); yb = min(tpb2$ymin); ye = max(tpb2$ymax)
    cols = pal_aaas()(10)
    if(length(unique(tp$Genotype)) > 10)
        cols = c(pal_aaas()(10),pal_simpsons()(15))[c(1:3,11:25,4:10)]
    px = ggplot(tp) +
        geom_rect(data=tpb, aes(xmin=xb,xmax=xe,ymin=ymin,ymax=ymax,fill=Treatment), alpha=.3) +
        geom_rect(data=tpb2, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=NA, color='black', size=.2) +
        geom_text(aes(x, y, color=Genotype), label=fa1, family=ff, size=6, vjust=-4)+
        geom_text(data=tpl, aes(x, y), label=fa2, color='black', alpha=.7, family='fas', size=5) +
        geom_point(data=tix, aes(x,y=ye+.1), shape=17, size=2.5) +
        geom_text(data=tiy, aes(x=xb-x.off, y=y+.1, label=Treatment), size=3, vjust=0) +
        geom_text(data=tfa, aes(x=xb-x.off, y=y-.1, label=fawt),size=4,family='fas',vjust=1) +
        scale_y_reverse(breaks=tiy$y, labels=tiy$Treatment, expand=c(0,0)) +
        scale_x_continuous(name='Hours after Treatment', breaks=tix$x, labels=tix$Timepoint, expand=expand_scale(mult=c(0,.01))) +
        expand_limits(y = c(yb-.2, ye)) +
        expand_limits(x = xb - 2*x.off) +
        scale_color_manual(values=cols) +
        scale_fill_manual(values=c('springgreen','steelblue1','firebrick1')) +
        otheme(legend.pos='none',legend.dir='v',legend.title=F,
               margin = c(0,0,0,mar.left),
               xtitle=xlab, xtext=T, ytext=F) +
        theme(legend.position=c(leg.pos.x,.85), legend.justification=c(1,1), legend.spacing=unit(0,'lines'), legend.key.size = unit(0,'lines')) +
        theme(panel.border = element_blank()) +
        guides(fill='none', col=guide_legend(nrow=leg.nrow, label.vjust=1, label.hjust=0))
    if(x3) px = px + geom_text(aes(x, y, color=Genotype), label='x3', size=2, hjust=-.6, vjust=0, show.legend=F)
    px
    #}}}
}

#{{{ TC
fi = file.path(dird, 'samples.xlsx')
gts = c('B73','Mo17','W22')
trs = c("Control",'Cold','Heat')
tms = c(0,.5,1,1.5,2,3,4,8,25)
ti = crossing(Genotype=gts,Treatment=trs,Timepoint=tms) %>%
    filter(Treatment == 'Control' | (Treatment!='Control' & Timepoint!=0)) %>%
    mutate(lost=ifelse(Genotype=='W22' & Timepoint %in% c(1.5,3,8), T, F)) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(Treatment=factor(Treatment, levels=trs)) %>%
    arrange(Treatment, Timepoint, Genotype)
tix1 = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x=ifelse(Timepoint==8, 6, Timepoint)) %>%
    mutate(x=ifelse(Timepoint==25, 10, x))
tiy1 = ti %>% distinct(Treatment) %>% mutate(y = 1:n())
tig1 = ti %>% distinct(Genotype) %>% mutate(ox=0, oy= (1:n() - (1+n())/2)/n())
tp1 = ti %>% inner_join(tix1, by='Timepoint') %>% rename(x0 = x) %>%
    inner_join(tiy1, by=c('Treatment')) %>% rename(y0 = y) %>%
    inner_join(tig1, by=c('Genotype')) %>%
    mutate(x = x0+ox, y = y0+oy)
#}}}
#{{{ HY
gts = c('B73','Mo17','W22','BxM', 'BxW','MxW')
trs = c("Control",'Cold','Heat')
tms = c(0,1,25)
ti = crossing(Genotype=gts,Treatment=trs,Timepoint=tms) %>%
    filter(Treatment == 'Control' | (Treatment!='Control' & Timepoint!=0)) %>%
    mutate(lost = F) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(Treatment=factor(Treatment, levels=trs)) %>%
    arrange(Treatment, Timepoint, Genotype)
ti = crossing(Genotype=gts,Treatment=trs,Timepoint=tms) %>%
    filter(Treatment == 'Control' | (Treatment!='Control' & Timepoint!=0)) %>%
    mutate(lost = F) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(Treatment=factor(Treatment, levels=trs)) %>%
    arrange(Treatment, Timepoint, Genotype)
tix2 = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x = Timepoint) %>%
    mutate(x=ifelse(Timepoint==25, 4, x))
tiy2 = ti %>% distinct(Treatment) %>% mutate(y = 1:n())
n=2; ix = (1:n - (1+n)/2)/(2*n)
n=3; iy = (1:n - (1+n)/2)/n
tig2 = ti %>% distinct(Genotype) %>% arrange(Genotype) %>%
    mutate(ox = rep(ix, each=3), oy = rep(iy, 2))
tp2 = ti %>% inner_join(tix2, by='Timepoint') %>% rename(x0 = x) %>%
    inner_join(tiy2, by=c('Treatment')) %>% rename(y0 = y) %>%
    inner_join(tig2, by=c('Genotype')) %>%
    mutate(x = x0+ox, y = y0+oy)
#}}}
#{{{ NM
fi = file.path(dird, 'samples.xlsx')
ti = read_xlsx(fi, 'NM')
gts = unique(ti$Genotype)
gts1 = c('B73','Mo17','W22')
gts = c('B73','Mo17','W22', gts[! gts %in% gts1])
trs = c("Control",'Cold')
tms = c(1,25)
ti = ti %>% select(SampleID,Genotype,Treatment,Timepoint) %>%
    mutate(lost = ifelse(SampleID %in% c("NM11"), T, F)) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(Treatment=factor(Treatment, levels=trs)) %>%
    arrange(Treatment, Timepoint, Genotype)
tix3 = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x=ifelse(Timepoint==25, 4, Timepoint))
tiy3 = ti %>% distinct(Treatment) %>% mutate(y = 1:n())
n=5; ix = (1:n - (1+n)/2)/n
n=5; iy = (1:n - (1+n)/2)/n
tig3 = ti %>% distinct(Genotype) %>% arrange(Genotype) %>%
    mutate(ox = rep(ix, each=5), oy = rep(iy, 5))
tp3 = ti %>% inner_join(tix3, by='Timepoint') %>% rename(x0 = x) %>%
    inner_join(tiy3, by=c('Treatment')) %>% rename(y0 = y) %>%
    inner_join(tig3, by=c('Genotype')) %>%
    mutate(x = x0+ox, y = y0+oy)
#}}}

p1 = plot_design(tp1,tix1,tiy1,box.wd=.25,leg.pos.x=0,mar.left=3.2)
fo = file.path(dirw, 'test.pdf')
#ggsave(p1, file=fo, width=8, height=3)
p2 = plot_design(tp2,tix2,tiy2,box.wd=.3,leg.pos.x=0,x.off=.2,mar.left=3.2,x3=T, leg.nrow=6)
fo = file.path(dirw, 'test2.pdf')
#ggsave(p2, file=fo, width=8, height=3)
p3 = plot_design(tp3,tix3,tiy3,box.wd=.5,leg.pos.x=0,x.off=.2,mar.left=11.5,leg.nrow=9,xlab=T)
fo = file.path(dirw, 'test3.pdf')
#ggsave(p3, file=fo, width=8, height=3)

fo = sprintf('%s/00.design.pdf', dirf)
ggarrange(p1, p2, p3, common.legend = F,
    nrow=3, ncol=1, heights=c(1,1,1),
    labels=c("Experiment 1",'Experiment 2','Experiment 3'), label.x=-.075) %>%
    ggexport(filename = fo, width=8, height=8)
#}}}


