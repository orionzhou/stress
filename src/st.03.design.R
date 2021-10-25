source('functions.R')
dirw = glue('{dird}/01_exp_design')

#{{{ create snk/nf sample list
fo = file.path(dirw, '01.meta.tsv')
write_tsv(th, fo)
#}}}

#{{{ BRB-seq design
fi = file.path(dirw, 'brbseq_layout.xlsx')
ti1 = read_xlsx(fi, 'batch1')
ti2 = read_xlsx(fi, 'batch2')
ti3 = read_xlsx(fi, 'batch3')

to = rbind(ti1, ti2, ti3)
to %>% dplyr::count(batch)
to %>% dplyr::count(SampleID) %>% dplyr::count(n)
to2 = to %>% mutate(Genotype='',Treatment=SampleID, Replicate=1) %>%
    mutate(directory=file.path('/scratch.global/zhoux379/stress/data',batch,'demul')) %>%
    mutate(SampleID=str_c(batch,coord,sep='-')) %>%
    select(SampleID, Tissue=batch, Genotype, Treatment, Replicate, directory,
        file_prefix=coord)

fo = file.path(dirw, '05.brb.local.tsv')
write_tsv(to2, fo)
#}}}

exps = c("TC", 'HY','NM')
tx = tibble(ExpID=exps, ExpTxt=sprintf("Experiment %d (%s)", 1:3, exps))
trs = c("Control",'Cold','Heat')
xlab = 'Hours since stress onset'

#{{{ process temperature data
load_fonts()
require(fs)
require(lubridate)
diri = file.path(dird, '03_temp_button')
ti = tibble(fp = dir_ls(file.path(diri, exps))) %>%
    filter(str_detect(fp, "\\.csv$")) %>%
    mutate(fn = basename(fp)) %>%
    mutate(fn = str_replace(fn, '\\.csv','')) %>%
    separate(fn, c("ExpID","stress","treatment",'rep'), sep="_", fill='right') %>%
    mutate(x = map(fp, read_csv, skip=14)) %>%
    select(ExpID,stress,treatment,rep, x) %>% unnest(x)

tp0 = tibble(ExpID = exps, datetime0=c('06/10/19 11:00:00', '06/20/19 11:00:00', '07/02/19 11:00:00')) %>%
    mutate(datetime0=mdy_hms(datetime0))

tp = ti %>% mutate(datetime = mdy_hms(`Date/Time`)) %>%
    select(ExpID, stress, treatment, rep, datetime, temp=Value) %>%
    inner_join(tp0, by='ExpID') %>%
    inner_join(tx, by='ExpID') %>%
    mutate(has = time_length(datetime-datetime0,'hour')) %>%
    mutate(stress = factor(stress, levels=stresses)) %>%
    mutate(treatment=as.double(str_replace(treatment, 'has', ''))) %>%
    mutate(treatment=ifelse(treatment==24, 25, treatment)) %>%
    mutate(treatment=factor(treatment, levels=sort(unique(treatment)))) %>%
    filter(is.na(rep) | rep %in% c("rep1",'rackA')) %>%
    group_by(ExpTxt, stress, treatment, has) %>%
    summarise(temp = mean(temp)) %>% ungroup() %>%
    filter(treatment == 25, has >= 0, has <= 25) %>% select(-treatment)

fo = file.path(diri, 'temperature.tsv')
write_tsv(tp, fo)
#}}}

#{{{ read & plot temperature data
fi = file.path(dird, '03_temp_button', 'temperature.tsv')
tp = read_tsv(fi) %>% mutate(stress=factor(stress, levels=trs)) %>%
    filter(ExpTxt == 'Experiment 1 (TC)')
tpa = tp %>% arrange(ExpTxt, stress, -has) %>%
    group_by(ExpTxt, stress) %>% slice(1) %>% ungroup()
tfa = tibble(stress=trs,
  faw = c('\uf6c4','\uf2dc','\uf185'),
  fat = c('\uf2c9','\uf2cb','\uf2c7')) %>%
    mutate(fawt = str_c(faw,sep='')) %>% inner_join(tpa, by='stress')
col0='corn silk'; colc='light sky blue'; colh='dark orange'
cols3 = pal_npg()(5)[c(3,2,1)]
cols3 = c('black',colc,colh)
p_temp = ggplot(tp, aes(has, temp)) +
    geom_rect(aes(xmin=14,xmax=22,ymin=-Inf,ymax=Inf), fill='#D3D3D3') +
    geom_line(aes(color=stress)) +
    geom_point(aes(color=stress), size=.5) +
    #geom_text(data=tfa, aes(x=25+.2, y=temp-.1, label=fawt),size=4,family='fas',vjust=1,hjust=0) +
    geom_text(data=tpa, aes(x=25+.3, y=temp+1, label=stress), hjust=.5, vjust=0, size=2.5) +
    #facet_wrap(~ExpTxt, ncol=1) +
    scale_x_continuous(name=xlab, breaks=c(0,1,14,22,25), expand=expansion(mult=c(.02,.06))) +
    scale_y_continuous(name=expression("Temperature " ( degree*C )), expand=expansion(mult=c(.02,.08))) +
    scale_color_manual(values=cols3) +
    otheme(legend.pos='none', margin=c(.1,.2,.1,.5),
           xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, xgrid=T) +
    theme(legend.spacing=unit(1,'lines'), legend.key.size = unit(1,'lines'))
fo = file.path(dirw, 'tempature.pdf')
#ggsave(fo, p_temp, width=6, height=6)
#}}}

#{{{ build design tibble
tms = c(0,.5,1,1.5,2,3,4,8,25)
ti1 = crossing(Genotype=gts3,Timepoint=tms) %>%
    filter(!(Genotype=='W22' & Timepoint %in% c(1.5,3,8))) %>%
    mutate(Genotype=factor(Genotype, levels=gts3)) %>%
    mutate(y = -as.integer(Genotype)) %>%
    mutate(ExpID='TC')
#
tms = c(0,1,25)
ti2 = crossing(Genotype=gts6,Timepoint=tms) %>%
    mutate(Genotype=factor(Genotype, levels=gts6)) %>%
    mutate(y = -10-as.integer(Genotype)) %>%
    mutate(ExpID='HY')
#
gts = c('B73','B97','CML322','...','PH207','Tx303','W22')
tms = c(1,25)
ti3 = crossing(Genotype=gts,Timepoint=tms) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(y = -20-as.integer(Genotype)) %>%
    mutate(ExpID='NM')
tp0 = rbind(ti1,ti2,ti3) %>% inner_join(tx, by='ExpID')
#}}}
ti = ti1 %>% inner_join(tx, by='ExpID')
#{{{ plot design
tpx = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x=ifelse(Timepoint==8, 6, Timepoint)) %>%
    mutate(x=ifelse(Timepoint==25, 10, x))
tpy = ti %>% distinct(y, Genotype)
tp = ti %>% inner_join(tpx, by='Timepoint')
tpa = tp %>% filter(Genotype!='...')
tpb = tp %>% filter(Genotype=='...')
tpr = tpa %>% mutate(txt=ifelse(ExpID=='HY', 'x3', 'x1'))
darkB = 6 + 6 * 4/17; darkE = 6 + 14 * 4/17
#
p = ggplot(tpa, aes(x, y)) +
    geom_rect(aes(xmin=darkB,xmax=darkE,ymin=-Inf,ymax=Inf), fill=alpha('#D3D3D3',.1)) +
    geom_point(aes(color=Genotype), size=5) +
    #geom_text(data=tpb, aes(x,y, label='...'), size=3) +
    geom_text(data=tpr, aes(x,y, label=txt),color='white',size=3,vjust=.5,hjust=.5) +
    scale_x_continuous(name=xlab, breaks=tpx$x, labels=tpx$Timepoint, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(breaks=tpy$y, labels=tpy$Genotype, expand=expansion(mult=c(.35,.35))) +
    scale_color_aaas() +#values=cols11) +
    otheme(legend.pos='none', margin=c(0,.3,0,.3),
           xtitle=T, ytitle=F, xtext=T, ytext=T, xtick=T, xgrid=T) +
    theme(axis.text.y = element_text(color='black')) +
    theme(plot.background=element_rect(fill=NA, color=NA))
#}}}
#{{{ overlay plot
col0='corn silk'; colc='light sky blue'; colh='dark orange'
p0 = p +
    ggtitle('Heat') + theme(plot.title=element_text(size=9, margin=margin(0,0,0,0))) +
    theme(axis.text.x=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(colh,.8),color=NA))
p1 = p +
    ggtitle('Cold') + theme(plot.title=element_text(size=9, margin=margin(0,0,00,0))) +
    theme(axis.text.x=element_text(color=NA)) +
    theme(axis.text.y=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(colc,.8),color=NA))
p2 = p +
    ggtitle('Control') + theme(plot.title=element_text(size=9, margin=margin(0,0,0,0))) +
    theme(axis.text.y=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(col0,.8),color=NA))
#
p = ggdraw() +
    draw_plot(p0, 0,.1,.9,.9)+
    draw_plot(p1,.05,.05,.9,.9)+
    draw_plot(p2,.1,0,.9,.9)
#}}}
pa=p
ti = ti2 %>% inner_join(tx, by='ExpID')
#{{{ plot design
tpx = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x=ifelse(Timepoint==8, 6, Timepoint)) %>%
    mutate(x=ifelse(Timepoint==25, 10, x))
tpy = ti %>% distinct(y, Genotype)
tp = ti %>% inner_join(tpx, by='Timepoint')
tpa = tp %>% filter(Genotype!='...')
tpb = tp %>% filter(Genotype=='...')
tpr = tpa %>% mutate(txt=ifelse(ExpID=='HY', 'x3', 'x1'))
darkB = 6 + 6 * 4/17; darkE = 6 + 14 * 4/17
#
p = ggplot(tpa, aes(x, y)) +
    geom_rect(aes(xmin=darkB,xmax=darkE,ymin=-Inf,ymax=Inf), fill=alpha('#D3D3D3',.1)) +
    geom_point(aes(color=Genotype), size=5) +
    #geom_text(data=tpb, aes(x,y, label='...'), size=3) +
    geom_text(data=tpr, aes(x,y, label=txt),color='white',size=3,vjust=.5,hjust=.5) +
    scale_x_continuous(name=xlab, breaks=tpx$x, labels=tpx$Timepoint, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(breaks=tpy$y, labels=tpy$Genotype, expand=expansion(mult=c(.15,.15))) +
    scale_color_aaas() +#values=cols11) +
    otheme(legend.pos='none', margin=c(0,.3,0,.3),
           xtitle=T, ytitle=F, xtext=T, ytext=T, xtick=T, xgrid=T) +
    theme(axis.text.y = element_text(color='black')) +
    theme(plot.background=element_rect(fill=NA, color=NA))
#}}}
#{{{ overlay plot
col0='corn silk'; colc='light sky blue'; colh='dark orange'
p0 = p +
    ggtitle('Heat') + theme(plot.title=element_text(size=9, margin=margin(0,0,0,0))) +
    theme(axis.text.x=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(colh,.8),color=NA))
p1 = p +
    ggtitle('Cold') + theme(plot.title=element_text(size=9, margin=margin(0,0,00,0))) +
    theme(axis.text.x=element_text(color=NA)) +
    theme(axis.text.y=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(colc,.8),color=NA))
p2 = p +
    ggtitle('Control') + theme(plot.title=element_text(size=9, margin=margin(0,0,0,0))) +
    theme(axis.text.y=element_text(color=NA)) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background=element_rect(fill=alpha(col0,.8),color=NA))
p0z = p0 + theme(axis.title.x=element_text(size=8, color=NA))
p1z = p1 + theme(axis.title.x=element_text(size=8, color=NA))
p2z = p2 + theme(axis.title.x=element_text(size=8))
#
p = ggdraw() +
    draw_plot(p0, 0,.1,.9,.9)+
    draw_plot(p1,.05,.05,.9,.9)+
    draw_plot(p2,.1,0,.9,.9)
pz = ggdraw() +
    draw_plot(p0z, 0,.1,.9,.9)+
    draw_plot(p1z,.05,.05,.9,.9)+
    draw_plot(p2z,.1,0,.9,.9)
#}}}
pb=p
pbz=pz
ti = ti3 %>% inner_join(tx, by='ExpID')
#{{{ plot design
tpx = ti %>% distinct(Timepoint) %>% arrange(Timepoint) %>%
    mutate(x=ifelse(Timepoint==8, 6, Timepoint)) %>%
    mutate(x=ifelse(Timepoint==25, 10, x))
tpy = ti %>% distinct(y, Genotype)
tp = ti %>% inner_join(tpx, by='Timepoint')
tpa = tp %>% filter(Genotype!='...')
tpb = tp %>% filter(Genotype=='...')
tpr = tpa %>% mutate(txt=ifelse(ExpID=='HY', 'x3', 'x1'))
darkB = 6 + 6 * 4/17; darkE = 6 + 14 * 4/17
#
p = ggplot(tpa, aes(x, y)) +
    geom_rect(aes(xmin=darkB,xmax=darkE,ymin=-Inf,ymax=Inf), fill=alpha('#D3D3D3',.1)) +
    geom_point(aes(color=Genotype), size=5) +
    geom_text(data=tpb, aes(x,y, label='...'), size=3) +
    geom_text(data=tpr, aes(x,y, label=txt),color='white',size=3,vjust=.5,hjust=.5) +
    scale_x_continuous(name=xlab, breaks=tpx$x, labels=tpx$Timepoint, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(breaks=tpy$y, labels=tpy$Genotype, expand=expansion(mult=c(.15,.15))) +
    scale_color_aaas() +#values=cols11) +
    otheme(legend.pos='none', margin=c(0,.3,0,.3),
           xtitle=T, ytitle=F, xtext=T, ytext=T, xtick=T, xgrid=T) +
    theme(axis.text.y = element_text(color='black')) +
    theme(plot.background=element_rect(fill=NA, color=NA))
#}}}
#{{{ overlay plot
col0='corn silk'; colc='light sky blue'; colh='dark orange'
p1 = p +
    ggtitle('Cold') + theme(plot.title=element_text(size=9, margin=margin(0,0,00,0))) +
    theme(axis.text.x=element_text(color=NA)) +
    theme(axis.title.x=element_text(color=NA)) +
    theme(panel.background=element_rect(fill=alpha(colc,.8),color=NA))
p2 = p +
    ggtitle('Control') + theme(plot.title=element_text(size=9, margin=margin(0,0,0,0))) +
    theme(axis.text.y=element_text(color=NA)) +
    theme(panel.background=element_rect(fill=alpha(col0,.8),color=NA))
#
p = ggdraw() +
    draw_plot(p1, 0,.1,.9,.9)+
    draw_plot(p2,.05,.05,.9,.9)
#}}}
pc=p
p = ggarrange(p_temp, pa, pb, pc, labels=LETTERS[1:4],
              ncol=1, heights=c(1,.7,.9,1), widths=c(1.3,1))

fo = glue('{dirf}/f1a.rds')
saveRDS(pbz, fo)
fo = glue('{dirf}/f1b.rds')
saveRDS(p_temp, fo)
fo = glue('{dirw}/design.pdf')
ggexport(p, filename=fo, width=5, height=7)


#{{{ RNA+ATAC design
gts3 = c("B73", "Mo17", 'B73xMo17')
tiss5 = c("embryo_imbibed_seed","coleoptile_tip",
          "seedling_leaf", "seedling_meristem","seedling_root")
tiss18 = c("blade_v12", "flag_leaf",
    "husk", "sheath", "auricle",
    "floret", "tassel_stem",
    "internode", "root_0DAP", "radicle_root",
    "silk", "tassel", "spikelets", "ear",
    "kernel_14DAP", "endosperm_14DAP", "embryo_27DAP", "endosperm_27DAP")
tiss23 = c(tiss5, tiss18)
tp = crossing(tissue=tiss23, genotype=gts3) %>%
    mutate(tissue=factor(tissue, levels=tiss23)) %>%
    mutate(genotype=factor(genotype,levels=gts3))
tpx = tp %>% distinct(genotype) %>% arrange(genotype) %>% mutate(x=1:n())
tpy = tp %>% distinct(tissue) %>% arrange(tissue) %>% mutate(y=n():1)
tp = tp %>% inner_join(tpx, by='genotype') %>%
    inner_join(tpy, by = 'tissue')
tpr = tp %>% mutate(txt='x5')

p = ggplot(tp, aes(x, y)) +
    geom_point(aes(color=genotype), size=5) +
    #geom_text(data=tpb, aes(x,y, label='...'), size=3) +
    geom_text(data=tpr, aes(x,y, label=txt),color='white',size=3,vjust=.5,hjust=.5) +
    geom_hline(yintercept=18.5,color='black',size=.5) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$genotype,
                       expand=expansion(add=c(.5,.9)), position='top') +
    scale_y_continuous(breaks=tpy$y, labels=tpy$tissue, expand=expansion(mult=c(.02,.02))) +
    scale_color_aaas() +
    annotate('text',x=4,y=21,size=2.5,hjust=.5,vjust=.5,label='Growth chamber') +
    annotate('text',x=4,y=10,size=2.5,hjust=.5,vjust=.5,label='Field') +
    otheme(legend.pos='none', margin=c(.3,.3,.3,.3),
           xtitle=F, ytitle=F, xtext=T, ytext=T, xtick=F, xgrid=F) +
    theme(axis.text.y = element_text(color='black')) +
    theme(plot.background=element_rect(fill=NA, color=NA))
fo = glue('{dirw}/atac.design.pdf')
ggexport(p, filename=fo, width=4, height=5)
#}}}


