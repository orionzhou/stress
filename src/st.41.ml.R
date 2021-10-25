source('functions.R')
require(progress)
#require(rstatix)
dirw = glue('{dird}/41_ml')
dirn = glue('{dird}/41_ml/00_nf')

#{{{ prepare training input: gene lists + motifs
fi = glue("{dird}/25_dreme/05.best.mtfs.rds")
r5 = readRDS(fi)
tk = r5$tk %>% filter(n_mtf >= 10)

fk = glue("{dirw}/05.tk.tsv")
write_tsv(tk %>% select(bid, n_mtf), fk)
diro = glue('{dirn}/03_motif_lists')
if(!dir.exists(diro)) dir.create(diro)
bkg = c(.2688,.2312,.2312,.2688)
tk %>% mutate(mtfs = map(mtfs, 'mtf')) %>%
    mutate(fo = glue("{diro}/{bid}.meme")) %>%
    mutate(l = map2(mtfs, fo, write_meme, bkg=bkg, overwrite=T))

# training scheme 1 (B), 2 (BMW_r) and 3 (BMW_nr)
tc123 = r5$tc %>% mutate(train=scope) %>% select(train,everything())
# training scheme 4 (variable)
tbc = r5$tc %>% filter(scope=='BMW') %>% select(bid,cid)
fi = glue("{dird}/17_cluster/50_modules/var1.rds")
tc4 = readRDS(fi) %>% mutate(train='var', scope='BMW') %>%
    inner_join(tbc, by='cid') %>%
    select(train, everything())

trains = c('B','BMW','BMW_nr', 'var')
tc = tc123 %>% bind_rows(tc4) %>%# filter(cid %in% tk$cid) %>%
    mutate(train = factor(train, levels=trains)) %>%
    arrange(train, cid) %>%
    mutate(tid = str_c('t', str_pad(1:n(), width=2,pad='0'))) %>%
    select(train,tid,scope,bid, everything()) %>% print(n=30)

r6 = list(tk=tk, tc=tc)
fo = glue("{dirw}/06.tk.tc.rds")
saveRDS(r6,fo)

diro = glue('{dirn}/05_gene_lists')
if(!dir.exists(diro)) dir.create(diro)
tc %>% mutate(fo = glue("{diro}/{tid}.tsv")) %>% mutate(j=map2(ts, fo, write_tsv))
#}}}

#{{{ motif meta plots - [obsolete]
fi = glue("{dirw}/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk
#{{{ function
fimo_locate <- function(mid,fm,fg) {
    #{{{
    cmd = glue("fimo.py prepare_ml {fg} {fm} tmp.bed --epi umr --nfea {mid} --fmt long")
    system(cmd)
    ti = read_tsv('tmp.bed', col_names=c("gid",'start','end','mid'))
    system("rm tmp.bed")
    ti
    #}}}
}
kmer_locate <- function(kmers,fg) {
    #{{{
    cmd = glue("kmer.py locate --fg \"{fg}\" {kmers} tmp.tsv")
    system(cmd)
    ti = read_tsv('tmp.tsv')
    system("rm tmp.tsv")
    ti
    #}}}
}
read_mtf_loc <- function(fi) {
    #{{{
    itvs = c(seq(0,4000,by=200), seq(4050,8050,by=200))
    tp = read_tsv(fi) %>%
        mutate(pos = (start+end)/2) %>%
        mutate(bin = cut(pos, breaks=itvs, include.lowest=T)) %>%
        mutate(bin = as.numeric(bin)) %>%
        distinct(fid,gid,bin)
    tp
    #}}}
}
#}}}

#{{{ prepare module gene lists to get  motif locations
tc0 = tc %>% filter(train=='B')
ts1 = tc0 %>% select(cid, cond, note, gids)
ts2 = tc0 %>% distinct(cond, gids_c) %>% mutate(cid=c('c29','c49')) %>%
    mutate(note='bg') %>% select(cid,cond,note,gids=gids_c)
ts = ts1 %>% bind_rows(ts2) %>% arrange(cid) %>%
    unnest(gids) %>% rename(gid=gids)

fo = glue("{dirw}/51.mod.genes.rds")
saveRDS(ts, fo)
to = ts %>% distinct(gid) %>% mutate(status=1)
fo = glue("{dirw}/51.mod.genes.tsv")
write_tsv(to, fo, na='')
#}}}

#{{{ read
fg = normalizePath(glue("{dirw}/51.mod.genes.tsv"))
fmd = glue("{dirw}/51.mod.genes.rds")
md = readRDS(fmd)
tc1 = tc %>% filter(train=='B') %>% select(bid,cid,cond,note)
rb = tk %>% inner_join(tc1, by='bid') %>%
    unnest(mtfs) %>%
    mutate(rate=pos/ng, rate.c=neg/ng_c) %>%
    mutate(fnote = glue("{bin} {epi}")) %>%
    select(bid,cid,cond,note,i,pval,rate,rate.c,mid,fid,fname,fnote)
rb1 = rb %>% filter(rate.c < .2)
#}}}
rb %>% count(bid,cid)
rb1 %>% count(bid,cid)

#{{{ explore
plot_tss_meta <- function(mid,fm,bid,cond,note,fid,fnote,fo,fg,md) {
    #{{{
    ti = fimo_locate(mid,fm,fg)
    mds = md %>% count(cid,cond,note) %>% rename(nt=n) %>%
        mutate(mod = glue("{cond}: {note} ({nt})")) %>%
        arrange(cid) %>% mutate(mod = as_factor(mod))
    itvs = c(seq(0,4000,by=200), seq(4050,8050,by=200))
    itvs = c(seq(0,4000,by=200))
    tp = ti %>%
        mutate(pos = (start+end)/2) %>%
        mutate(bin = cut(pos, breaks=itvs, include.lowest=T)) %>%
        mutate(bin = as.numeric(bin)) %>%
        distinct(gid,bin) %>%
        inner_join(md, by='gid') %>%
        inner_join(mds,  by='cid') %>%
        count(mod,nt, bin) %>% rename(nh=n) %>%
        mutate(prop = nh/nt)
    #
    tit = glue("{cond} {note} | {fnote}")
    tpx = tibble(x=c(.5,10.5,20.5,31.5,41.5),lab=c('-2kb','TSS','+2kb/-2kb','TTS','+2kb'))
    tpx = tibble(x=c(.5,10.5,20.5),lab=c('-2kb','TSS','+2kb'))
    #{{{
    p = ggplot(tp, aes(x=bin,y=prop)) +
        geom_line(aes(), size=.5, na.rm = F) +
        geom_point(aes(), size=1, na.rm = F) +
        scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
        scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
        scale_color_aaas(name='strand') +
        #scale_shape(labels=types) +
        #scale_linetype(labels=types) +
        facet_wrap(~mod, scale='free_x', ncol=5) +
        ggtitle(tit) +
        otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
               strip.style='white',margin = c(.3,.3,.3,.3),
               xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
        theme(plot.title=element_text(hjust=.5, size=10)) +
        guides(fill='none')
    #}}}
    #fo = glue("{dirw}/53_metaplots/{bid}_{fid}.pdf")
    ggsave(p, file=fo, width = 10, height = 5)
    #}}}
}
rb2 = rb %>%# filter(cid=='c10') %>%
    group_by(cid) %>% arrange(i) %>% slice(1:30) %>% ungroup() %>%
    mutate(fm = glue("{dirn}/03_motif_lists/{bid}.meme")) %>%
    mutate(fo = glue("{dirw}/53_metaplots/{bid}_{i}.pdf")) %>%
    mutate(x = pmap(list(mid,fm,bid,cond,note,fid,fnote,fo), plot_tss_meta,
                    fg=fg, md=md))
#}}}

#{{{ extract TFBS locations for picked motifs
fc = glue("{dirw}/config.xlsx")
cfg = read_xlsx(fc) %>% fill(bid, .direction='down')
tp1 = rb %>% inner_join(cfg, by=c('bid','fid')) %>%
    mutate(fm = glue("{dirn}/03_motif_lists/{bid}.meme")) %>%
    mutate(x = pmap(list(mid,fm), fimo_locate, fg=fg))
fo = glue("{dirw}/54.rds")
saveRDS(tp1, fo)
#}}}
fi = glue("{dirw}/54.rds")
tp1 = readRDS(fi)

tc2 = tc %>% distinct(cid) %>% mutate(i = as.integer(str_sub(cid,2,2))) %>%
    mutate(bg = ifelse(i<=2, 'c29', 'c49'))
tc2b = tc2 %>% distinct(bg, i) %>% rename(cid=bg)
tc2 = tc2 %>% bind_rows(tc2b)
tc2s = tc2 %>% group_by(i) %>% summarise(cids = list(cid)) %>% ungroup()
tc2 = tc2 %>% inner_join(tc2s, by='i') %>% select(cid, cids)

plot_meta_tss1 <- function(ti, tit,cids, md) {
    #{{{
    #{{{ prepare
    mds = md %>% filter(cid %in% cids) %>%
        count(cid,cond,note) %>% rename(nt=n) %>%
        mutate(note = ifelse(note=='bg', 'control', note)) %>%
        mutate(mod = glue("{cond}: {note} ({nt})")) %>%
        arrange(cid) %>% mutate(mod = as_factor(mod)) %>%
        mutate(modtype = ifelse(str_detect(cid,'9$'), 'bg','a'))
    itvs = c(seq(0,4000,by=200))
    tp = ti %>%
        mutate(pos = (start+end)/2) %>%
        mutate(bin = cut(pos, breaks=itvs, include.lowest=T)) %>%
        mutate(bin = as.numeric(bin)) %>%
        distinct(pnl,gid,bin) %>%
        mutate(pnl = as_factor(pnl)) %>%
        inner_join(md, by='gid') %>%
        inner_join(mds,  by='cid') %>%
        count(pnl,mod,modtype,nt, bin) %>% rename(nh=n) %>%
        complete(bin, nesting(pnl,mod,modtype,nt), fill=list(nh=0)) %>%
        mutate(prop = nh/nt)
    #}}}
    tpx = tibble(x=c(.5,10.5,20.5),lab=c('-2kb','TSS','+2kb'))
    #{{{
    cols6 = c(pal_npg()(length(cids)-1), 'black')
    ggplot(tp, aes(x=bin,y=prop, color=mod, shape=mod, linetype=modtype)) +
        geom_line(aes(), size=.5, na.rm = F) +
        geom_point(aes(), size=1, na.rm = F) +
        scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
        scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
        scale_color_manual(name='', values=cols6) +
        scale_shape(name='') +
        scale_linetype(name='') +
        facet_wrap(~pnl, scale='free_y', ncol=4) +
        otheme(legend.pos='none', legend.dir='v', legend.title=T,
               panel.spacing=.2, margin = c(0,9,0,.2),
               xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
        ggtitle(tit) +
        theme(plot.title=element_text(hjust=.5, size=10, face='bold', margin=margin(0,0,0,0))) +
        theme(legend.position=c(1,.5),legend.justification=c(0,.5)) +
        guides(fill='none',linetype='none')
    #}}}
    #}}}
}
tp2 = tp1 %>% filter(cid!='c30', !mid %in% c('b08_l446_1','b08_l465_2')) %>%
    mutate(tit = glue("{cond}: {note}")) %>%
    mutate(fname = ifelse(nchar(fname)>10, glue("{str_sub(fname,1,3)}...{str_sub(fname,-3,-1)}"), fname)) %>%
    mutate(pnl = glue("{fname} ({fnote})")) %>%
    select(tit,bid,cid,i,pnl,x) %>%
    inner_join(tc2, by='cid') %>%
    unnest(x) %>%
    group_by(bid,cid,tit,cids) %>% nest() %>% ungroup() %>% rename(x=data)
tp3 = tp2 %>%
    mutate(p = pmap(list(x,tit,cids), plot_meta_tss1, md=md))

ps = tp3$p
p = ggarrange(
              ggarrange(ps[[1]], NULL, ncol=2, widths=c(1,.23)),
              ggarrange(ps[[2]], NULL, ncol=2, widths=c(1,.6)),
              ggarrange(ps[[3]], NULL, ncol=2, widths=c(1,0)),
              ggarrange(ps[[4]], NULL, ncol=2, widths=c(1,.23)),
              ggarrange(ps[[5]], NULL, ncol=2, widths=c(1,.23)),
              ggarrange(ps[[6]], NULL, ncol=2, widths=c(1,.23)),
              ggarrange(ps[[7]], NULL, ncol=2, widths=c(1,0)),
              nrow=7, ncol=1, heights=c(2,2))
fo = glue("{dirw}/55.metaplot.pdf")
ggexport(p, filename=fo, width=8, height=10)
fo = glue("{dirf}/sf09.pdf")
#}}}

#{{{ create ML training config
fi = glue("{dirw}/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk
#{{{ params
bins = c(
         "TSS:-500","TSS:+500","TSS:-/+500",
         "TSS:-1k",'TSS:+1k','TSS:-/+1k',
         "TSS:-2k",'TSS:+2k','TSS:-/+2k'
         #"TTS:-500","TTS:+500","TTS:-/+500",
         #"TTS:-2k",'TTS:+2k','TTS:-/+2k',
         #"TSS:-/+2k,TTS:-500",
         #"TSS:-/+2k,TTS:-1k"
)
epis = c('raw','umr','acrL')
nfeas = c('top30', 'top50', 'top100', 'top200')
mods = c('zoops', 'anr')
bins1 = c("TSS:-/+2k")
epis1 = c('umr')
nfeas1 = c('top100')
mods1 = c('zoops')
#}}}

ctag = 'b1'; train = 'B'
ctag = 'b2'; train = 'BMW'
ctag = 'b3'; train = 'var'
ctag = 'b4'; train = 'BMW_nr'
#{{{
tc1 = tc %>% filter(train == !!train)
to1 = tc1 %>% select(tid,bid) %>%
  crossing(bin = bins, epi = epis1, nfea = nfeas1, mod = mods1)
to2 = tc1 %>% select(tid,bid) %>%
  crossing(bin = bins1, epi = epis, nfea = nfeas, mod = mods)
to = to1 %>% bind_rows(to2) %>%
    distinct(tid, bid, bin, epi, nfea, mod) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(tid, bid, bin, epi, nfea, mod) %>%
    mutate(did = sprintf("d%03d", 1:n())) %>%
    select(did, tid, bid, bin, epi, nfea, mod)
#
fo = glue("{dirn}/{ctag}.cfg.rds")
saveRDS(to, fo)
fo = glue("{dirn}/{ctag}.cfg.tsv")
write_tsv(to, fo)
#}}}
#}}}
# run pipeline to build models

######## MODEL EVALUATION ########
fi = glue("{dirw}/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk
tc0 = tc %>% mutate(cond_note=glue("{cond}: {note}")) %>%
    select(bid, cond_note) %>% mutate(cond_note=as_factor(cond_note))

#{{{ process trained models
#{{{ params
bins = c(
         "TSS:-500","TSS:+500","TSS:-/+500",
         "TSS:-1k",'TSS:+1k','TSS:-/+1k',
         "TSS:-2k",'TSS:+2k','TSS:-/+2k'
         #"TTS:-500","TTS:+500","TTS:-/+500",
         #"TTS:-2k",'TTS:+2k','TTS:-/+2k',
         #"TSS:-/+2k,TTS:-500",
         #"TSS:-/+2k,TTS:-1k"
)
epis = c('raw','umr')
nfeas = c('top30', 'top50', 'top100', 'top200')
mods = c('zoops', 'anr')
bins1 = c("TSS:-/+2k")
epis1 = c('umr')
nfeas1 = c('top100')
mods1 = c('zoops')
#}}}
collect_models <- function(tag, dirw, bin='TSS:-/+2k', epi=c('raw','umr'),
                           nfea='top100', mod='zoops') {
  #{{{
  fi = glue("{dirw}/00_nf/{tag}.cfg.rds")
  th = readRDS(fi)
  fi = glue("{dirw}/00_nf/{tag}/43.ml.rds")
  ml = readRDS(fi)
  metric = ml %>% select(did=sid, perm, metric) %>%
    unnest(metric) %>% spread(metric, estimate) %>%
    dplyr::rename(f1=f_meas, auroc=roc_auc, auprc=pr_auc) %>%
    inner_join(th, by='did')
  vis = ml %>% select(did=sid, perm, vis) %>%
    unnest(vis) %>%
    inner_join(th, by='did')
  best = th %>%
      filter(bin %in% !!bin, epi %in% !!epi, nfea %in% !!nfea, mod %in% !!mod) %>%
      inner_join(ml, by=c('did'='sid')) %>%
      filter(!is.na(fit)) %>% select(did,tid,bid,bin,epi,nfea,mod,fit,metric)
  list(metric=metric, best=best, vis=vis)
  #}}}
}
tags = c("b3")
tags = c("b1")
tags = c("b1",'b2','b4')

tm = tibble(tag = tags) %>%
  mutate(x = map(tag, collect_models, dirw=dirw)) %>%
  mutate(metric = map(x, 'metric')) %>%
  mutate(vis = map(x, 'vis')) %>%
  mutate(best = map(x,'best')) %>% select(-x)
fm = glue('{dirw}/11.models.rds')
saveRDS(tm, fm)
# write models
tb = tm %>% select(tag, best) %>% unnest(best) %>%
    select(tag,tid,bid,did,bin,epi,nfea,mod)
fb = glue('{dirw}/12.best.models.tsv')
write_tsv(tb, fb)
# write models
tb = tm %>% select(tag, best) %>% unnest(best)
tb %>% mutate(fo = glue("{dirw}/23_models/{tid}_{epi}.rds")) %>%
    mutate(map2(fit, fo, saveRDS))
#}}}

#{{{ eval model performance in B - f4b-c
fm = glue('{dirw}/11.models.rds')
tm = readRDS(fm)

tag='b3'; train='var'
tag='b2'; train='BMW'
tag='b1'; train='B'
tag='b4'; train='BMW_nr'
#{{{ prepare tc1
bins = c(
         "TSS:-500","TSS:+500","TSS:-/+500",
         "TSS:-1k",'TSS:+1k','TSS:-/+1k',
         "TSS:-2k",'TSS:+2k','TSS:-/+2k'
         #"TTS:-500","TTS:+500","TTS:-/+500",
         #"TTS:-2k",'TTS:+2k','TTS:-/+2k',
         #"TSS:-/+2k,TTS:-500",
         #"TSS:-/+2k,TTS:-1k"
)
bins = bins %>% str_replace("TSS:", '') %>% str_replace('[\\-][\\/][[\\+]]', '+/-')
epis = c('all','umr','acrL')
nfeas = c('top30', 'top50', 'top100', 'top200')
mods = c('binary', 'quantitative')
epi_map = c('raw'='all','umr'='umr','acrL'='acrL')
mod_map = c('zoops'='binary','anr'='quantitative')
tc1 = tc %>% filter(train==!!train) %>%
    filter(str_detect(note, "all")) %>%
    mutate(note = str_replace(note, 'all ', '')) %>%
    mutate(note = str_replace(note, '-', '')) %>%
    select(bid,cond,note) %>%
    crossing(bin = bins, epi = epis, nfea = nfeas, mod = mods) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(bid, bin, epi, nfea, mod) %>%
    mutate(cond_note = as_factor(str_c(cond, note, sep=' '))) %>%
    mutate(cond_note_nfea_mod = as_factor(str_c(cond,note,nfea,mod,sep=': '))) %>%
    mutate(bin_epi = as_factor(str_c(bin, epi, sep=': '))) %>%
    mutate(cond_mod = as_factor(str_c(cond, mod, sep=': '))) %>%
    mutate(cond_epi = as_factor(str_c(cond, epi, sep=': '))) %>%
    mutate(bin_epi_nfea_mod = as_factor(str_c(bin,epi,nfea,mod,sep=': '))) %>%
    mutate(bin_epi_nfea = as_factor(str_c(bin, epi, nfea, sep=': '))) %>%
    mutate(bin_epi_mod = as_factor(str_c(bin, epi, mod, sep=': '))) %>%
    mutate(epi_nfea_mod = as_factor(str_c(epi, nfea, mod, sep=': '))) %>%
    mutate(bin_nfea_mod = as_factor(str_c(bin, nfea, mod, sep=': ')))
#}}}
tm1 = tm %>% filter(tag == !!tag) %>% select(tag,metric) %>% unnest(metric) %>%
    mutate(epi=epi_map[epi], mod=mod_map[mod]) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    mutate(bin = str_replace(bin, "TSS:", '')) %>%
    mutate(bin = str_replace(bin, '[\\-][\\/][[\\+]]', '+/-')) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    inner_join(tc1, by=c('bid','bin','epi','nfea','mod'))

#{{{ sd1a, sd1b
to = tm1 %>% select(bin,epi,nfea,mod,cond_note,
    accuracy,precision,sens,spec,auroc,auprc,f1) %>%
    gather(metric, v, -bin,-epi,-nfea, -mod, -cond_note) %>%
    group_by(bin,epi,nfea,mod,cond_note,metric) %>%
    summarise(mean=mean(v), std=sd(v)) %>% ungroup() %>%
    mutate(txt = glue("{number(mean,accuracy=.001)} ({number(std,accuracy=.001)})")) %>%
    select(-mean,-std) %>% spread(metric,txt) %>%
    arrange(bin, epi, nfea, mod, cond_note) %>%
    select(TSS_Bin=bin, Epi_Filter=epi, NumFeatures=nfea, MotifEncoding=mod,
        GeneSet=cond_note,
        Accuracy=accuracy,Precision=precision,
        Sensitivity=sens, Specificity=spec,
        AUROC=auroc, AUPRC=auprc, F1=f1)

fo = glue("{dirf}/sd1a.tsv")
fo = glue("{dirf}/sd1b.tsv")
write_tsv(to, fo)
#}}}

#{{{ functions
plot_barplot_sig <- function(ti, tis, metric='auroc', x='bin', grp='epi',
    pnl='cond_mod', x.rotate=F, tit=F, cols=pal_npg()(8)) {
  #{{{
  tp = ti %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = get(x), pnl = get(pnl), grp=get(grp)) %>%
    group_by(cond, note, cond_note, bin, epi, nfea, mod, x, pnl, grp) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup()
  tps0 = tp %>% mutate(y = score + sd + .01) %>%
      group_by(pnl, x) %>% summarise(y = max(y)) %>% ungroup()
  ymax = max(tps0$y)
  tps = tis %>% mutate(x = get(x), pnl = get(pnl)) %>%
      mutate(xmin=as.integer(x) - .2) %>%
      mutate(xmax=as.integer(x) + .2) %>%
      inner_join(tps0, by=c('pnl','x'))
    #{{{
    pd = position_dodge(width=.8)
    p = ggplot(tp) +
      geom_col(aes(x=x,y=score,color=grp),fill=NA, position=pd, width=.7) +
      geom_errorbar(aes(x=x,y=score,ymin=score-sd, ymax=score+sd, col=grp), width=.2, size=.3, position=pd) +
      geom_segment(data=tps, aes(x=xmin,xend=xmax,y=y,yend=y), size=.5) +
      geom_text(data=tps, aes(x=x,y=y,label=sig), vjust=-.3, size=2.5) +
      scale_x_discrete(expand=expansion(add=c(.5,.5))) +
      scale_y_continuous(name=str_to_upper(metric), expand=expansion(add=c(0,.03))) +
      coord_cartesian(ylim= c(.5, ymax)) +
      scale_fill_manual(name='', values=cols) +
      scale_color_manual(name='', values=cols) +
      facet_wrap(pnl ~ ., scale='free_y', nrow=1) +
      otheme(legend.pos='none', legend.dir='h', legend.box='h',legend.title=F,
             legend.spacing.x=.1, legend.spacing.y=.1,
             margin=c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=F, panel.spacing=.1,
             xtick=T, ytick=T, ytitle=T, xtext=T, ytext=T) +
      theme(legend.position = c(.5,1), legend.justification = c(.5,1)) +
      theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y=element_text(size=7)) +
      guides(fill='none')
    #}}}
  if(x.rotate)
    p = p +
      theme(axis.text.x=element_text(angle=25, hjust=1, vjust=1, size=7.5))
  if(tit) {
    tit1 = str_to_upper(ifelse(metric=='f1','f-score',metric))
    p = p + ggtitle(tit)
  }
  p
  #}}}
}
get_bar_stats <- function(ti, metric='auroc', x='bin', grp='epi',
    pnl='cond_mod') {
    #{{{
    ti %>%
        mutate(score = get(metric)) %>%
        filter(!is.na(score)) %>%
        mutate(x = get(x), pnl = get(pnl), grp=get(grp)) %>%
        #group_by(cond, note, cond_note, bin, epi, nfea, mod, x, pnl, grp) %>%
        group_by(pnl, x, grp) %>%
        summarise(scores=list(score), std=sd(score), avg=mean(score)) %>%
        ungroup()
    #}}}
}
plot_barplot <- function(ti,tis, x.rotate=F, ylab='AUROC', wd=.7, tit=F, cols=pal_npg()(8)) {
  #{{{
  tps0 = tp %>% mutate(y = avg + std + .01) %>%
      group_by(pnl, x) %>% summarise(y = max(y)) %>% ungroup()
  tps = tis
  ymax = max(tps0$y, tps$ysig)
    #{{{
    pd = position_dodge(width=.8)
    p = ggplot(tp) +
      geom_col(aes(x=x,y=avg,fill=grp,color=grp),alpha=.2, position=pd, width=wd) +
      geom_errorbar(aes(x=x,y=avg,ymin=avg-std, ymax=avg+std, col=grp), width=.2, size=.3, position=pd) +
      geom_segment(data=tps, aes(x=xb,xend=xe,y=ysig,yend=ysig), size=.5) +
      geom_text(data=tps, aes(x=xm,y=ysig,label=sig), vjust=-.3, size=2.5) +
      scale_x_discrete(expand=expansion(add=c(.5,.5))) +
      scale_y_continuous(name=str_to_upper(ylab), expand=expansion(add=c(0,.03))) +
      coord_cartesian(ylim= c(.5, ymax)) +
      scale_fill_manual(name='', values=cols) +
      scale_color_manual(name='', values=cols) +
      facet_wrap(pnl ~ ., scale='free_y', nrow=1) +
      otheme(legend.pos='none', legend.dir='h', legend.box='h',legend.title=F,
             legend.spacing.x=.1, legend.spacing.y=.1,
             margin=c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=F, panel.spacing=.1,
             xtick=T, ytick=T, ytitle=T, xtext=T, ytext=T) +
      theme(legend.position = c(.5,1), legend.justification = c(.5,1)) +
      theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y=element_text(size=7)) +
      guides(fill='none')
    #}}}
  if(x.rotate)
    p = p +
      theme(axis.text.x=element_text(angle=25, hjust=1, vjust=1, size=7.5))
  if(tit) {
    tit1 = str_to_upper(ifelse(ylab=='f1','f-score',ylab))
    p = p + ggtitle(tit)
  }
  p
  #}}}
}
cols14 = c(brewer.pal(4,'Blues')[2:4],brewer.pal(4,'Greens')[2:4],
    brewer.pal(4,'Purples')[2:4],brewer.pal(4,'Greys')[2:4],
    brewer.pal(4,'Reds')[2:3])
#}}}

metric = 'auroc'
#{{{ f4b
tp0 = tm1 %>%
    filter(epi=='umr', nfea=='top100', mod=='binary') %>%
    filter(!str_detect(bin, "1k")) %>%
    mutate(bin=str_replace(bin, "TSS:", "")) %>%
    mutate(bin1 = str_replace_all(bin, "[\\+\\-\\/]", '')) %>%
    mutate(bin2=ifelse(str_detect(bin,"^\\+/\\-"),"+/-",
                       ifelse(str_detect(bin,"^\\+"), "+", "-"))) %>%
    mutate(bin1 = as_factor(bin1)) %>%
    mutate(bin2 = as_factor(bin2))
tp = get_bar_stats(tp0, x='bin1', grp='bin2', pnl='cond_note')
#{{{ sig
tps0 = tp %>% group_by(pnl,x) %>% summarise(ymax = max(avg+std)) %>% ungroup()
tps1 = tp %>% select(-avg,-std) %>% spread(grp, scores) %>%
    mutate(sig = map2_dbl(`-`, `+/-`, ttest_signif)) %>%
    inner_join(tps0, by=c('pnl','x')) %>%
    mutate(x = as.numeric(x)) %>%
    mutate(xm = x, xb = x - .3, xe = x +.3, ysig = ymax + .04) %>%
    select(pnl,xm,xb,xe,ysig,sig)
tps2 = tp %>% select(-avg,-std) %>% spread(grp, scores) %>%
    mutate(sig = map2_dbl(`+`, `+/-`, ttest_signif)) %>%
    inner_join(tps0, by=c('pnl','x')) %>%
    mutate(x = as.numeric(x)) %>%
    mutate(xm = x + .15, xb = x, xe = x +.3, ysig = ymax + .01) %>%
    select(pnl,xm,xb,xe,ysig,sig)
tps0b = tp %>% group_by(pnl) %>% summarise(ymax = max(avg+std)) %>% ungroup()
tps3 = tp %>% filter(grp=='+/-') %>% select(-avg,-std) %>%
    spread(x, scores) %>%
    mutate(sig = map2_dbl(`500`, `2k`, ttest_signif)) %>%
    inner_join(tps0b, by=c('pnl')) %>%
    mutate(xm = 1.75, xb = 1.25, xe = 2.25, ysig = ymax + .07) %>%
    select(pnl,xm,xb,xe,ysig,sig)
tps = rbind(tps1,tps2,tps3) %>%
    group_by(pnl) %>% mutate(sig = p.adjust(sig, method='BH')) %>% ungroup() %>%
    mutate(sig = map_chr(sig, map_signif))
#}}}
#
cols3 = pal_simpsons()(3)
pb = plot_barplot(tp, tps, x.rotate=F, ylab=metric, cols=cols3) +
  theme(legend.position = c(.37,.95), legend.justification = c(.5,.5)) +
  theme(legend.direction='horizontal')
#}}}
fo = glue("{dirw}/31.acc.b.pdf")
ggsave(pb, file=fo, width=8, height=4)
fo = glue("{dirw}/31.acc.b.rds")
saveRDS(pb, fo)

#{{{ #[obsolete] f4c
tp0 = tm1 %>% filter(bin=='TSS:-/+2k',epi=='umr',nfea=='top100',mod=='binary')
tp = get_bar_stats(tp0, x='cond_note', grp='cond_note', pnl='bin_epi_nfea_mod')
#{{{ sig
tps0 = tp %>% group_by(pnl,x) %>% summarise(ymax = max(avg+std)) %>% ungroup()
tps0 = tp %>% select(-avg,-std,-grp) %>%
    mutate(x = as.numeric(x)) %>%
    spread(x, scores)
tps1 = tps0 %>% mutate(sig = map2_chr(`1`,`2`, eval_signif)) %>%
    mutate(xm = 1.5, xb=1, xe=2, ysig=.9)
tps2 = tps0 %>% mutate(sig = map2_chr(`1`,`3`, eval_signif)) %>%
    mutate(xm = 2, xb=1, xe = 3, ysig=.94)
tps3 = tps0 %>% mutate(sig = map2_chr(`4`,`5`, eval_signif)) %>%
    mutate(xm = 4.5, xb=4, xe = 5, ysig=.78)
tps4 = tps0 %>% mutate(sig = map2_chr(`4`,`6`, eval_signif)) %>%
    mutate(xm=5, xb=4, xe=6, ysig=.82)
tps5 = tps0 %>% mutate(sig = map2_chr(`7`,`8`, eval_signif)) %>%
    mutate(xm=7.5, xb=7, xe=8, ysig=.85)
tps6 = tps0 %>% mutate(sig = map2_chr(`7`,`9`, eval_signif)) %>%
    mutate(xm=8, xb=7, xe=9, ysig=.89)
tps7 = tps0 %>% mutate(sig = map2_chr(`10`,`11`, eval_signif)) %>%
    mutate(xm=10.5, xb=10, xe=11, ysig=.79)
tps = rbind(tps1,tps2,tps3,tps4,tps5,tps6,tps7) %>%
    select(pnl,xm,xb,xe,ysig,sig)
#}}}
pb = plot_barplot(tp, tps, x.rotate=T, ylab=metric, cols=rep('black',20)) +
  o_margin(0,.2,.2,2.1) +
  theme(legend.position = 'none')
#}}}

#{{{ f4c
tp0 = tm1 %>% filter(bin=='+/-2k', nfea=='top100', mod=='binary',
                     epi %in% c("all","umr",'acrL'))
tp = get_bar_stats(tp0, x='cond_note', grp='epi', pnl='bin_nfea_mod')
#{{{ sig
tps0 = tp %>% group_by(pnl,x) %>% summarise(ymax = max(avg+std)) %>% ungroup()
tps1 = tp %>% select(-avg,-std) %>% spread(grp, scores) %>%
    mutate(sig = map2_dbl(all, umr, ttest_signif)) %>%
    inner_join(tps0, by=c('pnl','x')) %>%
    mutate(x = as.numeric(x)) %>%
    mutate(xm = x-.15, xb = x - .3, xe = x, ysig = ymax + .03) %>%
    select(pnl,xm,xb,xe,ysig,sig)
tps2 = tp %>% select(-avg,-std) %>% spread(grp, scores) %>%
    mutate(sig = map2_dbl(umr, acrL, ttest_signif)) %>%
    inner_join(tps0, by=c('pnl','x')) %>%
    mutate(x = as.numeric(x)) %>%
    mutate(xm = x+.15, xb = x, xe = x +.3, ysig = ymax + .01) %>%
    select(pnl,xm,xb,xe,ysig,sig)
tps = rbind(tps1,tps2) %>%
    group_by(pnl) %>% mutate(sig = p.adjust(sig, method='BH')) %>% ungroup() %>%
    mutate(sig = map_chr(sig, map_signif))
#}}}
cols3 = pal_aaas()(5)[c(1,3,2)]
cols3 = pal_jco()(3)
#pc = plot_barplot(tp, tps, x.rotate=T, ylab=metric, cols=cols3) +
pc = plot_barplot(tp, tps, x.rotate=F, ylab=metric, cols=cols3, wd=.4) +
  o_margin(0,.2,.2,2.1) +
  theme(legend.position = c(.5,1), legend.justification = c(.5,1)) +
  theme(legend.direction='horizontal')
#}}}
fo = glue("{dirw}/31.acc.c.pdf")
ggsave(pc, file=fo, width=8, height=4)
fo = glue("{dirw}/31.acc.c.rds")
saveRDS(pc, fo)

fo = glue("{dirw}/31.bar.sig.{tag}.pdf")
#p %>% ggexport(filename=fo, width=6, height=8)

#{{{ ## (removed) sf10
#{{{ functions
plot_model_heatmap <- function(ti, metric='auroc', x='cond_note', y='bin',
                               pnl='', x.rotate=F, tit=F) {
  #{{{ plot
  tp = ti %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = get(x), y = get(y), pnl = get(pnl)) %>%
    group_by(cond, note, bin, epi, nfea, mod, x, y, pnl) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup() %>%
    mutate(lab = glue("{number(score,accuracy=.01)}%+-%{number(sd,accuracy=.01)}")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}"))
    #{{{
    tp = tp %>% mutate(y = fct_rev(y))
    swit = (min(tp$score) + max(tp$score)) / 2
    p = ggplot(tp, aes(x=x,y=y)) +
      geom_tile(aes(fill=score), na.rm = F) +
      geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2, parse=T) +
      #geom_vline(xintercept=tpy$x, color='blue') +
      scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
      scale_y_discrete(expand=c(0,0)) +
      scale_fill_gradientn(name='metric',colors=cols100v) +
      #scale_fill_viridis(name='normalized eigengene value') +
      scale_color_manual(values=c('black','white')) +
      #facet_wrap(pnl ~ ., ncol=1) +
      otheme(legend.pos='none', legend.dir='v', legend.title=F,
             margin = c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=T, panel.spacing=.1,
             xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
      theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y = element_text(size=7)) +
      #theme(axis.text.y = element_markdown(size=7.5)) +
      guides(color='none')
    #}}}
  if(x.rotate)
    p = p +
      theme(axis.text.x=element_text(angle=25, hjust=0, vjust=0, size=7.5))
  if(tit) {
    tit = str_to_upper(ifelse(metric=='f1','f-score',metric))
    p = p + ggtitle(tit) +
      theme(plot.title=element_text(hjust=0, size=10))
  }
  p
  #}}}
}
plot_model_barplot <- function(ti, metric='auroc', x='bin', fill='epi',
    pnl='cond_mod', x.rotate=F, tit=F, cols=pal_npg()(8)) {
  #{{{
  ymax = 1.0 * max(ti[[metric]])
  tp = ti %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = get(x), pnl = get(pnl)) %>%
    group_by(cond, note, cond_note, bin, epi, nfea, mod, x, pnl) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup()
    #{{{
    pd = position_dodge(width=.7)
    p = ggplot(tp, aes(x=x,y=score)) +
      geom_col(aes(fill = get(fill)), position=pd, width=.7) +
      geom_errorbar(aes(ymin=score-sd, ymax=score+sd, col=get(fill)), width=.2, size=.3, position=pd) +
      scale_x_discrete(expand=expansion(mult=c(.03,.03))) +
      scale_y_continuous(name=str_to_upper(metric), expand=expansion(mult=c(0,.05))) +
      coord_cartesian(ylim= c(.5, ymax)) +
      scale_fill_manual(name='', values=cols) +
      scale_color_manual(values=rep('black',20)) +
      facet_wrap(pnl ~ ., scale='free_y', ncol=1) +
      otheme(legend.pos='none', legend.dir='h', legend.box='h',legend.title=F,
             legend.spacing.x=.1, legend.spacing.y=.1,
             margin=c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=F, panel.spacing=.1,
             xtick=T, ytick=T, ytitle=T, xtext=T, ytext=T) +
      theme(legend.position = c(.5,1), legend.justification = c(.5,1)) +
      theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y=element_text(size=7)) +
      guides(color='none')
    #}}}
  if(x.rotate)
    p = p +
      theme(axis.text.x=element_text(angle=25, hjust=1, vjust=1, size=7.5))
  if(tit) {
    tit1 = str_to_upper(ifelse(metric=='f1','f-score',metric))
    p = p + ggtitle(tit)
  }
  p
  #}}}
}
#}}}
#{{{ sf10 - barplot
metric = 'auroc'
bids4 = tm1 %>% filter(str_detect(note, '^all')) %>% pull(bid)
#{{{ pa
tpa = tm1 %>% filter(bid %in% bids4, epi=='umr', nfea=='top100', mod=='binary')
cols14 = c(brewer.pal(4,'Blues')[2:4],brewer.pal(4,'Greens')[2:4],
    brewer.pal(4,'Purples')[2:4],brewer.pal(4,'Greys')[2:4],
    brewer.pal(4,'Reds')[2:3])
pa = plot_model_barplot(tpa, metric=metric,, x='cond_note', fill="bin",
                        pnl='epi_nfea_mod', x.rotate=F, cols=cols14) +
  theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
  theme(legend.direction='horizontal') +
  guides(fill=guide_legend(nrow=3)) +
  o_margin(2.5,.5,.2,.2)
#}}}
#{{{ pb
tpb = tm1 %>% filter(bin=='TSS:-/+2k', nfea=='top100', mod=='binary')
pb = plot_model_barplot(tpb, metric=metric,, x='cond_note', fill='epi',
                        pnl='bin_nfea_mod', x.rotate=T) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.5,1), legend.justification = c(.5,1))
#}}}
#{{{ pc
tpc = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', nfea=='top100')
pc = plot_model_barplot(tpc, metric=metric,, x='cond_note', fill="mod",
                        pnl='bin_epi_nfea', x.rotate=T, cols=pal_aaas()(4)) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.5,1), legend.justification = c(.5,1))
#}}}
#{{{ pd
tpd = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', mod=='binary')
tpd2 = tpd %>% group_by(bid,cond_note) %>%
    summarise(ymax = mean(auroc) + sd(auroc) + .02) %>% ungroup()
fx = glue("{dirw}/05.tk.tsv")
tx = read_tsv(fx) %>% inner_join(tpd2) %>% mutate(txt=glue("({n_mtf})"))
pd = plot_model_barplot(tpd, metric=metric,, x='cond_note', fill="nfea",
                        pnl='bin_epi_mod', x.rotate=T, cols=pal_simpsons()(8)) +
  geom_text(data=tx, aes(x=cond_note,y=ymax,label=txt), vjust=0, size=2.5) +
  o_margin(0,.2,.2,2.1) +
  theme(legend.position = c(.5,1), legend.justification = c(.5,1))
#}}}
p = ggarrange(pa, pc, pd, nrow=3, ncol=1,
          labels = LETTERS[1:4], heights=c(1.5, 1, 1.3))
#}}}
fo = glue("{dirw}/30.auroc.{tag}.pdf")
p %>% ggexport(filename=fo, width=6, height=7)
fo = glue("{dirf}/sf10.pdf")
p %>% ggexport(filename=fo, width=6, height=7)
# {{{ F1 scores
#fo = glue("{dirw}/30.f1.{tag}.pdf")
#p %>% ggexport(filename=fo, width=6, height=9)
#fo = glue("{dirf}/sf10.pdf")
#p %>% ggexport(filename=fo, width=6, height=9)
#}}}
#}}}
#}}}

#{{{ eval feature importance - f4d
fm = glue('{dirw}/11.models.rds')
tm = readRDS(fm)

tp1 = tk %>% select(bid,mtfs) %>% unnest(mtfs) %>%
    mutate(fc = (pos/ng)/(neg/ng_c)) %>%
    select(bid,i,mid,pval,fc,fid,fname) %>%
    arrange(bid, desc(pval), mid) %>%
    group_by(bid) %>% mutate(i = 1:n()) %>% ungroup()
tp2 = tm %>% select(vis) %>% unnest(vis) %>%
    rename(mid=Variable,score=Importance) %>%
    group_by(did,bid,bin,epi,nfea,mod,mid) %>%
    summarise(avg=mean(score), std=sd(score), q25=quantile(score,.225),
              q50=quantile(score,.5), q75=quantile(score,.75)) %>%
    ungroup()
tp = tp1 %>% inner_join(tp2,by=c('bid','mid')) %>%
    inner_join(tc0, by='bid') %>% filter(bid <= 'b11') %>%
    filter(str_detect(cond_note, " all ")) %>%
    mutate(cond_note = str_replace(cond_note, ": all", '')) %>%
    mutate(cond_note = str_replace(cond_note, "-", ''))

get_cor <- function(xs, ys) cor(xs, ys, method='kendall')
tp0 = tp %>%
    filter(bin=='TSS:-/+2k', epi=='umr', nfea=='top200', mod=='zoops') %>%
    arrange(bid, desc(q50)) %>%
    group_by(bid) %>% mutate(j=1:n()) %>% ungroup() %>%
    mutate(top5 = j<=5)
tpt = tp0 %>% arrange(bid, desc(q50)) %>%
    group_by(bid) %>% dplyr::slice(1:5) %>% ungroup()
tpl = tp0 %>%
    group_by(bid,bin,epi,nfea,mod,cond_note) %>%
    summarise(imax=max(i), ymax=max(q75), n=n(),
              spc=cor(i, avg, method='spearman'),
              p.raw=cor.test(i,avg, method='spearman')$p.value) %>%
    ungroup() %>%
    mutate(txt = map_chr(p.raw, map_signif)) %>%
    mutate(txt = glue("{n} motifs\nrho = {number(spc,accuracy=.01)} ({txt})")) %>%
    mutate(ymax = max(tp0$q75) * .95) %>%
    print(n=20)

p = ggplot(tp0) +
    geom_pointrange(aes(x=i,y=q50,ymin=q25,ymax=q75,col=top5), size=.2) +
    geom_text_repel(data=tpt, aes(x=i,y=q50,label=fname), size=2) +
    geom_text(data=tpl, aes(x=0, y=ymax, label=txt), size=2.5, hjust=0, vjust=1) +
    scale_x_continuous(name="Level of Motif Enrichment (least enriched -> most enriched)", expand=expansion(mult=c(.03,.03))) +
    scale_y_continuous(name='Feature Importance', expand=expansion(mult=c(.01,.01)),
        breaks=c(0,.003,.006), labels=c('0','.003','.006')) +
    scale_color_manual(values=c('black','red')) +
    facet_wrap(~cond_note, nrow=2, scale='free_x') +
    otheme(legend.pos='none', xtext=T,xtick=T,xtitle=T,ytitle=T,ytext=T,ytick=T, strip.compact=F)
fo = glue("{dirw}/41.fea.imp.pdf")
ggsave(p, file=fo, width=6, height=6)
fo = glue("{dirw}/41.fea.imp.rds")
saveRDS(p, fo)
#}}}

#{{{ gather & save model predictions for all BMW genes
#{{{ write all gene list
fi = glue('{dird}/15_de/05.rds')
x =  readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12
md = deg12 %>% pluck('ds',1) %>% select(gid) %>%
    crossing(gt = gts32) %>% mutate(status = 1) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    select(gid, status)

fo = glue("{dirw}/21.all.tsv")
write_tsv(md, fo)
#}}}
# run j21 fimo.py prepare_ml
# run j21 ml_predict

fb = glue('{dirw}/12.best.models.tsv')
to = read_tsv(fb) %>% select(tag, tid, bid) %>%
    mutate(fi = glue("{dirw}/24_ml_out/{tag}_{tid}.tsv")) %>%
    mutate(pred = map(fi, read_tsv)) %>% select(-fi)

fo = glue('{dirw}/25.model.pred.rds')
saveRDS(to, fo)
#}}}

#{{{ evaluate model performance for different gene categories - f5b, f6, st5
fp = glue('{dirw}/25.model.pred.rds')
pd = readRDS(fp) %>%
    inner_join(tc %>% select(tid,train,bid,cid,cond,note), by=c('tid','bid')) %>%
    select(tag,train,tid,bid,cid,cond,note,pred) %>% unnest(pred) %>%
    separate(gid, c('gt','gid'), sep='_')
#{{{ prepare variable response
fg = glue('{dird}/15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2
#
fi = glue('{dird}/16_ase/20.rds')
res = readRDS(fi)
ddeg = res$ddeg
ddeg2 = ddeg %>% filter(condB != 'Control0') %>%
    mutate(time=as.integer(str_sub(cond,5,6))) %>%
    mutate(cond=str_to_lower(str_sub(cond,1,4))) %>%
    mutate(st = glue("d{deg}")) %>%
    mutate(st = factor(st, levels=levels(td2$st))) %>%
    select(cond,time,qry,tgt,gid,st,reg)
#}}}

#{{{ compare auroc in different subsets - f5b, st5
#{{{ prepare 
fi = glue("{dird}/17_cluster/50_modules/degB.rds")
x0 = readRDS(fi) %>% mutate(drc = rep(c('up','down'),2)) %>%
    select(cid,cond,drc,ts) %>% unnest(ts) %>%
    separate(gid, c('gt','gid'), sep='_') %>% group_by(cid,cond,drc) %>% nest() %>%
    rename(ts = data) %>% ungroup() %>% select(-drc) %>%
    crossing(drc = c('up','down')) %>%
    select(cid,cond,drc,ts)
fi = glue("{dird}/17_cluster/50_modules/dmodB2.rds")
x1 = tc %>% filter(scope=='B', str_detect(note, '^all')) %>% distinct(tid) %>%
    mutate(cond=rep(c('cold','heat'), each=2)) %>%
    mutate(drc=rep(c('up','down'), 2))
md = readRDS(fi)
xm = md %>% select(cid,cond,ts) %>% unnest(ts) %>%
    separate(gid, c('gt','gid'), sep='_') %>% group_by(cid,cond) %>% nest() %>%
    rename(ts = data) %>% ungroup() %>%
    crossing(drc = c('up','down')) %>%
    bind_rows(x0) %>%
    inner_join(x1, by=c('cond','drc'))
xp = pd %>% filter(tag %in% c("b1"), gt=='B73') %>%
    mutate(pred=factor(pred,levels=c(1,0))) %>%
    select(tid, gid, pred, prob)
x3 = xm %>% unnest(ts) %>% inner_join(xp, by=c('tid','gid')) %>%
    group_by(tid, cid, cond, drc) %>% nest() %>% ungroup()

perm=15
x4 = x3 %>% crossing(i = 1:perm) %>%
    mutate(x = map2(data, i, eval_pred, down_sample=T)) %>%
    select(-i) %>% unnest(x) %>% filter(metric=='roc_auc')

acc = x4 %>% select(cid,cond,drc,estimate) %>%
    group_by(cid,cond,drc) %>%
    summarise(vs=list(estimate), mean=mean(estimate), sd=sd(estimate)) %>%
    ungroup() %>%
    mutate(txt=glue("{number(mean,accuracy=.001)} ({number(sd,accuracy=.001)})")) %>%
    mutate(txt = str_replace_all(txt, "0\\.", "\\."))

fo = glue("{dirw}/27.coexp.acc.rds")
saveRDS(acc, file=fo)
#}}}

fi = glue("{dirw}/27.coexp.acc.rds")
acc = readRDS(fi)
acc.deg = acc %>% filter(cid %in% glue("c{1:4}")) %>%
    select(cid,cond,drc,txt) %>% spread(drc,txt)
acc.mod = acc %>% filter(! cid %in% glue("c{1:4}")) %>%
    select(cid,cond,drc,txt) %>% spread(drc,txt) %>%
    print(n=50)

#{{{ f5b
conds = c('cold','heat')
drcs = c('up','down')
f_cfg = glue('{dirw}/../17_cluster/config.xlsx')
cfg = read_xlsx(f_cfg) %>% filter(!is.na(idx)) %>%
    mutate(pnl = glue("{str_sub(cond,0,1)}{str_sub(drc,0,1)}{idx}")) %>%
    select(pnl, cond, drc, cid) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(drc = factor(drc, levels=drcs))
cfg2 = tibble(pnl='all', cond=rep(c('cold','heat'),each=2),
              drc = rep(c('up','down'), 2), cid=glue("c{1:4}")) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(drc = factor(drc, levels=drcs))

tp = acc %>% inner_join(cfg %>% bind_rows(cfg2), by=c('cond','drc','cid')) %>%
    rename(x=pnl) %>% mutate(pnl=glue("{cond} {drc}regulated")) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(drc = factor(drc, levels=drcs)) %>%
    arrange(cond,drc,pnl, x) %>% mutate(pnl = as_factor(pnl)) %>%
    print(n=35)
#{{{ sig
tps0 = tp %>% filter(x=='all') %>% select(pnl, vs0=vs)
tps = tp %>% filter(x!='all') %>% select(pnl, x, mean, sd, vs) %>%
    inner_join(tps0, by='pnl') %>%
    mutate(sig = map2_chr(vs0, vs, ttest_signif)) %>%
    mutate(sig = map_chr(sig, map_signif)) %>%
    mutate(y = mean + sd)
#}}}

pdg = position_dodge(width=.8)
tp1 = tp %>% filter(x=='all')
ymax = .95
p = ggplot(tp, aes(x=x, y=mean)) +
  geom_hline(data=tp1, aes(yintercept=mean), color=pal_npg()(2)[2], size=.5, linetype='dashed') +
  geom_col(aes(fill=x=='all'), col='black', size=.3, width=.5, alpha=.6) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, size=.3) +
  geom_text(data=tps, aes(x=x,y=y+.005,label=sig), size=2,hjust=.5,vjust=0) +
  scale_x_discrete(expand=expansion(mult=c(.05,.05))) +
  scale_y_continuous(name='AUROC', expand=expansion(mult=c(.02,.02))) +
  coord_cartesian(ylim=c(.5,ymax)) +
  scale_fill_manual(name='', values=pal_npg()(5)) +
  scale_color_manual(name='', values=pal_npg()(5)) +
  facet_wrap(pnl ~ ., scale='free', ncol=2, strip.position='top') +
  otheme(legend.pos='none',legend.dir='h',legend.box='h',
         legend.title=F, legend.vjust=-1.5,
         margin = c(.2,.5,.2,.2), panel.border=T,
         strip.style='white', strip.compact=F, panel.spacing=.1,
         ytick=T, ytitle=T, ytext=T, xtext=T) +
  theme(legend.justification=c(.5,-.8)) +
  #theme(axis.text.x = element_text(angle=0,size=8, hjust=1,vjust=1)) +
  theme(strip.text.y = element_text(angle=0,size=8))
#
fo = glue("{dirw}/28.coexp.acc.pdf")
p %>% ggexport(filename=fo, width=7, height=5)
fo = glue("{dirw}/28.coexp.acc.rds")
saveRDS(p, fo)
#}}}

#{{{ st5
require(gridExtra)
f_cfg = glue('{dirw}/../17_cluster/config.xlsx')
cfg = read_xlsx(f_cfg) %>% filter(!is.na(idx)) %>%
    mutate(pnl = glue("{str_sub(cond,0,1)}{str_sub(drc,0,1)}{idx}")) %>%
    select(pnl, cond, drc, cid) %>%
    mutate(cond = factor(cond)) %>%
    mutate(drc = factor(drc, levels=c('up','down')))

to1 = cfg %>% inner_join(acc.mod %>% select(-cond), by='cid') %>%
    select(cond, drc, pnl, up, down) %>% arrange(cond,drc) %>% print(n=40)
to0 = acc.deg %>% mutate(drc = rep(c("up",'down'), 2)) %>%
    mutate(pnl = glue("{str_sub(cond,0,1)}{str_sub(drc,0,1)}0")) %>%
    select(cond, drc, pnl, up, down) %>%
    mutate(cond=factor(cond, levels=levels(to1$cond))) %>%
    mutate(drc=factor(drc, levels=levels(to1$drc)))
tt = to0 %>% bind_rows(to1) %>% arrange(cond,drc, pnl) %>%
    mutate(cond = as.character(cond), drc = as.character(drc)) %>%
    group_by(cond,drc) %>% mutate(drc = ifelse(row_number()==1, drc, '')) %>% ungroup() %>%
    group_by(cond) %>% mutate(cond = ifelse(row_number()==1, cond, '')) %>% ungroup() %>%
    mutate(pnl = ifelse(str_detect(pnl, '0$'), glue("all {drc}-regulated"), pnl)) %>%
    rename(Stress=cond, Direction=drc, Cluster=pnl, `AUROC^up`=up,
        `AUROC^down`=down)

x = tt %>%
    kbl(format='latex', escape=T, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    #collapse_rows(1, latex_hline='major', valign='middle', longtable_clean_cut=F) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = glue("{dirf}/st5.rds")
saveRDS(x, file=fo)

tt2 = tableGrob(tt, rows=NULL,
                theme=ttheme_default(base_size=9,
                                     padding=unit(c(1,3), "mm"),
                                     colhead=list(fg_params = list(parse=T))))
p = ggarrange(tt2, nrow=1, ncol=1, labels = '', heights=c(1,1))
fo = glue("{dirw}/28.coexp.acc.pdf")
p %>% ggexport(filename=fo, width=5, height=8)
#}}}
#}}}


#{{{ overall performance - f6a
xm = tc %>% filter(train=='BMW', str_detect(note,'all')) %>%
    select(cid,ts) %>% unnest(ts) %>%
    separate(gid, c('gt','gid'), sep='_')
#
perm=10
x = pd %>% filter(tag %in% c("b1",'b4')) %>%
    mutate(pred=factor(pred,levels=c(1,0))) %>%
    inner_join(xm, by=c('cid','gt','gid')) %>%
    select(-tag,-tid) %>%
    group_by(train,cid,cond,note,gt) %>% nest() %>% rename(pred=data) %>%
    crossing(i = 1:perm) %>%
    mutate(x = map2(pred, i, eval_pred, down_sample=T)) %>%
    select(-i) %>% unnest(x) %>%
    group_by(train,cid,cond,note,gt,metric) %>%
    summarise(mean=mean(estimate), sd=sd(estimate)) %>% ungroup()

#{{{ sd1c
to = x %>% mutate(txt=glue("{number(mean,accuracy=.001)} ({number(sd,accuracy=.001)})")) %>%
    mutate(cond_note = glue("{cond}: {note}")) %>%
    select(train, cond_note, gt, metric, txt) %>%
    spread(metric, txt) %>%
    select(Model=train, GeneSet=cond_note, Genotype=gt,
           Accuracy=accuracy, Precision=precision,
           Sensitivy=sens, Specificity=spec,
           AUROC=roc_auc, AUPRC=pr_auc, F1=f_meas) %>% print(n=24)

fo = glue("{dirf}/sd1c.tsv")
write_tsv(to, fo)
#}}}

#{{{ bar plot - f6a
tp = x %>% filter(metric=='roc_auc') %>% rename(score=mean) %>%
    #mutate(gt = factor(gt, levels=rev(gts3))) %>%
    #filter(str_detect(note, 'all')) %>%
    mutate(note = str_replace(note, "all ", "")) %>%
    mutate(note = str_replace(note, "-", "")) %>%
    mutate(cond_note = str_c(cond,note,sep=": ")) %>%
    #mutate(train = factor(train, levels=c('BMW_nr','BMW','B'))) %>%
    mutate(train = factor(train, levels=c('BMW_nr','B'))) %>%
    mutate(cond_note = as_factor(cond_note))
tps = tp %>% distinct(cid,cond_note)
pdg = position_dodge(width=.8)
ymax = 1
p = ggplot(tp, aes(x=train, y=score)) +
  geom_col(aes(fill=gt), col='black', size=.3, position=pdg, width=.5, alpha=.6) +
  geom_errorbar(aes(ymin=score-sd, ymax=score+sd, color=gt), width=.2, size=.3, position=pdg) +
  scale_x_discrete(expand=expansion(mult=c(.05,.1))) +
  scale_y_continuous(name='AUROC', expand=expansion(mult=c(0,.02))) +
  #coord_cartesian(ylim=c(.5,ymax)) +
  coord_flip(ylim= c(.5, .93)) +
  scale_fill_manual(name='', values=pal_npg()(5)) +
  scale_color_manual(name='', values=pal_npg()(5)) +
  facet_wrap(cond_note ~ ., scale='free_y', ncol=1, strip.position='top') +
  otheme(legend.pos='top.center.out',legend.dir='h',legend.box='h',
         legend.title=F, legend.vjust=-1.5,
         margin = c(.2,.0,.2,.2), xgrid=T, panel.border=F, axis=T,
         strip.style='white', strip.compact=F, panel.spacing=.1,
         xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T) +
  theme(legend.justification=c(.5,-.8)) +
  theme(axis.line.y = element_line(size=.5)) +
  theme(strip.text.y = element_text(angle=0,size=8))
#}}}
fo = glue("{dirw}/32.auroc.pdf")
ggsave(p, file=fo, width=4, height=6)
fo = glue("{dirw}/32.auroc.rds")
saveRDS(p, fo) 
#}}}

#{{{ accuracy in predicting variable genotype response - f6b-c
#{{{ prepare
smap = c("+"=1,"-"=-1,"="=0)
smap2 = c('1'='+','-1'="-",'0'="=")
get_acc2 <- function(ti) {
  #{{{
  ti %>%
    mutate(qSt = smap[str_sub(st,3,3)], tSt = smap[str_sub(st,5,5)]) %>%
    group_by(st) %>%
    summarise(n.tot = n(), nc = sum(qSt==qPred & tSt ==tPred)) %>% ungroup() %>%
    mutate(acc = nc/n.tot,
           lab = glue("{number(acc,accuracy=.01)}\n{nc} / {n.tot}")) %>%
    select(st, acc, lab)
  #}}}
}
pd3 = pd %>% filter(tag %in% c("b1",'b4')) %>%
    select(-tag) %>% filter(str_detect(note, "all")) %>%
    mutate(drc = ifelse(str_detect(note, "all up"), "up", "down")) %>%
    mutate(pred=ifelse(drc=='down' & pred==1, -1, pred))
to1 = td2 %>%
  left_join(ddeg2, by=c('cond','time','qry','tgt','gid','st')) %>%
  #filter(!is.na(reg) & reg == 'cis') %>%
  inner_join(pd3, by=c('cond','qry'='gt','gid')) %>% rename(qPred=pred) %>%
  inner_join(pd3, by=c('train','cid','cond','drc','note','tgt'='gt','gid')) %>%
  rename(tPred=pred) %>%
  #group_by(bat,cond,drc,note,qry,tgt,time) %>% nest() %>% ungroup() %>%
  mutate(pred1 = smap2[as.character(qPred)]) %>%
  mutate(pred2 = smap2[as.character(tPred)]) %>%
  mutate(pred = as.character(glue("A{pred1}B{pred2}")))
to2 = to1 %>%
  count(train,cid,cond,drc,note,qry,tgt,time, st, pred,reg)
#}}}

#{{{ write accurately predicted variable response gene list
tv1 = to1 %>%
    filter(tgt=='B73', pred %in% c("A=B+",'A=B-'), st %in% c('dA=B+','dA=B-'))
tv2 = to1  %>%
    filter(tgt=='B73', pred %in% c("A+B=",'A-B='), st %in% c('dA+B=','dA-B='))
tv = tv1 %>% bind_rows(tv2) %>%
    mutate(st = str_sub(as.character(st), 2, 5)) %>% filter(st==pred) %>%
    select(train, cond,time,qry,tgt,gid,st,reg) %>%
    arrange(train, cond, time, tgt, qry, gid)
tv %>% filter(train=='BMW_nr') %>% distinct(cond,gid) %>% count(cond)

fo = glue("{dird}/71_share/28.variable.genes.tsv")
write_tsv(tv, fo)
#}}}

#{{{ dde acc bar plot - f6b
#{{{ prepare
tp = to2 %>% filter(tgt=='B73') %>%
  filter((drc=='up' & st %in% c("dA+B+",'nA+B+','dA+B=','dA=B+')) |
         (drc=='down' & st %in% c('dA-B-','nA-B-','dA-B=', 'dA=B-'))) %>%
  mutate(st = as.character(st)) %>%
  mutate(st = ifelse(str_detect(st, 'A\\+B\\+'),'A+B+', st)) %>%
  mutate(st = ifelse(str_detect(st, 'A\\-B\\-'),'A-B-', st)) %>%
  mutate(st = str_replace(st, 'd', '')) %>%
  mutate(st = str_replace(st, 'n', '')) %>%
  mutate(st = str_replace_all(st, '[\\+\\-]', '*')) %>%
  mutate(pred = str_replace_all(pred, '[\\+\\-]', '*')) %>%
  mutate(time = as_factor(time)) %>%
  mutate(note = as_factor(note)) %>%
  #mutate(pred = as_factor(pred)) %>%
  mutate(pnl = glue("{train} model")) %>%
  #mutate(pnl = glue("{cond} {note}: {train}")) %>%
  #group_by(train,cid,cond,drc,note,st,pred,pnl) %>%
  group_by(pnl,st,pred) %>%
  summarise(n=sum(n)) %>% ungroup() %>%
  rename(tag1=st, tag2=pred)
#tps = tp %>% distinct(drc,cid, train, pnl) %>% arrange(cid,train,drc)
stm = c('A=B='='Neither respond',"A*B="="Only Mo17/W22 respond","A=B*"="Only B73 respond","A*B*"='Both respond')
tp = tp %>% mutate(pnl = as_factor(pnl)) %>%
    mutate(tag1 = factor(stm[tag1], levels=stm)) %>%
    mutate(tag2 = factor(stm[tag2], levels=stm))
#preds1 = c('A=B=',"A+B=","A=B+","A+B+")
#preds2 = c('A=B=',"A-B=","A=B-","A-B-")
#tp1 = tp %>% filter(drc == 'up') %>% mutate(tag2=factor(tag2,levels=preds1))
#tp2 = tp %>% filter(drc == 'down') %>% mutate(tag2=factor(tag2,levels=preds2))
#}}}
plot_dde_bar <- function(ti) {
#{{{
tags1 = levels(ti$tag1)
tags2 = levels(ti$tag2)
sep = ifelse(T, " ", "\n")
col1 = 'black'; col2 = 'black'
tp = ti %>%
    arrange(tag1, desc(tag2)) %>%
    group_by(pnl, tag1) %>% mutate(n_tot = sum(n)) %>%
    mutate(prop = n/n_tot) %>%
    mutate(y = cumsum(prop)) %>%
    mutate(y = y - prop/2) %>%
    mutate(lab = glue("{number(n,accuracy=1)}{sep}({percent(prop,accuracy=1)})")) %>%
    ungroup() %>%
    mutate(tag2 = factor(tag2, levels=tags2))
tpx = tp %>% group_by(pnl, tag1) %>% summarise(n=number(sum(n))) %>% ungroup()
fillcols = brewer.pal("Accent",n=4)[c(2,1,4,3)]
tph = tp %>% filter(tag1==tag2)
tp = tp %>% mutate(tag1 = fct_rev(tag1))
p = ggplot(tp) +
    geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='fill', width=.8) +
    geom_tile(data=tph, aes(x=tag1,y=y,height=prop),width=.8,color='red',size=.5,fill=NA) +
    geom_text(aes(x=tag1,y=y,label=lab),size=2.5,lineheight=.8) +
    geom_text(data=tpx, aes(tag1,1.01,label=n), color='black',size=3, vjust=0) +
    scale_x_discrete(name='Observation', expand=expansion(mult=c(.25,.25))) +
    scale_y_continuous(name='Number/Proportion of Genes', expand=expansion(mult=c(.0,.05))) +
    scale_fill_manual(name="Prediction", values=fillcols) +
    facet_wrap(~pnl, nrow=2) +
    otheme(legend.pos='top.center.out',legend.dir='v',
           legend.title=T, legend.vjust=-.5, panel.spacing=.2,
           legend.spacing.x=.05, legend.spacing.y=.05,
           xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
           margin = c(2, .3, .3, .3)) +
    guides(fill=guide_legend(nrow=2,title.hjust=.5))
#}}}
}
#
#pa = plot_dde_bar(tp1) + theme(axis.title.x=element_blank()) + o_margin(.2,.2,0,.2)
#pb = plot_dde_bar(tp2) + o_margin(0,.2,.2,.2)
#p = ggarrange(pa, pb, nrow=2, ncol=1, labels = '', heights=c(1,1))
p = plot_dde_bar(tp)
#}}}
fo = glue("{dirw}/36.dde.bar.pdf")
ggsave(p, file=fo, width=4, height=7)
fo = glue("{dirw}/36.dde.bar.rds")
saveRDS(p, fo)

#{{{ dde cis/trans acc bar plot - f6c
#{{{ prepare
stmap = c("A=B="=0,"A+B+"=1,"A-B-"=1,"A+B="=2,"A-B="=2,"A=B+"=3,'A=B-'=3)
stm = c('A=B='='Neither respond',"A*B="="Only Mo17/W22 respond","A=B*"="Only B73 respond","A*B*"='Both respond')
tags = c('Both respond','Neither respond','Correct','Incorrect')
regs = c("cis","cis+trans",'conserved','unexpected','trans')
tp = to2 %>% filter(tgt=='B73') %>%
  filter((drc=='up' & st %in% c('dA+B=','dA=B+')) |
         (drc=='down' & st %in% c('dA-B=', 'dA=B-'))) %>%
  mutate(st = as.character(st)) %>%
  filter(!is.na(reg)) %>%
  mutate(st = str_replace(st, 'd', '')) %>%
  mutate(st = str_replace_all(st, '[\\+\\-]', '*')) %>%
  mutate(pred = str_replace_all(pred, '[\\+\\-]', '*')) %>%
  mutate(st = stm[st]) %>%
  mutate(tag = stm[pred]) %>%
  mutate(tag = ifelse(st==tag, 'Correct',
                      ifelse(pred %in% c("A*B*",'A=B='),tag,'Incorrect'))) %>%
  mutate(pnl = glue("{train} model")) %>%
  mutate(reg = factor(reg, levels=regs)) %>%
  mutate(tag = factor(tag, levels=tags)) %>%
  group_by(pnl,tag,reg) %>% summarise(n=sum(n)) %>% ungroup() %>%
  rename(tag1=tag, tag2=reg)
#}}}
plot_cis_bar <- function(ti) {
#{{{
fillcols = brewer.pal("Accent",n=4)[c(2,1,4,3)]
fillcols = brewer.pal("Set3",n=5)[c(1,2,3,5,4)]
tags1 = levels(ti$tag1)
tags2 = levels(ti$tag2)
sep = ifelse(T, " ", "\n")
col1 = 'black'; col2 = 'black'
tp = ti %>%
    arrange(tag1, desc(tag2)) %>%
    group_by(pnl, tag1) %>% mutate(n_tot = sum(n)) %>%
    mutate(prop = n/n_tot) %>%
    mutate(y = cumsum(prop)) %>%
    mutate(y = y - prop/2) %>%
    mutate(lab = glue("{number(n,accuracy=1)}{sep}({percent(prop,accuracy=1)})")) %>%
    ungroup() %>%
    mutate(tag2 = factor(tag2, levels=tags2))
tpx = tp %>% group_by(pnl, tag1) %>% summarise(n=number(sum(n))) %>% ungroup()
p = ggplot(tp) +
    geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='fill', width=.8) +
    geom_text(aes(x=tag1,y=y,label=lab),size=2.5,lineheight=.8) +
    geom_text(data=tpx, aes(tag1,1.01,label=n), color='black',size=3, vjust=0) +
    scale_x_discrete(name='Prediction of Variable Response Genes', expand=expansion(mult=c(.2,.2))) +
    scale_y_continuous(name='Number/Proportion of Genes', expand=expansion(mult=c(0,.05))) +
    scale_fill_manual(name="Mode", values=fillcols) +
    facet_wrap(~pnl, nrow=2) +
    otheme(legend.pos='top.center.out',legend.dir='h',
           legend.title=T, legend.vjust=-.7, panel.spacing=.2,
           legend.spacing.x=.05, legend.spacing.y=.05,
           xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
           margin = c(1.5, .3, .3, .3)) +
    guides(fill=guide_legend(nrow=2))
#}}}
}
#
p = plot_cis_bar(tp)
#}}}
fo = glue("{dirw}/38.dde.cis.pdf")
ggsave(p, file=fo, width=4, height=7)
fo = glue("{dirw}/38.dde.cis.rds")
saveRDS(p, fo)

#{{{ dde acc plot
tp = to2 %>% rename(score=acc) %>%
  mutate(bat = factor(bat, levels=bats)) %>%
  mutate(qry = as_factor(qry)) %>%
  mutate(tgt = as_factor(tgt)) %>%
  mutate(time = as_factor(time)) %>%
  mutate(note = as_factor(note)) %>%
  mutate(time_st = fct_cross(st, time, sep=':')) %>%
  mutate(bat_note_gt = str_c(bat, note, qry, tgt, sep=":")) %>%
  arrange(bat,note,qry,tgt,time,st) %>%
  mutate(bat_note_gt = as_factor(bat_note_gt))
tp = tp %>% mutate(bat_note_gt = factor(bat_note_gt, levels=rev(levels(tp$bat_note_gt))))
tpr1 = tibble(xmin=c(2,5), xmax=c(3,6), ymin=c(7,1), ymax=c(12,6))
tpr2 = tpr1 %>% mutate(xmin=xmin+18, xmax=xmax+18)
tpr = tpr1 %>% bind_rows(tpr2)
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=time_st,y=bat_note_gt)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    #geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(.05,.05)) +
    scale_fill_gradientn(name='Accuracy',colors=brewer_pal(type='seq',palette="Greens")(9)) +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
           margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    annotate(geom='rect',xmin=tpr$xmin-.5,xmax=tpr$xmax+.5,ymin=tpr$ymin-.5,ymax=tpr$ymax+.5, alpha=0, color='red', size=.4) +
    guides(color='none')
fo = glue("{dirw}/36.dde.acc.pdf", dirw)
p %>% ggexport(filename=fo, width=14, height=5)
#}}}
#}}}

#{{{ #[old] compare auroc in same subsets of genes - old st3
xm = tc %>% filter(train=='B') %>% select(cid,ts) %>% unnest(ts) %>%
    separate(gid, c('gt','gid'), sep='_') %>% group_by(cid) %>% nest() %>%
    rename(ts = data) %>% ungroup()
x1 = tc %>% filter(scope=='B', !str_detect(note, '^all')) %>% distinct(tid,cid)
x2 = x1 %>% mutate(tid = c(rep(c('t01','t04','t07'), each=2), 't10')) %>%
    bind_rows(x1) %>% inner_join(xm, by='cid') %>%
    arrange(cid)
xp = pd %>% filter(tag %in% c("b1"), gt=='B73') %>%
    mutate(pred=factor(pred,levels=c(1,0))) %>%
    select(tid, gid, pred, prob)
x3 = x2 %>% unnest(ts) %>% inner_join(xp, by=c('tid','gid')) %>%
    group_by(tid, cid) %>% nest() %>% ungroup()

perm=15
x4 = x3 %>% crossing(i = 1:perm) %>%
    mutate(x = map2(data, i, eval_pred, down_sample=T)) %>%
    select(-i) %>% unnest(x) %>% filter(metric=='roc_auc')

to1 = tc %>% filter(scope=='B', !str_detect(note, '^all')) %>%
    distinct(cid,cond,note,ng,ng_c) %>%
    mutate(cond=glue("{cond} {note} (pos:{ng} neg:{ng_c})")) %>%
    select(cid,cond)
to2 = tc %>% mutate(mdl = glue("{cond} {note}")) %>%
    select(tid,mdl)
to0 = x4 %>% select(tid, cid, estimate) %>%
    group_by(tid,cid) %>%
    summarise(vs=list(estimate), mean=mean(estimate), sd=sd(estimate)) %>%
    ungroup() %>%
    inner_join(to1, by='cid') %>%
    inner_join(to2, by='tid') %>%
    mutate(grp = ifelse(str_detect(mdl, "all"), "all model", "subset model")) %>%
    mutate(txt=glue("{mdl}\n{number(mean,accuracy=.001)} ({number(sd,accuracy=.001)})"))
to1 = to0 %>% select(cond,grp,txt) %>% spread(grp,txt)
to2 = to0 %>% select(cond,grp,vs) %>% spread(grp,vs) %>%
    mutate(sig = map2_chr(`all model`, `subset model`, eval_signif))
to = to1 %>% inner_join(to2 %>% select(cond, sig), by='cond')

x = to %>%
    mutate_all(linebreak, align = "c") %>%
    kbl(format='latex', escape=F, booktabs=T, linesep="",
        align='l', format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c('striped', "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirf, 'st3.rds')
saveRDS(x, file=fo)
#}}}
#{{{ #[old] accuracy in predicting U/D/N at 1h/25h time points
get_acc <- function(ti) {
  #{{{
  ti %>% group_by(st) %>%
    summarise(n.tot = n(), n0=sum(pred==0), n1=sum(pred==1)) %>% ungroup() %>%
    mutate(nc = ifelse(st %in% c("N",'zero'), n0, n1)) %>%
    mutate(acc = nc/n.tot,
           lab = glue("{number(acc,accuracy=.01)}\n{nc} / {n.tot}")) %>%
    select(st, acc, lab)
  #}}}
}
to1 = pd %>%
  inner_join(td1, by=c('cond','gt','gid')) %>%
  group_by(train,cid,cond,note,gt,time) %>% nest() %>% rename(ti=data) %>%
  mutate(acc = map(ti, get_acc)) %>%
  select(-ti) %>% unnest(acc)

#{{{ de acc plot
tp = to1 %>% rename(score=acc) %>%
  mutate(cid = as_factor(cid)) %>%
  mutate(gt = as_factor(gt)) %>%
  mutate(st = as_factor(st)) %>%
  mutate(time = as_factor(time)) %>%
  mutate(time_st = fct_cross(time, st, sep='_')) %>%
  mutate(cid_gt = fct_cross(cid, gt, sep=":"))
tp = tp %>% mutate(cid_gt = factor(cid_gt, levels=rev(levels(tp$cid_gt))))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=time_st,y=cid_gt)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    #geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='Accuracy',colors=brewer_pal(type='seq',palette="Greens")(9)) +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(~train, ncol=2) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
           margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color='none')
fo = glue("{dirw}/34.de.acc.pdf")
ggsave(p, file=fo, width=10, height=8)
#}}}
#}}}
#}}}

######## BELOW ARE DEPRECATED ########
#{{{ write out for cnn
x = r %>% filter(bat_mid=='cold_up:m16', bin=='+/-2k', nfea == 'top100')
x1 = x$to[[1]]
x2 = x$to[[2]] %>%
  #rename_at(starts_with("g"), function(x) {paste0("umr_", x)})
  rename_with(.fn = ~ paste0("umr_", .x))
xo = x1 %>% bind_cols(x2[,-1])

fo = file.path(dirw, 'z.tsv')
write_tsv(xo, fo)
#}}}

#{{{ [obsolete] evaluate model performance in M/W using downsampled datasets
eval_model <- function(mdl, ti) {
  #{{{
  pb$tick()
  features = names(mdl$fit$variable.importance)
  new.cols = features[!features %in% colnames(ti)]
  if (length(new.cols) > 0) {
    ti2 = matrix(0, nrow=nrow(ti), ncol=length(new.cols), dimnames=list(NULL, new.cols))
    ti = ti %>% bind_cols(as.data.frame(ti2))
  }
  x = predict(mdl, ti, type='class') %>% rename(pred=1) %>%
    bind_cols(predict(mdl, ti, type='prob')) %>%
    rename(prob = .pred_1) %>%
    bind_cols(truth=ti$status) %>%
    mutate(truth = factor(truth, levels=c(1,0)))
  metrics6 <- metric_set(sens,spec,precision,accuracy, f_meas, roc_auc, pr_auc)
  metrics6(x, truth=truth, estimate=pred, prob) %>%
    select(metric=.metric, estimate=.estimate)
  #}}}
}

gt2 = 'Zmays_Mo17'
gt2 = 'Zmays_W22'
tt = readRDS(glue("{dirw}/{gt2}/01.rds"))
to0 = tb %>%
  inner_join(tt, by=c('bat','note','bin','epi','nfea','mod')) %>%
  mutate(to = map(to, downsample))
pb = progress_bar$new(total=nrow(to0))
to = to0 %>%
  mutate(m2 = map2(fit, to, eval_model)) %>%
  select(bat,note,bin,epi,nfea,mod,m2) %>%
  unnest(m2) %>% spread(metric, estimate)

metric = 'roc_auc'
metric = 'f_meas'
fo = sprintf("%s/09.%s.%s.pdf", dirw, metric, gt2)
p = plot_model_perf(to, metric=metric)
p %>% ggexport(filename=fo, width=6, height=6)
#}}}
#{{{ prepare input for B/M/W, build B73 ML model [long time]
gt = 'Zmays_B73'
gt = 'Zmays_Mo17'
gt = 'Zmays_W22'
diro = file.path(dirw, gt)

#{{{ ts
fl = file.path(dird, '31_promoter', gt, '15.module.rds')
tl = readRDS(fl)
#
ts = tl %>% filter(ctag=='c1') %>% select(-ctag) %>%
  select(bat_mid, note, hit=tg, control=tg_c) %>%
  gather(status, tg, -bat_mid, -note) %>% unnest(tg) %>%
  mutate(status = factor(status, levels=c('hit', 'control'))) %>%
  mutate(status = 2-as.numeric(status)) %>%
  mutate(status = factor(status, levels=c(1,0))) %>%
  select(bat_mid,note,gid,status) %>%
  group_by(bat_mid, note) %>% nest() %>% rename(ts = data) %>% ungroup() %>%
  separate(bat_mid, c('bat','mid0'), sep=":") %>% select(-mid0)
#}}}

fi = file.path(dird, '32_kmer', gt, '05.mtf.loc.rds')
mloc = readRDS(fi)

bins = c("-500","+500","+/-500","-2k",'+2k','+/-2k')
epis = c('raw','umr')
nfeas = c('top30', 'top50', 'top100')
mods = c('zoops', 'anr')
to0 = mloc %>% filter(mtf >= 20) %>%
  filter(note == 'all') %>%
  inner_join(ts, by=c('bat','note')) %>%
  mutate(tl = map2(tl, ts, filter_kmer_loc)) %>%
  crossing(bin = bins, epi = epis, nfea = nfeas, mod = mods) %>%
  mutate(nfea = factor(nfea, levels=nfeas)) %>%
  mutate(bin = factor(bin, levels=bins)) %>%
  mutate(epi = factor(epi, levels=epis)) %>%
  mutate(mod = factor(mod, levels=mods)) %>%
  arrange(bin, epi, nfea, mod, bat, note)
pb = progress_bar$new(total=nrow(to0))
to = to0 %>%
  mutate(to = pmap(list(tl,bin,epi,nfea,mod,ts), prepare_ti)) %>%
  mutate(did = sprintf("d%03d", 1:n())) %>%
  select(did, bat, note, bin, epi, nfea, mod, to)

fo = file.path(diro, "01.rds")
saveRDS(to, fo)
fo = file.path(diro, "01.tsv")
write_tsv(to %>% select(-to), fo)

to %>% mutate(fo = sprintf("%s/01_ml_input/%s.tsv", diro, did)) %>%
  mutate(y = map2(to, fo, write_tsv))
#}}}
#{{{ test python script from shiulab/ML_workshop
tt = r %>%
  filter(bat_mid=='cold_up:m21',bin=='2k',epi=='umr') %>%
  pluck('ti',1) %>% mutate(gid=1:n()) %>%
  select(gid,Class=status,everything()) %>%
  mutate(across(starts_with("umr_"), str_replace, "y", '1')) %>%
  mutate(across(starts_with("umr_"), str_replace, 'n', '0'))
fo = file.path(dirw, 'test.tsv')
write_tsv(tt, fo)
#}}}
#{{{ test run
r1 = r %>%
  slice(1) %>% crossing(perm = 1:3) %>%
  mutate(ti2 = map2(ti, perm, downsample))
  #filter(bat_mid=='cold_up:m21',bin=='2k',epi=='umr+pos')
  #filter(bat_mid == 'cold_up:m21') %>% slice(13)
ti = r1 %>% pluck("ti2", 1)
x = r1 %>%  mutate(y = map(ti2, ml0, alg='rf',seed=5)) %>%
  select(-ti, -ti2) %>% unnest(y)
x = r1 %>%  mutate(y = map(ti2, ml1, alg='rf',fold=2, seed=5)) %>%
  unnest(y)

x %>% select(bat_mid,bin,epi,perm,metric) %>% unnest(metric) %>%
  spread(metric, estimate) %>% print(n=40)
x %>% select(bat_mid,bin,epi,metricB) %>% unnest(metricB) %>%
  spread(metric, estimate) %>% print(n=40)
x %>% select(bat_mid,bin,epi,pred) %>% unnest(pred) %>%
  count(bat_mid,bin,epi,pred, truth) %>% print(n=72)
#}}}
#{{{ VIP score
to = x %>% filter(bat_mid=='heat_up:m18', epi != 'pos') %>%
  arrange(bat_mid,bin,epi) %>%
  mutate(bin_epi = str_c(bin, epi, sep=':')) %>%
  mutate(bin_epi = factor(bin_epi, levels=tpl$bin_epi)) %>%
  filter(bin_epi == '+/-2k:umr') %>%
  select(bat_mid, bin, epi, bin_epi, vis) %>% unnest(vis) %>%
  separate(Variable, c('epi2','var'), sep="_", extra='merge')
tos = to %>%  group_by(bat_mid, bin_epi, var) %>%
  summarise(imp_avg=mean(Importance), imp_med=median(Importance),
    imp_sd=sd(Importance), imp_q25=quantile(Importance,.25),
    imp_q75 = quantile(Importance,.75)) %>%
  ungroup()

p = ggplot(to) +
    geom_violin(aes(x=var, y=Importance), trim=T, alpha=.8) +
    #geom_jitter(aes(x=bin_epi, y=estimate), color='darkgray') +
    geom_pointrange(data=tos, aes(x=var, y=imp_avg,ymin=imp_avg-imp_sd, ymax=imp_avg+imp_sd)) +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_continuous(expand=expansion(mult=c(.01,.01))) +
    coord_flip() +
    otheme(legend.pos='none', strip.style='white',
           xtext=T, xtick=T, ytext=T, ytick=T, ygrid=T) +
    theme(axis.text.x = element_text(angle = 20, hjust=.8, vjust=1.1)) +
    guides(color='none')
fo = file.path(dirw, 'vip2.pdf')
ggsave(p, filename=fo, width=8, height=15)
#}}}
#{{{ # obtain maximized F1 score
tp = x %>% pluck('pred', 1) %>% mutate(truth=ifelse(truth=='hit',1,0)) %>%
  mutate(truth=factor(truth, levels=c(1,0)))

f_meas2 <- function(thresh, pPred, truth) {
  #{{{
  tibble(pPred=pPred) %>%
    mutate(pred = ifelse(pPred >= thresh, 1, 0)) %>%
    mutate(pred = factor(pred, levels=c(1,0))) %>%
    f_meas(truth=truth, estimate=pred) %>%
    replace_na(list(.estimate=0)) %>%
    pluck('.estimate',1)
  #}}}
}
f_meas2(.5, tp$.pred_hit, tp$truth)

tp = x %>% select(bin,epi,pred) %>% unnest(pred) %>%
  mutate(truth=ifelse(truth=='hit',1,0)) %>%
  mutate(truth=factor(truth, levels=c(1,0))) %>%
  group_by(bin,epi) %>%
  summarise(pPred=list(.pred_hit), truth=list(truth)) %>%
  ungroup()

tp %>% mutate(o = map2(pPred, truth, maximize_f1)) %>% unnest(o)

maximize_f1 <- function(pPred, truth) {
  #{{{
  o = optim(.5, f_meas2, pPred=pPred, truth=truth,
            method='L-BFGS-B', lower=.01, upper=.99, control=list(fnscale=-1))
  tibble(thresh=o$par[[1]], f1=o$value[[1]])
  #}}}
}
#}}}

