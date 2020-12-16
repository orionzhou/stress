source('functions.R')
require(progress)
dirw = glue('{dird}/41_ml')
dirn = glue('{dird}/41_ml/00_nf')
#{{{ params
bins = c(
         "TSS:-500","TSS:+500","TSS:-/+500",
         "TSS:-2k",'TSS:+2k','TSS:-/+2k',
         "TTS:-500","TTS:+500","TTS:-/+500",
         "TTS:-2k",'TTS:+2k','TTS:-/+2k',
         "TSS:-/+2k,TTS:-500",
         "TSS:-/+2k,TTS:-1k"
)
epis = c('raw','umr')
nfeas = c('top30', 'top50', 'top100', 'top200')
mods = c('zoops', 'anr')
bins1 = c("TSS:-/+2k")
epis1 = c('umr')
nfeas1 = c('top100')
mods1 = c('anr')
#}}}

#{{{ prepare training input: gene lists + kmer lists
fi = glue("{dird}/25_dreme/05.best.mtfs.rds")
r5 = readRDS(fi)
tk = r5$tk %>% filter(n_mtf >= 10)

fk = glue("{dirw}/05.tk.tsv")
write_tsv(tk %>% select(cid, n_mtf), fk)
diro = glue('{dirn}/03_motif_lists')
if(!dir.exists(diro)) dir.create(diro)
tk %>% mutate(mtfs = map(mtfs, 'mtf')) %>%
    mutate(fo = glue("{diro}/{cid}.meme")) %>%
    mutate(l = map2(mtfs, fo, write_meme, overwrite=T))

tc1 = r5$tc %>% mutate(train='B') %>% select(train,everything())
#{{{ BMW training
f_cfg = glue('{dird}/25_dreme/config.xlsx')
cfg = read_xlsx(f_cfg) %>% fill(tag, .direction='down') %>% select(tag,ocid,cid)
#{{{ read in
diri = glue("{dird}/25_dreme/00_nf")
tag = 'degA'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc_1 = r$tc %>% mutate(tag = tag)
#
tag = 'dmodA'
fi = glue('{diri}/{tag}/08.tc.tl.rds')
r = readRDS(fi)
tc_2 = r$tc %>% mutate(tag = tag)
#}}}
tc2 = tc_1 %>% bind_rows(tc_2) %>% rename(ocid=cid) %>%
    inner_join(cfg, by=c('tag','ocid')) %>% arrange(cid) %>%
    select(cid,cond,note, ng,ng_c,gids,gids_c,ts) %>%
    mutate(train = 'BMW') %>% select(train, everything())
#}}}
fi = glue("{dird}/17_cluster/50_modules/var1.rds")
tc3 = readRDS(fi) %>%
    mutate(train = 'var') %>% select(train, everything())

trains = c('B','BMW','var')
tc = rbind(tc1,tc2,tc3) %>% filter(cid %in% tk$cid) %>%
    mutate(train = factor(train, levels=trains)) %>%
    arrange(train, cid) %>%
    mutate(tid = str_c('t', str_pad(1:n(), width=2,pad='0'))) %>%
    select(tid, everything())

r6 = list(tk=tk, tc=tc)
fo = glue("{dirw}/06.tk.tc.rds")
saveRDS(r6,fo)

diro = glue('{dirn}/05_gene_lists')
if(!dir.exists(diro)) dir.create(diro)
tc %>% mutate(fo = glue("{diro}/{tid}.tsv")) %>% mutate(j=map2(ts, fo, write_tsv))
#}}}

#{{{ motif meta plots - sf09
fi = glue("{dirw}/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk
#{{{ function
fimo_locate <- function(mid,fm,fg) {
    #{{{
    cmd = glue("fimo.py prepare {fg} {fm} tmp.bed --epi umr --nfea {mid} --fmt long")
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
tc1 = tc %>% filter(train=='B') %>% select(cid,cond,note)
rb = tk %>% inner_join(tc1, by='cid') %>%
    unnest(mtfs) %>%
    mutate(rate=pos/ng, rate.c=neg/ng_c) %>%
    mutate(fnote = glue("{cid}_{i} {fid} ({fname}): {opt} {bin} {epi}")) %>%
    select(cid,cond,note,i,pval,rate,rate.c,mid,fid,fname,fnote)
rb1 = rb %>% filter(rate.c < .2)
#}}}

#{{{ explore
plot_tss_meta <- function(mid,fm,cid,cond,note,fid,fnote,fg,md,dirw) {
    #{{{
    cmd = glue("fimo.py prepare_ml {fg} {fm} tmp.bed --epi umr --nfea {mid} --fmt long")
    system(cmd)
    ti = read_tsv('tmp.bed', col_names=c("gid",'start','end','mid'))
    system("rm tmp.bed")
    #
    mds = md %>% count(cid,cond,note) %>% rename(nt=n) %>%
        mutate(mod = glue("{cond}: {note} ({nt})")) %>%
        arrange(cid) %>% mutate(mod = as_factor(mod))
    itvs = c(seq(0,4000,by=200), seq(4050,8050,by=200))
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
        guides(fill=F)
    #}}}
    fo = glue("{dirw}/53_metaplots/{cid}_{fid}.pdf")
    ggsave(p, file=fo, width = 10, height = 5)
    #}}}
}
rb2 = rb1 %>% filter(cid=='c10') %>%
    mutate(fm = glue("{dirn}/03_motif_lists/{cid}.meme")) %>%
    mutate(x = pmap(list(mid,fm,cid,cond,note,fid,fnote), plot_tss_meta,
                    fg=fg, md=md, dirw=dirw))
#}}}

tp1 = rb1 %>% filter(cid=='c32', fid %in% c('f0013','f0736',"f0277")) %>%
    mutate(x = map(kmers, kmer_locate, fg=fg)) %>% unnest(x)

    #{{{ sf09
    mds = md %>% filter(cond=='heat', cid %in% c("c30",'c32','c49')) %>%
        count(cid,cond,note) %>% rename(nt=n) %>%
        mutate(mod = glue("{cond}: {note} ({nt})")) %>%
        arrange(cid) %>% mutate(mod = as_factor(mod)) %>%
        mutate(modtype = ifelse(str_detect(cid,'9$'), 'bg','a'))
    itvs = c(seq(0,4000,by=200), seq(4050,8050,by=200))
    tps = tp1 %>% distinct(fid,fname)
    tp = tp1 %>%
        mutate(pos = (start+end)/2) %>%
        mutate(bin = cut(pos, breaks=itvs, include.lowest=T)) %>%
        mutate(bin = as.numeric(bin)) %>%
        distinct(fname,sid,bin) %>%
        mutate(fname = as_factor(fname)) %>%
        rename(gid = sid) %>% inner_join(md, by='gid') %>%
        inner_join(mds,  by='cid') %>%
        count(fname,mod,modtype,nt, bin) %>% rename(nh=n) %>%
        mutate(prop = nh/nt)
    #
    tpx = tibble(x=c(.5,10.5,20.5,31.5,41.5),lab=c('-2kb','TSS','+2kb/-2kb','TTS','+2kb'))
    #{{{
    cols6 = c(pal_npg()(2), 'black')
    p = ggplot(tp, aes(x=bin,y=prop, color=mod, shape=mod, linetype=modtype)) +
        geom_line(aes(), size=.5, na.rm = F) +
        geom_point(aes(), size=1, na.rm = F) +
        scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
        scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
        scale_color_manual(name='', values=cols6) +
        scale_shape(name='') +
        scale_linetype(name='') +
        facet_wrap(~fname, scale='free_y', ncol=2) +
        otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
               panel.spacing=.2, margin = c(.3,.3,.3,.3),
               xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
        #theme(plot.title=element_text(hjust=.5, size=10)) +
        guides(fill=F,linetype=F)
    #}}}
    #}}}
fo = glue("{dirw}/54.metaplot.pdf")
ggsave(p, file=fo, width = 6, height = 4)
fo = glue("{dirf}/sf09.pdf")
ggsave(p, file=fo, width = 6, height = 4)
#}}}

#{{{ create ML training config
fi = glue("{dirw}/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk

#{{{ b1
ctag = 'b1'
tc1 = tc %>% filter(train == 'B')
to1 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins, epi = epis1, nfea = nfeas1, mod = mods1)
to2 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins1, epi = epis, nfea = nfeas, mod = mods)
to = to1 %>% bind_rows(to2) %>%
    distinct(tid, cid, bin, epi, nfea, mod) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(tid, cid, bin, epi, nfea, mod) %>%
    mutate(did = sprintf("d%03d", 1:n())) %>%
    select(did, tid, cid, bin, epi, nfea, mod)

fo = glue("{dirn}/{ctag}.cfg.rds")
saveRDS(to, fo)
fo = glue("{dirn}/{ctag}.cfg.tsv")
write_tsv(to, fo)
#}}}

#{{{ b2
ctag = 'b2'
tc1 = tc %>% filter(train == 'BMW')
to1 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins, epi = epis1, nfea = nfeas1, mod = mods1)
to2 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins1, epi = epis, nfea = nfeas, mod = mods)
to = to1 %>% bind_rows(to2) %>%
    distinct(tid, cid, bin, epi, nfea, mod) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(tid, cid, bin, epi, nfea, mod) %>%
    mutate(did = sprintf("d%03d", 1:n())) %>%
    select(did, tid, cid, bin, epi, nfea, mod)

fo = glue("{dirn}/{ctag}.cfg.rds")
saveRDS(to, fo)
fo = glue("{dirn}/{ctag}.cfg.tsv")
write_tsv(to, fo)
#}}}

#{{{ b3
ctag = 'b3'
tc1 = tc %>% filter(train == 'var')
to1 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins, epi = epis1, nfea = nfeas1, mod = mods1)
to2 = tc1 %>% select(tid,cid) %>%
  crossing(bin = bins1, epi = epis, nfea = nfeas, mod = mods)
to = to1 %>% bind_rows(to2) %>%
    distinct(tid, cid, bin, epi, nfea, mod) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(tid, cid, bin, epi, nfea, mod) %>%
    mutate(did = sprintf("d%03d", 1:n())) %>%
    select(did, tid, cid, bin, epi, nfea, mod)

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

#{{{ process trained models
collect_model_metrics <- function(tag, dirw) {
  #{{{
  fi = glue("{dirw}/00_nf/{tag}.cfg.rds")
  th = readRDS(fi)
  fi = glue("{dirw}/00_nf/{tag}.rds")
  ml = readRDS(fi)
  ml1 = ml %>% select(did=sid, perm, metric) %>%
    unnest(metric) %>% spread(metric, estimate) %>%
    dplyr::rename(f1=f_meas, auroc=roc_auc, auprc=pr_auc) %>%
    inner_join(th, by='did')
  ml1
  #}}}
}
get_best_model <- function(tag, dirw, bin='TSS:-/+2k', epi='umr',
                           nfea='top100', mod='anr') {
  #{{{
  fi = glue("{dirw}/00_nf/{tag}.cfg.rds")
  th = readRDS(fi)
  fi = glue("{dirw}/00_nf/{tag}.rds")
  ml = readRDS(fi)
  th %>% filter(bin==!!bin,epi==!!epi, nfea==!!nfea, mod==!!mod) %>%
      inner_join(ml, by=c('did'='sid')) %>%
      filter(!is.na(fit)) %>% select(did,tid,cid,fit,metric)
  #}}}
}
tags = c("b1",'b2')
tags = c("b3")

tm = tibble(tag = tags) %>%
  mutate(x = map(tag, collect_model_metrics, dirw=dirw))
fm = glue('{dirw}/11.model.metrics.rds')
saveRDS(tm, fm)

tb = tibble(tag = tags) %>%
  mutate(x = map(tag, get_best_model, dirw=dirw)) %>%
  unnest(x)
fb = glue('{dirw}/12.best.models.rds')
saveRDS(tb, fb)
fb = glue('{dirw}/12.best.models.tsv')
write_tsv(tb %>% select(tag,tid,cid), fb)
#}}}

#{{{ eval model performance in B - f4
fm = glue('{dirw}/11.model.metrics.rds')
tm = readRDS(fm)

#{{{ prepare tc1
bins = c(
         "TSS:-500","TSS:+500","TSS:-/+500",
         "TSS:-2k",'TSS:+2k','TSS:-/+2k',
         "TTS:-500","TTS:+500","TTS:-/+500",
         "TTS:-2k",'TTS:+2k','TTS:-/+2k',
         "TSS:-/+2k,TTS:-500",
         "TSS:-/+2k,TTS:-1k"
)
epis = c('all','umr')
nfeas = c('top30', 'top50', 'top100', 'top200')
mods = c('binary', 'quantitative')
epi_map = c('raw'='all','umr'='umr')
mod_map = c('zoops'='binary','anr'='quantitative')
tc1 = tc %>% filter(train=='B') %>%
    select(cid,cond,note) %>%
    crossing(bin = bins, epi = epis, nfea = nfeas, mod = mods) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(cid, bin, epi, nfea, mod) %>%
    mutate(cond_note = as_factor(str_c(cond, note, sep=': '))) %>%
    mutate(cond_note_nfea_mod = as_factor(str_c(cond,note,nfea,mod,sep=': '))) %>%
    mutate(bin_epi = as_factor(str_c(bin, epi, sep=': '))) %>%
    mutate(cond_mod = as_factor(str_c(cond, mod, sep=': '))) %>%
    mutate(cond_epi = as_factor(str_c(cond, epi, sep=': '))) %>%
    mutate(bin_epi_nfea = as_factor(str_c(bin, epi, nfea, sep=': '))) %>%
    mutate(bin_epi_mod = as_factor(str_c(bin, epi, mod, sep=': '))) %>%
    mutate(epi_nfea_mod = as_factor(str_c(epi, nfea, mod, sep=': '))) %>%
    mutate(bin_nfea_mod = as_factor(str_c(bin, nfea, mod, sep=': ')))
#}}}
tag = 'b1'
tm1 = tm %>% filter(tag == !!tag) %>% unnest(x) %>%
    mutate(epi=epi_map[epi], mod=mod_map[mod]) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    inner_join(tc1, by=c('cid','bin','epi','nfea','mod'))

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
      guides(color=F)
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
      otheme(legend.pos='none', legend.dir='v', legend.box='v',legend.title=F,
             legend.spacing.x=.1, legend.spacing.y=.1,
             margin=c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=F, panel.spacing=.1,
             xtick=T, ytick=T, ytitle=T, xtext=T, ytext=T) +
      theme(legend.position = c(.5,1), legend.justification = c(.5,1)) +
      theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y=element_text(size=7)) +
      guides(color=F)
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

#{{{ f4 - barplot
metric = 'auroc'
cids4 = tm1 %>% filter(str_detect(note, '^all')) %>% pull(cid)
tpa = tm1 %>% filter(cid %in% cids4, epi=='umr', nfea=='top100', mod=='quantitative')
cols14 = c(brewer.pal(4,'Blues')[2:4],brewer.pal(4,'Greens')[2:4],
    brewer.pal(4,'Purples')[2:4],brewer.pal(4,'Greys')[2:4],
    brewer.pal(4,'Reds')[2:3])
pa = plot_model_barplot(tpa, metric=metric,, x='cond_note', fill="bin",
                        pnl='epi_nfea_mod', x.rotate=F, cols=cols14) +
  theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
  theme(legend.direction='horizontal') +
  o_margin(2.5,.5,.2,.2)
#
tpb = tm1 %>% filter(bin=='TSS:-/+2k', nfea=='top100', mod=='binary')
pb = plot_model_barplot(tpb, metric=metric,, x='cond_note', fill='epi',
                        pnl='bin_nfea_mod', x.rotate=T) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
tpc = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', nfea=='top100')
pc = plot_model_barplot(tpc, metric=metric,, x='cond_note', fill="mod",
                        pnl='bin_epi_nfea', x.rotate=T, cols=pal_aaas()(4)) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
tpd = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', mod=='quantitative')
pd = plot_model_barplot(tpd, metric=metric,, x='cond_note', fill="nfea",
                        pnl='bin_epi_mod', x.rotate=T, cols=pal_simpsons()(8)) +
  o_margin(0,.2,.2,2.1) +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
#
p = ggarrange(pa, pb, pc, pd, nrow=4, ncol=1,
          labels = LETTERS[1:4], heights=c(1.5, 1, 1, 1.3))
#}}}
fo = glue("{dirw}/30.auroc.{tag}.pdf")
p %>% ggexport(filename=fo, width=6, height=9)
fo = glue("{dirf}/f4.pdf")
p %>% ggexport(filename=fo, width=6, height=9)

#{{{ sf09 - barplot
metric = 'f1'
cids4 = tm1 %>% filter(str_detect(note, '^all')) %>% pull(cid)
tpa = tm1 %>% filter(cid %in% cids4, epi=='umr', nfea=='top100', mod=='quantitative')
cols14 = c(brewer.pal(4,'Blues')[2:4],brewer.pal(4,'Greens')[2:4],
    brewer.pal(4,'Purples')[2:4],brewer.pal(4,'Greys')[2:4],
    brewer.pal(4,'Reds')[2:3])
pa = plot_model_barplot(tpa, metric=metric,, x='cond_note', fill="bin",
                        pnl='epi_nfea_mod', x.rotate=F, cols=cols14) +
  theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
  theme(legend.direction='horizontal') +
  o_margin(2.5,.5,.2,.2)
#
tpb = tm1 %>% filter(bin=='TSS:-/+2k', nfea=='top100', mod=='binary')
pb = plot_model_barplot(tpb, metric=metric,, x='cond_note', fill='epi',
                        pnl='bin_nfea_mod', x.rotate=T) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
tpc = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', nfea=='top100')
pc = plot_model_barplot(tpc, metric=metric,, x='cond_note', fill="mod",
                        pnl='bin_epi_nfea', x.rotate=T, cols=pal_aaas()(4)) +
  o_margin(0,.2,.2,2.1) + no_x_axis() +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
tpd = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', mod=='quantitative')
pd = plot_model_barplot(tpd, metric=metric,, x='cond_note', fill="nfea",
                        pnl='bin_epi_mod', x.rotate=T, cols=pal_simpsons()(8)) +
  o_margin(0,.2,.2,2.1) +
  theme(legend.position = c(.3,1), legend.justification = c(0,1))
#
p = ggarrange(pa, pb, pc, pd, nrow=4, ncol=1,
          labels = LETTERS[1:4], heights=c(1.5, 1, 1, 1.3))
#}}}
fo = glue("{dirw}/30.f1.{tag}.pdf")
p %>% ggexport(filename=fo, width=6, height=9)
fo = glue("{dirf}/sf10.pdf")
p %>% ggexport(filename=fo, width=6, height=9)

#{{{ [old] heatmap: sf09
pa = plot_model_heatmap(tm1, metric='f1', pnl='bat_note') +
  theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7))
pb = plot_model_heatmap(tm1, metric='auroc', pnl='bat_note') +
  theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7)) +
  theme(plot.margin = margin(.3, .3, .3, 0, "lines")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
fo = glue("{dirf}/sf09b.pdf")
ggarrange(pa, pb, nrow=1, ncol=2,
          labels = c(), widths=c(2,1.8)) %>%
ggexport(filename=fo, width=6, height=6)

pa = plot_model_heatmap(tm2, metric='f1', pnl='bat_note')
pb = plot_model_heatmap(tm2, metric='auroc', pnl='bat_note') +
  theme(plot.margin = margin(.3, .3, .3, 0, "lines")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
fo = glue("{dirf}/sf09b.pdf")
ggarrange(pa, pb, nrow=1, ncol=2,
          labels = c(), widths=c(2,1.4)) %>%
ggexport(filename=fo, width=6, height=7)
#}}}
#}}}

#{{{ gather & save model predictions for all B/M/W genes
#{{{ write all gene list
fi = glue('{dird}/15_de/05.rds')
x =  readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12
md = deg12 %>% pluck('ds',1) %>% select(gid) %>%
    crossing(gt = gts3) %>% mutate(status = 1) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    select(gid, status)

fo = glue("{dirw}/21.all.tsv")
write_tsv(md, fo)
#}}}
# run j21 kmer.py prepare_ml

fb = glue('{dirw}/12.best.models.rds')
tb = readRDS(fb)

# write models
tb %>% mutate(fo = glue("{dirw}/23_models/{tid}.rds")) %>%
    mutate(map2(fit, fo, saveRDS))
# run j21 ml_predict

to = tb %>% select(tag, tid, cid) %>%
    mutate(fi = glue("{dirw}/24_ml_out/{tag}_{tid}.tsv")) %>%
    mutate(pred = map(fi, read_tsv)) %>% select(-fi)

fo = glue('{dirw}/25.model.pred.rds')
saveRDS(to, fo)
#}}}

#{{{ evaluate model performance for different gene categories
fp = glue('{dirw}/25.model.pred.rds')
pd = readRDS(fp) %>%
    inner_join(tc %>% select(tid,train,cid,cond,note), by=c('tid','cid')) %>%
    select(tag,train,tid,cid,cond,note,pred) %>% unnest(pred) %>%
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

#{{{ overall performance - f6a
xm = tc %>% filter(train=='BMW') %>% select(cid,ts) %>% unnest(ts) %>%
    separate(gid, c('gt','gid'), sep='_')
#
perm=10
x = pd %>% mutate(pred=factor(pred,levels=c(1,0))) %>%
    inner_join(xm, by=c('cid','gt','gid')) %>%
    select(-tag,-tid) %>%
    group_by(train,cid,cond,note,gt) %>% nest() %>% rename(pred=data) %>%
    crossing(i = 1:perm) %>%
    mutate(x = map2(pred, i, eval_pred, downsample=T)) %>%
    select(-i) %>% unnest(x) %>%
    group_by(train,cid,cond,note,gt,metric) %>%
    summarise(mean=mean(estimate), sd=sd(estimate)) %>% ungroup()

#{{{ bar plot - f6a
tp = x %>% filter(metric=='roc_auc') %>% rename(score=mean) %>%
  #mutate(gt = factor(gt, levels=rev(gts3))) %>%
  #filter(str_detect(note, 'all')) %>%
  mutate(cond_note = str_c(cond,note,sep=": ")) %>%
  mutate(train = factor(train, levels=c('BMW','B'))) %>%
  mutate(cond_note = as_factor(cond_note))
tps = tp %>% distinct(cid,cond_note)
pdg = position_dodge(width=.8)
ymax = 1
p = ggplot(tp, aes(x=train, y=score)) +
  geom_col(aes(fill=gt), position=pdg, width=.7, alpha=.6) +
  geom_errorbar(aes(ymin=score-sd, ymax=score+sd, color=gt), width=.2, size=.3, position=pdg) +
  scale_x_discrete(expand=expansion(mult=c(.05,.1))) +
  scale_y_continuous(name='AUROC', expand=expansion(mult=c(0,0))) +
  #coord_cartesian(ylim=c(.5,ymax)) +
  coord_flip(ylim= c(.5, ymax)) +
  scale_fill_manual(name='', values=pal_npg()(5)) +
  scale_color_manual(name='', values=pal_npg()(5)) +
  facet_wrap(cond_note ~ ., scale='free_y', ncol=1, strip.position='top') +
  otheme(legend.pos='top.center.out',legend.dir='h',legend.box='h',
         legend.title=F, legend.vjust=-1.5,
         margin = c(.2,.2,.2,.2), xgrid=T, panel.border=F,
         strip.style='white', strip.compact=F, panel.spacing=.1,
         xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T) +
  #theme(axis.text.x = element_text(angle=0,size=8, hjust=1,vjust=1)) +
  theme(strip.text.y = element_text(angle=0,size=8))
#}}}
fo = glue("{dirw}/32.auroc.pdf")
ggsave(p, file=fo, width=4, height=8)
fo = glue("{dirf}/f6a.pdf")
#ggsave(p, file=fo, width=4, height=8)
saveRDS(p, glue("{dirf}/f6a.rds"))
#}}}

#{{{ accuracy in predicting variable genotype response - f6
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
pd3 = pd %>% select(-tag) %>% filter(str_detect(note, "all")) %>%
    mutate(drc = ifelse(str_detect(note, "all up"), "up", "down")) %>%
    mutate(pred=ifelse(drc=='down' & pred==1, -1, pred))
to2 = td2 %>%
  left_join(ddeg2, by=c('cond','time','qry','tgt','gid','st')) %>%
  #filter(!is.na(reg) & reg == 'cis') %>%
  inner_join(pd3, by=c('cond','qry'='gt','gid')) %>% rename(qPred=pred) %>%
  inner_join(pd3, by=c('train','cid','cond','drc','note','tgt'='gt','gid')) %>%
  rename(tPred=pred) %>%
  #group_by(bat,cond,drc,note,qry,tgt,time) %>% nest() %>% ungroup() %>%
  mutate(pred1 = smap2[as.character(qPred)]) %>%
  mutate(pred2 = smap2[as.character(tPred)]) %>%
  mutate(pred = as.character(glue("A{pred1}B{pred2}"))) %>%
  count(train,cid,cond,drc,note,qry,tgt,time, st, pred,reg)
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
    geom_tile(data=tph, aes(x=tag1,y=y,height=prop),width=.8,color='red',size=1,fill=NA) +
    geom_text(aes(x=tag1,y=y,label=lab),size=2.5,lineheight=.8) +
    geom_text(data=tpx, aes(tag1,1.01,label=n), color='black',size=3, vjust=0) +
    scale_x_discrete(name='Observation', expand=expansion(mult=c(.25,.25))) +
    scale_y_continuous(name='Number/Proportion of Genes', expand=expansion(mult=c(.02,.05))) +
    scale_fill_manual(name="Prediction", values=fillcols) +
    facet_wrap(~pnl, nrow=2) +
    otheme(legend.pos='top.center.out',legend.dir='v',
           legend.title=T, legend.vjust=-.5, panel.spacing=.2,
           legend.spacing.x=.05, legend.spacing.y=.05,
           xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
           margin = c(1.5, .3, .3, .3)) +
    guides(fill=guide_legend(nrow=2,title.hjust=.5))
#}}}
}

#pa = plot_dde_bar(tp1) + theme(axis.title.x=element_blank()) + o_margin(.2,.2,0,.2)
#pb = plot_dde_bar(tp2) + o_margin(0,.2,.2,.2)
#p = ggarrange(pa, pb, nrow=2, ncol=1, labels = '', heights=c(1,1))
p = plot_dde_bar(tp)
fo = glue("{dirw}/36.dde.bar.pdf")
ggsave(p, file=fo, width=4, height=7)
fo = glue("{dirf}/f6b.pdf")
ggsave(p, file=fo, width=4, height=7)
#p = ggarrange(pa, pb, nrow=2, ncol=1, heights=c(1,1))
saveRDS(p, glue("{dirf}/f6b.rds"))
#}}}

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
    scale_y_continuous(name='Number/Proportion of Genes', expand=expansion(mult=c(.02,.05))) +
    scale_fill_manual(name="Mode", values=fillcols) +
    facet_wrap(~pnl, nrow=2) +
    otheme(legend.pos='top.center.out',legend.dir='h',
           legend.title=T, legend.vjust=-.7, panel.spacing=.2,
           legend.spacing.x=.05, legend.spacing.y=.05,
           xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
           margin = c(1, .3, .3, .3)) +
    guides(fill=guide_legend(nrow=2))
#}}}
}

p = plot_cis_bar(tp)
fo = glue("{dirw}/38.dde.cis.pdf")
ggsave(p, file=fo, width=4, height=7)
fo = glue("{dirf}/f6c.pdf")
ggsave(p, file=fo, width=4, height=7)
saveRDS(p, glue("{dirf}/f6c.rds"))
#}}}

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
    guides(color=F)
fo = glue("{dirw}/36.dde.acc.pdf", dirw)
p %>% ggexport(filename=fo, width=14, height=5)
#}}}
#}}}

#{{{ [old] accuracy in predicting U/D/N at 1h/25h time points
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
    guides(color=F)
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
    guides(color=F)
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

