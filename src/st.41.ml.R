source('functions.R')
require(progress)
dirw = glue('{dird}/41_ml')

#{{{ prepare training input: gene lists + kmer lists
fi = glue("{dird}/25_dreme/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk

# gene lists
tc %>% mutate(fo=glue("{dirw}/01_gene_lists/{cid}.tsv")) %>%
    mutate(x = map2(ts, fo, write_tsv))

t1 = tk %>% select(mid,fid,fname,known, lid,pval) %>%
    inner_join(tl %>% select(lid, cid, bin, epi), by='lid')
t2 = t1 %>% arrange(cid, pval) %>%
    separate(bin,c('opt','bin'),sep=":", remove=F) %>%
    arrange(cid, pval) %>%
    group_by(cid, fid, fname) %>%
    slice(1) %>% ungroup() %>%
    arrange(cid, pval) %>%
    group_by(cid) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    inner_join(tk %>% select(mid,kmer,kmers), by='mid') %>%
    mutate(kmers = map_chr(kmers, str_c, collapse=',')) %>%
    select(cid, i, opt,bin,epi, pval,fid,fname,kmers)
t3 = t2 %>% group_by(cid) %>% nest() %>% rename(kmer=data) %>% ungroup()

fo = glue("{dirw}/03.best.mtfs.rds")
saveRDS(t3, fo)

t3 %>% mutate(fo = glue("{dirw}/03_motif_lists/{cid}.tsv")) %>%
    mutate(l = map2(kmer, fo, write_tsv))
#}}}

#{{{ motif meta plots
fi = glue("{dird}/25_dreme/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk

#{{{ prepare module genen lists to get  motif locations
ts1 = tc %>% select(cid, cond, note, gids)
ts2 = tc %>% distinct(cond, gids_c) %>% mutate(cid=c('c29','c49')) %>%
    mutate(note='bg') %>% select(cid,cond,note,gids=gids_c)
ts = ts1 %>% bind_rows(ts2) %>% arrange(cid) %>%
    unnest(gids) %>% rename(gid=gids)

fo = glue("{dirw}/51.mod.genes.rds")
saveRDS(ts, fo)
to = ts %>% distinct(gid) %>% mutate(status=1)
fo = glue("{dirw}/51.mod.genes.tsv")
write_tsv(to, fo, na='')
#}}}

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

fi = glue("{dirw}/03.best.mtfs.rds")
rb = readRDS(fi)
rb2 = rb %>% unnest(kmer) %>%
    mutate(fnote = glue("{cid}_{i} {fid} ({fname}): {opt} {bin} {epi}")) %>%
    select(cid,i,fid,fnote)

ti = tc %>% select(cid,cond,note) %>%
    mutate(fi = glue("{dirw}/52_mtf_loc/{cid}.tsv")) %>%
    mutate(x = map(fi, read_mtf_loc)) %>%
    select(-fi)

tss = ts %>% count(cid, cond, note) %>% rename(nt=n) %>%
    mutate(mod=glue("{cond}: {note} ({nt})")) %>%
    mutate(mod = as_factor(mod)) %>%
    mutate(fill = ifelse(note=='bg', 'grey', 'white')) %>%
    select(cid, nt, mod, fill)
ts2 = ts %>% inner_join(tss, by='cid') %>% select(nt, gid, mod)

cid = 'c31'
ti2 = ti %>% filter(cid == !!cid) %>% unnest(x) %>%
    inner_join(rb2, by=c('cid','fid'))
#
ti3 = ti2 %>% distinct(cid,cond,note, i,fid,fnote, gid, bin) %>%
    inner_join(ts2, by=c("gid")) %>%
    count(cid,cond,note, i,fid,fnote, mod,nt, bin) %>% rename(nh=n) %>%
    mutate(prop = nh/nt)

#{{{ meta plot of selected TFBS motifs
plot_tss_meta <- function(cid, i, tss, rb2, ti3, dirw) {
    #{{{
    tit1=tss %>% filter(cid==!!cid) %>% pull(mod)
tit2 = rb2 %>% filter(cid == !!cid, i==!!i) %>% pull(fnote)
tit = glue("{tit1} | {tit2}")
tp = ti3 %>% filter(cid == !!cid, i== !!i)
tpx = tibble(x=c(.5,10.5,20.5,31.5,41.5),lab=c('-2kb','TSS','+2kb/-2kb','TTS','+2kb'))
p = ggplot(tp, aes(x=bin,y=prop)) +
    geom_line(aes(), size=.5, na.rm = F) +
    geom_point(aes(), size=1, na.rm = F) +
    scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
    scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    #scale_shape(labels=types) +
    #scale_linetype(labels=types) +
    facet_wrap(~mod, scale='free_x', ncol=4) +
    ggtitle(tit) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           strip.style='white',margin = c(.3,.3,.3,.3),
           xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
    theme(plot.title=element_text(hjust=.5, size=10)) +
    guides(fill=F)
fo = glue("{dirw}/53_metaplots/{cid}_{i}.pdf")
ggsave(p, file=fo, width = 8, height = 6)
    #}}}
}
tibble(cid=!!cid, i=1:10) %>%
    mutate(l=map2(cid,i, plot_tss_meta, tss=tss,rb2=rb2,ti3=ti3,dirw=dirw))
#}}}
#}}}

#{{{ create ML training config
fi = glue("{dird}/25_dreme/03.mtf.grp.rds")
r3 = readRDS(fi)
tl = r3$tl; tc = r3$tc; tk = r3$tk
fi = glue("{dirw}/03.best.mtfs.rds")
tm = readRDS(fi)

#{{{ b1: 4 modules * 14 bin * 2 epis * 2 nfeas * 2 mods = 448
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

ctag = 'b1'
to1 = tc %>% select(cid) %>%
  crossing(bin = bins, epi = epis1, nfea = nfeas1, mod = mods1)
to2 = tc %>% select(cid) %>%
  crossing(bin = bins1, epi = epis, nfea = nfeas, mod = mods)
to = to1 %>% bind_rows(to2) %>%
    distinct(cid, bin, epi, nfea, mod) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(cid, bin, epi, nfea, mod) %>%
    mutate(did = sprintf("d%03d", 1:n())) %>%
    mutate(gt = 'Zmays_B73') %>%
    select(did, cid, bin, epi, nfea, mod, gt)

fo = glue("{dirw}/10_models/{ctag}.cfg.rds")
saveRDS(to, fo)
fo = glue("{dirw}/10_models/{ctag}.cfg.tsv")
write_tsv(to, fo)
#}}}

#{{{ b2: 4 modules * 7 bin * 1 epis * 1 nfeas * 2 mods = 56
ctag = 'b2'
to0 = tl %>% select(-ng0) %>% filter(gt == 'Zmays_B73') %>%
  inner_join(tm, by='cid') %>%
  mutate(n_mtf = map_int(kmer, nrow)) %>%
  filter(n_mtf >= 20)
#
bins = c(
         'TSS:-/+2k','TSS:-2k,TSS:+2k',
         "TSS:-2k,TSS:+2k,TTS:-500", 'TSS:-/+2k,TTS:-500'
)
epis = c('umr')
nfeas = c('top100')
mods = c('anr')
to = to0 %>%
  crossing(bin = bins, epi = epis, nfea = nfeas, mod = mods) %>%
  mutate(nfea = factor(nfea, levels=nfeas)) %>%
  mutate(bin = factor(bin, levels=bins)) %>%
  mutate(epi = factor(epi, levels=epis)) %>%
  mutate(mod = factor(mod, levels=mods)) %>%
  arrange(bnid, bin, epi, nfea, mod) %>%
  mutate(did = sprintf("d%03d", 1:n())) %>%
  select(did, bnid, bat, cond, drc, note, bin, epi, nfea, mod, gt)
#
fo = glue("{dirw}/10_models/{ctag}.cfg.rds")
saveRDS(to, fo)
fo = glue("{dirw}/10_models/{ctag}.cfg.tsv")
write_tsv(to, fo)
#}}}
#}}}
# run pipeline to build B models

######## MODEL EVALUATION ########
#{{{ process trained models
collect_model_metrics <- function(tag, dirw) {
  #{{{
  fi = glue("{dirw}/01_models/{tag}.cfg.rds")
  th = readRDS(fi)
  fi = glue("{dirw}/01_models/{tag}.rds")
  ml = readRDS(fi)
  ml1 = ml %>% select(did=sid, perm, metric) %>%
    unnest(metric) %>% spread(metric, estimate) %>%
    rename(f1=f_meas, auroc=roc_auc, auprc=pr_auc) %>%
    inner_join(th, by='did')
  ml1
  #}}}
}
get_model_fits <- function(tag, dirw) {
  #{{{
  fi = glue("{dirw}/01_models/{tag}.cfg.rds")
  th = readRDS(fi)
  fi = glue("{dirw}/01_models/{tag}.rds")
  ml = readRDS(fi)
  ml %>% filter(!is.na(fit)) %>% select(did=sid,perm,fit,metric) %>%
    inner_join(th, by='did')
  #}}}
}

tags = c("b1")
tm = tibble(tag = tags) %>%
  mutate(x = map(tag, collect_model_metrics, dirw=dirw))
fm = glue('{dirw}/11.model.metrics.rds')
saveRDS(tm, fm)

#}}}

#{{{ eval model performance in B
fm = glue('{dirw}/11.model.metrics.rds')
tm = readRDS(fm)
tm1 = tm %>% filter(tag == 'b2') %>% unnest(x) %>%
  mutate(note = as_factor(note)) %>%
  mutate(bat_note = as_factor(str_c(bat, note, sep=': '))) %>%
  mutate(bat_note_nfea_mod = as_factor(str_c(bat,note,nfea,mod,sep=': '))) %>%
  mutate(bin_epi = as_factor(str_c(bin, epi, sep=': '))) %>%
  mutate(bat_mod = as_factor(str_c(bat, mod, sep=': '))) %>%
  mutate(bat_epi = as_factor(str_c(bat, epi, sep=': '))) %>%
  mutate(epi_nfea = as_factor(str_c(epi, nfea, sep=': ')))

plot_model_heatmap <- function(ti, metric='auroc', pnl='bat_mod') {
  #{{{ plot
  tp = ti %>%# filter(tag == 'mod4_288_c2') %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = bin, y = epi_nfea, pnl = get(pnl)) %>%
    group_by(bat, note, bin, epi, nfea, mod, x, y, pnl) %>%
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
      facet_wrap(pnl ~ ., ncol=1) +
      otheme(legend.pos='none', legend.dir='v', legend.title=F,
             margin = c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=T, panel.spacing=.1,
             xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
      theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=.5, size=7)) +
      theme(axis.text.y = element_text(size=7)) +
      #theme(axis.text.y = element_markdown(size=7.5)) +
      guides(color=F)
    #}}}
  tit = str_to_upper(ifelse(metric=='f1','f-score',metric))
  p + ggtitle(tit)
  #}}}
}
plot_model_barplot <- function(ti, metric='auroc', x='bin', fill='epi',
                               pnl='bat_mod', x.rotate=F, tit=F) {
  #{{{
  ymax = 1.0 * max(ti[[metric]])
  tp = ti %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = get(x), pnl = get(pnl)) %>%
    group_by(bat, note, bin, epi, nfea, mod, x, pnl) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup()
    #{{{
    pd = position_dodge(width=.7)
    p = ggplot(tp, aes(x=x,y=score)) +
      geom_col(aes(fill = get(fill)), position=pd, width=.7) +
      geom_errorbar(aes(ymin=score-sd, ymax=score+sd, color=get(fill)), width=.2, size=.3, position=pd) +
      scale_x_discrete(expand=expansion(mult=c(.03,.03))) +
      scale_y_continuous(name=str_to_upper(metric), expand=expansion(mult=c(0,.05))) +
      coord_cartesian(ylim= c(.5, ymax)) +
      scale_fill_npg(name='') +
      scale_color_manual(values=c('black','black')) +
      facet_wrap(pnl ~ ., scale='free_y', ncol=1) +
      otheme(legend.pos='none', legend.dir='v', legend.box='v',legend.title=F,
             margin = c(.3,.3,.3,.3), ygrid=T,
             strip.style='white', strip.compact=F, panel.spacing=.1,
             xtick=T, ytick=T, ytitle=T, xtext=T, ytext=T) +
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

#{{{ heatmap: sf09
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

#{{{ f4 - barplot
tpa = tm1 %>% filter(nfea=='top100', mod=='zoops')
pa = plot_model_barplot(tpa, x='bin', fill='epi',
                        pnl='bat_note_nfea_mod', x.rotate=T) +
  o_margin(.2,.2,.2,.2)
tpb = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', nfea=='top100')
pb = plot_model_barplot(tpb, x='bat_note', fill="mod",
                        pnl='bin_epi', x.rotate=T) +
  o_margin(0,.2,.2,.2) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))
tpc = tm1 %>% filter(bin=='TSS:-/+2k', epi=='umr', mod=='anr')
pc = plot_model_barplot(tpc, x='bat_note', fill="nfea",
                        pnl='bin_epi', x.rotate=T) +
  o_margin(0,.2,.2,.2) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))
#
fo = glue("{dirf}/f4.pdf")
pbc = ggarrange(pb, pc, nrow=1, ncol=2, labels=LETTERS[2:3], widths=c(2,2))
ggarrange(pa, pbc, nrow=2, ncol=1,
          labels = c("A"), heights=c(2,1)) %>%
ggexport(filename=fo, width=5, height=6)
#}}}

#}}}

#{{{ gather & save model predictions for all B/M/W genes
#{{{ functions
get_ts_all <- function(gt, dird) {
  #{{{
  fp = file.path(dird, '31_promoter', gt, '01.bed')
  tp = read_tsv(fp, col_names=c('chrom','start','end','gid','score','srd'))
  ts = tp %>% select(gid) %>% mutate(status=1)
  #}}}
}
model_pred <- function(mdl, ti) {
  #{{{
  pb$tick()
  features = names(mdl$fit$variable.importance)
  new.cols = features[!features %in% colnames(ti)]
  if (length(new.cols) > 0) {
    ti2 = matrix(0, nrow=nrow(ti), ncol=length(new.cols), dimnames=list(NULL, new.cols))
    ti = ti %>% bind_cols(as.data.frame(ti2))
  }
  x1 = predict(mdl, ti, type='prob') %>% select(prob = .pred_1)
  x2 = predict(mdl, ti, type='class') %>% select(pred=1)
  tibble(gid=ti$gid, pred=x2$pred, prob=x1$prob)
  #}}}
}
#}}}

gts = str_c("Zmays", gts3, sep="_")
tt = tibble(gt=gts) %>%
  mutate(ts = map(gt, get_ts_all, dird=dird)) %>%
  mutate(floc = file.path(dird, '32_kmer', gt, '05.mtf.loc.rds')) %>%
  mutate(mloc = map(floc, readRDS)) %>% select(-floc) %>%
  unnest(mloc) %>%
  filter(note == 'all') %>%
  mutate(bin="+/-2k", epi='umr', nfea='top100', mod='anr')

pb = progress_bar$new(total=nrow(tt))
tt2 = tt %>%
  mutate(to = pmap(list(tl,bin,epi,nfea,mod,ts), prepare_ti, keep_gid=T)) %>%
  select(bat, note, bin, epi, nfea, mod, gt, to)

ti = tb %>% inner_join(tt2, by=c('bat','note','bin','epi','nfea','mod')) %>%
  slice(1:50)
pb = progress_bar$new(total=nrow(ti))
r = ti %>%
  mutate(pred = map2(fit, to, model_pred)) %>%
  select(-fit,-to) %>%
  unnest(pred) %>%
  select(did,perm,bat,note,bin,epi,nfea,mod,gt,gid,pred,prob)

fo = file.path(dirw, '10.model.pred.rds')
saveRDS(r, fo)
#}}}

#{{{ evaluate model performance for different gene categories
fp = file.path(dirw, '10.model.pred.rds')
pd = readRDS(fp) %>% mutate(gt=str_replace(gt,'Zmays_',''))
pd2 = pd %>% select(bat,note,gt,gid,pred,prob) %>%
  separate(bat, c('cond','drc'), sep='_', remove=F)

#{{{ overall performance
make_truth <- function(tg, tg_c) tg %>% mutate(status=1) %>% bind_rows(tg_c %>% mutate(status = 0))
xm = tibble(gt=gts3) %>%
  mutate(fm=glue("{dird}/31_promoter/{gt2}/15.module.rds")) %>%
  mutate(x = map(fm, readRDS)) %>%
  unnest(x) %>% filter(ctag=='c2') %>% select(-fm,-ctag) %>%
  mutate(truth = map2(tg, tg_c, make_truth)) %>%
  select(-tg, -tg_c) %>% unnest(truth) %>%
  mutate(status=factor(status, levels=c(1,0)))

perm=10
x = pd2 %>% inner_join(xm, by=c('bat','cond','drc','note','gt','gid')) %>%
  select(-bat_mid,-mid) %>%
  group_by(bat,cond,drc,note,ng0,ng,ng_c, gt) %>% nest() %>% rename(pred=data) %>%
  crossing(i = 1:perm) %>%
  mutate(x = map2(pred, i, eval_pred, downsample=T)) %>%
  select(-i) %>% unnest(x) %>%
  group_by(bat,note,gt,metric) %>%
  summarise(mean=mean(estimate), sd=sd(estimate)) %>% ungroup()

#{{{ bar plot - f6a
tp = x %>% filter(metric=='f_meas') %>% rename(score=mean) %>%
  mutate(gt = factor(gt, levels=rev(gts3))) %>%
  mutate(bat = factor(bat, levels=bats))
pd = position_dodge(width=.9)
ymax = .8
p = ggplot(tp, aes(x=gt,y=score)) +
  geom_col(aes(fill='0'), position=pd, width=.6, alpha=.8) +
  geom_errorbar(aes(ymin=score-sd, ymax=score+sd), width=.2, size=.3, position=pd) +
  scale_x_discrete(expand=expansion(mult=c(.2,.2))) +
  scale_y_continuous(name='AUROC', expand=expansion(mult=c(0,.05))) +
  coord_flip(ylim= c(.5, ymax)) +
  scale_fill_manual(name='', values=pal_npg()(2)[2]) +
  scale_color_manual(values=c('black','black')) +
  facet_wrap(bat ~ ., scale='free_y', ncol=1, strip.position='right') +
  otheme(legend.pos='none', legend.dir='v', legend.box='v',legend.title=F,
         margin = c(.2,.2,.2,.2), xgrid=T,
         strip.style='white', strip.compact=F, panel.spacing=.2,
         xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T) +
  theme(strip.text.y = element_text(angle=0,size=8))
ggsave(p, file=glue("{dirw}/09.bmw.auroc.pdf"), width=4.5, height=5)
saveRDS(p, glue("{dirf}/f6a.rds"))
#}}}
#}}}

# variable response
fg = file.path(dird, '15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2

#{{{ accuracy in predicting U/D/N at 1h/25h time points
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
to1 = pd2 %>%
  inner_join(td1, by=c('cond','gt','gid')) %>%
  group_by(bat,cond,drc,note,gt,time) %>% nest() %>% rename(ti=data) %>%
  mutate(acc = map(ti, get_acc)) %>%
  select(-ti) %>% unnest(acc)

#{{{ de acc plot
tp = to1 %>% rename(score=acc) %>%
  mutate(bat = factor(bat, levels=bats)) %>%
  mutate(gt = as_factor(gt)) %>%
  mutate(st = as_factor(st)) %>%
  mutate(time = as_factor(time)) %>%
  mutate(note = as_factor(note)) %>%
  mutate(time_st = fct_cross(time, st, sep='_')) %>%
  mutate(bat_note_gt = fct_cross(bat, note, gt, sep=":"))
tp = tp %>% mutate(bat_note_gt = factor(bat_note_gt, levels=rev(levels(tp$bat_note_gt))))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=time_st,y=bat_note_gt)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    #geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='Accuracy',colors=brewer_pal(type='seq',palette="Greens")(9)) +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
           margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/12.de.acc.pdf", dirw)
p %>% ggexport(filename=fo, width=6, height=6)
#}}}
#}}}

#{{{ accuracy in predicting variable genotype response
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
pd3 = pd2 %>% mutate(pred=as.integer(as.character(pred))) %>%
  mutate(pred=ifelse(drc=='down' & pred==1, -1, pred))
to2 = td2 %>%
  inner_join(pd3, by=c('cond','qry'='gt','gid')) %>% rename(qPred=pred) %>%
  inner_join(pd3, by=c('bat','cond','drc','note','tgt'='gt','gid')) %>% rename(tPred=pred) %>%
  #group_by(bat,cond,drc,note,qry,tgt,time) %>% nest() %>% ungroup() %>%
  mutate(pred1 = smap2[as.character(qPred)]) %>%
  mutate(pred2 = smap2[as.character(tPred)]) %>%
  mutate(pred = as.character(glue("A{pred1}B{pred2}"))) %>%
  count(bat,cond,drc,note,qry,tgt,time, st, pred)

#{{{ dde acc bar plot
tp = to2 %>% filter(tgt=='B73') %>%
  filter((drc=='up' & st %in% c("dA+B+",'nA+B+','dA+B=','dA=B+')) |
         (drc=='down' & st %in% c('dA-B-','nA-B-','dA-B=', 'dA=B-'))) %>%
  mutate(st = as.character(st)) %>%
  mutate(st = ifelse(str_detect(st, 'A\\+B\\+'),'A+B+', st)) %>%
  mutate(st = ifelse(str_detect(st, 'A\\-B\\-'),'A-B-', st)) %>%
  mutate(st = str_replace(st, 'd', '')) %>%
  mutate(st = str_replace(st, 'n', '')) %>%
  mutate(bat = factor(bat, levels=bats)) %>%
  mutate(time = as_factor(time)) %>%
  mutate(note = as_factor(note)) %>%
  mutate(pred = as_factor(pred)) %>%
  group_by(bat,cond,drc,note,st,pred) %>% summarise(n=sum(n)) %>% ungroup() %>%
  rename(tag1=st, tag2=pred, pnl=bat)
preds1 = c('A=B=',"A+B=","A=B+","A+B+")
preds2 = c('A=B=',"A-B=","A=B-","A-B-")
tp1 = tp %>% filter(drc == 'up') %>% mutate(tag2=factor(tag2,levels=preds1))
tp2 = tp %>% filter(drc == 'down') %>% mutate(tag2=factor(tag2,levels=preds2))

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
p = ggplot(tp) +
    geom_bar(aes(x=tag1,y=prop,fill=tag2), stat='identity', position='fill', width=.8) +
    geom_text(aes(x=tag1,y=y,label=lab),size=2.5,lineheight=.8) +
    geom_text(data=tpx, aes(tag1,1.01,label=n), color='black',size=3, vjust=0) +
    scale_x_discrete(name='Observed Status', expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number of Genes', expand=expansion(mult=c(0,.05))) +
    scale_fill_manual(name="Prediction", values=fillcols) +
    facet_wrap(~pnl, nrow=1) +
    otheme(legend.pos='right',legend.dir='v',
           legend.title=T, legend.vjust=-.3, panel.spacing=.2,
           strip.style='white',
           xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F,
           margin = c(.3, .3, .3, .3))
#}}}
}

pa = plot_dde_bar(tp1) + theme(axis.title.x=element_blank()) + o_margin(.2,.2,0,.2)
pb = plot_dde_bar(tp2) + o_margin(0,.2,.2,.2)
p = ggarrange(pa, pb, nrow=2, ncol=1, labels = LETTERS[1:2], heights=c(1,1))
ggsave(p, file=glue("{dirw}/13.dde.bar.pdf"), width=5, height=7)
p = ggarrange(pa, pb, nrow=2, ncol=1, heights=c(1,1))
saveRDS(p, glue("{dirf}/f6b.rds"))
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
fo = sprintf("%s/13.dde.acc.pdf", dirw)
p %>% ggexport(filename=fo, width=14, height=5)
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

