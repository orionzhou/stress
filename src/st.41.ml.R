source('functions.R')
require(progress)
dirw = file.path(dird, '41_ml')
#{{{ functions
filter_kmer_loc <- function(tl, ts) tl %>% filter(gid %in% ts$gid)
prepare_ti <- function(ti, bin, epi, nfea, mod, tst, keep_gid=F) {
  #{{{
  pb$tick()
  to = ti %>% mutate(pos = pos - 2000)
  if(bin == '-500') {
    to = to %>% filter(pos >= -500, pos <= 0)
  } else if(bin == '+500') {
    to = to %>% filter(pos >= 0, pos <= 500)
  } else if(bin == '+/-500') {
    to = to %>% filter(pos >= -500, pos <= 500)
  } else if(bin == '-2k') {
    to = to %>% filter(pos <= 0)
  } else if(bin == '+2k') {
    to = to %>% filter(pos >= 0)
  }
  if(epi == 'umr') {
    to = to %>% filter(umr>0)
  }
  if(nfea == 'top30') {
    to = to %>% filter(i <= 30)
  } else if(nfea == 'top50') {
    to = to %>% filter(i <= 50)
  } else if(nfea == 'top100') {
    to = to %>% filter(i <= 100)
  }
  if(mod == 'zoops') {
    to = to %>% distinct(gid, fid) %>% mutate(x = 1)
  } else if(mod == 'anr') {
    to = to %>% count(gid, fid) %>% rename(x = n)
  }
  fids = to %>% distinct(fid) %>% pull(fid)
  x = tst %>% crossing(fid=fids) %>%
    left_join(to, by=c('gid','fid')) %>%
    replace_na(list(x=0)) %>%
    spread(fid, x)
  if(!keep_gid)
    x %>% select(-gid)
  else
    x
  #}}}
}
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
# run pipeline to build B(M/W) models


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

#{{{ # obtain balanced F1 score and maximize F1 score during permutation
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

######## MODEL EVALUATION ########
#{{{ functions to read models
get_model_fits <- function(tag, dirw) {
  #{{{
  fi = glue("{dirw}/01_models/{tag}.tsv")
  th = read_tsv(fi)
  fi = glue("{dirw}/01_models/{tag}.rds")
  ml = readRDS(fi)
  ml %>% filter(!is.na(fit)) %>% select(did=sid,perm,fit,metric) %>%
    inner_join(th, by='did')
  #}}}
}
get_model_metrics <- function(tag, dirw) {
  #{{{
  fi = glue("{dirw}/01_models/{tag}.tsv")
  th = read_tsv(fi)
  fi = glue("{dirw}/01_models/{tag}.rds")
  ml = readRDS(fi)
  ml1 = ml %>% select(did=sid, perm, metric) %>%
    unnest(metric) %>% spread(metric, estimate) %>%
    rename(f1=f_meas, auroc=roc_auc, auprc=pr_auc) %>%
    inner_join(th, by='did')
  ml1
  #}}}
}

tm = tibble(tag = c('mod4_288_c2','mod7_252')) %>%
  mutate(x = map(tag, get_model_metrics, dirw=dirw)) %>%
  select(-mid0) %>% replace_na(list(mod='zoops'))
fm = file.path(dirw, '01.model.metrics.rds')
saveRDS(tm, fm)

fm = file.path(dirw, '01.model.metrics.rds')
tm = readRDS(fm)
  #{{{ prepare tm for plot
  bins = c("-500","+500","+/-500","-2k","+2k","+/-2k")
  nfeas = c("top30",'top50','top100')
  mods = c("zoops",'anr')
  epis = c("All Sequence", "UMR")
  notes = unique(tm$note)
  bat_notes = crossing(bat=bats, note=notes) %>%
    mutate(bat_note = str_c(bat, note, sep=': ')) %>%
    mutate(bat=factor(bat, levels=bats)) %>%
    mutate(note=factor(note, levels=notes)) %>%
    arrange(bat, note) %>% pull(bat_note)
  tm = tm %>%
    mutate(epi = ifelse(epi=='raw', 'All Sequence', 'UMR')) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(note = as_factor(note)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(mod = factor(mod, levels=mods)) %>%
    arrange(tag, bat, note, bin, epi, nfea, mod) %>%
    mutate(bat_note = str_c(bat, note, sep=': ')) %>%
    mutate(bat_note = factor(bat_note, levels=bat_notes)) %>%
    mutate(bat_mod = as_factor(str_c(bat, mod, sep=': '))) %>%
    mutate(bat_epi = as_factor(str_c(bat, epi, sep=': '))) %>%
    mutate(epi_nfea = as_factor(str_c(epi, nfea, sep=': ')))
  #}}}
tm1 = tm %>% filter(tag == 'mod4_288_c2')
tm2 = tm %>% filter(tag == 'mod7_252')

#}}}

#{{{ eval model performance in B
plot_model_heatmap <- function(ti, metric='f1', pnl='bat_mod') {
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
plot_model_barplot <- function(ti, metric='f1', x='bin', pnl='bat_mod', tit=F) {
  #{{{
  ymax = ifelse(x=='bat_note', .95, .8)
  tp = ti %>%# filter(tag == 'mod4_288_c2') %>%
    mutate(score = get(metric)) %>%
    filter(!is.na(score)) %>%
    mutate(x = get(x), pnl = get(pnl)) %>%
    group_by(bat, note, bin, epi, nfea, mod, x, pnl) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup()
    #{{{
    pd = position_dodge(width=.9)
    p = ggplot(tp, aes(x=x,y=score)) +
      geom_col(aes(fill = epi), position=pd, width=.6) +
      geom_errorbar(aes(ymin=score-sd, ymax=score+sd, color=epi), width=.2, size=.3, position=pd) +
      scale_x_discrete(expand=expansion(mult=c(.1,.1))) +
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
  if(x=='bat_note')
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
pa = plot_model_heatmap(tm1, metric='f1', pnl='bat_mod')
pb = plot_model_heatmap(tm1, metric='auroc', pnl='bat_mod') +
  theme(plot.margin = margin(.3, .3, .3, 0, "lines")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
fo = glue("{dirf}/sf09a.pdf")
ggarrange(pa, pb, nrow=1, ncol=2,
          labels = c(), widths=c(2,1.4)) %>%
ggexport(filename=fo, width=6, height=8)

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
tm3a = tm1 %>% filter(nfea=='top30', mod=='zoops')
#pa = plot_model_barplot(tm3a, metric='f1', pnl='bat_note')
pb = plot_model_barplot(tm3a, metric='auroc', pnl='bat_note') +
  o_margin(.2,.2,.2,.2)
#
tm3b = tm2 %>% filter(nfea=='top30', bin=='+/-2k', mod=='zoops')
#pc = plot_model_barplot(tm3b, metric='f1', x='bat_note', pnl='bin')
pd = plot_model_barplot(tm3b, metric='auroc', x='bat_note', pnl='bin') +
  o_margin(0,.2,.2,.2) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))
#
fo = glue("{dirf}/f4.pdf")
pab = ggarrange(pa, pb, nrow=1, ncol=2, widths=c(2,1.8))
pcd = ggarrange(pc, pd, nrow=1, ncol=2, widths=c(2,1.8))
ggarrange(pb, pd, nrow=2, ncol=1,
          labels = LETTERS[1:2], heights=c(2,1)) %>%
ggexport(filename=fo, width=4, height=6)
#}}}

#}}}

#{{{ evaluate model performance in M/W
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
tt = readRDS(file.path(dirw, gt2, "01.rds"))
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
fg = file.path(dird, '15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2
fp = file.path(dirw, '10.model.pred.rds')
pd = readRDS(fp) %>% mutate(gt=str_replace(gt,'Zmays_',''))
pd2 = pd %>% select(bat,note,gt,gid,pred) %>%
  separate(bat, c('cond','drc'), sep='_', remove=F)

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

#{{{ de acc bar plot
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

pa = plot_dde_bar(tp1)
pb = plot_dde_bar(tp2)
fo = glue("{dirf}/f6.pdf")
#sprintf("%s/13.dde.bar.pdf", dirw)
ggarrange(pa, pb, nrow=2, ncol=1,
          labels = LETTERS[1:2], heights=c(1,1)) %>%
ggexport(filename=fo, width=5, height=7)
#}}}

#{{{ de acc plot
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

######## BELOW ARE DEPRECATED ########
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
