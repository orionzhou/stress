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
#{{{ read models
tag = 'mod4_288_c2'
tag = 'mod7_252'
fi = glue("{dirw}/01_models/{tag}.tsv")
th = read_tsv(fi)
fi = glue("{dirw}/01_models/{tag}.rds")
ml = readRDS(fi)
ml1 = ml %>% select(did=sid, perm, metric) %>%
  unnest(metric) %>% spread(metric, estimate) %>% rename(f1=f_meas) %>%
  inner_join(th, by='did')
tb = ml %>% filter(!is.na(fit)) %>% select(did=sid,perm,fit,metric) %>%
  inner_join(th, by='did')
#}}}

#{{{ eval model performance in B
plot_model_perf <- function(ti, metric='f1') {
  #{{{ plot
  bins = c("-500","+500","+/-500","-2k","+2k","+/-2k")
  nfeas = c("top30",'top50','top100')
  tp = ti %>%# filter(mod=='zoops') %>%
    mutate(score = get(metric)) %>%
    group_by(bat, note, bin, epi, nfea, mod) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup() %>%
    #summarise(score=median(roc_auc)) %>% ungroup() %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(note = factor(note)) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(bat_note_mod = fct_cross(bat, note, mod, sep=': ')) %>%
    mutate(epi_nfea = fct_cross(nfea, epi, sep=':')) %>%
    mutate(y = str_c(bat, note, nfea, sep=":")) %>%
    mutate(x = str_c(bin, epi, sep=":")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}%+-%{number(sd,accuracy=.01)}")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}")) %>%
    filter(!is.na(score))
  swit = (min(tp$score) + max(tp$score)) / 2
  p = ggplot(tp, aes(x=bin,y=epi_nfea)) +
      geom_tile(aes(fill=score), na.rm = F) +
      geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2, parse=T) +
      #geom_vline(xintercept=tpy$x, color='blue') +
      scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
      scale_y_discrete(expand=c(0,0)) +
      scale_fill_gradientn(name='F1 score',colors=cols100v) +
      #scale_fill_viridis(name='normalized eigengene value') +
      scale_color_manual(values=c('black','white')) +
      facet_wrap(bat_note_mod ~ ., nrow=4) +
      otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
             margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
             xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
      theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
      theme(axis.text.y = element_text(size=7)) +
      #theme(axis.text.y = element_markdown(size=7.5)) +
      guides(color=F)
  p
  #}}}
}
plot_model_perf1 <- function(ti, metric='f1') {
  #{{{ plot
  bins = c("-500","+500","+/-500","-2k","+2k","+/-2k")
  nfeas = c("top30",'top50','top100')
  tp = ti %>%
    filter(nfea == 'top30', mod=='zoops') %>%
    mutate(score = get(metric)) %>%
    group_by(bat, note, bin, epi, nfea, mod) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup() %>%
    #summarise(score=median(roc_auc)) %>% ungroup() %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(note = factor(note)) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    group_by(bat, note, bin, epi, nfea, mod) %>%
    mutate(bat_note_mod = fct_cross(bat, note, mod, sep=': ')) %>%
    mutate(bat_note_epi = str_c(bat, note, epi, sep=': ')) %>%
    mutate(bat_note_epi = as_factor(bat_note_epi)) %>%
    mutate(epi_nfea = fct_cross(nfea, epi, sep=':')) %>%
    mutate(y = str_c(bat, note, nfea, sep=":")) %>%
    mutate(x = str_c(bin, epi, sep=":")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}%+-%{number(sd,accuracy=.01)}")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}")) %>%
    filter(!is.na(score))
  tp = tp %>% mutate(bat_note_epi = factor(bat_note_epi, levels=rev(levels(tp$bat_note_epi))))
  swit = (min(tp$score) + max(tp$score)) / 2
  p = ggplot(tp, aes(x=bin,y=bat_note_epi)) +
      geom_tile(aes(fill=score), na.rm = F) +
      geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2, parse=T) +
      #geom_vline(xintercept=tpy$x, color='blue') +
      scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
      scale_y_discrete(expand=c(0,0)) +
      scale_fill_gradientn(name='F1 score',colors=cols100v) +
      #scale_fill_viridis(name='normalized eigengene value') +
      scale_color_manual(values=c('black','white')) +
      #facet_wrap(bat_note_mod ~ ., nrow=4) +
      otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
             margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
             xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
      theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
      theme(axis.text.y = element_text(size=7)) +
      #theme(axis.text.y = element_markdown(size=7.5)) +
      guides(color=F)
  p
  #}}}
}
plot_model_perf2 <- function(ti, metric='f1') {
  #{{{ plot
  bins = c("-500","+500","+/-500","-2k","+2k","+/-2k")
  nfeas = c("top30",'top50','top100')
  tp = ti %>%
    filter(nfea == 'top30', bin=='+/-2k') %>%
    mutate(score = get(metric)) %>%
    group_by(bat, note, bin, epi, nfea) %>%
    summarise(sd = sd(score), score = mean(score)) %>% ungroup() %>%
    #summarise(score=median(roc_auc)) %>% ungroup() %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(note = factor(note)) %>%
    mutate(nfea = factor(nfea, levels=nfeas)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    group_by(bat, note, bin, epi, nfea) %>%
    mutate(bat_note = fct_cross(bat, note, sep=': ')) %>%
    mutate(bat_note_epi = str_c(bat, note, epi, sep=': ')) %>%
    mutate(bat_note_epi = as_factor(bat_note_epi)) %>%
    mutate(epi_nfea = fct_cross(nfea, epi, sep=':')) %>%
    mutate(y = str_c(bat, note, nfea, sep=":")) %>%
    mutate(x = str_c(bin, epi, sep=":")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}%+-%{number(sd,accuracy=.01)}")) %>%
    mutate(lab = glue("{number(score,accuracy=.01)}")) %>%
    filter(!is.na(score))
  swit = (min(tp$score) + max(tp$score)) / 2
  tp = mutate(tp, epi=factor(epi, levels=rev(epis)))
  p = ggplot(tp, aes(x=bat_note,y=epi)) +
      geom_tile(aes(fill=score), na.rm = F) +
      geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2, parse=T) +
      #geom_vline(xintercept=tpy$x, color='blue') +
      scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
      scale_y_discrete(expand=c(0,0)) +
      scale_fill_gradientn(name='F1 score',colors=cols100v) +
      #scale_fill_viridis(name='normalized eigengene value') +
      scale_color_manual(values=c('black','white')) +
      #facet_wrap(bat_note_mod ~ ., nrow=4) +
      otheme(legend.pos='none', legend.dir='v', legend.title=F, panel.spacing=.1,
             margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
             xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
      theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
      theme(axis.text.y = element_text(size=7)) +
      #theme(axis.text.y = element_markdown(size=7.5)) +
      guides(color=F)
  p
  #}}}
}

metric='f1'
fo = glue("{dirw}/08.{metric}.pdf")
p = plot_model_perf2(ml1, metric=metric)
p %>% ggexport(filename=fo, width=6, height=6)

p = plot_model_perf1(ml1, metric='f1')
saveRDS(p, glue("{dirf}/f.4a.rds"))
p = plot_model_perf1(ml1, metric='roc_auc')
saveRDS(p, glue("{dirf}/f.4b.rds"))

p = plot_model_perf2(ml1, metric='f1')
saveRDS(p, glue("{dirf}/f.4c.rds"))
p = plot_model_perf2(ml1, metric='roc_auc')
saveRDS(p, glue("{dirf}/f.4d.rds"))
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
  group_by(bat,cond,drc,note,qry,tgt,time) %>% nest() %>% ungroup() %>%
  rename(ti=data) %>%
  mutate(acc = map(ti, get_acc2)) %>%
  select(-ti) %>% unnest(acc)

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
#{{{ model performance [old]
readRDS2 <- function(fi) if(file.exists(fi)) readRDS(fi) else tibble()
x = tibble(i=1:54) %>% mutate(fi = sprintf("%s/10_rf/%d.rds", dirw, i)) %>%
  mutate(to = map(fi, readRDS2)) %>%
  unnest(to) %>% select(-fi,-i)

bins=c("-500","+500","+/-500","-2k","+2k","+/-2k")
epis=c("raw",'umr')
tpl = crossing(bin=bins,epi=epis) %>% mutate(bin_epi=str_c(bin,epi,sep=':')) %>%
  mutate(bin=factor(bin,levels=bins), epi=factor(epi,levels=epis)) %>%
  arrange(bin,epi)
tp = x %>% filter(epi != 'pos') %>%
  arrange(bat_mid,bin,epi) %>%
  mutate(bin_epi = str_c(bin, epi, sep=':')) %>%
  mutate(bin_epi = factor(bin_epi, levels=tpl$bin_epi)) %>%
  select(bat_mid, bin, epi, bin_epi, metric) %>% unnest(metric) %>% filter(metric=='f_meas')
tps = tp %>% group_by(bat_mid, bin_epi) %>%
  summarise(y = median(estimate), ymin = quantile(estimate,.25),
            ymax=quantile(estimate, .75)) %>% ungroup()

pv = ggplot(tp) +
    geom_violin(aes(x=bin_epi, y=estimate, fill=epi), trim=T, alpha=.8) +
    #geom_jitter(aes(x=bin_epi, y=estimate), color='darkgray') +
    geom_pointrange(data=tps, aes(x=bin_epi, y=y,ymin=ymin, ymax=ymax)) +
    scale_x_discrete(expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(limits=c(.5,1), expand=expansion(mult=c(0,.01))) +
    scale_fill_manual(values=pal_jco()(8)) +
    facet_wrap(~bat_mid, scales='free', ncol=1) +
    otheme(legend.pos='none', strip.style='white',
           xtext=T, xtick=T, ytext=T, ytick=T, ygrid=T) +
    theme(axis.text.x = element_text(angle = 20, hjust=.8, vjust=1.1)) +
    guides(color=F)
fo = file.path(dirw, 'tt.pdf')
ggsave(pv, filename=fo, width=8, height=8)

ti2 = x %>%
  mutate(nfeature = map_dbl(fit, pluck, 'fit','num.independent.variables')) %>%
  select(bat_mid, bin, epi, nfeature)
to = x %>% select(-metric) %>% rename(metric=metricB) %>%
  select(bat_mid, bin, epi, metric) %>% unnest(metric) %>%
  spread(metric, estimate) %>%
  inner_join(ti2, by=c('bat_mid','bin','epi')) %>%
  select(bat_mid, bin, epi, nfeature, sens, spec, accuracy, precision, f1=f_meas) %>%
  arrange(bat_mid, bin, epi)
fo = file.path(dirw, '11.tsv')
write_tsv(to, fo)
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
