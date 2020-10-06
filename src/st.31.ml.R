source('functions.R')
#require(vip)
dirw = file.path(dird, '31_ml')

#{{{ save for ML analysis
fln = file.path(dirw, '../21_seq/45.rds')
tln = readRDS(fln)

fi = file.path(dirw, '../23_mmd/05.mtf.enrich.rds')
#fi = file.path(dirw, '../23_mmd/07.mtf.grp.enrich.rds')
mer = readRDS(fi)
tf = mer %>% distinct(fid, fname) %>%
  mutate(fname=str_replace_all(fname, ".*\\|", "")) %>%
  mutate(fname=str_replace_all(fname, "\\[.*", "")) %>%
  mutate(fname=str_replace_all(fname, "-", ""))

bat_mids = c("cold_up:m16",'cold_up:m21','heat_up:m18')
ts = tln %>% filter(bat_mid %in% bat_mids) %>%
  select(bat_mid, hit=tg, control=tg_c) %>%
  gather(status, tg, -bat_mid) %>% unnest(tg) %>%
  mutate(status = factor(status, levels=c('hit', 'control'))) %>%
  mutate(status = 2-as.numeric(status)) %>%
  select(bat_mid,sid,status) %>%
  group_by(bat_mid) %>% nest() %>% rename(ts = data) %>% ungroup()

ta = mer %>% filter(bat_mid %in% bat_mids) %>%
  arrange(bat_mid, fid, pval.umr) %>%
  group_by(bat_mid, fid) %>%
  slice(1) %>%
  arrange(bat_mid, pval.umr) %>%
  group_by(bat_mid) %>% mutate(i = 1:n()) %>% ungroup() %>%
  unnest(tl) %>%
  mutate(pos = round((start+end)/2)) %>%
  select(bat_mid, fid, i, sid, kmer, pos, srd, umr) %>%
  #mutate(rPos = ifelse(umr==1, abs(pos-2000), NA)) %>%
  group_by(bat_mid) %>% nest() %>% rename(ti = data)

prepare_ti <- function(ti, bin, epi, nfea, tst) {
  #{{{
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
  to = to %>% distinct(sid, fid) %>% mutate(x = 1)
  fids = to %>% distinct(fid) %>% pull(fid)
  tst %>% crossing(fid=fids) %>%
    left_join(to, by=c('sid','fid')) %>%
    replace_na(list(x=0)) %>%
    spread(fid, x) %>% select(-sid)
  #}}}
}

bins = c("-500","+500","+/-500","-2k",'+2k','+/-2k')
epis = c('raw','umr')
nfeas = c('top30', 'top50', 'top100', 'all')
tops = crossing(bin = bins, epi = epis, nfea = nfeas) %>%
  mutate(nfea = factor(nfea, levels=nfeas)) %>%
  mutate(bin = factor(bin, levels=bins)) %>%
  mutate(epi = factor(epi, levels=epis)) %>%
  arrange(bin, epi, nfea)
ta2 = ta %>% crossing(tops) %>% inner_join(ts, by='bat_mid') %>%
  mutate(to = pmap(list(ti,bin,epi,nfea,ts), prepare_ti))
to = ta2 %>% mutate(did = sprintf("d%03d", 1:n())) %>%
  select(did, bat_mid, bin, epi, nfea, to)

to %>% mutate(fo = sprintf("%s/05_input/%s.tsv", dirw, did)) %>%
  mutate(y = map2(to, fo, write_tsv))

fo = file.path(dirw, "01.tsv")
to2 = to %>% select(-to) %>% crossing(perm=1:10)
write_tsv(to2, fo)
fo = file.path(dirw, "01.rds")
saveRDS(to, fo)
#}}}

fi = file.path(dirw, "01.rds")
r = readRDS(fi)

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

#{{{ obtain balanced F1 score and maximize F1 score during permutation
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

#{{{ model performance using different hyper-parameters
fi = file.path(dirw, '01.tsv')
th = read_tsv(fi)
fi = file.path(dirw, "03.rds")
ti = readRDS(fi) %>% select(sid, metric) %>%
  separate(sid, c("did", "perm"), sep='_') %>%
  mutate(perm = as.double(perm)) %>%
  unnest(metric) %>% spread(metric, estimate) %>%
  inner_join(th, by=c('did','perm'))

bins = c("-500","+500","+/-500","-2k","+2k","+/-2k")
nfeas = c("top30",'top50','top100','all')
tp = ti %>% group_by(bat_mid, bin, epi, nfea) %>%
  summarise(score=median(f_meas)) %>% ungroup() %>%
  mutate(y = str_c(bat_mid, nfea, sep=":")) %>%
  mutate(x = str_c(bin, epi, sep=":")) %>%
  mutate(nfea = factor(nfea, levels=nfeas)) %>%
  mutate(bin = factor(bin, levels=bins)) %>%
  mutate(epi = factor(epi, levels=epis)) %>%
  mutate(lab = number(score, accuracy=.01))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin,y=epi)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    #geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='F1 score',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    facet_grid(nfea ~ bat_mid, switch = "y") +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3), ygrid=T, strip.style='white',
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/04.pdf", dirw)
p %>% ggexport(filename=fo, width=7, height=5)
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
