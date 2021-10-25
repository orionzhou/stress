source('functions.R')
dirw = glue('{dird}/17_cluster')
#{{{ functions
make_status_tibble0 <- function(tg, tg_c) tg %>% mutate(status=1) %>%
    bind_rows(tg_c %>% mutate(status = 0)) %>%
    mutate(status=factor(status,levels=c(1,0)))
make_status_tibble <- function(gids, gids_c)
    tibble(gid=c(gids, gids_c),
           status=c(rep(1,length(gids)), rep(0, length(gids_c)))) %>%
    mutate(status=factor(status,levels=c(1,0)))
merge_top_ton <- function(top, ton, min_ng=50) {
#{{{
    ton1 = ton %>% select(cond,ng_c=ng,gids_c=gids)
    md = top %>% inner_join(ton1, by=c('cond')) %>%
        filter(ng >= min_ng) %>%
        select(cid, cond, note, ng, ng_c, gids, gids_c) %>%
        mutate(ts = map2(gids, gids_c, make_status_tibble))
    md
#}}}
}
get_nr_gids <- function(gids) {
    #{{{
    sample_w_seed <- function(xs, seed=31) {set.seed(seed); sample(xs, 1)}
    tx = tibble(ggid=gids) %>% separate(ggid, c('gt','gid'), sep='_', remove=F) %>%
        group_by(gid) %>% summarise(ggids = list(ggid)) %>% ungroup() %>%
        arrange(gid) %>% mutate(i = 1:n()) %>%
        mutate(ggid = map2_chr(ggids, i, sample_w_seed))
    tx %>% pull(ggid)
    #}}}
}
#}}}

#{{{ prepare TC matrix
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th1 = th %>% filter(Experiment=='TC') %>%
    select(sid=SampleID, gt=Genotype, cond=Treatment, time=Timepoint) %>%
    mutate(cond = str_to_lower(cond))
th1b = th1 %>% filter(time==0) %>% mutate(cond='cold')
th1c = th1 %>% filter(time==0) %>% mutate(cond='heat')
th2 = rbind(th1, th1b, th1c)# %>% filter(Timepoint != 8)

rr = tm %>% select(gid, sid=SampleID, cpm=CPM) %>%
    inner_join(th2, by ='sid') %>%
    mutate(has = sprintf("h%03d", time*10)) %>%
    select(-sid, -time) %>% spread(has, cpm) %>%
    mutate(h015 = if_else(is.na(h015), (h010+h020)/2, h015)) %>%
    mutate(h030 = if_else(is.na(h030), (h020+h040)/2, h015)) %>%
    mutate(h080 = if_else(is.na(h080), h040*.8+h250*.2, h080))
sum(is.na(rr))
rr %>% dplyr::count(gt,cond)

rd = rr %>% gather(time, rc, -gid, -gt, -cond) %>%
    spread(cond, rc) %>%
    mutate(rdc = cold-control, rdh = heat-control) %>%
    select(gid,gt,time,cold=rdc, heat=rdh) %>%
    gather(cond, rd, -gid,-gt,-time) %>%
    spread(time, rd) %>% arrange(gid, gt, cond)

rc = tm %>% select(gid, sid=SampleID, rc=ReadCount) %>%
    inner_join(th2, by ='sid') %>%
    mutate(has = sprintf("h%03d", time*10)) %>%
    select(-sid, -time) %>% spread(has, rc) %>%
    mutate(h015 = if_else(is.na(h015), (h010+h020)/2, h015)) %>%
    mutate(h030 = if_else(is.na(h030), (h020+h040)/2, h015)) %>%
    mutate(h080 = if_else(is.na(h080), h040*.8+h250*.2, h080))
rc0 = rc %>% filter(cond=='control') %>% select(-cond)
rc1 = rc %>% filter(cond!='control') %>%
    dplyr::rename(h505=h005,h510=h010,h515=h015,h520=h020,h530=h030,
        h540=h040,h580=h080,h750=h250) %>% select(-h000)
rc = rc0 %>% inner_join(rc1, by=c('gid','gt')) %>%
    select(gid,cond,gt, everything())

tc = list(raw=rr, diff=rd, wide_raw=rc)
fo = file.path(dirw, '01.tc.rds')
saveRDS(tc, fo)
#}}}
#run st.job.17.R

#{{{ read DE, TC, config
fi = glue('{dird}/15_de/05.rds')
x =  readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12
conds = c('cold','heat')
drcs = c('up','down')
deg = deg48 %>% filter(Genotype %in% gts3) %>%
    select(Genotype,Treatment,Timepoint,cond2,up,down) %>%
    gather(drc, gids, -Treatment,-Genotype,-Timepoint,-cond2) %>%
    spread(cond2, gids) %>%
    dplyr::rename(gids0 = time0, gids1 = timeM) %>%
    mutate(gids = map2(gids0, gids1, intersect)) %>%
    select(gt=Genotype,cond=Treatment,Timepoint,drc, gids) %>%
    unnest(gids) %>% mutate(gids = glue("{gt}_{gids}")) %>%
    mutate(cond = str_to_lower(cond)) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(drc = factor(drc, levels=drcs)) %>%
    group_by(cond, Timepoint, drc) %>% summarise(gids = list(gids)) %>%
    ungroup()
#
fi = glue('{dirw}/01.tc.rds')
tc = readRDS(fi)
tcr = tc$raw %>% mutate(gid=glue("{gt}_{gid}")) %>% select(-gt)
#}}}

#{{{ sf02-03
fi = glue("{dird}/06_tf_list/10.stress.tf.tsv")
tt = read_tsv(fi) %>% select(cond=stress,gid,name) %>%
    group_by(cond,gid) %>% slice(1) %>% ungroup() %>%
    mutate(name = str_replace(name, '[,/].*', '')) %>% filter(name!='')
ti = deg %>% filter(drc=='up') %>%
    unnest(gids) %>% select(cond,time=Timepoint,gid=gids) %>%
    separate(gid, c('gt','gid'), sep='_') %>% filter(gt=='B73') %>% select(-gt) %>%
    distinct(cond,time,gid) %>% arrange(cond,gid, time) %>%
    group_by(cond,gid) %>% summarise(time=str_c(time,collapse=',')) %>%
    ungroup()
xlab = 'Hours since stress onset'
plot_profile <- function(tp, opt='heat', xprop=1, ytitle=T, leg.pos='bottom.right') {
    #{{{
    col1 = ifelse(opt=='heat', 'red', 'blue')
    nc = ifelse(opt=='heat', 4, 5)
    sc = ifelse(opt=='heat', F, T)
    tpx = tp %>% distinct(time,x) %>% arrange(x)
    p = ggplot(tp) +
    geom_point(aes(x=time, y=value, color=cond), size=1) +
    geom_line(aes(x=time, y=value, group=cond, color=cond)) +
    scale_x_discrete(name=xlab,breaks=tpx$time, labels=tpx$x,
                     expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='Counts Per Million (CPM)', expand=expansion(mult=c(.02,.05))) +
    scale_color_manual(values=c("black",col1)) +
    facet_wrap(~pnl, ncol=nc, scale='free') +
    otheme(xtext=T, xtitle=T, ytext=T, ytick=T, ytitle=!!ytitle, strip.compact=sc,
           panel.border=F, axis=T,
           legend.pos=leg.pos, margin=c(.1,.2,.1,.2))
    if (xprop < 1) {
        ggarrange(p, NULL, nrow=1, ncol=2, widths=c(xprop, 1-xprop))
    } else {
        p
    }
    #}}}
}

#{{{ heat
cond = 'heat'
tt1 = tt %>% filter(cond=='heat', str_detect(name,'hsf') | str_detect(name,'HSF')) %>%
    mutate(name = str_to_upper(name))
ti2 = tt1 %>% left_join(ti, by=c('cond','gid')) %>% replace_na(list(time='unk')) %>%
    filter(time != 'unk' | name %in% glue("HSFTF{c(4,12,13,20)}")) %>%
    print(n=80)

tgrp = tibble(grp=c('1','25','1,25','unk'),
              grp0=LETTERS[1:4], grp1=c("1h_only",'25h_only','1h_and_25h', 'non-DE'))
tx = tibble(time=colnames(tc$raw)[-c(1:3)]) %>% mutate(x=str_sub(time,2)) %>%
    mutate(x = as.numeric(x)/10)
tp = ti2 %>% filter(cond==!!cond) %>% select(-cond) %>%
    inner_join(tc$raw %>% filter(gt=='B73'), by=c('gid')) %>% select(-gt) %>%
    filter(cond %in% c("control", !!cond)) %>% rename(grp=time) %>%
    gather(time, value, -grp, -gid, -name, -cond) %>%
    inner_join(tx, by='time') %>%
    inner_join(tgrp, by='grp') %>% mutate(pnl = name)

a = plot_profile(tp %>% filter(grp0=='A'))
b = plot_profile(tp %>% filter(grp0=='B'), xprop=.75, leg.pos='none', ytitle=F)
c = plot_profile(tp %>% filter(grp0=='C'), xprop=.5, leg.pos='none', ytitle=F)
d = plot_profile(tp %>% filter(grp0=='D'), xprop=1, leg.pos='none', ytitle=F)
fo = glue("{dirf}/sf02.pdf")
fo = glue("{dirw}/15.tf.{cond}.pdf")
ggarrange(a, b, c, d, nrow=4, ncol=1, labels=LETTERS[1:4],
          widths=c(2,2), heights=c(2,1,1,1)) %>%
ggexport(filename=fo, width=8, height=8)
#}}}

#{{{ cold
cond = 'cold'
tt1 = tt %>% filter(cond==!!cond) %>% mutate(name = str_to_upper(name))
ti2 = tt1 %>% inner_join(ti, by=c('cond','gid')) %>% replace_na(list(time='unk')) %>%
    group_by(name) %>% mutate(name = ifelse(row_number() == 1, name, glue("{name}_{row_number()}"))) %>% ungroup() %>%
    print(n=80)

tgrp = tibble(grp=c('1','25','1,25','unk'),
              grp0=LETTERS[1:4], grp1=c("1h_only",'25h_only','1h_and_25h', 'non-DE'))
tx = tibble(time=colnames(tc$raw)[-c(1:3)]) %>% mutate(x=str_sub(time,2)) %>%
    mutate(x = as.numeric(x)/10)
tp = ti2 %>% filter(cond==!!cond) %>% select(-cond) %>%
    inner_join(tc$raw %>% filter(gt=='B73'), by=c('gid')) %>% select(-gt) %>%
    filter(cond %in% c("control", !!cond)) %>% rename(grp=time) %>%
    mutate(cond=factor(cond, levels=c("control", !!cond))) %>%
    gather(time, value, -grp, -gid, -name, -cond) %>%
    inner_join(tgrp, by='grp') %>%
    inner_join(tx, by='time') %>%
    mutate(pnl = name)

a = plot_profile(tp %>% filter(grp0=='A'),opt='cold',xprop=.2,leg.pos='none',ytitle=F)
b = plot_profile(tp %>% filter(grp0=='B'),opt='cold',leg.pos='top.left')
c = plot_profile(tp %>% filter(grp0=='C'),opt='cold',xprop=.6,leg.pos='none',ytitle=F)
fo = glue("{dirw}/15.tf.{cond}.pdf")
ggarrange(a, b, c, nrow=3, ncol=1, labels=LETTERS[1:4],
          widths=c(2,2), heights=c(1,7,1)) %>%
ggexport(filename=fo, width=10, height=11)
#}}}
#}}}

######## create pos-neg gene lists ########
#{{{ degA & degA_nr - BMW DEGs
tag = 'degA'
top = deg %>% unnest(gids) %>%
    group_by(cond, drc) %>%
    summarise(ng=length(unique(gids)), gids=list(unique(gids))) %>% ungroup() %>%
    arrange(cond, drc) %>%
    mutate(cid = glue("c{1:n()}")) %>%
    mutate(note = glue("all {drc}-regulated")) %>%
    select(cid, cond, note, ng, gids)

#{{{ create background/negtavie list
ton = deg48 %>% filter(cond2=='timeM', Genotype %in% gts3) %>%
    mutate(gt = Genotype) %>%# str_c("Zmays", Genotype, sep='_')) %>%
    select(-up, -down, -Genotype) %>%
    mutate(cond = str_to_lower(Treatment)) %>%
    unnest(ds) %>%
    mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>%
    filter(nde) %>%
    count(gt, gid, cond) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>% select(-n_nde) %>%
    count(cond, gid) %>% rename(n_nde = n) %>%
    filter(n_nde == 3) %>% select(-n_nde) %>%
    crossing(gt = gts3) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>% summarise(gids=list(gid)) %>% ungroup() %>%
    mutate(ng = map_int(gids, length)) %>%
    mutate(cond = factor(cond, levels=conds)) %>% arrange(cond) %>%
    select(cond, ng, gids)
#}}}

md = merge_top_ton(top, ton)
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)

top2 = top %>% mutate(gids = map(gids, get_nr_gids)) %>%
    mutate(ng=map_int(gids, length))
ton2 = ton %>% mutate(gids = map(gids, get_nr_gids)) %>%
    mutate(ng=map_int(gids, length))
md2 = merge_top_ton(top2, ton2)
fo = glue("{dirw}/50_modules/{tag}_nr.rds")
saveRDS(md2, fo)
#}}}

#{{{ degB - B DEGs
tag = 'degB'
top = deg %>% unnest(gids) %>% filter(str_detect(gids, '^B')) %>%
    group_by(cond, drc) %>%
    summarise(ng=length(unique(gids)), gids=list(unique(gids))) %>% ungroup() %>%
    arrange(cond, drc) %>%
    mutate(cid = glue("c{1:n()}")) %>%
    mutate(note = glue("all {drc}-regulated")) %>%
    select(cid, cond, note, ng, gids)

#{{{ create background/negtavie list
tb = deg48 %>% filter(cond2=='timeM', Genotype == 'B73') %>%
    rename(gt = Genotype) %>%
    select(-up, -down) %>% unnest(ds)
ton = tb %>%
    mutate(cond = str_to_lower(Treatment)) %>%
    mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>%
    filter(nde) %>%
    count(gt, gid, cond) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>% select(-n_nde) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>% summarise(gids=list(unique(gid))) %>% ungroup() %>%
    mutate(ng = map_int(gids, length)) %>%
    select(cond, ng, gids)
#}}}

md = merge_top_ton(top, ton)
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}

#{{{ dmodA & dmodA_nr - BWM modules
#{{{ check ovlp w. DEG sets  - prepare t_deg
drcs = c('up','down')
t_deg = deg %>%
    mutate(st = 1) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    spread(Timepoint, st) %>% dplyr::rename(t1=`1`, t25=`25`) %>%
    replace_na(list(t1=0,t25=0)) %>% mutate(tag = '') %>%
    mutate(tag = ifelse(t1==1 & t25==0, '1h', tag)) %>%
    mutate(tag = ifelse(t1==0 & t25==1, '25h', tag)) %>%
    mutate(tag = ifelse(t1==1 & t25==1, '1h & 25h', tag)) %>%
    mutate(tag = glue("{drc}: {tag}")) %>%
    arrange(cond, drc, tag) %>% mutate(tag = as_factor(tag)) %>%
    select(cond, lab=tag, gid)
tags = levels(t_deg$lab)
#}}}

#{{{ cold
cond = 'cold'
gids = deg %>% filter(cond==!!cond) %>%
    select(gids) %>% unnest(gids) %>% distinct(gids) %>%
    filter(str_detect(gids, "^[BMW]")) %>% pull(gids)
#
ti1 = tc$diff %>% dplyr::filter(gt %in% gts3) %>%
    mutate(gid = glue("{gt}_{gid}")) %>% select(-gt) %>%
    filter(cond==!!cond, gid %in% gids) %>% select(-cond)
ti2 = filter_expr(ti1, wide=T, min_cpm=.5, num_sam_on=1, pct_sam_on=0, min_var_p=0, transform='norm')

w = run_wgcna(ti2, softPower=1, type='signed hybrid', corFnc='cor', TOM=F,
              TOMType='signed', hclust.opt='ward.D2')
#
rc = make_raw_modules(w$datExpr, w$diss, w$tree, minModuleSize=20, minGap=0,deepSplit=2)

mc = merge_modules(rc, w$datExpr, pre=str_sub(cond,1,1), cutHeight=.15,
    mms=c('c13'='c04'))
nc = 4
pa = plot_me(mc$tu, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
tcr1 = tcr %>% filter(cond %in% c(!!cond,'control'))
pb = plot_cluster(mc$tul, tcr1, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
#{{{ plot ovlp w. DEG sets
tp = mc$tu %>%
    select(cid, gids) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    left_join(t_deg %>% filter(cond==!!cond), by=c('gid')) %>%
    mutate(lab=as.character(lab)) %>%
    replace_na(list(lab='non-DE'))
tag2s = c(tags[1:3], 'non-DE', tags[4:6])
tp2 = tp %>% dplyr::count(cid,lab) %>% dplyr::rename(tag1=cid,tag2=lab) %>%
    mutate(tag1=as_factor(tag1), tag2=factor(tag2, levels=tag2s))
fills = brewer.pal(9, 'Pastel1')
fills = fills[c(1:3,4:6)]
fills = c(brewer.pal(3, 'Reds'), brewer.pal(3, 'Greens'))
#}}}
pc = multi_pie(tp2, ncol=nc, fills=fills, legend.pos='top.center.out',
               legend.dir='h', legend.vjust=-.2) + o_margin(2,.2,.2,.2)
fo = glue("{dirw}/61.{cond}.pdf")
ggarrange(pa, pb, pc, nrow=1, ncol=3, labels=LETTERS[1:3],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename=fo, width=12, height=7)
mc_cold = mc
me_cold = mc$tu %>% mutate(cond = cond)
#}}}

gids_tf = c('Zm00001d031736','Zm00001d021263','Zm00001d033987', 'Zm00001d032923', 'Zm00001d016255')
#{{{ heat
cond = 'heat'
gids = deg %>% filter(cond==!!cond) %>%
    select(gids) %>% unnest(gids) %>% distinct(gids) %>% pull(gids)
#
ti1 = tc$diff %>% filter(gt %in% gts3) %>%
    mutate(gid = glue("{gt}_{gid}")) %>% select(-gt) %>%
    filter(cond==!!cond, gid %in% gids) %>% select(-cond)
ti2 = filter_expr(ti1, wide=T, min_cpm=.5, num_sam_on=1, pct_sam_on=0, min_var_p=0, transform='norm')
#
w = run_wgcna(ti2, softPower=1, type='signed', corFnc='cor', TOM=F,
              TOMType='signed', hclust.opt='ward.D2')
rc = make_raw_modules(w$datExpr, w$diss, w$tree, minModuleSize=20, minGap=0, deepSplit=2)

mc = merge_modules(rc, w$datExpr, pre=str_sub(cond,1,1), cutHeight=.1,
    mms = c('h05'='h01', 'h10'='h02'))
nc = 4
pa = plot_me(mc$tu, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
tcr1 = tcr %>% filter(cond %in% c(!!cond,'control'))
pb = plot_cluster(mc$tul, tcr1, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
#{{{ plot ovlp w. DEG sets
tp = mc$tu %>%
    select(cid, gids) %>% unnest(gids) %>% dplyr::rename(gid=gids) %>%
    left_join(t_deg %>% filter(cond==!!cond), by=c('gid')) %>%
    mutate(lab=as.character(lab)) %>%
    replace_na(list(lab='non-DE'))
tag2s = c(tags[1:3], 'non-DE', tags[4:6])
tp2 = tp %>% dplyr::count(cid,lab) %>% dplyr::rename(tag1=cid,tag2=lab) %>%
    mutate(tag1=as_factor(tag1), tag2=factor(tag2, levels=tag2s))
fills = brewer.pal(9, 'Pastel1')
fills = fills[c(1:3,4:6)]
fills = c(brewer.pal(3, 'Reds'), brewer.pal(3, 'Greens'))
#}}}
pc = multi_pie(tp2, ncol=nc, fills=fills, legend.pos='top.center.out',
               legend.dir='h', legend.vjust=-.2) + o_margin(2,.2,.2,.2)
fo = glue("{dirw}/61.{cond}.pdf")
ggarrange(pa, pb, pc, nrow=1, ncol=3, labels=LETTERS[1:3],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename=fo, width=12, height=7)
me_heat = mc$tu %>% mutate(cond = cond)
mc_heat = mc
mc$tul %>% filter(str_replace(gid, '^.*_', '') %in% gids_tf)
#}}}

fo = glue("{dirw}/12.clustering.rds")
r = list(cold=mc_cold, heat=mc_heat)
saveRDS(r, fo)

f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg)
conds = c('cold','heat')
top = rbind(me_cold, me_heat) %>%
    arrange(cond, cid) %>%
    inner_join(cfg, by=c('cond','cid')) %>%
    filter(!is.na(note)) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    arrange(cond, cid) %>%
    select(cid, cond, note, ng=n, gids, me)

#{{{ create background/negtavie list
tb = deg48 %>% filter(cond2=='timeM', Genotype %in% gts3) %>%
    mutate(gt = Genotype) %>%# str_c("Zmays", Genotype, sep='_')) %>%
    select(-up, -down, -Genotype) %>% unnest(ds)
ton = tb %>%
    mutate(cond = str_to_lower(Treatment)) %>%
    mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>%
    filter(nde) %>%
    count(gt, gid, cond) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>% select(-n_nde) %>%
    count(cond, gid) %>% rename(n_nde = n) %>%
    filter(n_nde == 3) %>% select(-n_nde) %>%
    crossing(gt = gts3) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>% summarise(gids=list(gid)) %>% ungroup() %>%
    mutate(ng = map_int(gids, length)) %>%
    select(cond, ng, gids)
#}}}

md = merge_top_ton(top, ton)
tag = 'dmodA'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)

tag = 'dmodA'
fo = glue("{dirw}/50_modules/{tag}.rds")
md = readRDS(fo)
md2 = md %>%
    mutate(gids = map(gids, get_nr_gids)) %>%
    mutate(gids_c = map(gids_c, get_nr_gids)) %>%
    mutate(ng=map_int(gids, length)) %>%
    mutate(ng_c=map_int(gids_c, length)) %>%
    mutate(ts = map2(gids, gids_c, make_status_tibble))
fo = glue("{dirw}/50_modules/{tag}_nr.rds")
saveRDS(md2, fo)
#}}}

#{{{ dmodB - B modules
fi = glue("{dirw}/50_modules/dmodA.rds")
ti = readRDS(fi)
top = ti %>% select(cid,cond,note,ng,me,gids) %>% unnest(gids) %>%
    filter(str_detect(gids, '^B73')) %>%
    group_by(cid, cond, note) %>%
    summarise(ng=length(unique(gids)), gids=list(unique(gids))) %>% ungroup() %>%
    select(cid, cond, note, ng, gids)

#{{{ create background/negtavie list
tb = deg48 %>% filter(cond2=='timeM', Genotype == 'B73') %>%
    rename(gt = Genotype) %>%
    select(-up, -down) %>% unnest(ds)
ton = tb %>%
    mutate(cond = str_to_lower(Treatment)) %>%
    mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>%
    filter(nde) %>%
    count(gt, gid, cond) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>% select(-n_nde) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>% summarise(gids=list(unique(gid))) %>% ungroup() %>%
    mutate(ng = map_int(gids, length)) %>%
    select(cond, ng, gids)
#}}}

md = merge_top_ton(top, ton)
tag = 'dmodB'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}

#{{{ var1 - variable gene lists
#{{{
conds = c('cold','heat')
drcs = c('up','down')
deg0 = deg48 %>% filter(Genotype %in% gts3) %>%
    select(Genotype,Treatment,Timepoint,cond2,up,down) %>%
    gather(drc, gids, -Treatment,-Genotype,-Timepoint,-cond2) %>%
    spread(cond2, gids) %>%
    dplyr::rename(gids0 = time0, gids1 = timeM) %>%
    mutate(gids = map2(gids0, gids1, intersect)) %>%
    select(gt=Genotype,cond=Treatment,Timepoint,drc, gids) %>%
    unnest(gids) %>% dplyr::rename(gid=gids) %>%
    mutate(cond = str_to_lower(cond)) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(drc = factor(drc, levels=drcs)) %>% rename(time=Timepoint)
#}}}
x0 = deg0 %>% mutate(cond=glue("{cond}_{time}h_{drc}")) %>%
    select(gid,cond,gt)
x1 = x0 %>%
    count(cond, gid) %>% rename(n_de = n) %>%
    filter(!str_detect(cond, 'cold_1h'), n_de < 3)
xa = x1 %>% select(cond, gid) %>% crossing(gt = gts3) %>%
    left_join(x0 %>% mutate(st = 1), by=c('cond','gid','gt')) %>%
    replace_na(list(st=0))

make_status_tibble <- function(gids, gids_c) tibble(gid=c(gids, gids_c), status=c(rep(1,length(gids)), rep(0, length(gids_c)))) %>% mutate(status=factor(status,levels=c(1,0)))
md = xa %>% mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>%
    summarise(ng=sum(st==1), ng_c=sum(st==0),
              gids=list(gid[st==1]), gids_c=list(gid[st==0])) %>%
    ungroup() %>%
    separate(cond, c('cond','note','drc'), sep='_', extra='merge') %>%
    mutate(drc=factor(drc,levels=c('up','down'))) %>%
    arrange(cond,note,drc) %>%
    mutate(cid = c("c10",'c20','c30','c40','c30','c40')) %>%
    mutate(note=glue("{note}_{drc}")) %>%
    select(cid, cond, note, ng, ng_c, gids, gids_c) %>%
    mutate(ts = map2(gids, gids_c, make_status_tibble))

tag = 'var1'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}

#{{{ var2 - variable gene lists 2
fg = glue('{dird}/15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2

#{{{ old way
td3 = td2 %>% filter(st %in% c("dA+B=", "dA=B+", "dA-B=", 'dA=B-')) %>%
    mutate(drc = ifelse(str_detect(st, '-'), 'down', 'up')) %>%
    mutate(qgid = glue('{qry}_{gid}')) %>%
    mutate(tgid = glue('{tgt}_{gid}')) %>%
    mutate(gid1 = ifelse(str_detect(st,"A="), tgid, qgid)) %>%
    mutate(gid0 = ifelse(str_detect(st,"A="), qgid, tgid))
#
td4 = td3 %>% distinct(cond,drc,gid0,gid1) %>%
    gather(tag, gid, -cond,-drc) %>%
    distinct(cond,drc,gid,tag)
td4s = td4 %>%
    count(cond,drc,gid) %>% filter(n<=1)
td5 = td4 %>% inner_join(td4s, by=c("cond",'drc','gid')) %>%
    group_by(cond,drc,tag) %>% summarise(gids = list(gid)) %>% ungroup() %>%
    spread(tag, gids) %>%
    rename(gids=gid1,gids_c=gid0) %>%
    mutate(ng=map_int(gids, length)) %>%
    mutate(ng_c=map_int(gids_c, length))
#}}}

#{{{ new way
td3 = td2 %>% filter(st %in% c("dA+B=", "dA=B+", "dA-B=", 'dA=B-')) %>%
    mutate(drc = ifelse(str_detect(st, '-'), 'down', 'up')) %>%
    distinct(cond,drc,gid)
td1s = td1 %>% arrange(cond,gt,gid,st) %>%
    group_by(cond,gt,gid) %>% slice(1) %>% ungroup() %>% select(-time)
td4 = td3 %>% crossing(gt = gts3) %>%
    left_join(td1s, by = c('cond','gt','gid'))
td4s = td4 %>% mutate(drc1=str_to_upper(str_sub(drc,1,1))) %>%
    group_by(cond,drc,drc1,gid) %>%
    summarise(nd=sum(st==drc1), nu=sum(st!=drc1)) %>% ungroup() %>%
    filter(nd >0, nu>0) %>% select(cond,drc,drc1,gid)
td5 = td4 %>% inner_join(td4s, by=c('cond','drc','gid')) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond,drc) %>%
    summarise(gids = list(gid[st==drc1]), gids_c=list(gid[st!=drc1])) %>%
    ungroup() %>%
    mutate(ng=map_int(gids, length)) %>%
    mutate(ng_c=map_int(gids_c, length))
#}}}

make_status_tibble <- function(gids, gids_c) tibble(gid=c(gids, gids_c), status=c(rep(1,length(gids)), rep(0, length(gids_c)))) %>% mutate(status=factor(status,levels=c(1,0)))
md = td5 %>% mutate(note = glue("{drc}-regulated")) %>%
    mutate(drc = factor(drc, levels=c("up",'down'))) %>%
    arrange(cond,drc) %>%
    mutate(cid = c("c10",'c20','c30','c40')) %>%
    select(cid, cond, note, ng, ng_c, gids, gids_c) %>%
    mutate(ts = map2(gids, gids_c, make_status_tibble))

tag = 'var2'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}

#{{{ dmodB2 - all B modules
fi = glue("{dirw}/12.clustering.rds")
r = readRDS(fi)
me_cold = r$cold$tu %>% mutate(cond = 'cold')
me_heat = r$heat$tu %>% mutate(cond = 'heat')

conds = c('cold','heat')
top = rbind(me_cold, me_heat) %>%
    unnest(gids) %>% filter(str_detect(gids, '^B')) %>%
    group_by(cond, cid, me) %>%
    summarise(ng=length(unique(gids)), gids=list(unique(gids))) %>% ungroup() %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    mutate(note=cid) %>% arrange(cond, cid) %>%
    select(cid, cond, note, ng, gids, me)

#{{{ create background/negtavie list
tb = deg48 %>% filter(cond2=='timeM', Genotype == 'B73') %>%
    rename(gt = Genotype) %>%
    select(-up, -down) %>% unnest(ds)
ton = tb %>%
    mutate(cond = str_to_lower(Treatment)) %>%
    mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>%
    filter(nde) %>%
    count(gt, gid, cond) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>% select(-n_nde) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    group_by(cond) %>% summarise(gids=list(unique(gid))) %>% ungroup() %>%
    mutate(ng = map_int(gids, length)) %>%
    select(cond, ng, gids)
#}}}

md = merge_top_ton(top, ton, min_ng=1) %>% print(n=50)
tag = 'dmodB2'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}


#{{{ f5 - module eigengene w. model prediction acc.
fi = glue("{dirw}/12.clustering.rds")
r = readRDS(fi)
r1 = rbind(r$cold$tul, r$heat$tul) %>%
    group_by(cid) %>% summarise(ng=length(gid), gids=list(gid)) %>% ungroup()

f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg) %>% filter(!is.na(idx)) %>%
    mutate(pnl = glue("{str_sub(cond,0,1)}{str_sub(drc,0,1)}{idx}")) %>%
    select(pnl, cond, drc, cid) %>%
    mutate(cond = factor(cond)) %>%
    mutate(drc = factor(drc, levels=c('up','down')))
md1 = r1 %>% inner_join(cfg, by='cid') %>%
    #inner_join(acc %>% select(cid,acc=txt), by='cid') %>%
    select(cond,drc,pnl,cid,ng,gids) %>% arrange(cond, drc, pnl) %>%
    mutate(pnl = glue("{pnl}: {ng}")) %>% mutate(pnl = as_factor(pnl))

ti = r1 %>% inner_join(md1 %>% select(-gids), by='cid') %>% select(cid,cond,gids)
ti0 = ti %>% mutate(cond='control')
ti = ti %>% bind_rows(ti0) %>% unnest(gids) %>% rename(gid=gids)
tpa = extract_avg_expr(ti, tcr)

conds3 = c('cold','heat','control')
isum <- function(gids, sep="\n") {
    #{{{
    x = tibble(gid=gids) %>% unnest(gid) %>%
        separate(gid, c('gt','gid'), sep='_') %>%
        count(gt) %>% mutate(txt = glue("{gt}: {n}"))
    str_c(x$txt, collapse=sep)
    #}}}
}
tp = md1 %>%
    mutate(txt = map_chr(gids, isum, sep="\n")) %>%
    select(drc,pnl,cid,txt) %>%
    inner_join(tpa, by='cid') %>%
    mutate(cond = factor(cond, levels=conds3))
tp1 = tp %>% filter(str_detect(pnl, 'c'))
tp2 = tp %>% filter(str_detect(pnl, 'h'))

tit1 = 'Cold clusters: 7 up (cu1-cu7) + 9 down (cd1-cd9)'
tit2 = 'Heat clusters: 6 up (hu1-hu6) + 7 down (hd1-hd7)'
pa1 = plot_avg_expr(tp1, ncol=6, col.opt='c', tit=tit1, strip.compact=T)
pa2 = plot_avg_expr(tp2, ncol=6, col.opt='h', tit=tit2, strip.compact=T, xtitle=T)
p = pa1 + pa2 + plot_layout(nrow=2, heights=c(4.2,4))
fo = glue("{dirw}/22.all.expr.pdf")
ggsave(p, file=fo, width=7, height=8)
fo = glue("{dirw}/22.all.expr.rds")
saveRDS(p, fo)
#}}}


#{{{ share cluster gene IDs
to = md1 %>% mutate(gids = map_chr(gids, str_c, collapse = ','))
fo = glue("{dird}/71_share/19.coex.cluster.tsv")
write_tsv(to, fo)
#}}}

######## obsolete ########
#{{{ [obsolete] WGCNA-based modules (all genes)
gt = 'B73'
#{{{ check ovlp w. DEG sets  - prepare t_deg
drcs = c('up','down')
t_deg = deg48 %>% filter(Genotype==gt) %>%
    unnest(ds) %>% mutate(cond2=glue("{cond2}_{Timepoint}")) %>%
    mutate(cond=str_to_lower(Treatment)) %>%
    mutate(de = padj<.05 & abs(log2fc) >= log2(1.5)) %>%
    select(cond,gid,cond2,de) %>%
    spread(cond2,de) %>%
    mutate(de1 = timeM_1) %>%
    mutate(de25 = timeM_25) %>%
    mutate(tag = 'non-DE') %>%
    mutate(tag = ifelse(de1 & !de25, 'DE: 1h', tag)) %>%
    mutate(tag = ifelse(!de1 & de25, 'DE: 25h', tag)) %>%
    mutate(tag = ifelse(de1 & de25, 'DE: 1h & 25h', tag)) %>%
    arrange(cond, tag) %>% mutate(tag = as_factor(tag)) %>%
    select(cond, lab=tag, gid)
t_deg %>% count(cond, lab)
tags = levels(t_deg$lab)
#}}}

#{{{ cold
cond = 'cold'
ti1 = tc$diff %>% filter(Genotype==gt, str_to_lower(Treatment)==cond) %>%
    select(-Genotype, -Treatment)
ti2 = filter_expr(ti1, wide=T, min_cpm=1, num_sam_on=1, pct_sam_on=0, min_var_p=0, transform='asinh')
#sft = pick_soft_power(ti2)
#
w = run_wgcna(ti2, softPower=12, type='signed', corFnc='cor', TOM=T,
              TOMType='signed', hclust.opt='ward.D2')
#
rc = make_raw_modules(w$datExpr, w$diss, w$tree, minModuleSize=20,
                      deepSplit=3, minGap=0)

mc = merge_modules(rc, w$datExpr, cutHeight=.2, pre=str_sub(cond,1,1))
nc = 6
pa = plot_me(mc$me, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
#{{{ plot ovlp w. DEG sets
tp = mc$me %>%
    select(cid, gids) %>% unnest(gids) %>% rename(gid=gids) %>%
    left_join(t_deg %>% filter(cond==!!cond), by=c('gid')) %>%
    mutate(lab=as.character(lab)) %>%
    replace_na(list(lab='non-DE'))
tag2s = tags#c(tags[1:3], 'non-DE', tags[4:6])
tp2 = tp %>% count(cid,lab) %>% rename(tag1=cid,tag2=lab) %>%
    mutate(tag1=as_factor(tag1), tag2=factor(tag2, levels=tag2s))
fills = brewer.pal(9, 'Pastel1')[c(1:3,9)]
#}}}
pb = multi_pie(tp2, ncol=nc, fills=fills, legend.pos='top.center.out',
               legend.dir='h', legend.vjust=-.2) + o_margin(2,.2,.2,.2)
fo = glue("{dirw}/64.{cond}.pdf")
ggarrange(pa, pb, nrow=1, ncol=2, labels=LETTERS[1:2],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename=fo, width=10, height=8)
me_cold = mc$me %>% mutate(cond = cond)
#}}}

#{{{ heat
cond = 'heat'
ti1 = tc$diff %>% filter(Genotype==gt, str_to_lower(Treatment)==cond) %>%
    select(-Genotype, -Treatment)
ti2 = filter_expr(ti1, wide=T, min_cpm=1, num_sam_on=1, pct_sam_on=0, min_var_p=0, transform='asinh')
#sft = pick_soft_power(ti2)
#
w = run_wgcna(ti2, softPower=20, type='signed', corFnc='cor', TOM=T,
              TOMType='signed', hclust.opt='ward.D2')
#
rc = make_raw_modules(w$datExpr, w$diss, w$tree, minModuleSize=20,
                      deepSplit=2, minGap=0)

mc = merge_modules(rc, w$datExpr, cutHeight=.2, pre=str_sub(cond,1,1))
nc = 5
pa = plot_me(mc$me, ncol=nc, tit=glue("{cond} stress ({nrow(ti2)} genes)"))
#{{{ plot ovlp w. DEG sets
tp = mc$me %>%
    select(cid, gids) %>% unnest(gids) %>% rename(gid=gids) %>%
    left_join(t_deg %>% filter(cond==!!cond), by=c('gid')) %>%
    mutate(lab=as.character(lab)) %>%
    replace_na(list(lab='non-DE'))
tag2s = tags
tp2 = tp %>% count(cid,lab) %>% rename(tag1=cid,tag2=lab) %>%
    mutate(tag1=as_factor(tag1), tag2=factor(tag2, levels=tag2s))
fills = brewer.pal(9, 'Pastel1')[c(1:3,9)]
#}}}
pb = multi_pie(tp2, ncol=nc, fills=fills, legend.pos='top.center.out',
               legend.dir='h', legend.vjust=-.2) + o_margin(2,.2,.2,.2)
fo = glue("{dirw}/64.{cond}.pdf")
ggarrange(pa, pb, nrow=1, ncol=2, labels=LETTERS[1:2],
          widths=c(2,2,2), heights=c(2,2)) %>%
ggexport(filename=fo, width=10, height=8)
me_heat = mc$me %>% mutate(cond = cond)
#}}}

f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg, sheet='picked2')
conds = c('cold','heat')
top = rbind(me_cold, me_heat) %>%
    arrange(cond, cid) %>%
    inner_join(cfg, by=c('cond','cid')) %>%
    filter(!is.na(note)) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    arrange(cond, cid) %>%
    select(cid, cond, note, ng=n, gids)
#
ton = cfg %>% filter(!is.na(control)) %>% select(cond, cid) %>%
    inner_join(rbind(me_cold, me_heat), by=c('cond','cid')) %>%
    select(cond,gids) %>% unnest(gids) %>% distinct(cond,gids) %>%
    group_by(cond) %>% summarise(ng=n(), gids=list(gids)) %>% ungroup()
ton

#{{{ build module gene list pairs
min_ng = 50
ton1 = ton %>% select(cond,ng_c=ng,gids_c=gids)
make_status_tibble <- function(tg, tg_c) tg %>% mutate(status=1) %>% bind_rows(tg_c %>% mutate(status = 0)) %>% mutate(status=factor(status,levels=c(1,0)))
make_status_tibble <- function(gids, gids_c) tibble(gid=c(gids, gids_c), status=c(rep(1,length(gids)), rep(0, length(gids_c)))) %>% mutate(status=factor(status,levels=c(1,0)))
md = top %>% inner_join(ton1, by=c('cond')) %>%
    filter(ng >= min_ng) %>%
    select(cid, cond, note, ng, ng_c, gids, gids_c) %>%
    mutate(ts = map2(gids, gids_c, make_status_tibble))
#}}}

tag = 'wgcna'
fo = glue("{dirw}/50_modules/{tag}.rds")
saveRDS(md, fo)
#}}}
#{{{ [obsolete] old f3a - module eigengene, create module gene IDs
tag = 'dmodA'
fi = glue("{dirw}/50_modules/{tag}.rds")
md = readRDS(fi)

f_cfg = glue('{dirw}/config.xlsx')
cfg = read_xlsx(f_cfg)
md1 = md %>% inner_join(cfg %>% select(cid,idx), by='cid') %>%
    filter(!is.na(idx)) %>%
    select(idx,cid,cond,note, ng,me, gids) %>% arrange(idx) %>%
    mutate(pnl = glue("{cond}: {note} ({ng})")) %>%
    mutate(pnl=as_factor(pnl))

ti = md %>% select(cid,cond,gids)
ti0 = ti %>% mutate(cond='control')
ti = ti %>% bind_rows(ti0) %>% unnest(gids) %>% rename(gid=gids)
tpa = extract_avg_expr(ti, tcr)

conds3 = c('cold','heat','control')
isum <- function(gids, sep="\n") {
    #{{{
    x = tibble(gid=gids) %>% unnest(gid) %>%
        separate(gid, c('gt','gid'), sep='_') %>%
        count(gt) %>% mutate(txt = glue("{gt}: {n}"))
    str_c(x$txt, collapse=sep)
    #}}}
}
tp = md1 %>%
    mutate(txt = map_chr(gids, isum)) %>%
    select(cid,pnl,txt) %>%
    inner_join(tpa, by='cid') %>%
    mutate(cond = factor(cond, levels=conds3))
#{{{ plot
tpx = tp %>% distinct(x, xn) %>% arrange(xn)
#tpp = tp %>% distinct(pnl, n) %>% arrange(desc(n))
#tp = tp %>% mutate(pnl=factor(pnl,levels=tpp$pnl))
times = c(0,2,4,8,25)
cols3 = pal_npg()(6)[c(2,1,4)]
tpl1 = tp %>% arrange(cid, desc(y75)) %>%
    filter(cid %in% c("c01",'c05')) %>%
    group_by(cid) %>% slice(1) %>% ungroup() %>%
    mutate(x='h000') %>% select(cid,pnl,txt,x,y=y75)
tpl2 = tp %>% arrange(cid, y25) %>%
    filter(cid %in% c("c04",'c02')) %>%
    group_by(cid) %>% slice(1) %>% ungroup() %>%
    mutate(x='h000') %>% select(cid,pnl,txt,x,y=y25)
tpl3 = tp %>% arrange(cid, desc(y75)) %>%
    filter(str_detect(cid, 'h')) %>%
    group_by(cid) %>% slice(1) %>% ungroup() %>%
    mutate(x='h250') %>% select(cid,pnl,txt,x,y=y75)
pa = ggplot(tp, aes(x=x)) +
    geom_ribbon(aes(ymin=y25,ymax=y75,fill=cond,group=cond), alpha=.2) +
    geom_line(aes(y=y50, color=cond, group=cond)) +
    geom_text(data=tpl1, aes(x=x,y=y, label=txt), size=2, hjust=0, vjust=1) +
    geom_text(data=tpl2, aes(x=x,y=y, label=txt), size=2, hjust=0, vjust=0) +
    geom_text(data=tpl3, aes(x=x,y=y, label=txt), size=2, hjust=1, vjust=1) +
    scale_x_discrete(name="Hours", breaks=tpx$x, labels=tpx$xn, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
    facet_wrap(pnl~., ncol=2, scale='free_y') +
    scale_color_manual(values=cols3) +
    scale_fill_manual(values=cols3) +
    otheme(legend.pos='bottom.right', legend.dir='v',
           panel.spacing=.1, strip.compact=F,
           xtitle=T, xtext=T, xtick=T, margin=c(.2,.2,.2,.2))
#}}}
fo = glue("{dirw}/21.ME.pdf")
ggsave(pa, file=fo, width=4, height=7)
fo = glue("{dirf}/f3a.rds")
saveRDS(pa, fo)
#}}}
#{{{ # [old] wgcna-based DEG clustering
#{{{ test run wgcna
require(WGCNA)
enableWGCNAThreads()

cfg = read_xlsx(f_cfg)
i=9
cid=cfg$cid[i];cond=cfg$cond[i];drc=cfg$drc[i];
opt_deg=cfg$opt_deg[i];opt_clu=cfg$opt_clu[i];optQ=cfg$optQ[i];
softPower=cfg$softPower[i];deepSplit=cfg$deepSplit[i];
MEDissThres=cfg$MEDissThres[i]; minGap=cfg$minGap[i]
r = run_wgcna_pipe(cid,cond,drc,opt_deg,opt_clu,optQ,
    softPower,deepSplit,MEDissThres,minGap, gt_map, tc, deg, dirw)

r$clu %>% filter(str_detect(gid,gid0))
r$me %>% filter(ctag=='raw',clu==14) %>% pull(me)

cfg = read_xlsx(f_cfg)
res = cfg %>% dplyr::filter(cid == 'c09') %>%
    mutate(x=pmap(list(cid,cond,drc,opt_deg,opt_clu,optQ,
                       softPower,deepSplit,MEDissThres,minGap),
                  run_wgcna_pipe, tc=!!tc, deg=!!deg, dirw=!!dirw))

#}}}

#{{{ process wgcna results and save as module
fi = sprintf("%s/12.wgcna.rds", dirw)
res = readRDS(fi)# %>% rename(stress=cond) %>%
    mutate(ng=map_int(x, 'ng'), np=map_int(x, 'np')) %>%
    mutate(nm1=map_int(x, 'nm1'), nm2=map_int(x, 'nm2')) %>%
    mutate(me=map(x, 'me'), clu=map(x,'clu')) %>% select(-x)

#{{{ testing cpm_diff w. cpm_ratio / cpm_raw
tx = res %>% filter(cid >= 'c01', cid <= 'c12') %>%
    select(-clu) %>% unnest(me) %>%
    group_by(stress,drc,opt_deg,opt_clu,softPower,deepSplit,MEDissThres,minGap, ng, np) %>%
    nest() %>% rename(me = data)

plot_me <- function(stress,drc,opt_deg,opt_clu,
                    softPower,deepSplit,MEDissThres,minGap, ng, np, me, dirw) {
    #{{{
    tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
        number(ng), opt_deg, stress, drc, number(np), opt_clu)
    labx = ifelse(minGap==0, 'auto', minGap)
    lab = sprintf("deepSplit=%g minGap=%s", deepSplit, labx)
    #
    me0 = me %>% filter(ctag == 'merged') %>% select(-ctag) %>%
        arrange(cid, optQ, nm1, nm2, n, clu) %>%
        group_by(cid,optQ,nm1,nm2,n,clu) %>%
        mutate(val = scale_by_first(val)) %>%
        ungroup() %>%
        mutate(val = ifelse(val > 1, 1, val)) %>%
        mutate(val = ifelse(val < -1, -1, val))
    tpy0 = me0 %>% distinct(optQ, clu, n)
    tpy = me0 %>% select(optQ,clu,cond,val) %>%
        spread(cond, val) %>% rename(id = clu) %>%
        group_by(optQ) %>% nest() %>%
        mutate(clu = map(data, order_id_by_hclust)) %>%
        select(-data) %>% unnest(clu) %>%
        mutate(y = str_c(optQ, clu, sep="_")) %>%
        inner_join(tpy0, by=c("optQ",'clu')) %>%
        mutate(lab = sprintf("%s (%d)", clu, n))
    optQs = c('raw','diff','ratio')
    tp = me0 %>%
        mutate(y = str_c(optQ, clu, sep="_")) %>%
        mutate(y = factor(y, levels=rev(tpy$y))) %>%
        mutate(optQ = factor(optQ, levels=optQs)) %>%
        mutate(pan=sprintf("cpm [%s]: %s raw -> %d merged", optQ, nm1, nm2))
    pans = tp %>% distinct(optQ, pan) %>% arrange(optQ) %>% pull(pan)
    tp = tp %>% mutate(pan=factor(pan, levels=pans))
    p = ggplot(tp, aes(x=cond,y=y)) +
        geom_tile(aes(fill=val)) +
        #geom_text(aes(label=lab, color=fc>swit), hjust=.5, size=2.5) +
        scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
        scale_y_discrete(breaks=tpy$y, labels=tpy$lab, expand=c(0,0)) +
        scale_fill_gradientn(name='normalized eigengene value',colors=cols100,limits=c(-1,1)) +
        #scale_fill_viridis(name='normalized eigengene value') +
        #scale_color_manual(values=c('black','white')) +
        facet_wrap(~pan, nrow=1, scale='free_y') +
        otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
               margin = c(.2,.1,.2,.1),
               ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
        guides(color = F)
    p = annotate_figure(p, top=text_grob(tit, size=9.5, hjust=.5, color='maroon'),
        fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
    fo = sprintf("%s/17.%s.%s.pdf", dirw, stress, drc)
    p %>% ggexport(filename = fo, width = 10, height = 7)
    T
    #}}}
}

tx %>% mutate(r=pmap_lgl(list(stress,drc,opt_deg,opt_clu,
                              softPower,deepSplit,MEDissThres,minGap, ng,np,me),
    plot_me, dirw=!!dirw))
#}}}

process_wgcna <- function(opt_deg, opt_clu, ti, gt_map, tc, dirw) {
    #{{{
    gts_deg = gt_map[[opt_deg]]; gts_clu = gt_map[[opt_clu]]
    diro = sprintf("%s/20_%s_%s", dirw, opt_deg, opt_clu)
    if(!dir.exists(diro)) dir.create(diro)
    ti1 = ti %>% group_by(bat, stress, drc) %>%
        nest() %>%
        mutate(r = pmap(list(bat, stress, drc, data), process_wgcna_bat,
                        diro, opt_deg, opt_clu)) %>%
        mutate(toc=map(r, 'toc'), tom=map(r, 'tom')) %>%
        select(bat,stress,drc,toc,tom)
    ti1
    #}}}
}
process_wgcna_bat <- function(bat, stress, drc, ti, diro, opt_deg, opt_clu) {
#{{{ different module merging parameters
#{{{ read, filter & preprocess
me = ti %>% select(-clu) %>% unnest(me) %>%
    filter(ctag != 'raw' | (ctag=='raw' & MEDissThres==.1)) %>%
    select(cid,optQ,deepSplit,MEDissThres,minGap,
           np,ng,ctag,n,gids,clu,me)
mex = me %>% distinct(ctag, MEDissThres, clu) %>%
    count(ctag, MEDissThres) %>%
    mutate(tag = sprintf("%s_%s (%s)", ctag, MEDissThres, n)) %>%
    select(-n) %>% arrange(desc(ctag), MEDissThres) %>% mutate(x = 1:n())
#}}}
#
#{{{ process tc
clu = ti %>% select(MEDissThres, clu) %>% mutate(ctag = 'merged') %>%
    unnest(clu)
clu0 = clu %>% filter(MEDissThres==0.1) %>%
    mutate(col2=col1,clu2=clu1,ctag='raw')
tc = clu %>% bind_rows(clu0) %>% rename(clu = clu2, raw_clu = clu1) %>%
    count(ctag, MEDissThres, clu, raw_clu) %>%
    arrange(ctag, MEDissThres, clu, raw_clu) %>%
    group_by(ctag, MEDissThres, clu) %>%
    summarise(n_clu = n(), raw_clus = str_c(raw_clu,collapse='+'), n = sum(n), raw_clu=raw_clu[1]) %>%
    ungroup() %>%
    arrange(desc(ctag), MEDissThres) %>%
    group_by(raw_clus,n_clu,raw_clu,n) %>% nest() %>% ungroup() %>%
    mutate(ctag1 = map_chr(data, mf<-function(df) df$ctag[[1]])) %>%
    mutate(MEDissThres1 = map_dbl(data, mf<-function(df) df$MEDissThres[[1]])) %>%
    mutate(clu = map_dbl(data, mf<-function(df) df$clu[[1]])) %>%
    arrange(n_clu, raw_clu)  %>%
    mutate(mid = sprintf("m%02d", 1:n())) %>%
    select(mid, raw_clus, n_clu, n, ctag1, MEDissThres1, raw_clu, clu, data)
    #%>% print(n=40)
#}}}
#
#{{{ plot cluster merging
tp = tc %>%
    #rename(ctag=ctag1, MEDissThres=MEDissThres1) %>%
    select(mid, raw_clus, n_clu, n, data) %>%
    unnest(data) %>%# filter(! (ctag=='merged' & n_clu == 1)) %>%
    mutate(raw_clus = str_split(raw_clus, '[\\+]')) %>%
    unnest(raw_clus) %>% rename(raw_clu=raw_clus) %>%
    mutate(raw_clu = as.integer(raw_clu)) %>%
    mutate(fcol = ifelse(n_clu==1, 'm01', mid)) %>%
    mutate(lab = sprintf("%s (%d)", mid, n)) %>%
    inner_join(mex, by=c('ctag','MEDissThres'))
tpy = tp %>% filter(ctag=='raw') %>% distinct(raw_clu) %>%
    arrange(desc(raw_clu)) %>% mutate(y=1:n())
tp = tp %>% inner_join(tpy, by='raw_clu')
#
pc = ggplot(tp, aes(x, y)) +
    geom_tile(aes(fill=fcol), color='white', size=.5) +
    geom_text(aes(label=lab, color=fcol), size=2.5) +
    scale_x_continuous(breaks=mex$x, labels=mex$tag, expand=expansion(mult=c(0,0))) +
    scale_y_continuous(breaks=tpy$y, expand=c(0,0)) +
    scale_fill_manual(values=cols36) +
    scale_color_manual(values=brights36) +
    otheme(ygrid=F, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F,
           margin = c(.3,.5,.3,.5)) +
    guides(fill=F,color=F) +
    theme(axis.text.x = element_text(angle=15, hjust=.5, vjust=.8, size=8))
#fo = file.path(dirw, '19.1.pdf')
#ggsave(fo, p, width=6, height=6)
#}}}
#
#{{{ save module IDs,  gene IDs and MEs
toc = tp %>% distinct(ctag, MEDissThres, mid)
tom = tc %>% select(-n) %>% rename(ctag=ctag1,MEDissThres=MEDissThres1) %>%
    inner_join(me, by=c('ctag','MEDissThres','clu')) %>%
    select(mid, n, gids, me)
#}}}
#
#{{{ plot title and label
y = me[1,]
tit = sprintf("Using %s [%s] %s %s DEGs to cluster %s [%s] patterns",
    number(y$ng), opt_deg, stress, drc, number(y$np), opt_clu)
labx = ifelse(y$minGap==0, 'auto', y$minGap)
lab = sprintf("deepSplit=%g minGap=%s", y$deepSplit, labx)
#
tp = tc %>% select(-n) %>% rename(ctag=ctag1,MEDissThres=MEDissThres1) %>%
    mutate(fcol = ifelse(n_clu==1, 'm01', mid)) %>%
    inner_join(me, by=c('ctag','MEDissThres','clu')) %>%
    select(mid, n, clu, fcol, me) %>%
    unnest(me) %>% gather(cond, val, -mid,-n, -clu,-fcol) %>%
    group_by(mid,clu) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    mutate(val = ifelse(val > 1, 1, val)) %>%
    mutate(val = ifelse(val < -1, -1, val))
tpy = tp %>% distinct(mid,n) %>% arrange(desc(mid)) %>%
    mutate(lab=sprintf("%s (%d)", mid, n)) %>% mutate(y=1:n())
tpx = tp %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n())
tp = tp %>% inner_join(tpy[,c('mid','lab','y')], by='mid') %>%
    inner_join(tpx, by='cond')
#{{{ heatmap ph
p = ggplot(tp, aes(x=cond,y=y)) +
    geom_tile(aes(fill=val)) +
    #geom_text(aes(label=lab, color=fc>swit), hjust=.5, size=2.5) +
    scale_x_discrete(expand=expansion(mult=c(0,0))) +
    scale_y_continuous(breaks=tpy$y, labels=tpy$lab, expand=c(0,0)) +
    scale_fill_gradientn(name='normalized eigengene value',colors=cols100,limits=c(-1,1)) +
    #scale_fill_viridis(name='normalized eigengene value') +
    #scale_color_manual(values=c('black','white')) +
    #facet_wrap(~, nrow=2, dir='v', scale='free_y') +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           margin = c(1,.3,.3,.3),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=1, size=8)) +
    guides(color = F)
ph = annotate_figure(p, top=text_grob(tit, size=9, hjust=.5, color='maroon'),
    fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
#fo = sprintf("%s/18.1.pdf", dirw)
#p %>% ggexport(filename = fo, width = 5, height = 6)
#}}}
#{{{ lineplot pl
p = ggplot(tp, aes(x=x,y=val)) +
    geom_line(aes(color=fcol), size=.5) +
    geom_point(color='black', size=.5) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$cond, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(values=cols36) +
    facet_wrap(~lab, ncol=3, dir='v', scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    guides(color = F)
pl = annotate_figure(p, top=text_grob(tit, size=9, hjust=.5, color='maroon'),
    fig.lab=lab, fig.lab.pos='top.left', fig.lab.size=8)
#}}}
#}}}
#
fo = sprintf("%s/%s.pdf", diro, bat)
ggarrange(pc, pl,
    widths=c(1,1.7), heights=c(1), labels=LETTERS[1:6], nrow=1, ncol=2) %>%
ggexport(filename = fo, width = 10,  height = 8)
list(toc=toc, tom=tom)
#}}}
}
tr = res %>% rename(bat=batch) %>% filter(bat %in% bats) %>%
    group_by(opt_deg, opt_clu) %>% nest() %>% ungroup() %>%
    mutate(r = pmap(list(opt_deg, opt_clu, data), process_wgcna,
                    gt_map=!!gt_map, tc=!!tc, dirw=!!dirw))

to = tr %>% select(-data) %>% unnest(r)
fo = file.path(dirw, "25.wgcna.modules.rds")
saveRDS(to, fo)
#}}}

#{{{ plot MEs for selected merging parameters
#{{{ read modules and config
fi = file.path(dirw, "25.wgcna.modules.rds")
x = readRDS(fi)
fi = file.path(dirw, '01.tc.rds')
tc = readRDS(fi)
#
ff = file.path(dirw, 'config.xlsx')
tf = read_xlsx(ff, sheet='picked') %>%
    fill(opt_deg, .direction='down') %>%
    fill(opt_clu, .direction='down') %>%
    fill(bat, .direction='down') %>%
    group_by(opt_deg, opt_clu, bat) %>%
    mutate(i = 1:n()) %>% ungroup() %>% mutate(bat=factor(bat,levels=bats))
tfs = tf %>% count(opt_deg, opt_clu, bat) %>%
    mutate(y0 = c(0,cumsum(n)[1:(n()-1)]+(1:(n()-1))))
tf = tf %>% inner_join(tfs, by=c('opt_deg','opt_clu','bat')) %>%
    mutate(y=y0+i) %>% select(-n,-y0)
#
ti1 = tf %>% group_by(opt_deg, opt_clu,bat) %>% nest() %>% rename(picked=data)
ti = x %>% mutate(bat=factor(bat,levels=bats)) %>%
    inner_join(ti1, by=c("opt_deg",'opt_clu','bat')) %>%
    group_by(opt_deg,opt_clu) %>% nest() %>% ungroup()
#}}}

plot_me <- function(opt_deg, opt_clu, x, tc, dirw) {
    #{{{
    diro = sprintf("%s/20_%s_%s", dirw, opt_deg, opt_clu)
    if(!dir.exists(diro)) dir.create(diro)
    xc = x %>% select(bat, toc) %>% unnest(toc)
    xm = x %>% select(bat, tom) %>% unnest(tom)
    xp = x %>% select(bat, picked) %>% unnest(picked)
#
    #{{{ prop. of genes picked
    xm1 = xm %>% select(bat,gids) %>% unnest(gids) %>%
        distinct(bat, gids) %>% count(bat) %>% rename(ng_tot = n)
    xp1 = xp %>% select(bat, mid) %>% inner_join(xm, by=c('bat','mid')) %>%
        select(bat, gids) %>% unnest(gids) %>%
        distinct(bat, gids) %>% count(bat) %>% rename(ng_pick = n)
    bat_lab = xm1 %>% inner_join(xp1, by=c('bat')) %>% mutate(pg_pick=ng_pick/ng_tot) %>%
        mutate(lab = glue("{bat}\n{ng_pick}/{ng_tot}: {percent(pg_pick, accuracy=.1)}")) %>%
        select(bat, lab)
    #}}}
#{{{ heatmap of all modules
tp = xm %>% select(bat, mid, n, me) %>%
    unnest(me) %>% gather(cond, val, -bat,-mid,-n) %>%
    group_by(bat,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    mutate(bat = factor(bat, levels=bats))
tpx = tp %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n()) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tpy = tp %>% distinct(bat, mid, n) %>%
    mutate(lab=sprintf("%s (%d)", mid, n)) %>%
    mutate(y=str_c(bat, mid, sep=' ')) %>% arrange(bat, desc(mid)) %>% select(-n)
tp = tp %>% inner_join(tpx, by='cond') %>%
    inner_join(tpy, by=c('bat','mid')) %>%
    mutate(y = factor(y, levels=tpy$y))
p01 = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.2) +
    scale_x_continuous(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.01,.01))) +
    scale_y_discrete(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01))) +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    facet_wrap(~bat, ncol=4, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T) +
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '01.heatmap.all.pdf')
ggsave(p01, file=fo, width=10, height=8)
#}}}
#
#{{{ heatmap for picked modules
tp0 = xp %>% inner_join(xm, by=c('bat','mid')) %>%
    select(y,bat,i,mid,n,me,note) %>%
    unnest(me) %>% gather(cond, val, -bat,-i,-y,-mid,-n, -note) %>%
    group_by(bat,i,y,mid,n) %>%
    mutate(val = scale_by_first(val)) %>%
    ungroup() %>%
    #mutate(val = ifelse(val > 1, 1, val)) %>%
    #mutate(val = ifelse(val < -1, -1, val)) %>%
    mutate(lab=sprintf("%s (%d)", note, n)) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(i = factor(i))
tpx = tp0 %>% distinct(cond) %>% arrange(cond) %>% mutate(x=1:n()) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tp0 = tp0 %>% inner_join(tpx, by='cond')
tp = tp0
tpy = tp %>% distinct(y, bat, lab) %>% arrange(y)
tpys = tpy %>% group_by(bat) %>%
    summarise(ymin=min(y), ymax=max(y), ymid=(ymin+ymax)/2) %>% ungroup() %>%
    inner_join(bat_lab, by='bat')
p03 = ggplot(tp, aes(x=x,y=y)) +
    geom_tile(aes(fill=val),color='white',size=.3) +
    geom_segment(data=tpys, aes(x=.3,xend=.3,y=ymin,yend=ymax),color='dodgerblue',size=1) +
    geom_text(data=tpys, aes(x=-1,y=ymid,label=lab),size=3,angle=0,hjust=.5,vjust=.5) +
    scale_x_continuous(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.15,.01))) +
    scale_y_reverse(breaks=tpy$y, labels=tpy$lab, expand=expansion(mult=c(.01,.01)), position='right') +
    scale_fill_gradientn(name='module eigengene', colors=rev(cols100v)) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3),
           ygrid=F, xtick=T, ytick=T, xtitle=T, xtext=T, ytext=T)
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    #theme(strip.text.y = element_blank()) +
    #theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '03.heatmap.picked.pdf')
ggsave(p03, file=fo, width=6, height=6)
#}}}
#
#{{{ multi-line plot
# scale asinh expression
scale2 <- function(ti) {r=max(ti$val)-min(ti$val); ti %>% mutate(val=val/r)}
bat_lab2 = bat_lab %>% rename(blab = lab) %>%
    mutate(blab=str_replace(blab, "\n", " [")) %>%
    mutate(blab=str_c(blab, "]"))
bat_lab2 = bat_lab2 %>% mutate(blab=factor(blab, levels=bat_lab2$blab))
tp = tp0 %>% mutate(x=factor(x)) %>%
    inner_join(bat_lab2, by='bat')
tpx = tp %>% distinct(cond, x) %>% arrange(cond, x) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
get_itvs <- function(ymin, ymax, nl) tibble(yi=1:nl, y=seq(ymin+.05*(ymax-ymin), ymax-.05*(ymax-ymin), length=nl))
tply2 = tp %>% filter(cond=='h250') %>% arrange(blab, val) %>%
    group_by(blab) %>% mutate(yi = 1:n()) %>% ungroup() %>%
    select(blab, mid, yi)
tply = tp %>% select(blab, mid, lab, val) %>%
    group_by(blab) %>%
    summarise(ymin=min(val), ymax=max(val), nl = length(unique(mid))) %>% ungroup() %>%
    mutate(data = pmap(list(ymin,ymax,nl), get_itvs)) %>% unnest(data) %>%
    inner_join(tply2, by=c('blab','yi'))
tpl = tp %>% distinct(blab, mid,i, lab) %>%
    inner_join(tply, by=c('blab','mid'))
#
p06 = ggplot(tp) +
    geom_line(aes(x=x,y=val,group=y,color=as.character(i)), size=.5) +
    geom_point(tp, mapping=aes(x=x,y=val), color='gray35', size=1) +
    geom_text(tpl, mapping=aes(x=8.2,y=y,label=lab,color=as.character(i)), hjust=0, vjust=0, size=2.5) +
    #geom_text_repel(tpl, mapping=aes(x=8.1,y=i*.2,label=lab,color=as.character(i)), hjust=0, vjust=0, size=2.5, direction='y',nudge_x=.5, nudge_y=0, segment.size=.2) +
    scale_x_discrete(name='Hours after Treatment', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.02,.4))) +
    scale_y_continuous(name='Normalized Eigengene Expression',expand=expansion(mult=c(.01,.08))) +
    scale_color_aaas() +
    facet_wrap(~blab, ncol=1, scale='free') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3), strip.style='white',
           ygrid=T, xtick=T, ytick=F, xtitle=T, ytitle=T,xtext=T, ytext=T) +
    theme(strip.text = element_text(hjust=.5, size=9))
fo = file.path(diro, '06.multi.line.pdf')
ggsave(p06, file=fo, width=5, height=6)
#}}}
#
#{{{ write xm for sharing
to = xm %>% mutate(bat=factor(bat, levels=bats)) %>% arrange(bat,mid) %>%
    mutate(gids = map_chr(gids, str_c, collapse=',')) %>%
    mutate(me = map_chr(me, str_c, collapse=','))
fo = file.path(diro, '09.modules.tsv')
write_tsv(to, fo)
#}}}
    list(p01=p01,p03=p03,p06=p06,to=to)
    #}}}
}

#{{{ violin & line plot
# scale asinh expression
scale2 <- function(ti) {r=max(ti$val)-min(ti$val); ti %>% mutate(val=val/r)}
te = tc$diff %>% filter(Genotype %in% gt_map[[opt_clu]]) %>%
    mutate(gid = str_c(gid, Genotype, sep="_")) %>%
    mutate(stress=str_to_lower(Treatment)) %>%
    select(-Genotype, -Treatment) %>%
    gather(has, val, -stress, -gid) %>% mutate(val=asinh(val)) %>%
    group_by(stress, gid) %>% nest() %>% ungroup() %>%
    mutate(ndata = map(data, scale2)) %>%
    select(stress, gid, r=ndata) %>% unnest(r)
#
tp = tp0 %>% mutate(x=factor(x))
tpx = tp %>% distinct(cond, x) %>% arrange(cond, x) %>%
    mutate(xlab=as.double(str_sub(cond,2))/10)
tpl = tp %>% distinct(bat,i,mid,lab)
te1 = tp0 %>% distinct(bat,mid) %>% inner_join(xm, by=c('bat','mid')) %>%
    select(bat, mid, gids) %>% unnest(gids) %>%
    rename(gid=gids) %>%
    #mutate(gid = str_replace(gid, "_.*$", '')) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    separate(bat, c('stress','drc'), remove=F) %>% select(-drc) %>%
    inner_join(te, by=c('stress','gid')) %>%
    rename(cond=has) %>%
    inner_join(tpx, by='cond') %>%
    inner_join(tpl, by=c('bat','mid'))
#
p05 = ggviolin(te1, x='x', y='val', color='grey', fill='gray') +
    #geom_boxplot(te1, mapping=aes(x=x,y=val,group=x), width=.6, color='grey') +
    geom_line(tp, mapping=aes(x=x,y=val,group=y,color=bat), size=.5) +
    geom_point(tp, mapping=aes(x=x,y=val), color='gray35', size=1) +
    geom_text(tpl, mapping=aes(x=1,y=1.25,label=lab), color='royalblue', hjust=0, vjust=0, size=2.5) +
    scale_x_discrete(name='hour after stress', breaks=tpx$x, labels=tpx$xlab, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized eigengene value',expand=expansion(mult=c(.01,.08))) +
    scale_color_aaas() +
    facet_grid(i ~ bat, scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = F, margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=F, xtitle=T, xtext=T, ytext=F) +
    #theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(diro, '05.violin.picked.pdf')
ggsave(p05, file=fo, width=8, height=6)
#}}}

to = ti %>% slice(1:4) %>%
    mutate(r = pmap(list(opt_deg, opt_clu, data), plot_me, tc=!!tc, dirw=!!dirw)) %>%
    mutate(p01=map(r,'p01'), p03=map(r,'p03'), p06=map(r,'p06'), to=map(r,'to'))


saveRDS(to$p06[[1]], file.path(dirf, 'f.2a.rds'))
fo = file.path(dirw, "21.multi.lines.pdf")
ggarrange(to$p06[[1]], to$p06[[2]], to$p06[[3]],
    labels = c("B73", "Mo17", "W22"),
    nrow = 1, ncol = 3, widths=c(1,1,1,1), heights = c(1,1)) %>%
    ggexport(filename = fo, width = 13, height = 6)
fo = file.path(dirf, "sf3.pdf")
ggarrange(to$p06[[2]], to$p06[[3]],
    labels = c("Mo17", "W22"),
    nrow = 1, ncol = 2, widths=c(1,1,1,1), heights = c(1,1)) %>%
ggexport(filename = fo, width = 9, height = 6)
#}}}

#{{{ prepare picked module lists
fi = file.path(dird, "17_cluster/25.wgcna.modules.rds")
md0 = readRDS(fi)
#
ff = file.path(dird, '17_cluster/config.xlsx')
tf = read_xlsx(ff, sheet='picked') %>%
    fill(opt_deg, .direction='down') %>%
    fill(opt_clu, .direction='down') %>%
    fill(bat, .direction='down') %>% mutate(pick = T)
#{{{ create md
gdic = c(B="Zmays_B73",M='Zmays_Mo17',W='Zmays_W22')
tp = md0 %>% select(-toc) %>% unnest(tom) %>%
    select(-me) %>%
    left_join(tf, by=c("opt_deg","opt_clu","bat","mid")) %>%
    replace_na(list(pick=F)) %>%
    mutate(gt = gdic[opt_deg]) %>% select(-opt_deg,-opt_clu) %>%
    select(-stress,-drc) %>% rename(ng0 = n, gid = gids) %>% unnest(gid) %>%
    separate(gid, c("gid",'gt2'), sep='_') %>% select(-gt2) %>%
    select(gt, everything())
tp0 = tp %>% distinct(gt, bat, gid) %>%
    mutate(mid='m00', note='all', pick=T) %>%
    group_by(gt, bat) %>% mutate(ng0 = n()) %>% ungroup()
md = tp %>% bind_rows(tp0) %>%
    group_by(gt, bat, mid, ng0, pick, note) %>%
    summarise(gids = list(gid)) %>% ungroup() %>%
    mutate(bat = factor(bat, levels=bats))
mdp = md %>% filter(gt=='Zmays_B73', pick) %>%
    arrange(bat, mid) %>%
    #arrange(bat, note) %>%
    mutate(bnid = glue("bn{str_pad(1:n(), width=2, side='left',pad='0')}")) %>%
    select(bnid, bat, note)
md = md %>% left_join(mdp, by=c('bat','note')) %>%
    arrange(gt, bnid) %>%
    select(gt, bnid, bat,note,ng0,pick, gids)
#}}}

fo = file.path(dirw, '27.modules.rds')
saveRDS(md, fo)

to = md %>% filter(gt=='Zmays_B73', pick) %>%
    mutate(gids = map_chr(gids, str_c, collapse=',')) %>%
    select(bat, mid, ng0, note, gids)
fo = file.path(dirw, '27.modules.B73.tsv')
write_tsv(to, fo)
#}}}

#{{{ make B/M/W module lists with contols
#{{{ read in
fm = glue("{dird}/17_cluster/27.modules.rds")
md0 = readRDS(fm) %>% filter(!is.na(bnid)) %>% select(-pick)
#
fi = file.path(dirw, '../15_de/05.rds')
x = readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12
#}}}

#{{{ create background-CRE list
tb = deg48 %>% filter(cond2=='timeM', Genotype %in% gts3) %>%
    mutate(gt = str_c("Zmays", Genotype, sep='_')) %>%
    select(-up, -down, -Genotype) %>% unnest(ds)
tb1 = tb %>% mutate(nde = padj > .05 | abs(log2fc) <= log2(1.5)) %>% mutate(ctag='c1')
tb2 = tb %>% mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>% mutate(ctag='c2')
#tb = tb1 %>% bind_rows(tb2) %>%
tb = tb2 %>%
    filter(nde) %>%
    #count(gt, ctag, gid, Treatment) %>% rename(n_nde = n) %>%
    count(gt, gid, Treatment) %>% rename(n_nde = n) %>%
    filter(n_nde == 2) %>%
    select(-n_nde) %>% mutate(cond = str_to_lower(Treatment)) %>%
    select(-Treatment)
tc = tb %>%
    group_by(gt, cond) %>% nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(gt, cond, ng, tg = data)
#}}}

#{{{ prepare gids for each module
min_ng = 50
tl = md0 %>%
    unnest(gids) %>% rename(gid=gids) %>%
    group_by(gt, bnid, bat, note, ng0) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    filter(ng >= min_ng) %>%
    select(gt, bnid, bat, note, ng0, ng, tg=data)
#
tc1 = tc %>% select(gt, cond,  ng_c=ng,tg_c=tg)
make_status_tibble <- function(tg, tg_c) tg %>% mutate(status=1) %>% bind_rows(tg_c %>% mutate(status = 0)) %>% mutate(status=factor(status,levels=c(1,0)))
tl2 = tl %>% separate(bat, c("cond",'drc'), sep='_', remove=F) %>%
    inner_join(tc1, by=c('gt','cond')) %>%
    mutate(ts = map2(tg, tg_c, make_status_tibble)) %>%
    select(-ng0)
md = tl2
#}}}

fo = glue("{dirw}/50.modules.ctrl.rds")
saveRDS(md, fo)

md %>% mutate(fo=glue("{dirw}/{gt}/{bnid}.tsv")) %>%
    mutate(x = map2(ts, fo, write_tsv))
#}}}
#}}}
#{{{ [obsolete] test plot members within individual module
deg1 = deg$deg48 %>% filter(Genotype=='B73',Treatment=='Cold',cond2=='timeM') %>%
    select(Treatment,Timepoint,ds) %>% unnest(ds)
xm1 = xm %>% select(-me) %>% unnest(gids) %>% rename(gid=gids) %>%
    mutate(gid = str_replace(gid, "_.*$", ''))
gid0 = 'Zm00001d002065'
xm1 %>% filter(gid==gid0)
deg1 %>% filter(gid==gid0)

bat = 'cold_up'; mid = 'm20'
gids = xm %>% filter(bat == !!bat, mid == !!mid) %>% unnest(gids) %>%
    mutate(gid = str_replace(gids, "_.*$", '')) %>% pull(gid)

p = ggplot(tp, aes(x=x,y=val)) +
    #geom_line(aes(group=gid),size=.3, color='grey') +
    geom_boxplot(aes(group=x)) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$has, expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='normalized expression value',expand=expansion(mult=c(.05,.05))) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           panel.border = T, margin = c(.3,.3,.3,.3),
           ygrid=T, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(strip.text.y = element_blank()) +
    theme(strip.text.x = element_text(size=10))
fo = file.path(dirw, 't.pdf')
ggsave(p, file=fo, width=8, height=6)
#
#}}}

