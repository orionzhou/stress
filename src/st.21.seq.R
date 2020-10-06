source('functions.R')
dirw = file.path(dird, '21_seq')
get_coords <- function(tss, srd, opt) {
    #{{{
    offu=0; offd=0
    if (opt == '-500') {
        offu = 500; offd = 0
    } else if (opt == '+/-500') {
        offu = 500; offd = 500
    } else if (opt == '-1k') {
        offu = 1000; offd = 0
    } else if (opt == '+/-1k') {
        offu = 1000; offd = 1000
    } else if (opt == '-2k') {
        offu = 2000; offd = 0
    } else if (opt == '+/-2k') {
        offu = 2000; offd = 2000
    } else if (opt == '-10k') {
        offu = 10000; offd = 0
    } else if (opt == '+/-10k') {
        offu = 10000; offd = 10000
    } else {
        stop("unknown opt: ", opt)
    }
    start = tss; end = tss + 1
    if (srd == '+') {
        start = tss - offu; end = tss + offd
    } else if (srd == '-') {
        start = tss - offd; end = tss + offu
    } else {
        stop("unknown strand: ", srd)
    }
    list(start=start, end=end)
    #}}}
}

#{{{ create CRE regions
diri = file.path(dirg, 'Zmays_B73', '50_annotation')
fi = file.path(diri, "10.tsv")
ti = read_tsv(fi) %>% filter(etype == 'exon') %>% mutate(start=start-1) %>%
    filter(ttype=='mRNA') %>%
    group_by(gid, tid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
    ungroup() %>%
    arrange(chrom, start, end) %>%
    mutate(tss = ifelse(srd=='-', end, start)) %>%
    mutate(tss0 = ifelse(srd == '-', -tss, tss)) %>%
    arrange(gid, desc(tss0)) %>%
    group_by(gid) %>% slice(1) %>% ungroup() %>% select(-tss0, -tid)

tt = ti %>% crossing(bin = bins) %>%
    mutate(coord=pmap(list(tss,srd,bin), get_coords)) %>%
    mutate(pstart = map_dbl(coord, 'start')) %>%
    mutate(pend = map_dbl(coord, 'end')) %>%
    inner_join(gcfg$chrom[,c("chrom",'size')], by='chrom') %>%
    mutate(pstart = pmax(0, pstart), pend=pmin(pend, size)) %>%
    select(gid, chrom,start,end, srd, tss, bin, pstart, pend) %>%
    arrange(chrom,start,end)
tt %>% mutate(psize=pend-pstart) %>% count(psize)

tp = tt %>% select(chrom,start=pstart,end=pend,srd,gid,bin) %>%
    arrange(chrom,start)
tu = read_regions('umr')
tal = read_regions('acrL')
tae = read_regions('acrE')

tp0 = tp %>% mutate(epi = 'raw')
tp1 = intersect_s(tp, tu) %>% mutate(epi = 'umr')
tp2 = intersect_s(tp, tal) %>% mutate(epi = 'acrL')
tp3 = intersect_s(tp, tae) %>% mutate(epi = 'acrE')

bp_min1 = 50; bp_min = 100
tp = rbind(tp0, tp1, tp2, tp3) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(size=end-start) %>%
    filter(size >= bp_min1)
tps = tp %>% group_by(bin,epi,gid) %>% summarise(size=sum(size)) %>% ungroup() %>%
    filter(size >= bp_min) %>% select(-size)
tp = tp %>% inner_join(tps, by=c('bin','epi','gid')) %>%
    mutate(sid = sprintf("s%06d", 1:n()))
tp %>%
    group_by(bin,epi,gid) %>% summarise(size=sum(size)) %>% ungroup() %>%
    group_by(bin,epi) %>%
    summarise(ng = n(), mean=mean(size),
        q5 = quantile(size, .05), q25 = quantile(size, .25),
        q50 = quantile(size, .5), q75 = quantile(size, .75),
        q95 = quantile(size, .95)) %>% ungroup() %>% print(n=40)

to = tp %>% mutate(score=1) %>% select(chrom,start,end,sid,score,srd)
fo = file.path(dirw, '03.bed')
write_tsv(to, fo, col_names = F)
#bedtools getfasta -fi $ref/10.fasta -tab -bed 03.bed -fo 03.tsv -nameOnly -s
##fasta.py cleanid tmp.fas >03.cre.fas

fs = file.path(dirw, '03.tsv')
ts = read_tsv(fs, col_names=c('sid','seq'))
ts2 = ts %>% mutate(sid = str_replace(sid, '\\(.*\\)', ''))

tl = tp %>% inner_join(ts2, by='sid')
tl1 = tl %>% filter(srd == '+') %>% arrange(bin, epi, gid, chrom, start) %>%
    group_by(bin,epi,gid) %>%
    summarise(size=sum(size), nseg=n(), seq=str_c(seq, collapse='NNNNN')) %>%
    ungroup()
tl2 = tl %>% filter(srd == '-') %>% arrange(bin, epi, gid, chrom, desc(start)) %>%
    group_by(bin,epi,gid) %>%
    summarise(size=sum(size), nseg=n(), seq=str_c(seq, collapse='NNNNN')) %>%
    ungroup()
to = rbind(tl1, tl2) %>% arrange(bin, epi, gid) %>%
    mutate(sid = sprintf("s%06d", 1:n())) %>%
    select(sid, bin, epi, gid, size, nseg, seq)
to %>% count(bin, epi) %>% print(n=28)

to1 = to %>% select(-seq)
to2 = to %>% select(sid, seq)
fo = file.path(dirw, '05.cre.loc.tsv')
write_tsv(to1, fo, col_names = T)
fo = file.path(dirw, '05.cre.seq.tsv')
write_tsv(to2, fo, col_names = F)
#awk '{print ">"$1"\n"$2}' 05.cre.seq.tsv > 05.cre.fas
#}}}


#{{{ create cluster-CRE lists
fi = file.path(dirw, '05.cre.loc.tsv')
ti = read_tsv(fi, col_types='ccccii')

#{{{ read picked modules and gids
fi = file.path(dird, "17_cluster/25.modules.rds")
x = readRDS(fi)
#
ff = file.path(dird, '17_cluster/config.xlsx')
tf = read_xlsx(ff, sheet='picked') %>%
    fill(opt_deg, .direction='down') %>%
    fill(opt_clu, .direction='down') %>%
    fill(bat, .direction='down') %>% mutate(pick = T)
#
tp = x %>% select(-toc) %>% unnest(tom) %>%
    select(-me) %>%
    left_join(tf, by=c("opt_deg","opt_clu","bat","mid")) %>%
    replace_na(list(pick=F)) %>%
    filter(opt_deg == 'B', opt_clu == 'B') %>% select(-opt_deg, -opt_clu) %>%
    select(-stress,-drc) %>% rename(ng0 = n, gid = gids) %>% unnest(gid) %>%
    separate(gid, c("gid",'gt'), sep='_') %>% select(-gt)
#}}}

#{{{ prepare sids for each module
min_ng = 50
tl = ti %>% select(gid, bin, epi, size, sid) %>%
    inner_join(tp, by='gid') %>%
    group_by(bat, mid, ng0, note, bin, epi, pick) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(bat,mid,ng0,note, bin,epi, ng, pick, tg = data) %>%
    mutate(bat = factor(bat, levels=bats)) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    arrange(bat, mid, bin, epi) %>%
    mutate(bat_mid = fct_cross(bat, mid, sep=":")) %>%
    mutate(bin_epi = fct_cross(bin, epi, sep=":"))
y1 =  tl %>% distinct(bat,mid) %>% mutate(bat_mid=str_c(bat,mid,sep=":"))
y2 =  tl %>% distinct(bin,epi) %>% mutate(bin_epi=str_c(bin,epi,sep=":"))
tl = tl %>% mutate(bat_mid = factor(bat_mid, levels=y1$bat_mid)) %>%
    mutate(bin_epi = factor(bin_epi, levels=y2$bin_epi)) %>%
    mutate(lid = sprintf("l%04d", 1:n())) %>%
    filter(ng >= min_ng) %>%
    select(lid, bat,mid,ng0,note,bin,epi,bat_mid,bin_epi,ng,pick,tg)
#}}}

#{{{ plot #genes per cluster-CRE
tp = tl %>% filter(pick) %>%
    select(bat_mid, bin_epi, ng) %>%
    mutate(bat_mid = factor(bat_mid, levels=rev(y1$bat_mid)))
#
swit = (min(tl$ng) + max(tl$ng)) / 2
swit = 0
p = ggplot(tp, aes(x=bin_epi,y=bat_mid)) +
    geom_tile(aes(fill=ng)) +
    geom_text(aes(label=ng, color=ng>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='number genes in cluster',colors=cols100) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    #facet_wrap(~, nrow=2, dir='v', scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.3,1.9,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7.5)) +
    theme(axis.text.y = element_text(size=7.5)) +
    guides(color = F, fill=F)
#
fo = sprintf("%s/07.ngene.pdf", dirw)
p %>% ggexport(filename = fo, width = 8, height = 6)
#}}}

fo = file.path(dirw, '10.rds')
saveRDS(tl, fo)
#}}}

#{{{ create background-CRE list
fi = file.path(dirw, '../15_de/05.rds')
x = readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12

tb = deg48 %>% select(-up, -down) %>% unnest(ds) %>%
    filter(Genotype %in% gts3) %>%
    mutate(deg = padj < .05 & abs(log2fc) >= 1) %>%
    group_by(gid, Treatment, Genotype) %>%
    summarise(n_deg = sum(deg)) %>% ungroup()
tb2 = tb %>% group_by(gid, Treatment) %>%
    summarise(n_deg = sum(n_deg)) %>% ungroup()
tb3 = tb2 %>% filter(n_deg == 0) %>%
    select(-n_deg) %>% mutate(cond = str_to_lower(Treatment)) %>%
    select(-Treatment)

tc = ti %>% select(gid, bin, epi, size, sid) %>%
    inner_join(tb3, by='gid') %>%
    group_by(cond, bin, epi) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(cond, bin,epi, ng, tg = data) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    arrange(cond, bin, epi) %>%
    mutate(bin_epi = fct_cross(bin, epi, sep=":"))
y2 =  tc %>% distinct(bin,epi) %>% mutate(bin_epi=str_c(bin,epi,sep=":"))
tc = tc %>%
    mutate(bin_epi = factor(bin_epi, levels=y2$bin_epi)) %>%
    mutate(lid = sprintf("c%02d", 1:n())) %>%
    select(lid, cond, bin,epi,bin_epi,ng,tg)

tc1 = tc %>% mutate(fo = sprintf("%s/11_bg_lists/%s.txt", dirw, lid)) %>%
    mutate(sids = map(tg, 'sid')) %>%
    mutate(j = map2(sids, fo, write))

to = tc %>% select(-tg)
fo = file.path(dirw, '11.bg.tsv')
write_tsv(to, fo)
saveRDS(tc, fo)
#}}}

#{{{ create test-control CRE list pairs
fi = file.path(dirw, '10.rds')
tl0 = readRDS(fi)
fc = file.path(dirw, '11.bg.rds')
tc = readRDS(fc)

tc1 = tc %>% select(clid = lid, cond, bin_epi, ng_c=ng,tg_c=tg)
tl = tl0 %>% separate(bat, c("cond",'drc'), sep='_', remove=F) %>%
    select(-drc) %>% inner_join(tc1, by=c('cond','bin_epi'))

tl %>% filter(pick) %>%
    mutate(fo = sprintf("%s/15_lists/%s.txt", dirw, lid)) %>%
    mutate(sids = map(tg, 'sid')) %>%
    mutate(j = map2(sids, fo, write))

to = tl %>% select(-tg,-tg_c)
fo = file.path(dirw, '15.tsv')
write_tsv(to, fo)
fo = file.path(dirw, '15.rds')
saveRDS(tl, fo)
to = tl %>% filter(pick) %>% select(-tg, -pick, -tg_c)
fo = file.path(dirw, '15.picked.tsv')
write_tsv(to, fo)
#}}}


#---------- +/-2k TSS search -------#
#{{{ create +/-2kb promoter sequences
diri = file.path(dirg, 'Zmays_B73', '50_annotation')
fi = file.path(diri, "10.tsv")
ti = read_tsv(fi) %>% filter(etype == 'exon') %>% mutate(start=start-1) %>%
    filter(ttype=='mRNA') %>%
    group_by(gid, tid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
    ungroup() %>%
    arrange(chrom, start, end) %>%
    mutate(tss = ifelse(srd=='-', end, start)) %>%
    mutate(tss0 = ifelse(srd == '-', -tss, tss)) %>%
    arrange(gid, desc(tss0)) %>%
    group_by(gid) %>% slice(1) %>% ungroup() %>% select(-tss0, -tid) %>%
    filter(!str_detect(gid, "^(Zeam|Zema)"))

tt = ti %>% crossing(bin = "+/-2k") %>%
    mutate(coord=pmap(list(tss,srd,bin), get_coords)) %>%
    mutate(pstart = map_dbl(coord, 'start')) %>%
    mutate(pend = map_dbl(coord, 'end')) %>%
    inner_join(gcfg$chrom[,c("chrom",'size')], by='chrom') %>%
    mutate(pstart = pmax(0, pstart), pend=pmin(pend, size)) %>%
    select(gid, chrom,start,end, srd, tss, bin, pstart, pend) %>%
    arrange(chrom,start,end)
tt %>% mutate(psize=pend-pstart) %>% count(psize)

bp_min1 = 50; bp_min = 100
tp = tt %>% select(chrom,start=pstart,end=pend,srd,gid,bin) %>%
    arrange(chrom,start) %>%
    mutate(size=end-start) %>%
    filter(size >= bp_min1) %>%
    mutate(sid = sprintf("s%05d", 1:n()))

to = tp %>% mutate(score=1) %>% select(chrom,start,end,sid,score,srd)
fo = file.path(dirw, '41.bed')
write_tsv(to, fo, col_names = F)
#bedtools getfasta -fi $ref/10.fasta -bed 41.bed -fo tmp.fas -nameOnly -s
#fasta.py cleanid tmp.fas >41.fas

fo = file.path(dirw, '41.tsv')
write_tsv(tp, fo, col_names = T)

fo = file.path(dirw, '42.bed')
to = tp %>% mutate(qbeg=0, cid=1:n()) %>%
    select(chrom,start,end,srd,sid,qbeg,qend=size,cid)
write_tsv(to, fo, col_names = F)
# fasta.py size 41.fas > 41.sizes
# chain.py fromBed 42.bed $ref/15_intervals/01.chrom.sizes 41.sizes >tmp.chain
# chainSwap tmp.chain 42.chain
#}}}

#{{{ read picked modules and gids
fi = file.path(dird, "17_cluster/25.modules.rds")
x = readRDS(fi)
#
ff = file.path(dird, '17_cluster/config.xlsx')
tf = read_xlsx(ff, sheet='picked') %>%
    fill(opt_deg, .direction='down') %>%
    fill(opt_clu, .direction='down') %>%
    fill(bat, .direction='down') %>% mutate(pick = T)
#
tg = x %>% select(-toc) %>% unnest(tom) %>%
    select(-me) %>%
    left_join(tf, by=c("opt_deg","opt_clu","bat","mid")) %>%
    replace_na(list(pick=F)) %>%
    filter(opt_deg == 'B', opt_clu == 'B') %>% select(-opt_deg, -opt_clu) %>%
    select(-stress,-drc) %>% rename(ng0 = n, gid = gids) %>% unnest(gid) %>%
    separate(gid, c("gid",'gt'), sep='_') %>% select(-gt)
#}}}

#{{{ create background-CRE list
fi = file.path(dirw, '../15_de/05.rds')
x = readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12

tb = deg48 %>% select(-up, -down) %>% unnest(ds) %>%
    filter(Genotype %in% gts3) %>%
    mutate(deg = padj < .05 & abs(log2fc) >= 1) %>%
    group_by(gid, Treatment, Genotype) %>%
    summarise(n_deg = sum(deg)) %>% ungroup()
tb2 = tb %>% group_by(gid, Treatment) %>%
    summarise(n_deg = sum(n_deg)) %>% ungroup()
tb3 = tb2 %>% filter(n_deg == 0) %>%
    select(-n_deg) %>% mutate(cond = str_to_lower(Treatment)) %>%
    select(-Treatment)

tc = ti %>% select(gid, size, sid) %>%
    inner_join(tb3, by='gid') %>%
    group_by(cond) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(cond, ng, tg = data)
#}}}

#{{{ prepare sids for each module
fi = file.path(dirw, '41.tsv')
ti = read_tsv(fi)
min_ng = 50
tl = ti %>% select(gid, bin, size, sid) %>%
    inner_join(tg, by='gid') %>%
    group_by(bat, mid, ng0, note, bin, pick) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(bat,mid,ng0,note, bin, ng, pick, tg = data) %>%
    mutate(bat=factor(bat, levels=bats)) %>%
    #mutate(bat_mid = fct_cross(bat, mid, sep=":")) %>%
    mutate(bat_mid=str_c(bat,mid,sep=":")) %>%
    filter(ng >= min_ng) %>%
    select(bat_mid,pick,note,bat,mid,bin,ng0,ng,tg)
y1 = tl %>% distinct(bat,mid,bat_mid) %>% mutate(bat=factor(bat, levels=bats)) %>%
    arrange(bat, mid)
tl = tl %>% mutate(bat_mid = factor(bat_mid, levels=y1$bat_mid)) %>%
    arrange(bat_mid) %>%
    select(-bin)

tc1 = tc %>% select(cond, ng_c=ng,tg_c=tg)
tl2 = tl %>% separate(bat, c("cond",'drc'), sep='_', remove=F) %>%
    select(-drc) %>% inner_join(tc1, by=c('cond'))

fo = file.path(dirw, '45.rds')
saveRDS(tl2, fo)
#}}}
