source('functions.R')
dirw = glue('{dird}/21_seq')
setwd(dirw)
xref = read_xref()

#{{{ prepare promoter db for B/M/W/
fi = glue("{dird}/21_seq/regions.xlsx")
tr = read_xlsx(fi) %>% filter(offu==2000,offd==2000) %>% mutate(bin=as_factor(bin))
bins = levels(tr$bin)

tssB = get_tss_tts('Zmays_B73') %>% mutate(gid=glue("B73_{gid}"))
#{{{ M
tssM = get_tss_tts('Zmays_Mo17')
tssM = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd) %>%
    mutate(gid=glue("Mo17_{gid}"))
#}}}
#{{{ W
tssW = get_tss_tts('Zmays_W22')
tssW = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssW, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd) %>%
    mutate(gid=glue("W22_{gid}"))
#}}}

#{{{ prepare tl
tss = rbind(tssB, tssM, tssW)
chrom_size = read_tsv('00.sizes', col_names=c('chrom','size'))
tt = tss %>% crossing(tr) %>%
    mutate(pos = ifelse(str_detect(bin,"^TTS"), tts, tss)) %>%
    mutate(start = pos, end = pos + 1) %>%
    mutate(start = ifelse(srd=='-', pos-offd, pos-offu)) %>%
    mutate(end = ifelse(srd=='-', pos+offu, pos+offd)) %>%
    inner_join(chrom_size, by='chrom') %>%
    mutate(start=pmax(0, start), end=pmin(end, size)) %>%
    select(chrom, start, end, srd, gid, bin) %>%
    arrange(chrom,start,end) %>%
    mutate(sid = sprintf("s%05d", 1:n()))
tt %>% mutate(size=end-start) %>% count(bin,size)
tl1 = tt %>% filter(srd == '+') %>% arrange(gid, chrom, start) %>%
    group_by(gid) %>% mutate(i = 1:n()) %>% ungroup()
tl2 = tt %>% filter(srd == '-') %>% arrange(gid, chrom, desc(start)) %>%
    group_by(gid) %>% mutate(i = 1:n()) %>% ungroup()
tl = tl1 %>% bind_rows(tl2) %>% arrange(gid, i)
#}}}

#{{{ obtain segment sequence ts
to = tt %>% mutate(score='.') %>% select(chrom,start,end,sid,score,srd)
write_tsv(to, '01.bed', col_names = F)
system(glue("bedtools getfasta -fi 00.fasta -tab -bed 01.bed -fo 01.tsv -nameOnly -s"))
ts = read_tsv('01.tsv', col_names=c('sid','seq')) %>%
    mutate(sid = str_replace(sid, '\\(.*\\)', ''))
#}}}

#{{{ make tl/to
gap = str_c(rep("N",50), collapse='')
tl0 = tl %>% inner_join(ts, by='sid') %>% arrange(gid, i) %>%
    group_by(gid) %>%
    summarise(size=sum(nchar(seq)), nseg=n(), seq=str_c(seq, collapse=gap)) %>%
    ungroup() %>% mutate(size2=map_int(seq, nchar)) %>%
    select(gid, nseg, size, size2, seq)
tl0 %>% count(nseg, size,size2)
#}}}

to = tl0 %>% select(gid, seq)
write_tsv(to, '02.tsv', col_names = F)
system("bioawk -t '{print \">\"$1\"\\n\"$2}' 02.tsv > 02.fas")
system("fasta.py extract 02.fas s1")
system("fasta.py size 02.fas > 02.sizes")
system("fasta-get-markov -dna 02.fas 02.bg")

tb0 = tl %>% distinct(gid) %>% mutate(cid = 1:n())
tb = tl %>% inner_join(tb0, by='gid') %>%
    mutate(start2=0, end2=end-start) %>%
    mutate(start2 = ifelse(i==2, start2+4050, start2)) %>%
    mutate(end2=ifelse(i==2, end2+4050, end2)) %>%
    arrange(chrom,start,end) %>%
    select(chrom,start,end, srd, gid,start2,end2, cid)
write_tsv(tb, '10.bed', col_names = F)
#
system(glue("chain.py fromBed 10.bed 00.sizes 02.sizes > 10.reverse.chain"))
system("chainSwap 10.reverse.chain 10.forward.chain")
#}}}

# create UMR/acrE/acrL BED



#{{{ [obsolete] create CRE regions
fi = glue("{dirw}/regions.xlsx")
tr = read_xlsx(fi) %>% mutate(bin=as_factor(bin))
bins = levels(tr$bin)
epis = c('raw','umr','acrL','acrE')
tssB = get_tss_tts('Zmays_B73')
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

#{{{ prepare location BED
tt = tssB %>% crossing(tr) %>%
    mutate(pos = ifelse(str_detect(bin,"^TTS"), tts, tss)) %>%
    mutate(start = pos, end = pos + 1) %>%
    mutate(start = ifelse(srd=='-', pos-offd, pos-offu)) %>%
    mutate(end = ifelse(srd=='-', pos+offu, pos+offd)) %>%
    inner_join(gcfg$chrom[,c("chrom",'size')], by='chrom') %>%
    mutate(start = pmax(0, start), pend=pmin(end, size)) %>%
    select(chrom, start, end, srd, gid, bin) %>%
    arrange(chrom,start,end)
tt %>% mutate(size=end-start) %>% count(bin,size)
#
tp = tt
tu = read_regions('umr')
tal = read_regions('acrL')
tae = read_regions('acrE')
#
tp0 = tp %>% mutate(epi = 'raw')
tp1 = intersect_s(tp, tu) %>% mutate(epi = 'umr')
tp2 = intersect_s(tp, tal) %>% mutate(epi = 'acrL')
tp3 = intersect_s(tp, tae) %>% mutate(epi = 'acrE')
#
bp_min1 = 50; bp_min = 100
tp = rbind(tp0, tp1, tp2, tp3) %>%
    mutate(bin = factor(bin, levels=bins)) %>%
    mutate(epi = factor(epi, levels=epis)) %>%
    mutate(size=end-start) %>%
    filter(size >= bp_min1)
tps = tp %>% group_by(bin,epi,gid) %>%
    summarise(size=sum(size)) %>% ungroup() %>%
    filter(size >= bp_min) %>% select(-size)
tp = tp %>% inner_join(tps, by=c('bin','epi','gid')) %>%
    mutate(sid = sprintf("s%07d", 1:n()))
tp %>%
    group_by(bin,epi,gid) %>% summarise(size=sum(size)) %>% ungroup() %>%
    group_by(bin,epi) %>%
    summarise(ng = n(), mean=mean(size),
        q5 = quantile(size, .05), q25 = quantile(size, .25),
        q50 = quantile(size, .5), q75 = quantile(size, .75),
        q95 = quantile(size, .95)) %>% ungroup() %>% print(n=40)
#}}}

to = tp %>% mutate(score=1) %>% select(chrom,start,end,sid,score,srd)
fo = glue("{dirw}/03.bed")
write_tsv(to, fo, col_names = F)
system("bedtools getfasta -fi $ref/10.fasta -tab -bed 03.bed -fo 03.tsv -nameOnly -s")
####system("fasta.py cleanid tmp.fas >05.cre.fas")

fs = file.path(dirw, '03.tsv')
ts = read_tsv(fs, col_names=c('sid','seq')) %>%
    mutate(sid = str_replace(sid, '\\(.*\\)', ''))

#{{{ make tl/to
tl = tp %>% inner_join(ts, by='sid')
gap = str_c(rep("N",10), collapse='')
tl1 = tl %>% filter(srd == '+') %>% arrange(bin, epi, gid, chrom, start) %>%
    group_by(bin,epi,gid) %>%
    summarise(size=sum(size), nseg=n(), seq=str_c(seq, collapse=gap)) %>%
    ungroup()
tl2 = tl %>% filter(srd == '-') %>% arrange(bin,epi,gid,chrom,desc(start)) %>%
    group_by(bin,epi,gid) %>%
    summarise(size=sum(size), nseg=n(), seq=str_c(seq, collapse=gap)) %>%
    ungroup()
to = rbind(tl1, tl2) %>% arrange(bin, epi, gid) %>%
    mutate(sid = sprintf("s%06d", 1:n())) %>%
    select(sid, bin, epi, gid, size, nseg, seq)
to %>% count(bin, epi) %>% spread(epi, n) %>% print(n=36)
#}}}

to1 = to %>% select(-seq)
fo = file.path(dirw, '05.cre.loc.rds')
saveRDS(to1, fo)
ts = to %>% select(sid, seq)
fo = file.path(dirw, '05.cre.seq.tsv')
write_tsv(ts, fo, col_names = F)

system("bioawk -t '{print \">\"$1\"\\n\"$2}' 05.cre.seq.tsv > 05.cre.fas")
system("fasta.py extract 05.cre.fas s000001")
system("fasta-get-markov -dna 05.cre.fas 05.cre.bg")
#}}}

######## create pos/neg seq list pairs for motif scanning/mining #########
fs = glue("{dirw}/05.cre.seq.tsv")
ts = read_tsv(fs, col_names=c('sid','seq'))
#
fi = glue('{dirw}/05.cre.loc.rds')
tl = readRDS(fi)
tls0 = tl %>% distinct(bin, epi) %>% arrange(bin, epi)
diro = glue("{dirw}/20_seq_pairs")

##### prepare seq lst pairs
#{{{ prepare all-DEG seq list pairs
tag = 'deg'
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)# %>% filter(gt=='Zmays_B73') %>%
mdp = md %>% select(cid, gid=gids)
mdn = md %>% select(cond,gid=gids_c) %>% distinct(cond,gid)
tls = tls0
#}}}

#{{{ prepare DEG-based module seq list pairs
tag = 'dmod'
fi = glue("{dird}/17_cluster/50_modules/{tag}.rds")
md = readRDS(fi)# %>% filter(gt=='Zmays_B73') %>%
mdp = md %>% select(cid, gid=gids)
mdn = md %>% select(cond,gid=gids_c) %>% distinct(cond,gid)
tls = tls0
#}}}

#{{{ prepare WGCNA-based seq list pairs
tag = 'wgcna'
fi = glue("{dird}/17_cluster/65.modules.ctrl.rds")
i_tibble <- function(gids) tibble(gid=gids)
mdp = md %>% select(cid, gid=gids)
mdn = md %>% select(cond,gid=gids_c) %>% distinct(cond,gid)
#tls = tls0 %>% filter(str_detect(bin, "\\+/\\-"),!epi %in% c('acrL','acrE'))
tls = tls0 %>% filter(!epi %in% c('acrL','acrE'))
#}}}

#{{{ make seq lists and pairs
min_ng = 50
mdp1 = mdp %>% crossing(tls) %>% unnest(gid) %>%
    inner_join(tl, by=c('gid','bin','epi')) %>%
    arrange(cid, bin, epi, gid) %>%
    group_by(cid, bin, epi) %>% nest() %>% ungroup() %>%
    rename(tg = data) %>%
    mutate(ng = map_int(tg, xf <- function(x) length(unique(x$gid)))) %>%
    mutate(lid = sprintf("l%04d", 1:n()))
mdn1 = mdn %>% crossing(tls) %>% unnest(gid) %>%
    inner_join(tl, by=c('gid','bin','epi')) %>%
    arrange(cond, bin, epi, gid) %>%
    group_by(cond, bin, epi) %>% nest() %>% ungroup() %>%
    rename(tg = data) %>%
    mutate(ng = map_int(tg, xf <- function(x) length(unique(x$gid)))) %>%
    rename(tg_c=tg, ng_c=ng) %>%
    mutate(clid = sprintf("cl%02d", 1:n()))
to = md %>% select(cid,cond,ng0=ng,ng0_c=ng_c) %>%
    inner_join(mdp1, by=c('cid')) %>%
    inner_join(mdn1, by=c('cond','bin','epi')) %>%
    select(lid,cid,cond,bin,epi,clid,ng0,ng,ng0_c,ng_c,tg,tg_c) %>%
    filter(ng >= min_ng)
tol = to %>% select(lid, tg) %>%
    bind_rows(to %>% select(lid=clid, tg=tg_c) %>% distinct(lid, tg))
#}}}
#{{{ write seq lists and fas
fo = glue("{diro}/{tag}.tsv")
write_tsv(to %>% select(-tg,-tg_c), fo, na='')
fo = glue("{diro}/{tag}.rds")
saveRDS(to, fo)

tos = tol %>% unnest(tg) %>% select(lid, sid) %>%
    inner_join(ts, by='sid') %>%
    group_by(lid) %>% nest() %>% rename(ts = data) %>% ungroup()
diro1 = glue("{diro}/{tag}_seqlsts")
if(!dir.exists(diro1)) dir.create(diro1)
write_seqlst <- function(ts, fo)
    write_delim(ts %>% mutate(sid=str_c(">",sid)), fo, delim="\n", col_names=F)
tos %>% mutate(fo = glue("{diro1}/{lid}.fas")) %>%
    mutate(j = map2(ts, fo, write_seqlst))
#}}}

#{{{ # plot #genes per cluster-CRE
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
p %>% ggexport(filename = fo, width = 10, height = 6)
#}}}
