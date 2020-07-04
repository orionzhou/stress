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

fo = file.path(dirw, '02.cre.tsv')
write_tsv(tp, fo, col_names = T)

to = tp %>% mutate(score=1) %>% select(chrom,start,end,sid,score,srd)
fo = file.path(dirw, '03.cre.bed')
write_tsv(to, fo, col_names = F)
#bedtools getfasta -fi $ref/10_genome.fna -bed 03.cre.bed -fo tmp.fas -nameOnly -s
#fasta.py cleanid tmp.fas >03.cre.fas
#}}}

#{{{ create cluster-CRE lists
fi = file.path(dirw, '02.cre.tsv')
ti = read_tsv(fi)

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

tl1 = tl %>% mutate(fo = sprintf("%s/10_lists/%s.txt", dirw, lid)) %>%
    mutate(sids = map(tg, 'sid')) %>%
    mutate(j = map2(sids, fo, write))

to = tl %>% select(-tg)
fo = file.path(dirw, '10.tsv')
write_tsv(to, fo)
fo = file.path(dirw, '10.rds')
saveRDS(tl, fo)
#}}}




