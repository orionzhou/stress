source('functions.R')
dirw = file.path(dird, '31_promoter')
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
xref = read_xref()

gt = 'Zmays_B73'
gt = 'Zmays_Mo17'
gt = 'Zmays_W22'
#{{{ prepare promoter db for B/M/W/
diro = file.path(dirw, gt)
if(!dir.exists(diro)) system(sprintf("mkdir -p %s", diro))

#{{{ B
tss = get_tss(gt)
#}}}
#{{{ M
tssM = get_tss(gt)
tss = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,pos,srd)
#}}}
#{{{ W
tssW = get_tss(gt)
tss = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssW, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,pos,srd)
#}}}

tt = tss %>% crossing(bin="+/-2k") %>%
    mutate(coord=pmap(list(pos,srd,bin), get_coords)) %>%
    mutate(pstart = map_dbl(coord, 'start')) %>%
    mutate(pend = map_dbl(coord, 'end')) %>%
    inner_join(read_chrom_size(gt), by='chrom') %>%
    mutate(pstart = pmax(0, pstart), pend=pmin(pend, size)) %>%
    select(chrom, start=pstart,end=pend, srd, gid, bin) %>%
    arrange(chrom,start,end)
#
bp_min1 = 50; bp_min = 100
tp = tt %>%
    mutate(size=end-start) %>%
    filter(size >= bp_min1)

to = tp %>% mutate(score='.') %>% select(chrom,start,end,gid,score,srd)
fo = file.path(diro, '01.bed')
write_tsv(to, fo, col_names = F)
# run readme.sh
#}}}

#{{{ prepare module list
#{{{ read in
fi = file.path(dird, "17_cluster/27.modules.rds")
md0 = readRDS(fi)
#
fi = file.path(dirw, '../15_de/05.rds')
x = readRDS(fi)
deg48 = x$deg48; deg12 = x$deg12
#}}}

gt = 'Zmays_B73'
gt = 'Zmays_Mo17'
gt = 'Zmays_W22'

md1 = md0 %>% filter(gt==!!gt, pick) %>% select(-gt) %>%
    unnest(gids) %>% rename(gid=gids)

fp = file.path(dirw, gt, '01.bed')
tp = read_tsv(fp, col_names=c("chrom",'start','end','gid','score','srd')) %>%
    select(gid,chrom,start,end,srd)
#{{{ create background-CRE list
tb = deg48 %>% filter(cond2=='timeM', Genotype==str_replace(gt,"^Zmays_",'')) %>%
    select(-up, -down) %>% unnest(ds)
tb1 = tb %>% mutate(nde = padj > .05 | abs(log2fc) <= log2(1.5)) %>% mutate(ctag='c1')
tb2 = tb %>% mutate(nde = padj > .05 & abs(log2fc) <= log2(1.5)) %>% mutate(ctag='c2')
tb = tb1 %>% bind_rows(tb2) %>%
    group_by(ctag, gid, Treatment) %>%
    summarise(n_nde = sum(nde)) %>% ungroup() %>%
    filter(n_nde == 2) %>%
    select(-n_nde) %>% mutate(cond = str_to_lower(Treatment)) %>%
    select(-Treatment)
tb %>% count(ctag, cond)
#
tc = tp %>% select(gid, chrom, start, end, srd) %>%
    inner_join(tb, by='gid') %>%
    group_by(ctag, cond) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(ctag, cond, ng, tg = data)
#}}}

#{{{ prepare gids for each module
min_ng = 50
tl = md1 %>%
    inner_join(tp, by='gid') %>%
    group_by(bat, mid, ng0, note, pick) %>%
    nest() %>% ungroup() %>%
    mutate(ng = map_int(data, xf <- function(x) length(unique(x$gid)))) %>%
    select(bat,mid,ng0,note,  ng, pick, tg = data) %>%
    mutate(bat=factor(bat, levels=bats)) %>%
    #mutate(bat_mid = fct_cross(bat, mid, sep=":")) %>%
    mutate(bat_mid=str_c(bat,mid,sep=":")) %>%
    filter(ng >= min_ng) %>%
    select(bat_mid,pick,note,bat,mid,ng0,ng,tg)
y1 = tl %>% distinct(bat,mid,bat_mid) %>% mutate(bat=factor(bat, levels=bats)) %>%
    arrange(bat, mid)
tl = tl %>% mutate(bat_mid = factor(bat_mid, levels=y1$bat_mid)) %>%
    arrange(bat_mid)
#
tc1 = tc %>% select(ctag, cond, ng_c=ng,tg_c=tg)
tl2 = tl %>% separate(bat, c("cond",'drc'), sep='_', remove=F) %>%
    select(-pick) %>% inner_join(tc1, by=c('cond'))
md = tl2

fo = file.path(dirw, gt, '15.module.rds')
saveRDS(md, fo)
#}}}
#}}}

#{{{ # check liftover TSSs vs. ATGs in M and W
fi = file.path(dirw, 'Zmays_B73/01.bed')
ti = read_tsv(fi, col_names=c('chrom','start','end','gid','score','srd'))

to = ti %>% select(chrom,start,end,gid)
fo = file.path(dirw, '01.promoter.bed')
write_tsv(to, fo, col_names=F)

to = ti %>% mutate(start = start+2000, end = start+1) %>% select(chrom,start,end,gid)
fo = file.path(dirw, '01.tss.bed')
write_tsv(to, fo, col_names=F)

fi = file.path(dirw, '05.tss.BtoM.rds')
fi = file.path(dirw, '05.tss.BtoW.rds')
x = readRDS(fi)

x1 = x %>% filter(n_block>0) %>% select(id,coord) %>% unnest(coord) %>%
    select(gid1=id,chrom2,beg2,end2)
x %>% count(size2, aln, n_block)
x2 = x %>% filter(n_block>0) %>% select(id,vnt) %>% unnest(vnt)
x2 %>% print(width=Inf)

tssB = get_tss('Zmays_B73')
tssM = get_tss('Zmays_Mo17')
tssW = get_tss('Zmays_W22')

xref = read_xref()
tx = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type)
tx = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type)

tx2 = tx %>% inner_join(x1, by='gid1') %>%
    #inner_join(tssW, by=c('gid2'='gid')) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    mutate(dist = ifelse(chrom==chrom2, abs(end2-pos), -1)) %>%
    mutate(opt = ifelse(dist==-1, '?', ifelse(dist <= 5e3, 'good', 'bad')))
tx2 %>% count(type, opt)
skim(tx2 %>% filter(opt=='good') %>% pull(dist))
#}}}


