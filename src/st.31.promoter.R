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
fi = glue("{dird}/21_seq/regions.xlsx")
tr = read_xlsx(fi) %>% filter(offu==2000,offd==2000) %>% mutate(bin=as_factor(bin))
bins = levels(tr$bin)
setwd(diro)

#{{{ B
tss = get_tss_tts(gt)
#}}}
#{{{ M
tssM = get_tss_tts(gt)
tss = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd)
#}}}
#{{{ W
tssW = get_tss_tts(gt)
tss = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssW, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd)
#}}}

#{{{ prepare tl
chrom_size = read_chrom_size(gt)
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
system(glue("bedtools getfasta -fi $genome/data/{gt}/10.fasta -tab -bed 01.bed -fo 01.tsv -nameOnly -s"))
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

tb0 = tl %>% distinct(gid) %>% mutate(cid = 1:n())
tb = tl %>% inner_join(tb0, by='gid') %>%
    mutate(start2=0, end2=end-start) %>%
    mutate(start2 = ifelse(i==2, start2+4050, start2)) %>%
    mutate(end2=ifelse(i==2, end2+4050, end2)) %>%
    arrange(chrom,start,end) %>%
    select(chrom,start,end, srd, gid,start2,end2, cid)
write_tsv(tb, '10.bed', col_names = F)
#
system(glue("chain.py fromBed 10.bed $genome/data/{gt}/15_intervals/01.chrom.sizes 02.sizes > 10.reverse.chain"))
system("chainSwap 10.reverse.chain 10.forward.chain")
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

# many genes have 0-bp UTR5
x = md1 %>% inner_join(utr5, by='gid') %>%
    group_by(bat,note) %>% summarise(n_utr0=sum(size.utr5==0), nt=n()) %>%
    ungroup() %>% mutate(pct=n_utr0/nt) %>% print(n=30)
