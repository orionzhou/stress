source('functions.R')
dirw = glue('{dird}/31_promoter')
xref = read_xref()

tssB = get_tss_tts('Zmays_B73') %>% mutate(gt="B73")
#{{{ M
tssM = get_tss_tts('Zmays_Mo17')
tssM = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd) %>%
    mutate(gt="Mo17")
#}}}
#{{{ W
tssW = get_tss_tts('Zmays_W22')
tssW = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssW, by=c('gid2'='gid')) %>%
    select(gid=gid1,chrom,tss,tts,srd) %>%
    mutate(gt="W22")
#}}}

tss = rbind(tssB, tssM, tssW)

#{{{ 
chrom_size = read_tsv(glue('{dirw}/../21_seq/00.sizes'), col_names=c('chrom','size'))
ti = tss %>% filter(gt=='B73') %>%
    mutate(start = ifelse(srd == '-', tts, tss)) %>%
    mutate(end = ifelse(srd == '-', tss, tts)) %>%
    mutate(start = start - 2000, end = end + 2000) %>%
    inner_join(chrom_size, by='chrom') %>%
    mutate(start=pmax(0, start), end=pmin(end, size)) %>%
    select(chrom, start, end, srd, tss, tts, gid) %>%
    arrange(chrom,start,end)

to = ti %>% select(chrom,start,end,gid)
fo = glue('{dirw}/01.B.bed')
write_tsv(to, fo, col_names=F)


refine_indel <- function(ti, b1,e1,b2,e2) {
    #{{{
    if( is.null(ti) ) {
        NULL
    } else {
    ti %>%
        mutate(beg1 = pmax(beg1, b1)) %>%
        mutate(beg2 = pmax(beg2, b2)) %>%
        mutate(end1 = pmin(end1, e1)) %>%
        mutate(end2 = pmin(end2, e2)) %>%
        mutate(del=end1-beg1, ins=end2-beg2)
    }
    #}}}
}
sum_crossmap <- function(fm, fv) {
    #{{{
    #{{{ map
    ti = read_tsv(fm, col_names=c('chrom1','beg1','end1','id','opt','chrom2','beg2','end2','id2')) %>%
        mutate(size1 = end1-beg1, size2 = end2-beg2) %>%
        select(-id2)
    #
    ti1 = ti %>% filter(opt %in% c("Unmap")) %>%
        select(id,c1=chrom1,b1=beg1,e1=end1,size1)
    ti2 = ti %>% filter(!opt %in% c("Unmap"))
    #
    tm = ti2 %>% mutate(size1=end1-beg1, aln=end2-beg2) %>%
        group_by(id,chrom1,beg1,end1,size1,chrom2) %>%
        summarise(beg2=min(beg2), end2=max(end2), aln=sum(aln)) %>%
        ungroup() %>%
        mutate(size2=end2-beg2) %>%
        arrange(id, desc(aln)) %>%
        group_by(id) %>% slice(1) %>% ungroup()
    #}}}
    #{{{
    tv = read_tsv(fv, col_names=c('chrom0','beg0','end0','id',
                                  'chrom1','beg1','end1','srd',
                                  'chrom2','beg2','end2','cid','vtype',
                                  'ref','alt','bp')) %>%
        select(-chrom0,-beg0,-end0,-bp)
    tvs = tv %>% filter(vtype=='snp') %>%
        select(id,chrom1,pos1=end1,chrom2,pos2=end2,ref,alt) %>%
        group_by(id,chrom1,chrom2) %>% nest() %>% rename(snp=data)
    tvi = tv %>% filter(vtype=='indel') %>% select(-vtype) %>%
        group_by(id,chrom1,chrom2) %>% nest() %>% rename(indel = data)
    #}}}
    isum1 <- function(x) ifelse(is.null(x), 0, sum(x$ins))
    isum2 <- function(x) ifelse(is.null(x), 0, sum(x$del))
    inrow <- function(x) ifelse(is.null(x), 0, nrow(x))
    to = tm %>% left_join(tvs, by=c('id','chrom1','chrom2')) %>%
        left_join(tvi, by=c('id','chrom1','chrom2')) %>%
        mutate(indel=pmap(list(indel,beg1,end1,beg2,end2), refine_indel)) %>%
        mutate(ins = map_dbl(indel, isum1)) %>%
        mutate(del = map_dbl(indel, isum2)) %>%
        mutate(n_snp = map_dbl(snp, inrow)) %>%
        select(id,size1,size2,aln,ins,del,n_snp,
               c1=chrom1,b1=beg1,e1=end1,c2=chrom2,b2=beg2,e2=end2,snp,indel)
    to %>%# filter(c1!='B99') %>%
        mutate(x = size1-del-aln, y =size2-ins-aln) %>%
        select(size1,del,x, size2, ins,y) %>%
        filter(x!=0 | y!=0) %>%
        print(n=10)
    to
    #}}}
}

fm=glue("{dirw}/11.BtoM.map.bed")
fv=glue("{dirw}/11.BtoM.vnt.bed")
to = sum_crossmap(fm, fv)
#}}}




##### obsolete #####
#{{{ prepare promoter db for B/M/W (individually)
gt = 'Zmays_B73'
gt = 'Zmays_Mo17'
gt = 'Zmays_W22'
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
