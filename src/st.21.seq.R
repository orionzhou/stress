source('functions.R')
dirw = glue('{dird}/21_seq')
setwd(dirw)
xref = read_xref()

#{{{ prepare promoter db for B/M/W/
fi = glue("{dird}/21_seq/regions.xlsx")
tr = read_xlsx(fi) %>% filter(offu==2000,offd==2000) %>% mutate(bin=as_factor(bin))
bins = levels(tr$bin)
gts = gts32_ph207
get_ortho_tss <- function(gt, xref) {
#{{{
    if (gt == "B73") {
        tss = get_tss_tts('Zmays_B73') %>% mutate(gid=glue("B73_{gid}"))
    } else {
    tss = get_tss_tts(glue("Zmays_{gt}"))
    tss = xref %>% filter(qry==gt,tgt=='B73') %>%
        select(gid1,gid2,type) %>%
        inner_join(tss, by=c('gid2'='gid')) %>%
        select(gid=gid1,chrom,tss,tts,srd) %>%
        mutate(gid=glue("{gt}_{gid}"))
    }
    tss
#}}}
}

t_ch = tibble(gt=gts) %>%
    mutate(fi=glue("{dirg}/Zmays_{gt}/15_intervals/01.chrom.sizes")) %>%
    mutate(ti = map(fi, read_tsv, col_names=c("chrom",'size'))) %>%
    select(gt,ti) %>% unnest(ti)

tss = tibble(gt=gts) %>% mutate(tss = map(gt, get_ortho_tss, xref=xref)) %>%
    select(gt, tss) %>% unnest(tss)

#{{{ prepare tl
fz = glue("{dirw}/01.sizes")
to = t_ch %>% mutate(chrom=glue("{gt}_{chrom}")) %>% select(chrom, size)
write_tsv(to, fz, col_names=F)

tt = tss %>% crossing(tr) %>%
    mutate(pos = ifelse(str_detect(bin,"^TTS"), tts, tss)) %>%
    mutate(start = pos, end = pos + 1) %>%
    mutate(start = ifelse(srd=='-', pos-offd, pos-offu)) %>%
    mutate(end = ifelse(srd=='-', pos+offu, pos+offd)) %>%
    inner_join(t_ch, by=c('gt','chrom')) %>%
    mutate(start=pmax(0, start), end=pmin(end, size)) %>%
    select(gt, chrom, start, end, srd, gid, bin) %>%
    arrange(gt, chrom,start,end) %>%
    mutate(sid = sprintf("s%06d", 1:n()))
tt %>% mutate(size=end-start) %>% count(bin,size)
tl1 = tt %>% filter(srd == '+') %>% arrange(gid, chrom, start) %>%
    group_by(gid) %>% mutate(i = 1:n()) %>% ungroup()
tl2 = tt %>% filter(srd == '-') %>% arrange(gid, chrom, desc(start)) %>%
    group_by(gid) %>% mutate(i = 1:n()) %>% ungroup()
tl = tl1 %>% bind_rows(tl2) %>% arrange(gid, i)
#}}}

#{{{ obtain segment sequence ts
to = tt %>% mutate(score='.') %>% select(gt,chrom,start,end,sid,score,srd) %>%
    group_by(gt) %>% nest() %>% ungroup() %>%
    mutate(fo = glue("{dirw}/01_tss/{gt}.bed")) %>%
    mutate(l = map2(data, fo, write_tsv, col_names=F))

# cat $genome/data2/gts32_ph207.txt | parallel echo bedtools getfasta -fi $genome/data/Zmays_{}/10.fasta -tab -bed {}.bed -fo {}.tsv -nameOnly -s

ts = to %>% select(gt) %>% mutate(fs = glue("{dirw}/01_tss/{gt}.tsv")) %>%
    mutate(x = map(fs, read_tsv, col_names=c('sid','seq'))) %>%
    select(gt, x) %>% unnest(x) %>%
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
system("fasta.py size 02.fas 02.sizes")
system("fasta-get-markov -dna 02.fas 02.bg")

tb0 = tl %>% distinct(gid) %>% mutate(cid = 1:n())
tb = tl %>% inner_join(tb0, by='gid') %>%
    mutate(chrom = glue("{gt}_{chrom}")) %>%
    mutate(start2=0, end2=end-start) %>%
    mutate(start2 = ifelse(i==2, start2+4050, start2)) %>%
    mutate(end2=ifelse(i==2, end2+4050, end2)) %>%
    arrange(chrom,start,end) %>%
    select(chrom,start,end, srd, gid,start2,end2, cid)
write_tsv(tb, '10.bed', col_names = F)
#
system(glue("chain.py fromBed 10.bed 01.sizes 02.sizes 10.reverse.chain"))
system("chainSwap 10.reverse.chain 10.forward.chain")
#}}}

# create UMR/ACL BED


