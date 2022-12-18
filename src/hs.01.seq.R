source('functions.R')
dirw = glue('{dird}/73_zach_hsf')
setwd(dirw)

gts = c('Zmays_B73v4', 'Sviridis_A10v2')
idsM = read_tsv("00.gid.maize.txt", col_names=F)$X1
idsS = read_tsv("00.gid.setaria.txt", col_names=F)$X1
th = tibble(gt = gts, gid = list(idsM, idsS)) %>% unnest(gid)

#{{{ prepare promoter db
fi = glue("{dird}/21_seq/regions.xlsx")
tr = read_xlsx(fi) %>% filter(offu==2000,offd==2000) %>% mutate(bin=as_factor(bin))
bins = levels(tr$bin)

t_ch = tibble(gt=gts) %>%
    mutate(fi=glue("{dirg}/{gt}/15_intervals/01.chrom.sizes")) %>%
    mutate(ti = map(fi, read_tsv, col_names=c("chrom",'size'))) %>%
    select(gt,ti) %>% unnest(ti)

tss = tibble(gt=gts) %>% mutate(tss = map(gt, get_tss_tts)) %>%
    select(gt, tss) %>% unnest(tss)

tss.hsf = tss %>% inner_join(th, by=c('gt','gid'))
write_tsv(tss.hsf, '62.hsf.tss.tsv', col_names = T)

#{{{ prepare tl
fz = glue("{dirw}/chrom.sizes")
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
    mutate(sid = sprintf("s%08d", 1:n()))
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

# cat ../genomes.txt | parallel echo bedtools getfasta -fi $genome/data/{}/21_dbs/gatk/db.fasta -tab -bed {}.bed -fo {}.tsv -nameOnly -s

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
#system("fasta.py extract 02.fa s1")
system("fasta.py size 02.fa 02.sizes")
#system("fasta-get-markov -dna 02.fa 02.bg")

tb0 = tl %>% distinct(gid) %>% mutate(cid = 1:n())
tb1 = tl %>% inner_join(tb0, by='gid') %>%
    mutate(chrom = glue("{gt}_{chrom}")) %>%
    mutate(start2=0, end2=end-start, len=end2-start2)
tb2a = tb1 %>% filter(i==1)
tb2as = tb2a %>% select(gid, off=len)
tb2b = tb1 %>% filter(i==2) %>% inner_join(tb2as, by='gid') %>%
    mutate(start2 = start2+off+50) %>%
    mutate(end2 = end2+off+50) %>% select(-off)
tb2 = tb2a %>% bind_rows(tb2b) %>%
    arrange(chrom,start,end) %>%
    select(chrom,start,end, srd, gid,start2,end2, cid)
write_tsv(tb2, '10.bed', col_names = F)
#
system(glue("chain.py fromBed 10.bed chrom.sizes 02.sizes 10.reverse.chain"))
system("chainSwap 10.reverse.chain 10.forward.chain")
#}}}


fi = glue("{dird}/25_dreme/05.best.mtfs.rds")
r = readRDS(fi)
tk = r$tk %>% filter(bid=='b07') %>% select(-bid, -n_mtf) %>% unnest(mtfs)
# bash readme.sh
tk1 = tk %>% select(i, motif_id=mid, motif=fname, known, conseq, pval)
fo = glue("{dirw}/67.motifs.tsv")
write_tsv(tk1, fo)

fi = glue("{dirw}/63.motif.tsv")
tx = read_tsv(fi, col_names=c('mid','beg','end')) %>%
    separate('mid', c('mid', 'gid'), sep='%') %>%
    mutate(pos = round((beg+end)/2)) %>%
    select(mid, gid, pos)
to = tx %>% filter(pos <= 4000) %>% mutate(pos = pos - 2000) %>%
    inner_join(tk %>% select(mid, fname), by='mid') %>%
    select(gid, pos, motif_id=mid, motif=fname) %>%
    arrange(gid, pos)

fo = glue("{dirw}/68.motif.loc.tsv")
write_tsv(to, fo)
