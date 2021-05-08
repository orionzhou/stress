source('functions.R')
dirw = glue('{dird}/31_promoter')
setwd(dirw)
xref = read_xref()
gts = c("Mo17",'B73','W22')
txdbs = tibble(gt = gts) %>% mutate(txdb=map(gt, load_txdb, primary=T))

#{{{ obtian TSS/TTS/promoter/terminator coords
tssB = get_tss_tts('Zmays_B73') %>% mutate(gid2=gid,gt="B73")
#{{{ M
tssM = get_tss_tts('Zmays_Mo17')
tssM = xref %>% filter(qry=='Mo17',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssM, by=c('gid2'='gid')) %>%
    select(gid=gid1,gid2,chrom,tss,tts,srd) %>%
    mutate(gt="Mo17")
#}}}
#{{{ W
tssW = get_tss_tts('Zmays_W22')
tssW = xref %>% filter(qry=='W22',tgt=='B73') %>%
    select(gid1,gid2,type) %>%
    inner_join(tssW, by=c('gid2'='gid')) %>%
    select(gid=gid1,gid2,chrom,tss,tts,srd) %>%
    mutate(gt="W22")
#}}}
tss = rbind(tssB, tssM, tssW)

chrom_size = read_tsv(glue('{dirw}/../21_seq/00.sizes'), col_names=c('chrom','size'))
to = tss %>%
    mutate(gb = ifelse(srd == '-', tts, tss)) %>%
    mutate(ge = ifelse(srd == '-', tss, tts)) %>%
    mutate(pb = ifelse(srd == '-', tss, tss-2000)) %>%
    mutate(pe = ifelse(srd == '-', tss+2000, tss)) %>%
    mutate(pb1 = ifelse(srd == '-', tss, tss-500)) %>%
    mutate(pe1 = ifelse(srd == '-', tss+500, tss)) %>%
    mutate(tb = ifelse(srd == '-', tts-2000, tts)) %>%
    mutate(te = ifelse(srd == '-', tss, tts+2000)) %>%
    inner_join(chrom_size, by='chrom') %>%
    mutate(pb=pmax(0, pb), tb=pmax(0,tb),
           pe=pmin(pe, size), te=pmin(te, size)) %>%
    select(gt, gid, gid2, chrom, gb, ge, srd, tss, tts, pb, pe, pb1, pe1, tb, te) %>%
    arrange(gt,chrom,gb,ge)

fo = glue('{dirw}/00.coord.tsv')
write_tsv(to, fo)
#}}}

#{{{ find syntenic promoter/gene region to M/W
fi = glue('{dirw}/00.coord.tsv')
tss = read_tsv(fi)

# lift 500 / 2k promoter region
ti1 = tss %>% filter(gt == 'B73') %>% select(chrom,pb=pb1,pe=pe1,gid) %>%
    arrange(chrom,pb,pe)
fo = glue('{dirw}/11.promoter.B.bed')
write_tsv(ti1, fo, col_names=F)

gt = 'W22'
gt = 'Mo17'
fi = glue('{dirw}/11.promoter.B.bed')
fo = glue('{dirw}/13.{gt}.rds')
system(glue("liftover.R {fi} {fo} --aln $wgc/data/raw/Zmays_{gt}-Zmays_B73/aln.bed --vnt $wgc/data/raw/Zmays_{gt}-Zmays_B73/vnt.bed"))

# liftover entire genic region +/- 3kb
ti1 = tss %>% filter(gt == 'B73') %>%
    mutate(gb = gb - 3000, ge = ge + 3000) %>% select(chrom,gb,ge,gid) %>%
    arrange(chrom,gb, ge)
fo = glue('{dirw}/15.gene.B.bed')
write_tsv(ti1, fo, col_names=F)

gt = 'Mo17'
gt = 'W22'
fi = glue('{dirw}/15.gene.B.bed')
fo = glue('{dirw}/16.{gt}.rds')
system(glue("liftover.R {fi} {fo} --aln $wgc/data/raw/Zmays_{gt}-Zmays_B73/aln.bed --vnt $wgc/data/raw/Zmays_{gt}-Zmays_B73/vnt.bed"))
#}}}

# variable response
#{{{ read
#{{{ DEG
fg = glue('{dird}/15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2
#
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
td = td4 %>% inner_join(td4s, by=c('cond','drc','gid')) %>%
    spread(gt,st) %>% mutate(st=glue('{B73}{Mo17}{W22}')) %>%
    select(cond,drc,gid,st)
#}}}
#
tpro = tibble(qry=gts[c(1,3)], tgt=gts[2]) %>%
    mutate(fa = glue("{dirw}/13.{qry}.rds")) %>%
    mutate(ta = map(fa, readRDS)) %>% select(-fa)
tpro1 = tpro %>%
    unnest(ta) %>% select(qry,gid=id,size) %>%
    spread(qry,size) %>% rename(sizeM=Mo17,sizeW=W22)
#
fi = glue("{dird}/41_ml/06.tk.tc.rds")
r6 = readRDS(fi)
tc = r6$tc; tk = r6$tk
tc1 = tc %>% filter(train=='B') %>% filter(str_detect(note,'all')) %>%
    mutate(drc = ifelse(str_detect(note, 'up'), 'up', 'down')) %>%
    select(cid,cond,drc)

fi = glue('{dird}/16_ase/20.rds')
res = readRDS(fi)
ddeg = res$ddeg
ddeg2 = ddeg %>% filter(condB != 'Control0') %>%
    mutate(time=as.integer(str_sub(cond,5,6))) %>%
    mutate(cond=str_to_lower(str_sub(cond,1,4))) %>%
    mutate(st = glue("d{deg}")) %>%
    mutate(st = factor(st, levels=levels(td2$st))) %>%
    select(cond,time,qry,tgt,gid,st,reg)
#}}}

#{{{ box/scatter plot
tags1 = c("dA+B=", "dA=B+", "dA-B=", "dA=B-")
tags2 = c("nA+B+",'nA-B-','nA=B=')
    #mutate(st1 = ifelse(st %in% c("nA+B+","nA-B-"), 'n', st1)) %>%
conds = c("cold_25h","heat_1h", "heat_25h")
x = td2 %>%
    filter(qry==gt) %>%
    left_join(ddeg2, by=c('cond','time','qry','tgt','gid','st')) %>%
    filter(!is.na(reg) & reg == 'cis') %>%
    mutate(cond = glue("{cond}_{time}h")) %>%
    filter(cond %in% conds) %>% mutate(cond = factor(cond, levels=conds)) %>%
    filter(st %in% c(tags1,tags2))
#
tp = x %>% inner_join(ti, by=c('gid'='id')) %>%
    mutate(dst = ifelse(size >= 495, 1, 0))
#tp %>% count(st, dst) %>% spread(dst,n) %>% mutate(p=`1`/(`1`+`0`))
tps = tp %>% count(cond, st) %>% mutate(lab=number(n,accuracy=1))
#
p = ggplot(tp, aes(x=st, y=size, color=st)) +
    geom_boxplot(size=.5, outlier.shape = NA) +
    geom_jitter(width=.3, size=.3, alpha=.2) +
    geom_text(data=tps, aes(x=st,y=500*1.01,label=lab), size=2.5,vjust=0) +
    scale_y_continuous(expand=expansion(add=c(5,20))) +
    facet_wrap(~cond, ncol=3) +
    scale_color_npg() +
    otheme(ytitle=T,ytext=T,ytick=T, xtick=T, xtext=T,
           legend.pos='none', panel.spacing=.2, strip.compact=F) +
    theme(axis.text.x = element_text(size=7.5, angle=0, hjust=.5, vjust=.5))
fo = glue("{dirw}/21.bp.syn.{gt}.2.pdf")
ggsave(p, file=fo, width=10, height=6)
#}}}

##### synteny plots #####
fi = glue('{dirw}/00.coord.tsv')
tss = read_tsv(fi)
#
cmps = tibble(qry=gts[c(1,3)], tgt=gts[2]) %>%
    mutate(fa = glue("{dirw}/16.{qry}.rds")) %>%
    mutate(ta = map(fa, readRDS)) %>% select(-fa)
x = td %>%# filter(st == 'UNN') %>%
    #left_join(ddeg2, by=c('cond','time','qry','tgt','gid','st')) %>%
    #filter(!is.na(reg) & reg == 'cis') %>%
    inner_join(tpro1, by='gid') %>%
    inner_join(tc1, by=c('cond','drc')) %>%
    arrange(sizeM) %>% print(n=50) #%>% filter(sizeM < 300, sizeM > 0)

#{{{ functions
plot_kmer <- function(gt, gid, cid, srd, pb, pe) {
    #{{{
    if(is.na(gid)) {
        return(ggplot()+otheme(panel.border=F,xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0)))
    }
    #{{{ prepare
    ggid = glue("{gt}_{gid}")
    ti = tibble(gid=ggid,st=1)
    write_tsv(ti, 'tmp.tsv', col_names=F)
    fm = glue("{dird}/41_ml/00_nf/03_motif_lists/{cid}.tsv")
    cmd = glue("kmer.py prepare_ml tmp.tsv {fm} ",
               "--bin 'TSS:-/+2k,TTS:-/+2k' ",
               "--epi umr --nfea top100 --mod anr ",
               "--fmt long out.tsv")
    system(cmd)
    tp = read_tsv('out.tsv')
    if (nrow(tp)==0) 
        return(ggplot()+otheme(panel.border=F,xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0)))
    tp = tp %>% mutate(pos = round(start+end)) %>%
        select(fid,pos,srd)
    if(srd == '-') tp = tp %>% mutate(pos = pe - pos)
    system("rm tmp.tsv out.tsv")
    #}}}
    ggplot(tp) +
        #geom_dotplot(aes(x=pos), binwidth=(pe-pb)/100, fill="white", stroke=.1) +
        geom_jitter(aes(x=pos,y=0), height=1, width=50, size=.2) +
        scale_x_continuous(limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(expand=expansion(mult=c(.1,.1))) +
        otheme(panel.border=F, xtick=F,xtext=F,margin=c(0,0,0,0))
    #}}}
}
txdb_gene <- function(txdb, gid, gb=0) {
    #{{{
    if(is.na(gid)) return(NA)
    ti = AnnotationDbi::select(txdb, keys=gid, columns=columns(txdb), keytype='GENEID') %>%
        rename(gid=GENEID, tid=TXNAME, chrom=TXCHROM, tb=TXSTART, te=TXEND,
            srd=TXSTRAND, eb=EXONSTART, ee=EXONEND, cb=CDSSTART, ce=CDSEND)
    ti1 = ti %>% distinct(gid,tid,beg=tb,end=te,srd) %>% mutate(type='rna')
    ti2 = ti %>% filter(!is.na(ee)) %>%
        select(gid,tid,beg=eb,end=ee,srd) %>% mutate(type='exon')
    ti3 = ti %>% filter(!is.na(cb)) %>%
        select(gid,tid,beg=cb,end=ce,srd) %>% mutate(type='cds')
    rbind(ti1,ti2,ti3) %>%
        mutate(beg=beg-gb, end=end-gb)
    #}}}
}
plot_genes <- function(tg, pb, pe, n_arrows=10, ht.exon=.2, ht.cds=.4, x=T) {
    #{{{
    if(is.na(tg)) {
        p = ggplot() +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
        return(p)
    }
    #{{{ prepare
    if (! 'y' %in% colnames(tg)) {
        tg0 = tg %>% filter(type=='rna') %>% arrange(gid) %>% mutate(y = 1:n())
        tg = tg %>% inner_join(tg0 %>% dplyr::select(tid,y), by='tid')
    }
    ir = function(xs) tibble(b=xs[-length(xs)], e = xs[-1])
    tgr = tg %>% filter(type=='rna') %>% mutate(pos=(beg+end)/2)
    tgx = tgr %>% mutate(size=end-beg) %>% arrange(desc(size))
    brks = c(tgx$beg[[1]], tgx$end[[1]])
    labs = c("TSS", "TTS")
    if (tgx$srd[[1]] == '-')  labs = c("TTS", "TSS")
    tga = tgr %>% mutate(x=map2(beg,end,seq,length.out=n_arrows)) %>%
        mutate(x1 = map(x, ir)) %>%
        dplyr::select(gid,tid,srd,y,x1) %>% unnest(x1)
    tga1 = tga %>% filter(srd=='+')
    tga2 = tga %>% filter(srd=='-')
    tge = tg %>% filter(type=='exon')
    tgc = tg %>% filter(type=='cds')
    #}}}
    #{{{
    col.exon = 'royalblue'; col.cds = 'royalblue'
    col.syn = 'grey'
    wd.arrow1 = .1; wd.arrow2 = .08
    arrow11 = arrow(length=unit(wd.arrow1,'cm'), angle=30, ends='last',type="open")
    arrow12 = arrow(length=unit(wd.arrow1,'cm'), angle=30, ends='first',type="open")
    arrow21 = arrow(length=unit(wd.arrow2,'cm'), angle=30, ends='last',type="open")
    arrow22 = arrow(length=unit(wd.arrow2,'cm'), angle=30, ends='first',type="open")
    #}}}
    ggplot() +
        geom_segment(data=tgr,aes(x=beg,xend=end,y=y,yend=y),size=.5) +
        geom_segment(data=tga1,aes(x=e-5,xend=e,y=y,yend=y),
                     color='black', size=1, arrow=arrow11) +
        geom_segment(data=tga2,aes(x=b,xend=b+5,y=y,yend=y),
                     color='black', size=1, arrow=arrow12) +
        geom_rect(data=tge,aes(xmin=beg,xmax=end,ymin=y-ht.exon,ymax=y+ht.exon),fill=col.exon,color=NA,alpha=1) +
        geom_rect(data=tgc,aes(xmin=beg,xmax=end,ymin=y-ht.cds,ymax=y+ht.cds),fill=col.cds,color=NA,alpha=1) +
        geom_segment(data=tga1,aes(x=e-2,xend=e-1,y=y,yend=y),
                     color='white', size=.3, arrow=arrow21) +
        geom_segment(data=tga2,aes(x=b+1,xend=b+2,y=y,yend=y),
                     color='white', size=.3, arrow=arrow22) +
        geom_text(data=tgr,aes(x=pos,y=y+ht.cds*1.1,label=tid), vjust=0, size=2) +
        scale_x_continuous(breaks=brks, labels=labs,
                           limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        #scale_y_continuous(expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(expand=expansion(add=c(0,.4),mult=c(0,0))) +
        otheme(panel.border=F, xtick=T,xtext=T,margin=c(0,0,0,0)) +
        theme(axis.text.x = element_text(size=6))
    #}}}
}
prepare_syn <- function(ta,qgid,tgid,qtss,ttss) {
    #{{{
    if (is.na(qgid))
        return(tibble(srd='*', aln=list(tibble())))
    ta1 = ta %>% filter(id==tgid)
    o1 = ttss - ta1$b1; o2 = qtss - ta1$b2
    srd = ta1$srd[[1]]
    aln = ta1$aln[[1]] %>%
        mutate(rb1=rb1-o1, re1=re1-o1, rb2=rb2-o2, re2=re2-o2)
    tibble(srd=srd, aln=list(aln))
    #}}}
}
plot_syn <- function(ta, srd, pb, pe, qtop=F,size=.2) {
    #{{{
    if(is.na(ta) || nrow(ta)==0) {
        p = ggplot() +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
        return(p)
    }
    #{{{
    ta = ta %>% filter(rb1<pe,re1>pb, rb2<pe, re2>pb) %>%
        mutate(ob1 = pmax(0, pb-rb1)) %>%
        mutate(oe1 = pmax(0, re1-pe)) %>%
        mutate(rb1 = rb1 + ob1, re1 = re1 - oe1) %>%
        #mutate(rb2 = rb2 + ob1, re2 = re2 - oe1) %>%
        mutate(ob2 = pmax(0, pb-rb2)) %>%
        mutate(oe2 = pmax(0, re2-pe)) %>%
        #mutate(rb1 = rb1 + ob2, re1 = re1 - oe2) %>%
        mutate(rb2 = rb2 + ob2, re2 = re2 - oe2)
    if(qtop) {
        tay = tibble(xt=c('rb1','re1','rb2','re2'), y=c(-1,-1,1,1))
    } else {
        tay = tibble(xt=c('rb1','re1','rb2','re2'), y=c(1,1,-1,-1))
    }
    if(srd == "-") {
        tay = tay %>% mutate(xo=c(1,2,1,2))
    } else {
        tay = tay %>% mutate(xo=c(1,2,2,1))
    }
    tp = ta %>% mutate(i=1:n()) %>% gather(xt, x, -i) %>%
        inner_join(tay, by='xt') %>%
        arrange(i,y,xo)
    #}}}
    col.syn='coral'
    ggplot() +
        geom_polygon(data=tp, aes(x=x,y=y,group=i), fill=col.syn, alpha=.7,
                     color=ifelse(size>0,'black',NA),size=size) +
        scale_x_continuous(limits=c(pb,pe),expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(limits=c(-1,1), expand=expansion(mult=c(.0,.0))) +
        otheme(panel.border=F, xtick=F,xtext=F,panel.spacing=0,margin=c(0,0,0,0))
    #}}}
}
plot_combo <- function(gid,cid,cond,drc,st,fo, tss, cmps, gts, off=2000) {
    #{{{
    #drc = ifelse(str_detect(st, "\\+"), "+", "-")
    tit = glue("{gid} {cond} {drc} (B|M|W) {st}")
    cat(tit,'\n')
    #tss %>% filter(gid==!!gid) %>% print(width=Inf)
    cfg = tibble(gt=gts) %>% mutate(gid=!!gid) %>%
        mutate(gt=factor(gt,levels=gts)) %>% arrange(gt) %>%
        left_join(tss, by=c('gt','gid')) %>%
        select(gt,gid,gid2,chrom,gb,ge,tss,srd) %>%
        mutate(pb=gb-tss-off, pe=ge-tss+off) %>%
        mutate(pb = min(pb,na.rm=T)) %>%
        mutate(pe = max(pe,na.rm=T)) %>%
        mutate(pk = pmap(list(gt,gid,cid,srd,pb,pe), plot_kmer)) %>%
        inner_join(txdbs, by='gt') %>%
        mutate(tg = pmap(list(txdb, gid2, tss), txdb_gene))
    #
    o = cfg %>%
        mutate(pg = pmap(list(tg,pb,pe), plot_genes, n_arrows=10, ht.exon=.2, ht.cds=.4))
    #
    cmp1 = cmps %>%
        inner_join(cfg %>% select(qry=gt,qgid=gid2,qtss=tss),by='qry') %>%
        inner_join(cfg %>% select(tgt=gt,tgid=gid2,ttss=tss,pb,pe),by='tgt')
    oc = cmp1 %>%
        mutate(x = pmap(list(ta,qgid,tgid,qtss,ttss), prepare_syn)) %>%
        unnest(x) %>% mutate(qtop=c(T,F)) %>%
        mutate(pc = pmap(list(aln,srd,pb,pe,qtop), plot_syn, size=0))
    #
    p = o$pg[[1]] + o$pk[[1]] +
        oc$pc[[1]] +
    o$pg[[2]] + o$pk[[2]] +
        oc$pc[[2]] +
    o$pg[[3]] + o$pk[[3]] +
        plot_layout(ncol=1, heights=c(1,1,1.5, 1,1,1.5, 1,1)) +
        plot_annotation(title=tit,theme=theme(plot.title=element_text(size=8)))
    ggsave(p, filename=fo, width=5, height=3)
    p
    #}}}
}
#}}}
x %>% print(n=50)
i = which(x$gid == gid)
gid=x$gid[[i]];cid=x$cid[[i]];cond=x$cond[[i]];drc=x$drc[[i]];st=x$st[[i]]
j = x %>% slice(i) %>%
    mutate(fo = glue("{dirw}/50_syn_plots/{str_sub(gid,10)}.pdf")) %>%
    mutate(p = pmap(list(gid,cid,cond,drc,st,fo), plot_combo,
                    tss=tss, cmps=cmps, gts=gts, off=2000))# %>%
    #mutate(x = map2(p, fo, ~ggsave(.x,filename=.y, width=5, height=3)))

to = x %>% filter(st=='dA=B+') %>% filter(size < 300, size>0) %>%
    select(cond,gid,qry,st,size,ins,del,n_snp) %>% print(width=Inf,n=30)
fo = glue("{dirw}/tmp.tsv")
write_tsv(to,fo)

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
#{{{ read kmer locations [long time]
read_mtf_loc1 <- function(fi) {
    #{{{
    read_tsv(fi) %>%
        mutate(pos = round((start+end)/2)) %>%
        #separate(gid, c('gt','gid'), sep='_') %>%
        distinct(fid,gid,pos)
    #}}}
}
# save for promoter visulization
km = tc %>% filter(train=='BMW', str_detect(note, '^all')) %>%
    select(cid,cond,note) %>%
    mutate(fi = glue("{dird}/41_ml/52_mtf_loc_umr/{cid}.tsv")) %>%
    mutate(x = map(fi, read_mtf_loc1)) %>%
    select(-fi)
#}}}

