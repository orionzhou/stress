source('functions.R')
dirw = glue('{dird}/45_nam')
tgl = gcfg$gene %>% select(gid,chrom,start,end)
#{{{ functions
read_vnt <- function(chrom,beg,end,
                     fv='~/projects/wgc/data/05_maize_merged/08.vcf.gz',
                     gts=gts31) {
    #{{{
    cmd = glue("bcftools view -r {chrom}:{beg}-{end} {fv} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\n' > tmp.tsv")
    system(cmd)
    if (file.info('tmp.tsv')$size == 0) return(tibble())
    gt_map = c('0/0'=0, '1/1'=1, './.'=NA)
    tv = read_tsv('tmp.tsv', col_names=c('chrom','pos','ref','alt',gts),
        col_types=str_c(c('c','i','c','c', rep('c',length(gts))), sep='', collapse='')) %>%
        mutate(B73 = '0/0') %>%
        gather(sid, gt, -chrom, -pos, -ref, -alt) %>%
        mutate(gt = gt_map[gt])
    system("rm tmp.tsv")
    tv
    #}}}
}
#}}}

#{{{ read variable DEG
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
    select(cond,drc,gid,st) %>% filter(cond=='cold')
#}}}

#{{{ prepare NM matrix
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm

th1 = th %>% filter(Experiment=='NM') %>%
    select(sid=SampleID, gt=Genotype, cond=Treatment, time=Timepoint) %>%
    mutate(cond = str_to_lower(cond))

rr = tm %>% select(gid, sid=SampleID, cpm=CPM) %>%
    inner_join(th1, by ='sid') %>%
    select(-sid) %>% spread(time, cpm)
rd = rr %>% gather(time, rc, -gid, -gt, -cond) %>%
    spread(cond, rc) %>%
    mutate(log2fc = log2(cold/control)) %>%
    arrange(gid, gt)

rd2 = td %>% inner_join(rd, by='gid') %>%
    arrange(drc,time,gid,gt) %>%
    group_by(drc,time,st,gid) %>% summarise(log2fcs = list(log2fc)) %>%
    ungroup()

require(mclust)
xs=rd2$log2fcs[[i]]
bimod <- function(i, xs, drc, gts25, gts3) {
    #{{{
    #cat(i,'\n')
    pb$tick()
    names(xs) = gts25
    xs1 = xs[!is.na(xs)]
    if (length(xs1) < 10) return(tibble())
    if(length(which(xs1==Inf)) >0) xs1[which(xs1==Inf)] = max(xs1[xs1!=Inf],na.rm=T)
    if(length(which(xs1==-Inf)) >0) xs1[which(xs1==-Inf)] = min(xs1[xs1!=-Inf],na.rm=T)
    mod = densityMclust(xs1, verbose=F)
    #summary(mod)
    n=length(mod$parameters$pro)
    pros=mod$parameters$pro
    means=mod$parameters$mean
    classes = mod$classification
    if (is.null(names(means))) names(means) = classes[1]
    if (is.null(names(classes))) names(classes) = names(xs1)
    clu = tibble(clu=as.integer(names(means)), avg=as.numeric(means), prop=pros)
    if (drc == 'up') {
        clu = clu %>% mutate(nclu = ifelse(avg < 1, 0, 1))
    } else if (drc == 'down') {
        clu = clu %>% mutate(nclu = ifelse(avg <= -1, -1, 0))
    }
    t_lfc = tibble(gt=names(xs1), log2fc=as.numeric(xs1))
    t_clu = tibble(gt=names(classes), clu=as.integer(classes)) %>%
        inner_join(clu %>% select(clu,nclu), by='clu') %>%
        select(-clu) %>% rename(clu=nclu) %>%
        inner_join(t_lfc, by='gt')
    clu = t_clu %>% group_by(clu) %>%
        summarise(avg = mean(log2fc), prop=n()/nrow(t_clu)) %>% ungroup()
    tibble(n=nrow(clu), bic=BIC(mod), clu=list(clu), t_clu=list(t_clu))
    #}}}
}
pb <- progress_bar$new(total=nrow(rd2))
rd3 = rd2 %>% mutate(i=1:n()) %>%
    mutate(x=pmap(list(i,log2fcs,drc), bimod,gts25=gts25,gts3=gts3)) %>%
    select(-i,-log2fcs) %>% unnest(x)

rd3 %>% count(drc,time, n, bic>50) %>% print(n=40)
rd4 = rd3 %>% filter(n>1, bic>50)

fo = glue('{dirw}/01.mclust.rds')
saveRDS(rd3, fo)
#}}}

fi = glue('{dirw}/01.mclust.rds')
ti = readRDS(fi)
ti %>% count(drc,time, n, bic>50) %>% print(n=40)
rd = ti %>% filter(n>1, bic>50)

gid = 'Zm00001d046561'
gid = 'Zm00001d044762'
rd %>% filter(gid==!!gid)
rd %>% filter(gid==!!gid) %>% pluck('clu',1)
rd %>% filter(gid==!!gid) %>% pluck('t_clu',1) %>% print(n=25)

#{{{ prepare variable response
fg = glue('{dird}/15_de/09.gene.status.rds')
x = readRDS(fg)
td1=x$td1; td2=x$td2
#
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
#{{{
i=2; vnt=x$vnt[[i]]; t_clu=x$t_clu[[i]]
fet <- function(x00,x01,x10,x11) fisher.test(matrix(c(x00,x10,x01,x11), nrow=2))$p.value
vnt_scan <- function(t_clu, chrom,start,end) {
    #{{{
    pb$tick()
    vnt = read_vnt(chrom,start,end)
    if(nrow(vnt)==0) return(tibble())
    y = vnt %>% filter(!is.na(gt)) %>% inner_join(t_clu, by=c('sid'='gt'))
    y1 = y %>% group_by(chrom,pos,ref,alt) %>% nest() %>% rename(v=data) %>% ungroup()
    ys = y %>% group_by(chrom,pos,ref,alt) %>%
        summarise(x00=sum(gt==0 & clu==0), x01=sum(gt==0 & clu==1),
        x10 = sum(gt==1 & clu==0), x11 = sum(gt==1 & clu==1)) %>%
        ungroup() %>%
        mutate(p.raw=pmap_dbl(list(x00,x01,x10,x11), fet)) %>%
        mutate(p.adj=p.adjust(p.raw)) %>%
        inner_join(y1, by=c('chrom','pos','ref','alt'))
    ys
    #}}}
}

pb <- progress_bar$new(total=nrow(rd))
x0 = rd %>% #filter(gid == !!gid) %>%
    inner_join(tgl, by='gid') %>%
    mutate(start = start - 2000, end = end + 2000) %>%
    mutate(x = pmap(list(t_clu,chrom,start,end), vnt_scan)) %>%
    select(drc,time,st,gid,x) %>% unnest(x)

fo = glue("{dirw}/05.test.rds")
saveRDS(x0, fo)
x0 = readRDS(fo)

x1 = x0 %>% mutate(refl=nchar(ref), altl=nchar(alt),
                   vtype = ifelse(refl==1&altl==1, 'snp', 'idl'),
                   size=glue("{refl}->{altl}"), x=str_c(x00,x01,x10,x11,sep=','))
x1s = x1 %>% group_by(drc,time,st,gid) %>%
    summarise(nsnp=sum(vtype=='snp'),nidl=sum(vtype=='indel'),
        nsig=sum(p.raw<.05)) %>%
    ungroup()
x1s %>% filter(nsig> 0)
x2 = x1 %>% filter(p.raw < .05) %>% arrange(p.raw)
x3 = x2 %>% select(time,st,gid,chrom,pos, size,x,p.raw,v)
x3 %>% select(-v) %>% print(width=Inf,n=20)

i = 292
gid = x3 %>% pluck('gid', i)
x3 %>% select(-v) %>% filter(gid==!!gid) %>% print(width=Inf)

#{{{
tp = x2 %>% slice(i) %>%
    mutate(tit=glue("{time}h {gid} [B|M|W]={st} {chrom}:{pos}:{ref}->{alt}")) %>%
    unnest(v) %>% select(gid,tit,p.raw,sid,gt,clu,log2fc) %>%
    arrange(clu) %>%
    mutate(gtc = ifelse(sid %in% gts3, sid, 'ZZZ'))
#
p = ggplot(tp, aes(x=gt, y=log2fc)) +
    geom_boxplot(aes(group=gt), size=.5, width=.7) +
    geom_jitter(aes(color=gtc), width=.3, size=2, alpha=.9) +
    #geom_text(data=tps, aes(x=st,y=500*1.01,label=lab), size=2.5,vjust=0) +
    scale_x_continuous(breaks=c(0,1),labels=c('ref','alt'),expand=expansion(mult=c(.2,.2))) +
    scale_y_continuous(name='log2FoldChange', expand=expansion(mult=c(.02,.02))) +
    scale_color_aaas() +
    ggtitle(tp$tit[1]) +
    otheme(ytitle=T,ytext=T,ytick=T, xtick=T, xtext=T,
           legend.pos='top.right', panel.spacing=.2, strip.compact=F) +
    theme(plot.title=element_text(size=8))
#}}}
fo = glue("{dirw}/{tp$gid[1]}.pdf")
ggsave(p, file=fo, width=4, height=4)

#}}}

ddeg3 = ddeg2 %>% filter(cond=='cold',tgt=='B73') %>% select(time,gid,reg)
j = x1s %>%
    mutate(time=as.integer(time)) %>%
    #filter(time==25) %>%
    mutate(sig = ifelse(nsig > 0,'sig','insig')) %>%
    left_join(ddeg3, by=c("time",'gid'))
j %>% count(reg, sig) %>% spread(sig,n) %>% mutate(prop=sig/(sig+insig))
