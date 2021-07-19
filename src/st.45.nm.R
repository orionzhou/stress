source('functions.R')
dirw = glue('{dird}/45_nam')
setwd(dirw)
loci = read_loci()
tgl = gcfg$gene %>% select(gid,chrom,start,end,srd)
tga = read_loci() %>% select(gid,symbol,note)
txdb = load_txdb('B73', primary=F)
#{{{ functions
require(cluster)
tstat <- function(x, drc) {
    #{{{
    thresh = ifelse(drc == 'down', -1, 1)
    xs = x %>% mutate(clu=ifelse(log2fc < thresh, 0, 1)) %>%
        group_by(clu) %>% summarise(vals=list(log2fc)) %>% ungroup() %>% pull(vals)
    if (length(xs) <= 1) {
        0
    } else if( min(length(xs[[1]]), length(xs[[2]])) <= 1) {
        0
    } else {
        abs(t.test(xs[[1]], xs[[2]])$statistic[[1]])
    }
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
    distinct(cond,time,drc,gid)
td1s = td1 %>% arrange(cond,time,gt,gid,st) %>%
    group_by(cond,time,gt,gid) %>% dplyr::slice(1) %>% ungroup()
td4 = td3 %>% crossing(gt = gts3) %>%
    left_join(td1s, by = c('cond','time','gt','gid'))
td4s = td4 %>% mutate(drc1=str_to_upper(str_sub(drc,1,1))) %>%
    group_by(cond,time,drc,drc1,gid) %>%
    summarise(nd=sum(st==drc1), nu=sum(st!=drc1)) %>% ungroup() %>%
    filter(nd >0, nu>0) %>% select(cond,time,drc,drc1,gid)
td = td4 %>% inner_join(td4s, by=c('cond','time','drc','gid')) %>%
    spread(gt,st) %>% mutate(st=glue('{B73}{Mo17}{W22}')) %>%
    select(cond,time,drc,gid,st) %>% filter(cond=='cold')
#}}}

#{{{ read in NM matrix and run mclust - sf10a
#{{{ read in
yid = 'zm.rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm
#
th1 = th %>% filter(Experiment=='NM') %>%
    select(sid=SampleID, gt=Genotype, cond=Treatment, time=Timepoint) %>%
    mutate(cond = str_to_lower(cond)) %>%
    mutate(gt = ifelse(gt=='MS71', 'Ms71', gt))
#
rr = tm %>% select(gid, sid=SampleID, cpm=CPM) %>%
    inner_join(th1, by ='sid') %>%
    select(-sid) %>% spread(time, cpm)
rd = rr %>% gather(time, rc, -gid, -gt, -cond) %>%
    mutate(time = as.numeric(time)) %>%
    spread(cond, rc) %>%
    mutate(log2fc = log2(cold/control)) %>%
    arrange(gid, gt)
#}}}

#{{{ HY & NM consistency - sf10a
fi = glue('{dird}/15_de/01.de.rds')
ti = readRDS(fi)

fc1 = ti %>% filter(Genotype %in% gts3, Treatment=='Cold', Timepoint==25, TimepointB == 25) %>%
    select(gt=Genotype,ds) %>% unnest(ds) %>%
    select(gid, gt, log2fc.hy=log2fc)
fc2 = rd %>% filter(gt %in% gts3, time==25) %>%
    select(gid,gt,log2fc.nm=log2fc)
tp = fc1 %>% inner_join(fc2, by=c('gid','gt')) %>%
    filter(!is.nan(log2fc.nm), is.finite(log2fc.nm)) %>%
    rename(x=log2fc.hy, y=log2fc.nm)
lfc.fix = tp %>% rename(lfc.hy=x,lfc.nm=y) %>% mutate(time=25)
#
tpl = tp %>% group_by(gt) %>% nest() %>% ungroup() %>%
    mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
    mutate(tidied = map(fit, glance))%>%
    unnest(tidied) %>%
    mutate(lab = glue("adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}")) %>%
    print(width=Inf)
fcm = 5
tp = tp %>%
    mutate(x = pmax(pmin(x,fcm), -fcm)) %>%
    mutate(y = pmax(pmin(y,fcm), -fcm))
#{{{ plot
p = ggplot(tp, aes(x=x,y=y)) +
    geom_hex(bins=80) +
    geom_richtext(data=tpl, aes(x=-5,y=5,label=lab,hjust=0,vjust=1), size=2.5) +
    scale_x_continuous(name='mean log2fc in experiment 1', limits=c(-fcm,fcm), expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='log2fc in experiment 3', limits=c(-fcm,fcm), expand=expansion(mult=c(.02,.02))) +
    scale_fill_viridis(name='density', direction=-1) +
    #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    facet_wrap(gt~., nrow=3) +
    otheme(legend.pos='none', legend.title=T, legend.dir='v',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
#}}}
fo = glue("{dirw}/08.fc.pdf")
ggsave(p, file=fo, width=3, height=7)
fo = glue("{dirf}/sf10a.rds")
saveRDS(p, fo)
#}}}

#{{{ mclust
# replace lfc.nm with lfc.hy
#x = rd %>% left_join(lfc.fix, by=c('gid','gt','time')) %>%
    #mutate(log2fc = ifelse(is.na(lfc.hy), log2fc, lfc.hy)) %>%
    #select(-lfc.hy,-lfc.nm)
rd2 = td %>% inner_join(rd, by=c('time','gid')) %>%
    arrange(drc,time,gid,gt) %>%
    group_by(drc,time,st,gid) %>% summarise(log2fcs = list(log2fc)) %>%
    ungroup()
#
require(mclust)
#xs=rd2$log2fcs[[i]]
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
    #if (drc == 'up') {
        #clu = clu %>% mutate(nclu = ifelse(avg < 1, 0, 1))
    #} else if (drc == 'down') {
        #clu = clu %>% mutate(nclu = ifelse(avg <= -1, -1, 0))
    #}
    t_lfc = tibble(gt=names(xs1), log2fc=as.numeric(xs1))
    t_clu = tibble(gt=names(classes), clu=as.integer(classes)) %>%
        #inner_join(clu %>% select(clu,nclu), by='clu') %>%
        #select(-clu) %>% rename(clu=nclu) %>%
        inner_join(t_lfc, by='gt')
    clu = t_clu %>% group_by(clu) %>%
        summarise(avg = mean(log2fc), prop=n()/nrow(t_clu)) %>% ungroup()
    ebic = mod$BIC[,'E']
    vbic = mod$BIC[,'V']
    t_bic = tibble(nc=as.integer(names(ebic)), bic=as.numeric(ebic), vbic=as.numeric(vbic))
    tibble(nc=nrow(clu), t_bic=list(t_bic), clu=list(clu), t_clu=list(t_clu))
    #}}}
}
pb <- progress_bar$new(total=nrow(rd2))
rd3 = rd2 %>% mutate(i=1:n()) %>%
    mutate(x=pmap(list(i,log2fcs,drc), bimod,gts25=gts25,gts3=gts3)) %>%
    select(-i,-log2fcs) %>% unnest(x)
#}}}

#{{{ refine mclust results
refine_bimod0 <- function(clu, t_clu) {
    #{{{
    pb$tick()
    assign_clu <- function(x) tibble(clu=c(-1,0,1), d=c(x+1,x,x-1)) %>%
        mutate(d = abs(d)) %>% arrange(d) %>% pluck('clu',1)
    clu = clu %>% mutate(nclu = map_dbl(avg, assign_clu))
    t_clu = t_clu %>%
        inner_join(clu %>% select(clu,nclu), by='clu') %>%
        select(-clu) %>% rename(clu=nclu)
    clu = t_clu %>% group_by(clu) %>%
        summarise(avg = mean(log2fc), prop=n()/nrow(t_clu)) %>% ungroup()
    dn = -1 %in% clu$clu
    nd = 0 %in% clu$clu
    up = 1 %in% clu$clu
    tibble(nc1=nrow(clu), nd=!!nd, dn=!!dn, up=!!up,
           clu1=list(clu), t_clu1=list(t_clu))
    #}}}
}
assign_clu <- function(x, drc) {
    #{{{
    cent = ifelse(drc=='down', -1, 1)
    tibble(clu=c(0,1), d=c(x,x-cent)) %>%
        mutate(d = abs(d)) %>% arrange(d) %>% pluck('clu',1)
    #}}}
}
refine_bimod <- function(drc, clu, t_clu) {
    #{{{
    pb$tick()
    clu = clu %>% mutate(nclu = map2_dbl(avg, drc, assign_clu))
    t_clu = t_clu %>%
        inner_join(clu %>% select(clu,nclu), by='clu') %>%
        select(-clu) %>% rename(clu=nclu)
    n0 = sum(t_clu$clu==0)
    n1 = sum(t_clu$clu==1) 
    tibble(n0=!!n0, n1=!!n1, t_clu1=list(t_clu))
    #}}}
}
rd4 = rd3 %>% filter(nc>=2)
pb <- progress_bar$new(total=nrow(rd4))
rd4 = rd4 %>%
    mutate(x = pmap(list(drc,clu,t_clu), refine_bimod)) %>% unnest(x)
rd5 = rd4 %>% filter(time==25, n0>0,n1>0)
#}}}

r = list(rd3=rd3, rd4=rd4, rd5=rd5)
fo = glue('{dirw}/01.mclust.rds')
saveRDS(r, fo)
#}}}

fi = glue('{dirw}/01.mclust.rds')
r = readRDS(fi)
rd3 = r$rd3; rd5 = r$rd5

#{{{ uni- to bi- model composition - sf10b
clus = c('uni-modal','bi-modal', 'multi-modal (>=3)')
tp = rd3 %>% mutate(clu = ifelse(nc==1, clus[1],
                         ifelse(nc==2, clus[2], clus[3]))) %>%
    filter(time==25) %>%
    mutate(time = glue("cold_{time}h")) %>%
    count(drc, clu) %>%
    mutate(drc = factor(drc, levels=c('up','down'))) %>%
    mutate(clu = factor(clu, levels=clus)) %>%
    rename(tag1=drc,tag2=clu)
#
pa = cmp_proportion1(tp, ytext=T, xangle=0, oneline=F, legend.title='') +
    o_margin(.3,.3,.3,.3) +
    theme(legend.position='none')
fo = glue("{dirw}/10.modal.dist.pdf")
ggsave(pa, file=fo, width=3, height=5)
fo = glue("{dirf}/sf10b.rds")
saveRDS(pa, fo)
#}}}

#{{{ example of uni- and bi- model genes - sf10c
plot_lfc_cluster <- function(ti, tit, drc, gts3) {
    #{{{
    tp = ti %>% mutate(clu=factor(clu)) %>%
        mutate(gt0 = ifelse(gt %in% gts3, gt, 'other')) %>%
        mutate(gt0 = factor(gt0, levels=c(gts3,'other'))) %>%
        mutate(tit = !!tit)
    bw = diff(range(tp$log2fc)) / 15
    off = ifelse(drc == 'down', -1, 1)
    cols4 = c(pal_aaas()(length(unique(tp$gt0))-1),'white')
    ggplot(tp, aes(x=log2fc, fill=gt0)) +
        geom_dotplot(method='histodot', binwidth=bw, stackgroups=T, dotsize=1) +
        geom_vline(xintercept=off, size=.5, color='gray', linetype='dashed') +
        scale_fill_manual(values=cols4) +
        scale_color_manual(values=cols4) +
        facet_wrap(tit ~., scale='free', nrow=1) +
        otheme(xtext=T, xtick=T, xtitle=T,
            legend.pos='top.right') +
        o_margin(0,.3,0,.3)
    #}}}
}

#{{{ explore
tp0 = rd5 %>%
    #filter(drc=='down') %>% slice(31:39)
    filter(drc=='up') %>% slice(181:189)
    #filter(nc1 == 1, drc=='up') %>% slice(11:19)
tp0 = rd3 %>%
    filter(nc == 1, time==25, drc=='up') %>% slice(51:59) %>%
    rename(t_clu1=t_clu)
tp = tp0 %>%
    mutate(tit=glue("{st} {gid}")) %>%
    mutate(p = pmap(list(t_clu1, tit,drc), plot_lfc_cluster, gts3=gts3))
ps = tp$p
fo = glue("{dirw}/10.modal.pdf")
ggarrange(plotlist=ps, nrow=3, ncol=3) %>%
    ggexport(filename=fo, width=7, height=7)
#}}}

#final
gids = c('Zm00001d011276', 'Zm00001d012312',
    'Zm00001d007630','Zm00001d040011',
    'Zm00001d038577','Zm00001d045565')
gids=gids[c(2,5)]
tp1 = rd3 %>% filter(gid %in% gids,time==25,drc=='up') %>% select(drc,gid,st,t_clu1=t_clu)
gids = c('Zm00001d017130','Zm00001d007890',
         'Zm00001d032811', 'Zm00001d003460')
gids=gids[c(2,4)]
tp2 = rd5 %>% filter(gid %in% gids) %>% select(drc,gid,st,t_clu1)
tp = tp1 %>% mutate(tag='unimodal') %>%
    bind_rows(tp2 %>% mutate(tag='bimodal')) %>%
    mutate(tit=glue("{tag} {gid} ({st})")) %>%
    mutate(p = pmap(list(t_clu1, tit,drc), plot_lfc_cluster, gts3=gts3))
ps = tp$p
fo = glue("{dirw}/10.modal.pdf")
pb = ggarrange(plotlist=ps, nrow=2, ncol=2)
pb %>% ggexport(filename=fo, width=5, height=5)
fo = glue("{dirf}/sf10c.rds")
saveRDS(pb, fo)
#}}}

#{{{ read variable response
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
#{{{ local association test
#{{{
i=37
drc=x$drc[i]; t_clu=x$t_clu1[[i]]
chrom=x$chrom[i]; start=x$start[i]; end=x$end[[i]]
fet <- function(x00,x01,x10,x11) fisher.test(matrix(c(x00,x10,x01,x11), nrow=2))$p.value
read_vnt <- function(chrom,beg,end,
                     fv='~/projects/wgc/data/05_maize_v4/08.vcf.gz',
                     gts=gts31_ph207) {
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
vnt_scan <- function(t_clu, chrom,start,end) {
    #{{{
    #pb$tick()
    pheno = t_clu %>% select(sid=gt,clu)
    vnt = read_vnt(chrom,start,end)
    if(nrow(vnt)==0) return(tibble())
    vnt0 = vnt %>% distinct(chrom,pos,ref,alt) %>%
        mutate(snp=glue("v{1:n()}")) %>%
        mutate(snp = as_factor(snp))
    y = vnt %>% mutate(st = ifelse(is.na(gt), '0 0',
        ifelse(gt==0, "A A", "T T"))) %>%
        inner_join(vnt0, by=c('chrom','pos','ref','alt')) %>%
        select(chrom,snp,pos, sid,st)
    t_ped = y %>% select(-chrom,-pos) %>%
        spread(snp,st) %>% left_join(pheno, by='sid') %>%
        replace_na(list(clu=-9)) %>% select(sid,clu,everything())
    t_map = y %>% distinct(chrom,snp,pos)
    #
    pre = glue("vs{sample(1e5,1)}")
    f_ped = glue("{pre}.ped")
    write.table(t_ped, f_ped, sep=" ", row.names=F, col.names=F, quote=F)
    f_map = glue("{pre}.map")
    write_delim(t_map, f_map, delim=" ", col_names=F, quote_escape=F)
    #
    system(glue("plink --file {pre} --no-fid --no-parents --no-sex --1 --allow-extra-chr --blocks --blocks-min-maf 0 --blocks-strong-highci 1 --blocks-inform-frac 0.8"))
    #system(glue("plink --file {pre} --no-fid --no-parents --no-sex --1 --allow-extra-chr --indep-pairwise 50 5 0.99"))
    #system(glue("plink --file {pre} --no-fid --no-parents --no-sex --1 --allow-extra-chr --assoc --adjust --extract plink.prune.in"))
    system(glue("plink --file {pre} --no-fid --no-parents --no-sex --1 --allow-extra-chr --assoc --adjust"))
    system(glue("rm {pre}*"))
    #{{{ haplo block
    if (file.info("plink.blocks")$size == 0) {
        tb = tibble()
    } else {
        fi = "plink.blocks.det"
        tb = read_delim(fi, delim=' ', trim_ws=T, col_names=T) %>%
            select(snps=SNPS) %>% mutate(bid=1:n()) %>%
            mutate(snps = map(snps, str_split, pattern="[\\|]")) %>%
            mutate(snps = map(snps, 1)) %>% unnest(snps) %>% rename(snp=snps) %>%
            inner_join(vnt0, by='snp') %>%
            group_by(bid) %>%
            summarise(start=min(pos), end=max(pos),
                      snp1=min(snp), snp2=max(snp)) %>%
            ungroup()
    }
    #}}}
    fi = "plink.assoc.adjusted"
    tr = read_delim(fi, delim=' ', trim_ws=T, col_names=T) %>%
        select(snp=SNP, p.raw=UNADJ, p.bon=BONF, p.fdr=FDR_BH)
    p = vnt0 %>% inner_join(tr, by='snp') %>%
        select(chrom,pos,ref,alt,snp,everything())
    nsig = sum(p$p.bon < 0.05)
    nsig.fdr = sum(p$p.fdr < 0.05)
    tibble(vnt=list(vnt), p=list(p), block=list(tb), nsig=nsig, nsig.fdr=nsig.fdr)
    #}}}
}

x = rd5 %>%# slice(1:2) %>% #filter(gid == !!gid) %>%
    inner_join(tgl, by='gid') %>%
    mutate(start = start - 2000, end = end + 2000)
#pb <- progress_bar$new(total=nrow(x))
la = x %>%
    mutate(data = pmap(list(t_clu1, chrom,start,end), vnt_scan)) %>%
    select(drc,time,st,gid,data) %>% unnest(data)

fo = glue("{dirw}/51.assoc.rds")
saveRDS(la, fo)
#}}}

fi = glue("{dirw}/51.assoc.rds")
la = readRDS(fi)
la %>% count(drc, nsig > 0)

# add reg
ddeg3 = ddeg2 %>% filter(cond=='cold',tgt=='B73') %>% select(time,gid,reg)
ti2 = la %>%
    mutate(time=as.integer(time)) %>%
    mutate(sig = ifelse(nsig > 0,'sig','insig')) %>%
    left_join(ddeg3, by=c("time",'gid')) %>%
    filter(!is.na(reg))
ti2 %>% count(reg, sig) %>% spread(sig,n) %>% mutate(prop=sig/(sig+insig))

#{{{ plot cis/trans prop - 
regs = c("bkg","cis","cis+trans",'trans','conserved','unexpected')
sigs = c('sig','insig')
tp = ti2 %>% count(reg, sig) %>% rename(tag1=reg, tag2=sig) %>%
    mutate(tag1=factor(tag1, levels=regs)) %>%
    mutate(tag2 = factor(tag2, levels=sigs))
p = cmp_proportion1(tp, ytext=T, xangle=15, lab.size=2) +
    o_margin(.3,.3,.3,.3) +
    theme(legend.position='none')
#}}}
fo = glue("{dirw}/17.reg.pdf")
ggsave(p, file=fo, width=4, height=4)

# add symbol & note
ti3 = ti2 %>% inner_join(tga, by='gid')
j = ti3

j %>% filter(nsig>0,reg %in% c('cis','cis+trans')) %>%
    select(gid,p,nsig,reg,symbol) %>% print(n=40,width=Inf)
#}}}

#{{{ ##(obsolete) sf13
#{{{ read
fi = glue('{dirw}/01.mclust.rds')
r = readRDS(fi)
rd3 = r$rd3; rd5 = r$rd5
rd5b = rd5 %>% select(drc,time,gid,clu=t_clu1)
#
fi = glue("{dirw}/05.test.rds")
ti = readRDS(fi)
#
fi = glue("{dird}/25_dreme/05.best.mtfs.rds")
r = readRDS(fi)
tk = r$tk %>% filter(bid=='b01') %>% select(mtfs) %>% unnest(mtfs)
fm = glue("{dird}/41_ml/00_nf/03_motif_lists/b01.meme")
nfea = 'top30'
#}}}
#{{{ functions
plot_vcf <- function(vnt) {
    #{{{ snp haplotype plot
    cols_vcf = c(pal_aaas()(2), 'white')
    tp = vnt %>% replace_na(list(gt=9)) %>% mutate(gt = factor(gt)) %>%
        mutate(x = as_integer(factor(pos)))
    ggplot(tp) +
        geom_tile(aes(x=x, y=sid, fill=gt), alpha=1) +
        scale_x_continuous(expand=expansion(mult=c(.0,.0))) +
        scale_y_discrete(expand=expansion(mult=c(.0,.0))) +
        scale_fill_manual(values = cols_vcf) +
        otheme(legend.pos='none', panel.border=F, ytext=T,ytick=T,xtext=F)
    #}}}
}
fimo_locate <- function(gid,fm,nfea) {
    #{{{
    tg = tibble(gid = glue("B73_{gid}"), status = 1)
    write_tsv(tg, 'tmp.tsv')
    cmd = glue("fimo.py prepare_ml tmp.tsv {fm} tmp.bed --epi raw --nfea {nfea} --fmt long")
    system(cmd)
    ti = read_tsv('tmp.bed', col_names=c("gid",'start','end','mid'))
    system("rm tmp.tsv tmp.bed")
    ti
    #}}}
}
get_mtf_loc <- function(gid, fm, nfea, tk, srd, pb) {
    #{{{
    x = fimo_locate(gids[1], fm, nfea)
    x1 = x %>% inner_join(tk %>% select(mid,fid,fname,known,conseq), by='mid') %>%
        mutate(i = 1:n()) %>%
        mutate(fname = ifelse(known, fname, fid)) %>%
        select(gid,start,end,i,known,fname,conseq)
    x1m = bed_merge(x1) %>%
        inner_join(x1 %>% select(i,known,fname,conseq), by='i') %>%
        arrange(chr,beg,end, desc(known), fname) %>%
        group_by(chr,beg,end) %>%
        summarise(txt = str_c(fname,collapse=',')) %>% ungroup() %>%
        mutate(pos = (beg+end)/2) %>%
        select(gid=chr,pos,txt)
    if (srd=='-') {
        x1m %>% mutate(pos=pe-pos)
    } else {
        x1m %>% mutate(pos=pb+pos)
    }
    #}}}
}
plot_box <- function(clu,assoc,vnt, gts3) {
    #{{{ boxplot
    tp1 = assoc %>% arrange(p.bon) %>% slice(1) %>%
        select(chrom,pos,ref,alt)
    xtit = tp1 %>%
        mutate(ref = ifelse(nchar(ref)>=5, glue('[{nchar(ref)}]'), ref)) %>%
        mutate(xtit=glue("{chrom}:{pos}:{ref}/{alt}")) %>%
        pluck('xtit',1)
    tp = tp1 %>%
        inner_join(vnt, by=c('chrom','pos','ref','alt')) %>%
        inner_join(clu, by=c('sid'='gt')) %>%
        select(sid,gt,clu,log2fc) %>%
        arrange(clu) %>%
        mutate(gtc = ifelse(sid %in% gts3, sid, 'other')) %>%
        mutate(gtc = factor(gtc, levels=c(gts3,'other')))
    #
    ggplot(tp, aes(x=gt, y=log2fc)) +
        geom_boxplot(aes(group=gt), size=.5, width=.7, outlier.shape=NA) +
        geom_jitter(aes(color=gtc), width=.3, size=2, alpha=.9) +
        #geom_text(data=tps, aes(x=st,y=500*1.01,label=lab), size=2.5,vjust=0) +
        scale_x_continuous(name=xtit,breaks=c(0,1),labels=c('ref','alt'),expand=expansion(mult=c(.2,.2))) +
        scale_y_continuous(name='log2FoldChange', expand=expansion(mult=c(.02,.02))) +
        scale_color_aaas() +
        #ggtitle(tp$tit[1]) +
        otheme(ytitle=T,ytext=T,ytick=T, xtitle=T,xtick=T, xtext=T,
               legend.pos='top.center.out',legend.dir='h', legend.vjust=-.4,
               panel.spacing=.2, strip.compact=F) +
        o_margin(2,.3,.3,.1) +
        guides(color=guide_legend(nrow=2)) +
        theme(axis.title.x=element_text(size=7.5)) +
        theme(plot.title=element_text(size=8))
    #}}}
}
#}}}

gids = c('Zm00001d006065','Zm00001d047307','Zm00001d003252','Zm00001d024425',
    'Zm00001d052653','Zm00001d018794','Zm00001d025016','Zm00001d020970',
    'Zm00001d041920','Zm00001d021891')
gids = c('Zm00001d003252')
gids = c('Zm00001d024425')
gids = c('Zm00001d017130')
gids = c('Zm00001d020970') # sf13b
gids = c('Zm00001d052653') # sf13a
plot_gene_vnt <- function(gid,pb,pe,srd,vnt,assoc,blk,tit,clu,gts,fo,
                          txdb,fm,nfea,tk) {
    #{{{
    tg = txdb_gene(txdb, gid, gb=0)
    pg = plot_genes(tg, pb, pe, label.top=T)
    tm = get_mtf_loc(gid, fm, nfea, tk, srd, pb)
    #{{{ prepare
    assoc2 = assoc %>% mutate(sig=p.bon < 0.05) %>%
        select(chrom,pos,ref,alt,sig)
    vnt0 = vnt %>% distinct(chrom,pos,ref,alt) %>%
        mutate(snp=glue("v{1:n()}")) %>%
        mutate(snp = as_factor(snp))
    tp0 = vnt %>%
        inner_join(vnt0, by=c('chrom','pos','ref','alt')) %>%
        mutate(i = as_integer(snp))
    tp0w = tp0 %>% select(snp,sid,gt) %>% spread(snp,gt)
    sids = hc_order_row(tp0w, cor.opt='gower')
    tp = tp0 %>% replace_na(list(gt=9)) %>% mutate(gt = factor(gt)) %>%
        mutate(sid = factor(sid, levels=sids[sids %in% clu$gt]))
    tpy = tp %>% distinct(sid) %>%# filter(sid %in% gts) %>%
        #mutate(sid = factor(sid, levels=gts[gts %in% tp$sid])) %>%
        mutate(y = -1 -as_integer(sid) * .5)
    tp = tp %>% inner_join(tpy, by='sid') %>%
        mutate(x = pb + (pe-pb)/max(tp$i) * i)
    tpx = tp %>% distinct(i, snp, chrom, pos, ref, alt, x) %>%
        left_join(assoc2, by=c("chrom",'pos','ref','alt')) %>%
        replace_na(list(sig=F))
    xmin=min(tpx$pos,tpx$x); xmax = max(tpx$pos,tpx$x)
    tpx0 = tpx %>% select(snp,pos,x)
    #}}}
    #{{{ block
    tpb = blk %>%
        inner_join(tpx0, by=c('snp1'='snp')) %>%
        rename(pos1=pos,x1=x) %>%
        inner_join(tpx0, by=c('snp2'='snp')) %>%
        rename(pos2=pos,x2=x) %>%
        select(snp1,snp2,pos1,pos2,x1,x2) %>%
        mutate(hcol = rep(1:5, as.integer(nrow(blk)/5)+1)[1:nrow(blk)]) %>%
        mutate(hcol = factor(hcol))
    tpt = clu %>% rename(sid=gt) %>% inner_join(tpy, by='sid') %>%
        mutate(txt = str_pad(number(log2fc, accuracy=1), width=2, pad=' '))
    #}}}
    #{{{ plot
    cols_vcf = c(pal_simpsons()(2), 'white')
    seg.sizes = c(.1,.5); blk.size = .7
    seg.cols = c('gray','orange red')
    blk.col = pal_lancet()(3)[1]
    pa = pg +
        geom_point(data=tm, aes(x=pos,y=.3), size=1, color='red') +
        geom_text(data=tm, aes(x=pos,y=.3,label=txt), size=1.5,vjust=.5,angle=30) +
        geom_segment(data=tpx, aes(x=pos,y=0,xend=pos,yend=-.5, color=sig, size=sig), lineend='round') +
        geom_segment(data=tpx, aes(x=pos,y=-.5,xend=x,yend=-1.1, color=sig, size=sig), lineend='round') +
        geom_segment(data=tpb, aes(x=pos1,xend=pos2,y=0,yend=0),color=blk.col,size=blk.size) +
        geom_segment(data=tpb, aes(x=x1,xend=x2,y=-1.1,yend=-1.1),color=blk.col,size=blk.size) +
        geom_text(data=tpt, aes(x=xmax, y=y, label=txt),size=2,hjust=0) +
        geom_tile(data=tp, aes(x=x, y=y, fill=gt), alpha=1) +
        scale_x_continuous(label=label_number(accuracy=1),
                           position='top',
                           limits=c(pb,pe),expand=expansion(mult=c(0,.03))) +
        scale_y_continuous(breaks=tpy$y,labels=tpy$sid,
                           expand=expansion(add=c(0,.4),mult=c(0,0))) +
        scale_color_manual(values = seg.cols) +
        scale_fill_manual(values = cols_vcf) +
        scale_size_manual(values = seg.sizes) +
        otheme(legend.pos='none', panel.border=F,
               ytext=T,ytick=T,xtext=T,xtick=T,
               margin=c(0,0,0,0)) +
        theme(axis.text.x = element_text(size=6)) +
        theme(axis.text.y = element_text(size=6))
    #}}}
    pb = plot_box(clu,assoc,vnt, gts3) + ggtitle(tit)
    p = ggarrange(pa, pb, nrow=1, ncol=2, widths=c(3,1))
    p %>% ggexport(filename=fo, width=8, height=3)
    p
    #}}}
}
tp0 = ti %>% filter(gid %in% gids) %>% rename(assoc = p) %>%
    inner_join(rd5b, by=c('drc','time','gid')) %>%
    inner_join(tgl, by='gid') %>%
    mutate(pb=start-2000, pe=end+2000) %>%
    mutate(tit=glue("[B|M|W]={st}")) %>%
    mutate(fo = glue("{dirw}/{gid}.pdf"))
i=1; gid=tp0$gid[i];pb=tp0$pb[i];pe=tp0$pe[i];srd=tp0$srd[i]
vnt=tp0$vnt[[i]];assoc=tp0$assoc[[i]]; blk=tp0$block[[i]]
tit=tp0$tit[i];clu=tp0$clu[[i]];gts=gts25;fo=tp0$fo[i]
tp = tp0 %>% mutate(gts = list(gts25)) %>%
    mutate(p = pmap(list(gid, pb, pe, srd, vnt,assoc,block, tit,clu,gts,fo),
                    plot_gene_vnt, txdb=!!txdb,
                    fm=!!fm, nfea=!!nfea, tk=!!tk)) %>%
    print(width=Inf)

fo = glue("{dirf}/f7c.rds")
saveRDS(tp$p[[1]], fo)

fo = glue("{dirf}/sf13a.rds")
saveRDS(tp$p[[1]], fo)
fo = glue("{dirf}/sf13b.rds")
saveRDS(tp$p[[1]], fo)
#}}}


#{{{ eval model pred in gts26 cold_25h panel
#{{{ write gene list
x0 = rd5 %>% select(drc,gid,t_clu1) %>% unnest(t_clu1) %>%
    mutate(gt = ifelse(gt=='MS71', 'Ms71', gt)) %>%
    mutate(status = clu) %>%
    mutate(gid = glue("{gt}_{gid}")) %>%
    select(drc, gid, status)

fo = glue("{dirw}/31.genes.rds")
saveRDS(x0, fo)
to = x0 %>% distinct(gid) %>% mutate(status=1)
fo = glue("{dirw}/31.genes.tsv")
write_tsv(to, fo)
#}}}

#{{{ extract models
fi = glue('{dird}/41_ml/06.tk.tc.rds')
r6 = readRDS(fi)
fb = glue('{dird}/41_ml/12.best.models.tsv')
tb = read_tsv(fb)
tm = tb %>% inner_join(r6$tc, by=c('bid','tid')) %>%
    select(-gids, -gids_c, -ts) %>%
    filter(cond=='cold', str_starts(note, 'all'))

fo = glue("{dirw}/32.models.rds")
saveRDS(tm, fo)
to = tm %>% select(tag,tid,bid,epi)
fo = glue("{dirw}/32.models.tsv")
write_tsv(to, fo)
#}}}

# run j39 fimo.py prepare_ml
# run j39 ml_predict

fi = glue("{dirw}/32.models.rds")
tm = readRDS(fi)
to = tm %>%# select(tag, tid, bid) %>%
    mutate(fi = glue("{dirw}/37_ml_out/{tag}_{tid}_{epi}.tsv")) %>%
    mutate(pred = map(fi, read_tsv)) %>% select(-fi) %>%
    mutate(drc = ifelse(str_detect(note, 'down'), 'down', 'up')) %>%
    select(tid,bid,epi,train,drc,pred) %>% unnest(pred)

fg = glue("{dirw}/31.genes.rds")
tg = readRDS(fg)

pred = tg %>% inner_join(to, by=c('drc','gid')) %>%
    separate(gid, c('gt','gid'), sep='_') %>% select(-prob) %>%
    inner_join(rd5 %>% select(drc,st,gid,t_clu1) %>% unnest(t_clu1),
        by=c('drc','gid','gt')) %>%
    mutate(mdl=glue("{train}_{epi}")) %>%
    select(gid,drc,st, gt, log2fc,clu,status, mdl, pred) %>%
    spread(mdl, pred)
pred %>% print(n=30,width=Inf)

to3 = pred %>% rename(pred=B_umr) %>%
    count(gid, drc, status, pred) %>%
    mutate(p = glue("n{status}{pred}")) %>%
    select(-status, -pred) %>%
    spread(p, n) %>%
    replace_na(list(n00=0,n01=0,n10=0,n11=0)) %>%
    mutate(acc = n00/(n00+n01)/2 + n11/(n10+n11)/2)
to3 %>% filter(acc >= .8) %>% print(n=40)

gid='Zm00001d017363'
pred %>% filter(gid==!!gid) %>% select(-gid,-BMW_nr_raw,-BMW_nr_umr) %>% print(n=30)

fo = glue('{dirw}/25.model.pred.rds')
saveRDS(pred, fo)
#}}}

#{{{ identify signif. associated motifs
fi = glue('{dirw}/32.models.rds')
ti = readRDS(fi) %>%
    mutate(fi=glue("{dirw}/36_ml_in/{bid}_{epi}.tsv")) %>%
    mutate(ti = map(fi, read_tsv)) %>%
    select(train,cond,note,epi,ti)
fi = glue('{dirw}/25.model.pred.rds')
pred = readRDS(fi)

spread2 <- function(ti) ti %>% select(-status) %>% separate(gid,c('gt','gid'), sep='_') %>% gather(mid, pav, -gid, -gt)
ti2 = ti %>% mutate(drc=ifelse(str_detect(note,'down'), 'down', 'up')) %>%
    filter(train=='B') %>%
    mutate(x = map(ti, spread2)) %>%
    select(epi,drc,x) %>% unnest(x)

#{{{ model performance
tk = pred %>% select(drc,gid, gt, status,raw=B_raw, umr=B_umr) %>%
    gather(epi, pred, -drc,-gid,-status,-gt) %>%
    group_by(gt, epi) %>%
    summarise(n_tot=n(), n1=sum(status==pred)/n_tot) %>% ungroup() %>%
    select(-n_tot) %>% spread(epi, n1) %>% print(n=40)

tk = pred %>% select(drc,gid,status,gt,raw=B_raw, umr=B_umr) %>%
    gather(epi, pred, -drc,-gid,-status,-gt) %>%
    count(epi, gid, drc, status, pred) %>%
    mutate(p = glue("n{status}{pred}")) %>% select(-status, -pred) %>%
    spread(p, n) %>%
    replace_na(list(n00=0,n01=0,n10=0,n11=0)) %>%
    mutate(acc = n00/(n00+n01)/2 + n11/(n10+n11)/2) %>% filter(!is.na(acc)) %>%
    arrange(epi, desc(acc)) %>% print(n=40)

tk %>% count(epi, acc > 0.9)
tk %>% count(epi, acc > 0.8)
#}}}

#{{{ motif association
pred1 = pred %>% group_by(drc, gid) %>% summarise(gts = list(gt)) %>% ungroup()
tj = ti2 %>% inner_join(pred, by=c('drc','gid','gt')) %>%
    count(epi, gid, drc, mid, status, pav) %>%
    mutate(p = glue("n{status}{pav}")) %>%
    select(-status, -pav) %>%
    spread(p, n) %>%
    replace_na(list(n00=0,n01=0,n10=0,n11=0)) %>%
    mutate(acc = n00/(n00+n01)/2 + n11/(n10+n11)/2) %>% filter(!is.na(acc))

get_lfc  <- function(gid, pav, pred) {
    #{{{
    pred %>% filter(gid == !!gid) %>%
        select(gt,clu,log2fc,B_raw,B_umr) %>%
        mutate(txt = ifelse(is.na(log2fc), 'NA', number(log2fc,accuracy=.01))) %>%
        left_join(pav, by='gt') %>%
        mutate(tag_de = ifelse(clu==1, "x", "")) %>%
        mutate(tag_mtf = ifelse(pav==1, "x", ""))
    #}}}
}
tj2 = tj %>% arrange(epi, gid, drc, desc(acc)) %>%
    group_by(epi, gid, drc) %>% dplyr::slice(1:1) %>% ungroup()
mtf_pav = ti2 %>% inner_join(tj2 %>% select(epi,gid,drc,mid), by=c('epi','drc','gid','mid')) %>%
    group_by(epi,gid,drc,mid) %>% nest() %>% rename(pav=data) %>% ungroup()
tj2 = tj2 %>% inner_join(mtf_pav, by=c('epi','drc','gid','mid')) %>%
    inner_join(pred1, by=c('gid','drc')) %>%
    mutate(lfc = map2(gid, pav, get_lfc, pred=pred)) %>% select(-pav) %>%
    arrange(epi, desc(acc)) %>% print(n=20)
tj2 %>% mutate(n_gt=map_dbl(gts,length)) %>% count(n_gt-n00-n01-n10-n11)
tj2 %>% count(epi, acc > 0.9)
tj2 %>% count(epi, acc > 0.8)

fo = glue("{dirw}/61.mtf.assoc.rds")
saveRDS(tj2, fo)
#}}}
#}}}

#{{{ final synteny plots - f7c-d & sf09
#{{{ read & functions
require(Biostrings)
get_fasta <- function(gid, gts, beg=0, end=4000, db=glue("{dird}/21_seq/02.fas")) {
    #{{{
    pre = glue("tmp.gf.{sample(10000, 1)}")
    tt = tibble(gid=!!gid, gt=gts) %>% mutate(sid=glue("{gt}_{gid}")) %>%
        mutate(beg=!!beg, end=!!end) %>% select(sid, beg, end)
    write_tsv(tt, file=glue("{pre}.bed"), col_names=F)
    system(glue("fasta.py extract {db} {pre}.bed > {pre}.fa"))
    y = readDNAStringSet(glue("{pre}.fa"))
    snames = tibble(sid = names(y)) %>%
        separate(sid, c("sid",'other'), sep='-', extra='merge') %>%
        separate(sid, c('gt', 'gid'), sep='_') %>%
        pull(gt)
    names(y) = snames
    #
    system(glue("rm {pre}.*"))
    y
    #}}}
}
msa_tree <- function(seqs) {
    #{{{
    pre = glue("tmp.msa.{sample(10000, 1)}")
    writeXStringSet(seqs, glue('{pre}.fa'))
    system(glue("muscle -in {pre}.fa -out {pre}.2.fa -tree2 {pre}.2.nwk"))
    msa = readDNAMultipleAlignment(filepath=glue("{pre}.2.fa"), format='fasta')
    tree = read.tree(glue("{pre}.2.nwk"))
    labs = with(subset(fortify(tree), isTip), label[order(y, decreasing=T)])
    nl = length(labs)
    ty = tibble(lab = labs, y = nl:1) %>%
        separate(lab, c('gt','extra'), sep='_', remove=F) %>%
        select(gt, lab, y)
    #
    system(glue("rm {pre}.*"))
    list(msa=msa, tree=tree, ty=ty)
    #}}}
}
get_mtf_loc <- function(mid, seqs, mtf_dir=glue("{dird}/41_ml/00_nf/03_motif_lists")) {
    #{{{
    pre = glue("tmp.mtf.{sample(10000, 1)}")
    writeXStringSet(seqs, glue('{pre}.fa'))
    bid = stringr::str_split(mid,'_')[[1]][1]
    system(glue("fimo.py locate --motif {mid} {mtf_dir}/{bid}.meme {pre}.fa {pre}.txt"))
    mtf = read_tsv(glue("{pre}.txt"), col_names=c("sid",'beg','end')) %>%
        separate('sid',c('mid','sid'), sep='%') %>%
        separate('sid',c('gt','extra'), sep='_') %>%
        mutate(pos = (beg+end)/2) %>%
        select(gt, beg, end, pos)
    #
    system(glue("rm {pre}.*"))
    mtf
    #}}}
}
get_umr_loc <- function(gid, umr=glue("{dird}/21_seq/15.UMR.bed")) {
    #{{{
    pre = glue("tmp.umr.{sample(10000, 1)}")
    system(glue("grep {gid} {umr} > {pre}.bed"))
    umr = read_tsv(glue("{pre}.bed"), col_names=c('gid','start','end')) %>%
        mutate(start=start+1) %>%
        separate("gid", c("gt",'gid'), sep='_') %>%
        select(gt, start, end)
    system(glue("rm {pre}.*"))
    umr
    #}}}
}
fi = glue("{dirw}/61.mtf.assoc.rds")
ta = readRDS(fi)
#}}}

#{{{ extract seqs & motifs - 62.rds
ta2 = ta %>% filter(epi=='raw', acc >= .8)
ta2 = ta2 %>%
    #dplyr::slice(1:2) %>%
    mutate(seqs = map2(gid, gts, get_fasta)) %>%
    mutate(mtf = map2(mid,seqs, get_mtf_loc))

fo = glue("{dirw}/62.seqs.mtf.rds")
saveRDS(ta2, fo)
#}}}

#{{{ characterize motif variation - 63.rds
fi = glue('{dirw}/../41_ml/06.tk.tc.rds')
r06 = readRDS(fi)
loci = read_loci()
#
fi = glue("{dirw}/62.seqs.mtf.rds")
ta2 = readRDS(fi)
ti = ta2
j=3; mtf=ti$mtf[[j]]; seqs=ti$seqs[[j]]; gts=ti$gts[[j]]; lfc=ti$lfc[[j]]
mtf_aln  <- function(n1, n2, tb, te, seqs) {
    #{{{
    cat(n1, n2, '\n')
    pw = pairwiseAlignment(seqs[[n1]],seqs[[n2]], type='local')
    toff = start(pattern(pw)); qoff = start(subject(pw))
    y1 = as.character(pattern(pw))
    y2 = as.character(subject(pw))
    h0 = as.character(compareStrings(y2,y1))
    h = str_replace_all(h0, "[ATCGN?]", '*')
    h2=unlist(strsplit(h, split = ""))
    j = tibble(v = rle(h2)$values, aSize = rle(h2)$lengths) %>%
        mutate(aEnd=cumsum(aSize), aBeg = aEnd-aSize) %>%
        mutate(tSize = ifelse(v=='+', 0, aSize)) %>%
        mutate(qSize = ifelse(v=='-', 0, aSize)) %>%
        mutate(tEnd=cumsum(tSize), tBeg=tEnd-tSize) %>%
        mutate(qEnd=cumsum(qSize), qBeg=qEnd-qSize) %>%
        select(v, aBeg,aEnd,aSize,tBeg,tEnd,tSize,qBeg,qEnd,qSize) %>%
        mutate(tBeg = toff+tBeg-1, tEnd=tBeg+tSize) %>%
        mutate(qBeg = qoff+qBeg-1, qEnd=qBeg+qSize)
    get_tpos <- function(v, aSize, tBeg, tEnd)
        ifelse(v=='+', list(rep(tBeg, aSize)), list((tBeg+1):tEnd))
    get_qpos <- function(v, aSize, qBeg, qEnd)
        ifelse(v=='-', list(rep(qBeg, aSize)), list((qBeg+1):qEnd))
    j2 = j %>%
        mutate(tpos = pmap(list(v, aSize, tBeg, tEnd), get_tpos)) %>%
        mutate(qpos = pmap(list(v, aSize, qBeg, qEnd), get_qpos))
    tposs = unlist(j2$tpos)
    qposs = unlist(j2$qpos)
    mm = str_locate_all(h0, "\\?")[[1]] %>% as_tibble() %>%
        mutate(qPos = qposs[start], tPos = tposs[start]) %>%
        select(aPos=start,tPos,qPos)
    ab = which(tposs == tb)[1] 
    ae = which(tposs == te)[1] 
    tseq = str_sub(y1, ab, ae); qseq=str_sub(y2,ab,ae)
    mm1 = mm %>% filter(aPos >= ab, aPos <= ae)
    if(nrow(mm1) > 0) {
        for (i in 1:nrow(mm1)) {
            mp = mm1$aPos[i] - ab + 1
            str_sub(qseq,mp,mp) <- str_to_lower(str_sub(qseq,mp,mp))
        }
    }
    qb=qposs[ab]; qe=qposs[ae]
    ins=str_count(tseq,'-')
    del=str_count(qseq,'-')
    tibble(tseq=tseq,qseq=qseq,qb=qb,qe=qe, mm=nrow(mm1),ins=ins,del=del)
    #}}}
}
make_msa <- function(seqs) {
    #{{{
    pre = glue("tmp.msa.{sample(10000, 1)}")
    writeXStringSet(seqs, glue('{pre}.fa'))
    system(glue("muscle -in {pre}.fa -out {pre}.2.fa"))
    msa = readDNAMultipleAlignment(filepath=glue("{pre}.2.fa"), format='fasta')
    msa = DNAStringSet(msa)
    #
    system(glue("rm {pre}.*"))
    msa
    #}}}
}
mark_mm  <- function(tseq, qseq) {
    #{{{ mark mismatch and indel in query seqs
    h0 = as.character(compareStrings(qseq,tseq))
    h2 = unlist(strsplit(h0, split = ""))
    th2 = tibble(p=1:length(h2), op = h2) %>%
        filter(op %in% c("?", "+"))
    if (nrow(th2) > 0) {
        for (p in th2$p) {
            str_sub(qseq,p,p) <- str_to_lower(str_sub(qseq,p,p))
        }
    }
    qseq
    #}}}
}
check_mtf <- function(j, mtf, seqs, gts) {
    #{{{
    cat(j, '\n')
    mtf0 = mtf %>% mutate(dst = abs(pos-2000)) %>% arrange(dst)
    i=1; gt0 = mtf0$gt[[i]]; b0 = mtf0$beg[[i]]; e0 = mtf0$end[[i]]
    seq0 = seqs[[gt0]][b0:e0]
    aln = tibble(gt0=!!gt0, gt1=gts[! gts %in% gt0]) %>%
        filter(gt1 %in% names(seqs)) %>%
        mutate(tb=b0, te=e0) %>%
        mutate(x = pmap(list(gt0,gt1,tb,te), mtf_aln, seqs=seqs)) %>%
        unnest(x)
    #
    aln_na = aln %>% filter(is.na(tseq))
    aln2 = aln %>% filter(!is.na(tseq))
    if (length(unique(aln2$tseq)) > 1) {
        ss = aln2 %>% select(gt=gt1, seq=qseq) %>%
            bind_rows(aln2 %>% dplyr::slice(1) %>% select(gt=gt0, seq=tseq)) %>%
            mutate(seq=str_replace(seq,'-',''))
        ss0 = ss$seq; names(ss0) = ss$gt
        ss1 = DNAStringSet(ss0)
        ss2 = make_msa(ss1)
        ss3 = as.character(ss2)
        aln2 = aln2 %>% mutate(tseq = ss3[gt0], qseq=ss3[gt1])
    }
    aln2 = aln2 %>% mutate(qseq = map2_chr(tseq, qseq, mark_mm))
    aln2 %>% bind_rows(aln_na) %>% arrange(gt1)
    #}}}
}
var_type <- function(aln) {
    #{{{
    aln = aln %>% replace_na(list(mm=0, ins=0, del=1))
    mm = sum(aln$mm) > 0
    indel = sum(aln$ins) + sum(aln$del) > 0
    ifelse(mm, ifelse(indel, 'snp+indel', 'snp'), ifelse(indel, 'indel', 'none'))
    #}}}
}
grab_mtf_gene <- function(mid, gid, lfc, kmer, loci) {
    #{{{ #get motif and gene info (motif + gene ID + alias + annotation)
    bid = stringr::str_split(mid,'_')[[1]][1]
    x = kmer %>% filter(bid==!!bid) %>% pluck('mtfs', 1) %>%
        filter(mid==!!mid)
    m.pwm = x$mtf[[1]]; m.name = ifelse(x$known, x$fname, 'unknown motif')
    #
    g.name = loci %>% filter(gid==!!gid) %>% pluck('symbol2', 1)
    g.note = loci %>% filter(gid==!!gid) %>% pluck('note', 1)
    # acc
    hy = lfc %>% rename(status=clu)
    acc.raw = sum(hy$status==hy$B_raw, na.rm=T) / sum(!is.na(hy$B_raw))
    acc.umr = sum(hy$status==hy$B_umr, na.rm=T) / sum(!is.na(hy$B_umr))
    #
    #tit = glue("{m.name}; {g.name}; model accuracy = {number(acc,accuracy=.01)}")
    tibble(mname=m.name, pwm=list(m.pwm), gname=g.name, gnote=g.note,
           acc.raw=acc.raw, acc.umr=acc.umr)
    #}}}
}
ta3 = ta2 %>% mutate(j=row_number()) %>%
    #dplyr::slice(6) %>%
    mutate(aln = pmap(list(j, mtf, seqs, gts), check_mtf)) %>%
    select(-epi,-j) %>% mutate(vt = map_chr(aln, var_type)) %>%
    mutate(x = pmap(list(mid,gid,lfc), grab_mtf_gene, kmer=r06$tk, loci=loci)) %>%
    unnest(x)

fi = glue("{dirw}/51.assoc.rds")
la = readRDS(fi)
la1 = la %>% select(gid,drc,nsig)
ta4 = ta3 %>% left_join(la1, by=c('gid','drc'))

fo = glue("{dirw}/63.mtf.stat.rds")
saveRDS(ta4, fo)
fo = glue("{dirw}/63.mtf.stat.tsv")
to = ta4 %>% select(gid,gname,gnote,drc,mid,mname,vt,n00,n01,n10,n11,
                    acc.mtf=acc,acc.raw,acc.umr,plink.nsig=nsig)
write_tsv(to, fo, na='')
#}}}

#{{{ st6
fi = glue("{dirw}/63.mtf.stat.tsv")
ti = read_tsv(fi) %>% replace_na(list(gname='', gnote='')) %>%
    mutate(gnote=str_sub(gnote, 0, 30)) %>%
    select(GeneID=gid, Alias=gname, Description=gnote,
        Direction=drc, MotifID=mid, Motif=mname)

x = ti %>%
    kbl(format='latex', escape=T, longtable=T, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirf, 'st6.rds')
saveRDS(x, file=fo)
#}}}

#{{{ prepare for synteny plot - 65.rds
fi = glue("{dirw}/63.mtf.stat.rds")
ti = readRDS(fi)
tk = ti %>% filter(n00+n01 >= 3, n10+n11>=3, acc >.9) %>% arrange(gid) %>%
    print(n=40, width=Inf) %>%
    mutate(x = map(seqs, msa_tree)) %>%
    mutate(msa = map(x, 'msa'), tree = map(x, 'tree'), ty = map(x, 'ty')) %>%
    select(-x) %>%
    mutate(umr = map(gid, get_umr_loc))

fo = glue("{dirw}/65.plot.data.rds")
saveRDS(tk, fo)
#}}}

#{{{ read & functions
require(ape)
require(ggtree)
require(tidytree)
fx = glue("{dirg}2/syntelog/xref.maize.v4.rds")
tg = readRDS(fx)
#
fi = glue("{dirw}/65.plot.data.rds")
tk = readRDS(fi)

make_syn  <- function(name1, name2, seqs) {
    #{{{ extract synteny blocks and mismatch positions from alignment
    require(Biostrings)
    #x = unmasked(msa)
    #s1 = RemoveGaps(x[name1])[[1]]
    #s2 = RemoveGaps(x[name2])[[1]]
    pw = pairwiseAlignment(seqs[[name1]], seqs[[name2]], type='local')
    toff = start(pattern(pw)); qoff = start(subject(pw))
    y1 = as.character(pattern(pw))
    y2 = as.character(subject(pw))
    h0 = as.character(compareStrings(y2,y1))
    h = str_replace_all(h0, "[ATCGN?]", '*')
    h2=unlist(strsplit(h, split = ""))
    j = tibble(v = rle(h2)$values, aSize = rle(h2)$lengths) %>%
        mutate(aEnd=cumsum(aSize), aBeg = aEnd-aSize) %>%
        mutate(tSize = ifelse(v=='+', 0, aSize)) %>%
        mutate(qSize = ifelse(v=='-', 0, aSize)) %>%
        mutate(tEnd=cumsum(tSize), tBeg=tEnd-tSize) %>%
        mutate(qEnd=cumsum(qSize), qBeg=qEnd-qSize) %>%
        select(v, aBeg,aEnd,aSize,tBeg,tEnd,tSize,qBeg,qEnd,qSize) %>%
        mutate(tBeg = toff+tBeg-1, tEnd=tBeg+tSize) %>%
        mutate(qBeg = qoff+qBeg-1, qEnd=qBeg+qSize)
    get_tpos <- function(v, aSize, tBeg, tEnd)
        ifelse(v=='+', list(rep(tBeg, aSize)), list((tBeg+1):tEnd))
    get_qpos <- function(v, aSize, qBeg, qEnd)
        ifelse(v=='-', list(rep(qBeg, aSize)), list((qBeg+1):qEnd))
    j2 = j %>%
        mutate(tpos = pmap(list(v, aSize, tBeg, tEnd), get_tpos)) %>%
        mutate(qpos = pmap(list(v, aSize, qBeg, qEnd), get_qpos))
    tposs = unlist(j2$tpos)
    qposs = unlist(j2$qpos)
    mm = str_locate_all(h0, "\\?")[[1]] %>% as_tibble() %>%
        mutate(qPos = qposs[start], tPos = tposs[start]) %>%
        select(aPos=start,tPos,qPos)
    list(aln = j, mm = mm)
    #}}}
}
plot_syn <- function(ty, seqs) {
    #{{{
    #{{{ prepare syn tible
    labs = ty %>% arrange(desc(y)) %>% pull(lab)
    to = tibble(tgt=labs[1:nrow(ty)-1], qry=labs[2:nrow(ty)]) %>%
        mutate(x = map2(tgt,qry, make_syn, seqs=seqs)) %>%
        mutate(syn = map(x, 'aln'), mm = map(x, 'mm')) %>%
        select(tgt, qry, syn, mm)
    #
    tp = to %>% select(tgt,qry,syn) %>% unnest(syn) %>%
        filter(v == '*') %>%
        inner_join(ty %>% select(lab,y), by=c('tgt'='lab')) %>% rename(y1=y) %>%
        inner_join(ty %>% select(lab,y), by=c('qry'='lab')) %>% rename(y2=y) %>%
        mutate(tBeg=tBeg + 1, qBeg = qBeg + 1) %>%
        mutate(i = 1:n())
    tp1 = tp %>% select(i, y=y1, tBeg, tEnd) %>%
        gather(type, pos, -i, -y) %>% mutate(y=y-.1)
    tp2 = tp %>% select(i, y=y2, qBeg, qEnd) %>%
        gather(type, pos, -i, -y) %>% mutate(y=y+.1)
    coordmap = c("tBeg"=1,'tEnd'=2,'qEnd'=3,'qBeg'=4)
    tp = tp1 %>% rbind(tp2) %>% mutate(i2 = coordmap[type]) %>% arrange(i, i2)
    #}}}
    #
    p_syn = ggplot(tp) +
        geom_polygon(aes(x=pos,y=y,group=i), fill='royalblue', alpha=.2,
                     size=0,color=NA) +
        coord_cartesian(xlim = c(0,4000)) +
        scale_x_continuous(breaks=c(0,2000,4000), labels=c('-2k','TSS','+2k'),
                           expand=expansion(mult=c(.01,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.01))) +
        otheme(xtext=T,xtick=T,ytext=T,ytick=T, panel.border=F,
            margin=c(.2,.5,.2,0))
    p_syn
    #}}}
}
add_gene_track <- function(gid, tg, ty, p, ht.exon=.05, ht.cds=.1) {
    #{{{
    tg1 = tg %>% filter(gid==!!gid) %>% select(gt, gene) %>% unnest(gene) %>%
        mutate(gt = str_replace(gt, '^Zmays_', '')) %>%
        group_by(gt, srd) %>%
        mutate(tss = ifelse(srd=='-', max(end), min(start))) %>% ungroup() %>%
        mutate(rb = ifelse(srd=='-', tss-end+1, start-tss+1)) %>%
        mutate(re = ifelse(srd=='-', tss-start+1, end-tss+1)) %>%
        select(-start, -end) %>% rename(beg=rb, end=re) %>%
        inner_join(ty, by='gt')
    tg1 = tg1 %>% mutate(beg = beg+2000, end = end+2000)
    tgr = tg1 %>% group_by(gt, y) %>%
        summarise(beg=min(beg), end=max(end)) %>% ungroup()
    tge = tg1 %>% filter(etype=='exon')
    tgc = tg1 %>% filter(etype=='CDS')
    col.intron='grey'; col.exon='grey'; col.cds='royalblue'
    arw = arrow(length=unit(.1,'cm'), angle=30, ends='last',type="open")
    p +
        geom_segment(data=tgr,aes(x=beg,xend=end,y=y,yend=y),col=col.intron,size=.5) +
        geom_rect(data=tge,aes(xmin=beg,xmax=end,ymin=y-ht.exon,ymax=y+ht.exon),fill=col.exon,color=NA,alpha=1) +
        geom_rect(data=tgc,aes(xmin=beg,xmax=end,ymin=y-ht.cds,ymax=y+ht.cds),fill=col.cds,color=NA,alpha=1) +
        geom_segment(data=tgr,aes(x=2000,xend=2100,y=y+ht.cds*1.1,yend=y+ht.cds*1.1),
                     color='black', size=.2, arrow=arw) +
        geom_rect(xmin=1999,xmax=2001,ymin=-Inf,ymax=Inf, fill='yellow', alpha=.2)
    #}}}
}
mtf_title <- function(gid, gname,gnote, pwm,mname, acc,nr=1, font.size=3) {
    #{{{ #plot title (motif + gene ID + alias + annotation)
    glab = ifelse(is.na(gname), gid, glue("{gid} ({gname})"))
    glab = ifelse(is.na(gnote), glue("{glab} (unknown)"), glue("{glab}: {gnote}"))
    #
    tit = glue("{mname}; {glab} model accuracy = {number(acc,accuracy=.01)}")
    if (nr == 2) {
        tit = glue("{mname}; {glab}\n model accuracy = {number(acc,accuracy=.01)}")
    }
    p01 = view_motifs(pwm, method="PCC",min.mean.ic=0,min.overlap=6) +
        otheme(panel.border=F, margin=c(.2,.2,0,2))
    p02 = ggplot() +
        annotate('text', x=0, y=0, label=tit, size=font.size, hjust=.5, vjust=.5) +
        #scale_x_continuous(limits=c(0,3),expand=expansion(mult=c(.01,.01))) +
        #scale_y_continuous(limits=c(0,1), expand=expansion(add=c(.2,0))) +
        otheme(panel.border=F, margin=c(.2,.2,0,.2))
    p0 = ggarrange(p01, p02, nrow=1, ncol=2, widths=c(1,3))
    p0
    #}}}
}
plot_msa <- function(aln, ty) {
    #{{{
    tseq = aln %>% filter(!is.na(tseq)) %>% pluck('tseq', 1)
    if( is.null(tseq) ) {
        ggplot() + geom_blank()
    } else {
    gap = str_c(rep('-', nchar(tseq)), collapse='')
    aln = aln %>% replace_na(list(tseq=gap, qseq=gap))
    ts0 = aln %>% dplyr::slice(1) %>% select(gt=gt0, seq=tseq)
    ts1 = aln %>% select(gt=gt1, seq=qseq)
    str_explode <- function(s) tibble(nt=unlist(str_split(s,''))) %>% mutate(i=1:n())
    ts = ts0 %>% bind_rows(ts1) %>%
        mutate(x = map(seq, str_explode)) %>% select(gt, x) %>% unnest(x) %>%
        mutate(high = nt %in% letters[1:26]) %>%
        inner_join(ty, by='gt')
    cols2 = c('white', pal_simpsons()(1))
    ggplot(ts) +
        geom_tile(aes(x=i, y=y, fill=high)) +
        geom_text(aes(x=i, y=y, label=nt), size=2) +
        scale_x_continuous(expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(0,.02))) +
        scale_fill_manual(values=cols2) +
        otheme(legend.pos='none', panel.border=F, margin=c(.2,.05,.2,0))
    }
    #}}}
}
plot_lfc <- function(gid, ty, lfc) {
#{{{ log2fc
    lfc = lfc %>% right_join(ty, by='gt')
    ymin0 = min(lfc$log2fc, na.rm=T); ymax0 = max(lfc$log2fc, na.rm=T)
    lfc.mean = (ymin0 + ymax0) / 2
    lfc = lfc %>% mutate(txt.col = ifelse(is.na(log2fc), ymax0, log2fc))
    t_na = lfc %>% filter(is.na(log2fc)) %>% mutate(pos=(ymin0+ymax0)/2)
    ymin = ymin0; ymax = ymax0
    ymin = ifelse(ymin < 0, ceiling(ymin), floor(ymin))
    ymax = ifelse(ymax < 0, ceiling(ymax), floor(ymax))
    p.a = ggplot(lfc) +
        #geom_hline(yintercept=0, color='gray', alpha=.4) +
        geom_tile(aes(x=0, y=y, fill=log2fc), width=.4, height=.4) +
        #geom_text(aes(x=y, y=ymax0, label=B_umr), hjust=0, size=2) +
        geom_text(aes(x=0, y=y, label=txt, col=txt.col>lfc.mean), hjust=.5,size=2) +
        scale_x_continuous(name='log2fc', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.02))) +
        scale_color_manual(values=c('white','black')) +
        scale_fill_gradientn(colors=rev(cols100v), na.value="white") +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.b = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=tag_de), hjust=.5,size=2.5) +
        scale_x_continuous(name='DE', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.c = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=tag_mtf), hjust=.5,size=2.5) +
        scale_x_continuous(name='mtf', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.d = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=B_raw), size=2) +
        scale_x_continuous(name='model', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        scale_fill_manual(values=c('white','black')) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,.2,0,0)) +
        theme(axis.title.x = element_text(size=7))
    ggarrange(p.a, p.b, p.c, p.d, nrow=1, ncol=4, widths=c(3,1.5,1.5,3.5), align='v')
#}}}
}
combo_plot <- function(gid,gname,gnote,pwm,mname,acc, ty, seqs, mtf, tree,aln, lfc, fp, tg) {
    #{{{
    p0 = mtf_title(gid, gname,gnote, pwm,mname, acc, font.size=3)
    # tree plot
    ty2 = ty %>% select(taxa=lab, gt)
    p_tree = ggtree(tree) %<+% ty2 +
        #geom_tiplab(aes(label=gt), hjust=0, align=T) +
        scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
        scale_y_continuous(expand=expansion(mult=c(.01,.05))) +
        otheme(panel.border=F, margin=c(.2,.2,.2,0))
    #
    p1 = plot_syn(ty, seqs)
    pg = add_gene_track(gid,tg,ty,p1)
    #umr2 = umr %>% inner_join(ty, by='gt')
    #pg = pg + geom_segment(data=umr2, aes(x=start,xend=end,y=y-.1,yend=y-.1), color='green', size=.5, alpha=1)
    mtf2 = mtf %>% inner_join(ty, by='gt')
    pg = pg + geom_point(data=mtf2,aes(x=pos, y=y-.1),col='red',size=1)
    #
    p_lfc = plot_lfc(gid, ty, lfc)
    #
    p_msa = plot_msa(aln, ty)
    pb = ggarrange(p_tree, pg, p_msa, p_lfc, nrow=1, ncol=4, widths=c(1,5,.5,1))
    #pb = ggarrange(p_tree, pg, p_lfc, nrow=1, ncol=3, widths=c(1,5,1))
    ggarrange(p0, pb, nrow=2, heights=c(1,20)) %>%
        ggexport(filename=fp, width=8, height=7)
    #}}}
}
combo_plot1 <- function(gid,gname,gnote,pwm,mname,acc, ty,seqs,mtf,tree,aln,lfc,tg) {
    #{{{
    p0 = mtf_title(gid, gname,gnote, pwm,mname, acc, nr=2, font.size=3)
    # tree plot
    ty2 = ty %>% select(taxa=lab, gt)
    p_tree = ggtree(tree) %<+% ty2 +
        #geom_tiplab(aes(label=gt), hjust=0, align=T) +
        scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
        scale_y_continuous(expand=expansion(mult=c(.01,.05))) +
        otheme(panel.border=F, margin=c(.2,.2,.2,0))
    #
    p1 = plot_syn(ty, seqs)
    pg = add_gene_track(gid,tg,ty,p1)
    #umr2 = umr %>% inner_join(ty, by='gt')
    #pg = pg + geom_segment(data=umr2, aes(x=start,xend=end,y=y-.1,yend=y-.1), color='green', size=.5, alpha=1)
    mtf2 = mtf %>% inner_join(ty, by='gt')
    pg = pg + geom_point(data=mtf2,aes(x=pos, y=y-.1),col='red',size=1)
    #
    p_lfc = plot_lfc(gid, ty, lfc)
    #
    p_msa = plot_msa(aln, ty)
    #pb = ggarrange(p_tree, pg, p_msa, p_lfc, nrow=1, ncol=4, widths=c(1,5,.5,1))
    pb = ggarrange(p_tree, pg, p_lfc, nrow=1, ncol=3, widths=c(.7,5,1))
    ggarrange(p0, pb, nrow=2, heights=c(1,15))
    #}}}
}
combo_plot2 <- function(gid,gname,gnote,pwm,mname,acc, ty,seqs,mtf,tree,aln,lfc,tg) {
    #{{{
    p0 = mtf_title(gid, gname,gnote, pwm,mname, acc, nr=1, font.size=3)
    # tree plot
    ty2 = ty %>% select(taxa=lab, gt)
    p_tree = ggtree(tree) %<+% ty2 +
        #geom_tiplab(aes(label=gt), hjust=0, align=T) +
        scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
        scale_y_continuous(expand=expansion(mult=c(.01,.05))) +
        otheme(panel.border=F, margin=c(.2,.2,.2,0))
    #
    p1 = plot_syn(ty, seqs)
    pg = add_gene_track(gid,tg,ty,p1)
    #umr2 = umr %>% inner_join(ty, by='gt')
    #pg = pg + geom_segment(data=umr2, aes(x=start,xend=end,y=y-.1,yend=y-.1), color='green', size=.5, alpha=1)
    mtf2 = mtf %>% inner_join(ty, by='gt')
    pg = pg + geom_point(data=mtf2,aes(x=pos, y=y-.1),col='red',size=1)
    #
    p_lfc = plot_lfc(gid, ty, lfc)
    #
    p_msa = plot_msa(aln, ty)
    pb = ggarrange(p_tree, pg, p_msa, p_lfc, nrow=1, ncol=4, widths=c(1,5,.5,1))
    #pb = ggarrange(p_tree, pg, p_lfc, nrow=1, ncol=3, widths=c(.7,5,1))
    ggarrange(p0, pb, nrow=2, heights=c(1,20))
    #}}}
}
#}}}

i=10
gid=tk$gid[[i]]; gname=tk$gname[[i]]; gnote=tk$gnote[[i]]; pwm=tk$pwm[[i]];
mid=tk$mid[[i]]; mname=tk$mname[[i]]; acc=tk$acc.raw[[i]]; ty=tk$ty[[i]];
aln=tk$aln[[i]]; lfc=tk$lfc[[i]];
seqs=tk$seqs[[i]]; tree=tk$tree[[i]]; mtf=tk$mtf[[i]]; umr=tk$umr[[i]]
fp = glue("{dirw}/z{i}-{gid}.pdf")
combo_plot(gid,gname,gnote,pwm,mname,acc, ty,seqs,mtf,tree,aln,lfc,fp, tg=tg)

tk %>% mutate(i = row_number()) %>%# dplyr::slice(11:n()) %>%
    mutate(fp = glue("{dirw}/z{i}-{gid}.pdf")) %>% print(width=Inf) %>%
    mutate(j = pmap(list(gid,gname,gnote,pwm,mname,acc.raw,ty,seqs,mtf,tree,aln,lfc,fp), combo_plot, tg=tg))

tp = tk %>% mutate(i = row_number()) %>% filter(i %in% c(7,17)) %>%
    mutate(fp = glue("{dirw}/z{i}-{gid}.pdf")) %>% print(width=Inf) %>%
    mutate(p = pmap(list(gid,gname,gnote,pwm,mname,acc.raw,ty,seqs,mtf,tree,aln,lfc), combo_plot1, tg=tg))
#tp %>% mutate(j = pmap(list(p, filename=fp), ggexport, width=8, height=7))
saveRDS(tp$p[[2]], glue("{diro}/69.case.c.rds"))
saveRDS(tp$p[[1]], glue("{diro}/69.case.d.rds"))

tp = tk %>% mutate(i = row_number()) %>% filter(i %in% c(6,18)) %>%
    mutate(fp = glue("{dirw}/z{i}-{gid}.pdf")) %>% print(width=Inf) %>%
    mutate(p = pmap(list(gid,gname,gnote,pwm,mname,acc.raw,ty,seqs,mtf,tree,aln,lfc), combo_plot2, tg=tg))
p1 = tp$p[[2]]; p2 = tp$p[[1]]
p12 = ggarrange(p1, p2, nrow=2, ncol=1, heights=c(1,1), labels=LETTERS[1:2])
fo = glue("{dirw}/69.sup.pdf")
ggexport(p12, filename=fo, width=8, height=8)
fo = glue("{dirf}/sf09.pdf")
ggexport(p12, filename=fo, width=8, height=8)
#}}}

#{{{ f7a-b
#{{{ functions
read_mtf <- function(f_bed) {
    #{{{
    read_tsv(f_bed, col_names=c("sid",'beg','end')) %>%
        separate('sid',c('mid','sid'), sep='%') %>%
        mutate(pos = (beg+end)/2) %>%
        select(gt=sid, beg, end, pos)
    #}}}
}
mtf_aln  <- function(n1, n2, tb, te, seqs) {
    #{{{
    cat(n1, n2, '\n')
    pw = pairwiseAlignment(seqs[[n1]],seqs[[n2]], type='local-global')
    pw = pairwiseAlignment(seqs[[n1]],seqs[[n2]], type='global-local')
    toff = start(pattern(pw)); qoff = start(subject(pw))
    y1 = as.character(pattern(pw))
    y2 = as.character(subject(pw))
    h0 = as.character(compareStrings(y2,y1))
    h = str_replace_all(h0, "[ATCGN?]", '*')
    h2=unlist(strsplit(h, split = ""))
    j = tibble(v = rle(h2)$values, aSize = rle(h2)$lengths) %>%
        mutate(aEnd=cumsum(aSize), aBeg = aEnd-aSize) %>%
        mutate(tSize = ifelse(v=='+', 0, aSize)) %>%
        mutate(qSize = ifelse(v=='-', 0, aSize)) %>%
        mutate(tEnd=cumsum(tSize), tBeg=tEnd-tSize) %>%
        mutate(qEnd=cumsum(qSize), qBeg=qEnd-qSize) %>%
        select(v, aBeg,aEnd,aSize,tBeg,tEnd,tSize,qBeg,qEnd,qSize) %>%
        mutate(tBeg = toff+tBeg-1, tEnd=tBeg+tSize) %>%
        mutate(qBeg = qoff+qBeg-1, qEnd=qBeg+qSize)
    get_tpos <- function(v, aSize, tBeg, tEnd)
        ifelse(v=='+', list(rep(tBeg, aSize)), list((tBeg+1):tEnd))
    get_qpos <- function(v, aSize, qBeg, qEnd)
        ifelse(v=='-', list(rep(qBeg, aSize)), list((qBeg+1):qEnd))
    j2 = j %>%
        mutate(tpos = pmap(list(v, aSize, tBeg, tEnd), get_tpos)) %>%
        mutate(qpos = pmap(list(v, aSize, qBeg, qEnd), get_qpos))
    tposs = unlist(j2$tpos)
    qposs = unlist(j2$qpos)
    mm = str_locate_all(h0, "\\?")[[1]] %>% as_tibble() %>%
        mutate(qPos = qposs[start], tPos = tposs[start]) %>%
        select(aPos=start,tPos,qPos)
    ab = which(tposs == tb)[1] 
    ae = which(tposs == te)[1] 
    tseq = str_sub(y1, ab, ae); qseq=str_sub(y2,ab,ae)
    mm1 = mm %>% filter(aPos >= ab, aPos <= ae)
    if(nrow(mm1) > 0) {
        for (i in 1:nrow(mm1)) {
            mp = mm1$aPos[i] - ab + 1
            str_sub(qseq,mp,mp) <- str_to_lower(str_sub(qseq,mp,mp))
        }
    }
    qb=qposs[ab]; qe=qposs[ae]
    ins=str_count(tseq,'-')
    del=str_count(qseq,'-')
    tibble(tseq=tseq,qseq=qseq,qb=qb,qe=qe, mm=nrow(mm1),ins=ins,del=del)
    #}}}
}
make_msa <- function(seqs) {
    #{{{
    pre = glue("tmp.msa.{sample(10000, 1)}")
    writeXStringSet(seqs, glue('{pre}.fa'))
    system(glue("muscle -in {pre}.fa -out {pre}.2.fa"))
    msa = readDNAMultipleAlignment(filepath=glue("{pre}.2.fa"), format='fasta')
    msa = DNAStringSet(msa)
    #
    system(glue("rm {pre}.*"))
    msa
    #}}}
}
mark_mm  <- function(tseq, qseq) {
    #{{{ mark mismatch and indel in query seqs
    h0 = as.character(compareStrings(qseq,tseq))
    h2 = unlist(strsplit(h0, split = ""))
    th2 = tibble(p=1:length(h2), op = h2) %>%
        filter(op %in% c("?", "+"))
    if (nrow(th2) > 0) {
        for (p in th2$p) {
            str_sub(qseq,p,p) <- str_to_lower(str_sub(qseq,p,p))
        }
    }
    qseq
    #}}}
}
check_mtf <- function(mtf, seqs, gts) {
    #{{{
    #cat(j, '\n')
    mtf0 = mtf %>% mutate(dst = abs(pos-2000)) %>% arrange(dst)
    i=1; gt0 = mtf0$gt[[i]]; b0 = mtf0$beg[[i]]; e0 = mtf0$end[[i]]
    seq0 = seqs[[gt0]][b0:e0]
    aln = tibble(gt0=!!gt0, gt1=gts[! gts %in% gt0]) %>%
        filter(gt1 %in% names(seqs)) %>%
        mutate(tb=b0, te=e0) %>%
        mutate(x = pmap(list(gt0,gt1,tb,te), mtf_aln, seqs=seqs)) %>%
        unnest(x)
    #
    aln_na = aln %>% filter(is.na(tseq))
    aln2 = aln %>% filter(!is.na(tseq))
    if (length(unique(aln2$tseq)) > 1) {
        ss = aln2 %>% select(gt=gt1, seq=qseq) %>%
            bind_rows(aln2 %>% dplyr::slice(1) %>% select(gt=gt0, seq=tseq)) %>%
            mutate(seq=str_replace(seq,'-',''))
        ss0 = ss$seq; names(ss0) = ss$gt
        ss1 = DNAStringSet(ss0)
        ss2 = make_msa(ss1)
        ss3 = as.character(ss2)
        aln2 = aln2 %>% mutate(tseq = ss3[gt0], qseq=ss3[gt1])
    }
    aln2 = aln2 %>% mutate(qseq = map2_chr(tseq, qseq, mark_mm))
    aln2 %>% bind_rows(aln_na) %>% arrange(gt1)
    #}}}
}
var_type <- function(aln) {
    #{{{
    aln = aln %>% replace_na(list(mm=0, ins=0, del=1))
    mm = sum(aln$mm) > 0
    indel = sum(aln$ins) + sum(aln$del) > 0
    ifelse(mm, ifelse(indel, 'snp+indel', 'snp'), ifelse(indel, 'indel', 'none'))
    #}}}
}
anno_gene <- function(gid, loci) {
    #{{{ #get gene info (gene ID + alias + annotation)
    gname = loci %>% filter(gid==!!gid) %>% pluck('symbol2', 1)
    gnote = loci %>% filter(gid==!!gid) %>% pluck('note', 1)
    glab = ifelse(is.na(gname), gid, glue("{gid} ({gname})"))
    glab = ifelse(is.na(gnote), glue("{glab} (unknown)"), glue("{glab}: {gnote}"))
    glab
    #}}}
}
grab_mtf_gene <- function(mid, gid, lfc, kmer, loci) {
    #{{{ #get motif and gene info (motif + gene ID + alias + annotation)
    bid = stringr::str_split(mid,'_')[[1]][1]
    x = kmer %>% filter(bid==!!bid) %>% pluck('mtfs', 1) %>%
        filter(mid==!!mid)
    m.pwm = x$mtf[[1]]; m.name = ifelse(x$known, x$fname, 'unknown motif')
    #
    g.name = loci %>% filter(gid==!!gid) %>% pluck('symbol2', 1)
    g.note = loci %>% filter(gid==!!gid) %>% pluck('note', 1)
    # acc
    hy = lfc %>% rename(status=clu)
    acc.raw = sum(hy$status==hy$B_raw, na.rm=T) / sum(!is.na(hy$B_raw))
    acc.umr = sum(hy$status==hy$B_umr, na.rm=T) / sum(!is.na(hy$B_umr))
    #
    #tit = glue("{m.name}; {g.name}; model accuracy = {number(acc,accuracy=.01)}")
    tibble(mname=m.name, pwm=list(m.pwm), gname=g.name, gnote=g.note,
           acc.raw=acc.raw, acc.umr=acc.umr)
    #}}}
}
plot_title <- function(glab, pwm,mname, acc=NA, nr=1, font.size=3) {
    #{{{ #plot title (motif + gene ID + alias + annotation)
    tit = glue("{mname}; {glab}")
    sep = ifelse(nr == 2, "\n", " ")
    if (!is.na(acc)) {
        tit = glue("{tit}{sep}model accuracy = {number(acc,accuracy=.01)}")
    }
    p01 = view_motifs(pwm, method="PCC",min.mean.ic=0,min.overlap=6) +
        otheme(panel.border=F, margin=c(.2,.2,0,2))
    p02 = ggplot() +
        annotate('text', x=0, y=0, label=tit, size=font.size, hjust=.5, vjust=.5) +
        #scale_x_continuous(limits=c(0,3),expand=expansion(mult=c(.01,.01))) +
        #scale_y_continuous(limits=c(0,1), expand=expansion(add=c(.2,0))) +
        otheme(panel.border=F, margin=c(.2,.2,0,.2))
    p0 = ggarrange(p01, p02, nrow=1, ncol=2, widths=c(1,3))
    p0
    #}}}
}
plot_mtf_msa <- function(aln, ty) {
    #{{{
    if( is.null(aln) ) {
        ggplot() + geom_blank()
    }
    size = aln$te[[1]] - aln$tb[[1]] + 1
    #tseq = aln %>% filter(!is.na(tseq)) %>% pluck('tseq', 1)
    gap = str_c(rep('-', size), collapse='')
    aln = aln %>% replace_na(list(tseq=gap, qseq=gap))
    ts0 = aln %>% dplyr::slice(1) %>% select(gt=gt0, seq=tseq)
    ts1 = aln %>% select(gt=gt1, seq=qseq)
    str_explode <- function(s) tibble(nt=unlist(str_split(s,''))) %>% mutate(i=1:n())
    ts = ts0 %>% bind_rows(ts1) %>%
        mutate(x = map(seq, str_explode)) %>% select(gt, x) %>% unnest(x) %>%
        mutate(high = nt %in% letters[1:26]) %>%
        inner_join(ty, by='gt')
    cols2 = c('white', pal_simpsons()(1))
    ggplot(ts) +
        geom_tile(aes(x=i, y=y, fill=high)) +
        geom_text(aes(x=i, y=y, label=nt), size=2) +
        scale_x_continuous(expand=expansion(mult=c(.01,.01))) +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(0,.02))) +
        scale_fill_manual(values=cols2) +
        otheme(legend.pos='none', panel.border=F, margin=c(.2,.05,.2,0))
    #}}}
}
plot_lfc0 <- function(gid, ty, lfc) {
#{{{ log2fc
    lfc = lfc %>% right_join(ty, by='gt')
    ymin0 = min(lfc$log2fc, na.rm=T); ymax0 = max(lfc$log2fc, na.rm=T)
    lfc.mean = (ymin0 + ymax0) / 2
    lfc = lfc %>% mutate(txt.col = ifelse(is.na(log2fc), ymax0, log2fc))
    t_na = lfc %>% filter(is.na(log2fc)) %>% mutate(pos=(ymin0+ymax0)/2)
    ymin = ymin0; ymax = ymax0
    ymin = ifelse(ymin < 0, ceiling(ymin), floor(ymin))
    ymax = ifelse(ymax < 0, ceiling(ymax), floor(ymax))
    p.a = ggplot(lfc) +
        #geom_hline(yintercept=0, color='gray', alpha=.4) +
        geom_tile(aes(x=0, y=y, fill=log2fc), width=.4, height=.4) +
        #geom_text(aes(x=y, y=ymax0, label=B_umr), hjust=0, size=2) +
        geom_text(aes(x=0, y=y, label=txt, col=txt.col>lfc.mean), hjust=.5,size=2) +
        scale_x_continuous(name='log2fc', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.02))) +
        scale_color_manual(values=c('white','black')) +
        scale_fill_gradientn(colors=rev(cols100v), na.value="white") +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.b = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=tag_de), hjust=.5,size=2.5) +
        scale_x_continuous(name='DE', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.c = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=tag_mtf), hjust=.5,size=2.5) +
        scale_x_continuous(name='mtf', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.d = ggplot(lfc) +
        geom_text(aes(x=0, y=y, label=B_raw), size=2) +
        scale_x_continuous(name='model', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.02,.03))) +
        scale_fill_manual(values=c('white','black')) +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,.2,0,0)) +
        theme(axis.title.x = element_text(size=7))
    ggarrange(p.a, p.b, p.c, p.d, nrow=1, ncol=4, widths=c(3,1.5,1.5,3.5), align='v')
#}}}
}
plot_lfc <- function(gid, ty, lfc) {
#{{{ log2fc
    lfc = lfc %>% right_join(ty, by='gt') %>%
        mutate(txt = number(log2fc, accuracy=.01))
    ymin0 = min(lfc$log2fc, na.rm=T); ymax0 = max(lfc$log2fc, na.rm=T)
    lfc.mean = (ymin0 + ymax0) / 2
    lfc = lfc %>% mutate(txt.col = ifelse(is.na(log2fc), ymax0, log2fc))
    t_na = lfc %>% filter(is.na(log2fc)) %>% mutate(pos=(ymin0+ymax0)/2)
    ymin = ymin0; ymax = ymax0
    ymin = ifelse(ymin < 0, ceiling(ymin), floor(ymin))
    ymax = ifelse(ymax < 0, ceiling(ymax), floor(ymax))
    p.a = ggplot(lfc) +
        #geom_hline(yintercept=0, color='gray', alpha=.4) +
        geom_tile(aes(x=0, y=y, fill=log2fc), width=.4, height=.4) +
        #geom_text(aes(x=y, y=ymax0, label=B_umr), hjust=0, size=2) +
        geom_text(aes(x=0, y=y, label=txt, col=txt.col>lfc.mean), hjust=.5,size=2) +
        scale_x_continuous(name='log2fc', expand=expansion(mult=c(.02,.02)), position='top') +
        scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.02))) +
        scale_color_manual(values=c('white','black')) +
        scale_fill_gradientn(colors=rev(cols100v), na.value="white") +
        otheme(xtitle=T, legend.pos='none', panel.border=F, margin=c(0,0,0,0)) +
        theme(axis.title.x = element_text(size=7))
    p.a
#}}}
}
#}}}
require(universalmotif)
require(Biostrings)
fm = glue("{dirw}/70_cases/hsf.meme")
pwm = read_meme(fm)

# f7a
gid = 'Zm00001d048532'
fi1 = glue("{dirw}/70_cases/hsf1.rds")
r1 = readRDS(fi1); r=r1
p_tree=r$p.tree; p_aln=r$p.aln; seqs=r$seqs; ty=r$ty; msa=r$msa; tree=r$tree
p_tree=p_tree+o_margin(1.2,.1,.5,.1)
p_title = plot_title(anno_gene(gid, loci), pwm, 'HSF')
mtf = read_mtf(glue("{dirw}/70_cases/hsf1.bed"))
aln = check_mtf(mtf, seqs, ty$gt)
p_mtf = plot_mtf_msa(aln, ty) + o_margin(.5,.1,0,.1)
lfc = read_xlsx(glue("{dirw}/70_cases/lfc.hsf.xlsx"), sheet='hsf1')
p_lfc = plot_lfc(gid, ty, lfc) + o_margin(.2,.4,.2,0)
a = ggarrange(p_title,
              ggarrange(p_tree, p_aln, p_mtf, p_lfc,
                        nrow=1, ncol=4, widths=c(.5,5,.5,.3)),
              nrow=2, heights=c(1,5))

# f7b
gid = 'Zm00001d039670'
fi2 = glue("{dirw}/70_cases/hsf2.rds")
r2 = readRDS(fi2); r=r2
p_tree=r$p.tree; p_aln=r$p.aln; seqs=r$seqs; ty=r$ty; msa=r$msa; tree=r$tree
p_tree=p_tree+o_margin(1.2,.1,.5,.1)
p_title = plot_title(anno_gene(gid, loci), pwm, 'HSF')
mtf = read_mtf(glue("{dirw}/70_cases/hsf2.bed"))
aln = check_mtf(mtf, seqs, ty$gt)
p_mtf = plot_mtf_msa(aln, ty) + o_margin(.5,.1,0,.1)
lfc = read_xlsx(glue("{dirw}/70_cases/lfc.hsf.xlsx"), sheet='hsf2')
p_lfc = plot_lfc(gid, ty, lfc) + o_margin(.2,.4,.2,0)
b = ggarrange(p_title,
              ggarrange(p_tree, p_aln, p_mtf, p_lfc,
                        nrow=1, ncol=4, widths=c(.5,5,.5,.3)),
              nrow=2, heights=c(1,5))


p3 = readRDS(glue("{dird}/45_nam/69.c.rds"))
p4 = readRDS(glue("{dird}/45_nam/69.d.rds"))
ab = ggarrange(a, b, nrow=1, ncol=2, widths=c(1,1), labels=LETTERS[1:2])
cd = ggarrange(p3, p4, nrow=1, ncol=2, widths=c(1,1), labels=LETTERS[3:4])
fo = glue("{dirf}/f7.pdf")
ggarrange(ab, cd, nrow=2, ncol=1, heights=c(2,8)) %>%
    ggexport(filename=fo, width=11, height=8)
#}}}





