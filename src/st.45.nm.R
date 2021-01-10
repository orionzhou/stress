source('functions.R')
dirw = glue('{dird}/45_nam')
setwd(dirw)
tgl = gcfg$gene %>% select(gid,chrom,start,end,srd)
tga = read_loci() %>% select(gid,symbol,note)
txdb = load_txdb('Zmays_B73', primary=T)
#{{{ functions
require(cluster)
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
    group_by(cond,time,gt,gid) %>% slice(1) %>% ungroup()
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

#{{{ read in NM matrix and run mclust
#{{{ read in
yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm
#
th1 = th %>% filter(Experiment=='NM') %>%
    select(sid=SampleID, gt=Genotype, cond=Treatment, time=Timepoint) %>%
    mutate(cond = str_to_lower(cond))
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

#{{{ HY & NM consistency - sf15a
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
    scale_x_continuous(name='mean log2fc in HY experiment', limits=c(-fcm,fcm), expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='log2fc in NM experiment', limits=c(-fcm,fcm), expand=expansion(mult=c(.02,.02))) +
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
fo = glue("{dirf}/sf15a.rds")
saveRDS(p, fo)
#}}}

#{{{ mclust
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
refine_bimod <- function(clu, t_clu) {
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
rd4 = rd3 %>% filter(nc>=2)
pb <- progress_bar$new(total=nrow(rd4))
rd4 = rd4 %>%
    mutate(x = map2(clu,t_clu, refine_bimod)) %>% unnest(x)
rd5 = rd4 %>% filter(time==25, (drc=='down' & dn) | (drc=='up' & up), nc1>=2)
#}}}

r = list(rd3=rd3, rd4=rd4, rd5=rd5)
fo = glue('{dirw}/01.mclust.rds')
saveRDS(r, fo)
#}}}

fi = glue('{dirw}/01.mclust.rds')
r = readRDS(fi)
rd3 = r$rd3; rd5 = r$rd5

#{{{ uni- to bi- model composition - f7a
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
fo = glue("{dirf}/f7a.rds")
saveRDS(pa, fo)
#}}}

#{{{ example of uni- and bi- model genes - f7b
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
fo = glue("{dirf}/f7b.rds")
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
vnt_scan <- function(t_clu,drc, chrom,start,end) {
    #{{{
    #pb$tick()
    if (drc == 'down') {
        pheno = t_clu %>% mutate(clu=ifelse(clu==-1, 1, 0))
    } else {
        pheno = t_clu %>% mutate(clu=ifelse(clu==1, 1, 0))
    }
    pheno = pheno %>% select(sid=gt,clu)
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
fet <- function(x00,x01,x10,x11) fisher.test(matrix(c(x00,x10,x01,x11), nrow=2))$p.value
vnt_scan0 <- function(t_clu, chrom,start,end) {
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

x = rd5 %>%# slice(1:2) %>% #filter(gid == !!gid) %>%
    inner_join(tgl, by='gid') %>%
    mutate(start = start - 2000, end = end + 2000)
#pb <- progress_bar$new(total=nrow(x))
x = x %>%
    mutate(data = pmap(list(t_clu1,drc, chrom,start,end), vnt_scan)) %>%
    select(drc,time,st,gid,data) %>% unnest(data)

fo = glue("{dirw}/05.test.rds")
saveRDS(x, fo)
#}}}

fi = glue("{dirw}/05.test.rds")
ti = readRDS(fi)
ti %>% count(drc, nsig > 0)

# add reg
ddeg3 = ddeg2 %>% filter(cond=='cold',tgt=='B73') %>% select(time,gid,reg)
ti2 = ti %>%
    mutate(time=as.integer(time)) %>%
    mutate(sig = ifelse(nsig > 0,'sig','insig')) %>%
    left_join(ddeg3, by=c("time",'gid')) %>%
    filter(!is.na(reg))
ti2 %>% count(reg, sig) %>% spread(sig,n) %>% mutate(prop=sig/(sig+insig))

#{{{ plot cis/trans prop - sf15b
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
fo = glue("{dirf}/sf15b.rds")
saveRDS(p, fo)

# add symbol & note
ti3 = ti2 %>% inner_join(tga, by='gid')
j = ti3

j %>% filter(nsig>0,reg %in% c('cis','cis+trans')) %>%
    select(gid,p,nsig,reg,symbol) %>% print(n=40,width=Inf)
#}}}

#{{{ final plot
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
gids = c('Zm00001d020970') # sf16b
gids = c('Zm00001d052653') # sf16a
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

fo = glue("{dirf}/sf16a.rds")
saveRDS(tp$p[[1]], fo)
fo = glue("{dirf}/sf16b.rds")
saveRDS(tp$p[[1]], fo)
#}}}



