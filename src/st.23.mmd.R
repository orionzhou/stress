source('functions.R')
require(universalmotif)
#require(dynamicTreeCut)
dirw = file.path(dird, '23_mmd')
#
fl = file.path(dirw, '../21_seq/15.rds')
tl0 = readRDS(fl)
#tl = tl0 %>% filter(pick) %>% select(lid,ng0,ng,tg,clid,ng_c,tg_c)
fln = file.path(dirw, '../21_seq/45.rds')
tln = readRDS(fln)

#{{{ process kmer location and save [long time]
fi = file.path(dirr, 'x.bed')
ti = read_tsv(fi, col_names=c('sid','start','end','kmer','idx','srd','kmer2','idx2','umr','acrL','acrE','cds','utr5','utr3','intron'))
identical(ti$kmer,ti$kmer2)
identical(ti$idx,ti$idx2)
kl = ti %>% select(-idx2,-kmer2) %>% group_by(kmer) %>% nest() %>% rename(loc=data)

fo = file.path(dirw, '02_kmer_loc.rds')
saveRDS(kl, fo)
#}}}

#{{{ read DREME kmers
rename_mid <- function(lid, mid) sprintf("%s_%03d", lid, as.numeric(str_replace(mid,'DREME-','')))
fi = sprintf("%s/dreme_kmer.tsv", dirr)
tk = read_tsv(fi) %>%
    mutate(mid = map2_chr(lid, mid, rename_mid))
tk1 = tk %>% filter(!is.na(re))# %>% inner_join(r0, by=c('lid','mid'))
tk2 = tk %>% filter(is.na(re))
make_kmer_motif <- function(kmer) create_motif(kmer, name=kmer, type='PPM', alphabet='DNA')
kmers = tk2 %>% rename(kmer=seq) %>% distinct(kmer) %>%
    mutate(mtf = map(kmer, make_kmer_motif))
#}}}

#{{{ collapse DREME motifs w. known cisbp motifs
#{{{ read in & merge
fi = sprintf("%s/dreme.rds", dirr)
rd = readRDS(fi)
diri = '~/projects/cre/data/01_tfbs'
fi = file.path(diri, '10.fam.rds')
rk = readRDS(fi)
mtfs_known = rk$fam$pwm
mtfs_known = rk$mtf$pwm

mtfs_dreme = rd$mtf
mtfs = c(mtfs_known, mtfs_dreme)
tm1 = rd %>% mutate(ctag='dreme') %>% rename(pwm=mtf) %>%
    mutate(conseq=map_chr(pwm, 'consensus'), name=conseq) %>%
    select(ctag, mid, name, conseq, pwm)
tm2 = rk$fam %>% mutate(ctag='cisbp', mid=fid, name=fname) %>%
    select(ctag, mid, name, conseq, pwm)
tm2 = rk$mtf %>% mutate(ctag='cisbp', mid, name=str_replace(gene, " .*", "")) %>%
    select(ctag, mid, name, conseq, pwm)
tm = tm1 %>% bind_rows(tm2)
#}}}

ti = tm
#{{{ # 1-pass clustering
cmp0 = compare_motifs(mtfs, method="PCC", min.mean.ic=.0, min.overlap=5, score.strat="a.mean")
dst0 = as.dist(1 - cmp0)
hcu0 = hclust(dst0, method='average')

x = cutree(hcu0, h=.1)
tx0 = tibble(mid = hcu0$labels, grp = as.integer(x))
tx0$grp[tx0$grp==0] = max(tx0$grp) + 1:sum(tx0$grp==0)
tx = tx0 %>%
    inner_join(ti, by='mid') %>%
    mutate(icscore = map_dbl(pwm, 'icscore')) %>%
    arrange(grp, desc(icscore)) %>%
    select(ctag, mid, name, icscore, grp, conseq, pwm) %>% print(n=30)
#
txs = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mids = str_c(mid, collapse=' '),
        names = str_c(name, collapse=' ')) %>%
    ungroup()
tx = tx %>% inner_join(txs, by='grp')
tx %>% count(grp) %>% arrange(desc(grp))
#tx %>% filter(grp==4) %>% print(n=30)
#}}}
# tx txs

#{{{ # 2-pass clustering
tx2 = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mid=mid[1], pwm=pwm[1]) %>% ungroup()

cmp1 = compare_motifs(tx2$pwm, method="PCC", min.mean.ic=.0,
                      min.overlap=6, score.strat="a.mean")
dst = as.dist(1 - cmp1)
hcu = hclust(dst, method='average')

x = cutree(hcu, h=.15)
tx3 = tibble(mid = hcu$labels, grp2 = as.integer(x)) %>%
    inner_join(tx2, by='mid') %>% select(grp, grp2)
#
tr = tx %>% inner_join(tx3, by='grp') %>%
    select(-grp) %>% rename(grp=grp2) %>%
    select(-mids, -names) %>%
    arrange(grp, ctag, desc(icscore))
trs = tr %>%
    group_by(grp) %>%
    summarise(n_mid = n(), n_cisbp = sum(ctag=='cisbp'),
              n_dreme = sum(ctag=='dreme'),
              mids = str_c(mid, collapse=' '),
              pwms=list(pwm),
              mid = mid[1], conseq = conseq[1], pwm=pwm[1],
              fname = name[1],
              fname_l = str_c(name, collapse=' ')) %>%
    ungroup()
tr %>% count(grp) %>% arrange(desc(grp))
trs %>% count(n_cisbp)
trs %>% filter(n_cisbp>0)
trs %>% filter(n_cisbp==0) %>% count(n_dreme)
#}}}
# tr trs

#{{{ similarity among motif clusters
mtfs = trs$pwm
cmp = compare_motifs(mtfs, method="PCC", min.mean.ic=0,
                     min.overlap = 6, score.strat="a.mean")
dst = as.dist(1 - cmp)
hcu = hclust(dst, method='average')

#{{{ distance heatmap
mids_sorted = hcu$labels[hcu$order]
tp = as_tibble(cmp) %>% mutate(mid = rownames(cmp)) %>%
    gather(mid2, score, -mid) %>%
    mutate(mid = factor(mid, levels=mids_sorted)) %>%
    mutate(mid2 = factor(mid2, levels=mids_sorted))
#
p = ggplot(tp) +
    geom_tile(aes(x=mid, y=mid2, fill=score)) +
    scale_x_discrete(position='bottom', expand=c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(name='PCC', na.value='grey50', colors=cols100) +
    #scale_fill_viridis(name=leg) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.5,.5,.5,.5))
fo = file.path(dirw, '10.heatmap.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

#{{{ ggtree for all motifs
tree = ape::as.phylo(hcu0)
tpt = tr %>% select(taxa=mid, grp,gene) %>%
    mutate(lgd=str_c(grp, str_sub(gene,1,20), sep=' '))
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=lgd), size=1) +
    scale_x_continuous(expand=expansion(mult=c(0,.2))) +
    scale_y_continuous(expand=expansion(mult=c(.001,.001))) +
    scale_color_aaas() +
    theme_tree2() +
    theme(plot.margin = margin(.3,.5,.3,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '10.tree.all.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 60)
#}}}

#{{{ ggtree for motif centroids
tree = ape::as.phylo(hcu)
tpt = trs %>% select(taxa=mid, grp,grpname) %>%
    mutate(lgd=str_c(grp, str_sub(grpname,1,20), sep=' '))
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=lgd), size=2) +
    scale_x_continuous(expand=expansion(mult=c(0,.2))) +
    scale_y_continuous(expand=expansion(mult=c(.01,.01))) +
    scale_color_aaas() +
    theme_tree2() +
    theme(plot.margin = margin(.3,.5,.3,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '10.tree.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 20)
#}}}
#}}}

#{{{ rename motif clusters
rename_mtf <- function(mtf, id) { mtf['name'] = id; mtf }
mtf = tr %>% mutate(grp = sprintf("g%04d", grp)) %>% rename(fid=grp)
grp = trs %>% mutate(grp = sprintf("g%04d", grp)) %>% rename(fid = grp) %>%
    mutate(pwm = map2(pwm, fid, rename_mtf))

fo = file.path(dirw, '01.motif.grp.rds')
res = list(mtf = mtf, grp = grp)
saveRDS(res, fo)
#}}}


#{{{ [obsolete] group motifs and save
dst = as.dist(1 - cmp)
hcu = hclust(dst, method='average')
#hcu$height = round(hcu$height,6)
#x = cutreeHybrid(hcu, as.matrix(dst), deepSplit=1, minClusterSize=1, pamStage=T, maxPamDist=.15)
x = cutree(hcu, h=.15)
tx0 = tibble(mid = hcu$labels, grp = sprintf("g%04d", as.integer(x)))
tx = tx0 %>%
    inner_join(tm, by='mid') %>%
    mutate(icscore = map_dbl(pwm, 'icscore')) %>%
    arrange(grp, desc(icscore)) %>%
    print(n=30)
grp = tx %>% arrange(grp, ctag, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(n_cisbp = sum(ctag=='cisbp'),
              n_dreme = sum(ctag=='dreme'),
              mids = str_c(mid, collapse=' '),
        mid = ifelse(n_cisbp>0, mid[n_cisbp>0][1], mid[n_dreme>0][1]),
        name = ifelse(n_cisbp>0, name[n_cisbp>0][1], name[n_dreme>0][1])
        ) %>%
    ungroup()
grp %>% count(n_cisbp)
grp %>% filter(n_cisbp==0) %>% count(n_dreme)
mtf = tx
#tx2 %>% count(n_dreme>0, n_cisbp) %>% print(n=40)

#{{{ ggtree for all motifs
tree = ape::as.phylo(hcu)
tpt = tm %>% select(taxa=mid, name)
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=name), size=1) +
    scale_x_continuous(expand=expansion(mult=c(0,.1))) +
    scale_y_continuous(expand=expansion(mult=c(.001,.001))) +
    scale_color_aaas() +
    theme_tree2() +
    theme(plot.margin = margin(.3,.5,.3,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '08.tree.all.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 60)
#}}}

fo = file.path(dirw, '01.motif.grp.rds')
res = list(mtf = mtf, grp = grp)
saveRDS(res, fo)
#}}}

find_known_mtf <- function(kmer, mtfs) {
    #{{{
    cmp = compare_motifs(c(list(kmer), mtfs), compare.to=1, method='PCC', min.mean.ic=0, score.strat='a.mean', min.overlap=6)
    if(length(cmp) == 0) {
        tibble()
    } else {
        cmp %>% as_tibble() %>%
        arrange(desc(score)) %>%
        select(tgt=target, PCC=score, logP=logPval)
    }
    #}}}
}
#}}}

#{{{ assess significance and store resutls
#{{{ read in
fi = file.path(dirw, '01.motif.grp.rds')
tg = readRDS(fi)
mtf = tg$mtf; grp=tg$grp
#
tls = tl0 %>% distinct(bat_mid, ng0, note) %>% filter(!is.na(note))
tln2 = tln %>% filter(pick) %>%
    select(bat_mid, ng0, ng, tg, ng_c, tg_c) %>%
    mutate(sids = map(tg, 'sid'), sids_c = map(tg_c,'sid')) %>%
    select(-tg, -tg_c, -ng0)
#{{{ read DREME results
rename_mid <- function(lid, mid) sprintf("%s_%03d", lid, as.numeric(str_replace(mid,'DREME-','')))
fi = sprintf("%s/dreme_kmer.tsv", dirr)
tk = read_tsv(fi) %>%
    mutate(mid = map2_chr(lid, mid, rename_mid))
tk0 = tk %>% filter(is.na(re)) %>% select(-re,-seq_rc) %>%
    group_by(mid) %>% summarise(kmers = list(seq)) %>% ungroup()
tk2 = tk %>% filter(is.na(re)) %>% select(mid, kmer=seq)
#}}}
filter_kmer_loc <- function(sids, sids_c, loc) {
    #{{{
    t1 = loc %>% filter(sid %in% sids) %>% mutate(type = 'cold/heat-responsive')
    t2 = loc %>% filter(sid %in% sids_c) %>% mutate(type = 'control')
    t1 %>% bind_rows(t2)
    #}}}
}
loc_freq <- function(tl, ng, ng_c) {
    #{{{
    to1 = tibble(type=c('cold/heat-responsive','control'), nc=c(ng,ng_c))
    to = tl %>% mutate(pos = round((start + 1 + end)/2)) %>%
        mutate(bin = cut(pos, breaks=seq(0,4000,by=200))) %>%
        mutate(bin = as.numeric(bin)) %>%
        #distinct(sid,type,srd,bin,umr)
        distinct(sid,type,bin,umr)
    tc1 = to %>% count(type, bin) %>% mutate(epi = 'raw')
    tc2 = to %>% filter(umr > 0) %>% count(type, bin) %>% mutate(epi = 'umr')
    #tc3 = to %>% filter(acrL > 0) %>% count(type, srd, bin) %>% mutate(epi = 'acrL')
    #tc4 = to %>% filter(acrE > 0) %>% count(type, srd, bin) %>% mutate(epi = 'acrE')
    rbind(tc1,tc2) %>% complete(type, bin, epi, fill=list(n=0)) %>%
        inner_join(to1, by='type') %>%
        mutate(p=n/nc)
    #}}}
}
loc_ftype <- function(tl) {
    #{{{
    tl %>% mutate(ftype = 'intergenic') %>%
        mutate(ftype = ifelse(intron > 0, 'intron', ftype)) %>%
        mutate(ftype = ifelse(utr3 > 0, 'utr3', ftype)) %>%
        mutate(ftype = ifelse(utr5 > 0, 'utr5', ftype)) %>%
        mutate(ftype = ifelse(cds > 0, 'cds', ftype)) %>%
        filter(umr > 0) %>%
        count(type, ftype) %>%
        group_by(type) %>% mutate(p = n / sum(n)) %>% ungroup()
    #}}}
}
assess_signif <- function(tl, ng, ng_c) {
    #{{{
    nh = tl %>% filter(type!='control') %>% distinct(sid) %>% pull(sid) %>% length()
    nh_c = tl %>% filter(type=='control') %>% distinct(sid) %>% length()
    nhu = tl %>% filter(type!='control',umr>0) %>% distinct(sid) %>% length()
    nhu_c = tl %>% filter(type=='control',umr>0) %>% distinct(sid) %>% length()
    fc = (nh/ng) / (nh_c/ng_c)
    fcu = (nhu/ng) / (nhu_c/ng_c)
    pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)
    pval.u = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F)
    tibble(nh=nh,nh_c=nh_c,nhu=nhu,nhu_c=nhu_c,fc=fc,fcu=fcu,pval=pval,pval.u=pval.u)
    #}}}
}

fi = file.path(dirw, '02_kmer_loc.rds')
kl = readRDS(fi)
#}}}

#{{{ first pass and save
x = tk %>% filter(!is.na(re)) %>% select(lid, mid) %>%
    inner_join(tk0, by='mid')
x1 = tl0 %>% filter(pick) %>% select(bat_mid, lid)
x2 = x %>% inner_join(x1, by='lid') %>%
    inner_join(tln2, by='bat_mid') %>%
    select(mid, kmers, bat_mid, ng, ng_c, sids, sids_c)

x3 = x2 %>% unnest(kmers) %>% rename(kmer=kmers) %>%
    inner_join(kl, by='kmer') %>%
    mutate(tl = pmap(list(sids, sids_c, loc), filter_kmer_loc)) %>%
    select(mid, kmer, bat_mid, ng, ng_c, tl) %>% unnest(tl) %>%
    group_by(mid, bat_mid, ng, ng_c) %>%
    nest() %>% rename(tl=data) %>% ungroup() %>%
    mutate(tp = pmap(list(tl, ng, ng_c), loc_freq)) %>%
    mutate(tf = map(tl, loc_ftype))

x4a = x3 %>% select(mid,tl) %>% unnest(tl) %>%
    distinct(mid,sid,type,umr) %>%
    group_by(mid) %>%
    summarise(nh = length(unique(sid[type!='control'])),
              nh_c = length(unique(sid[type=='control'])),
              nhu = length(unique(sid[type!='control' & umr > 0])),
              nhu_c = length(unique(sid[type=='control' & umr >0]))) %>%
    ungroup()
x4 = x3 %>% inner_join(x4a, by='mid') %>%
    mutate(fc = (nh/ng)/(nh_c/ng_c), fcu = (nhu/ng)/(nhu_c/ng_c)) %>%
    mutate(pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)) %>%
    mutate(pval.umr = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F)

mer = x4 %>% inner_join(mtf[,c('mid','fid')], by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid')
fo = file.path(dirw, '05.mtf.enrich.rds')
saveRDS(mer, fo)
#}}}

#{{{ heatmap
tp = x5 %>% filter(nh_c >= 25) %>%
    mutate(score = -log10(pval.umr)) %>%# filter(!is.infinite(score)) %>%
    filter(score >= 10) %>%
    mutate(score2 = -log10(pval)) %>%# filter(!is.infinite(score2)) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} ({ng0}) {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    arrange(ylab, fid, fname, desc(score)) %>%
    group_by(ylab, fid, fname) %>% slice(1) %>% ungroup() %>%
    mutate(lab = sprintf("%.0f %.0f %.1f %.1f", score2, score, fc, fcu)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
#fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
fids = rev(tpm$gid)
tp = tp %>% mutate(fid=factor(fid, levels=fids))
tps = tp %>% distinct(fid, name_s) %>% arrange(fid)
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/11.mtf.grp.enrich.pdf", dirw)
p %>% ggexport(filename=fo, width=10, height=80)
#}}}

fi = file.path(dirw, '05.mtf.enrich.rds')
mer = readRDS(fi)

#{{{ second pass to unify motifs within groups
x6 = mer %>% inner_join(mtf[,c('mid','fid')], by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid') %>%
    filter(nh_c >= 25, pval.umr < 5) %>%
    arrange(fid, fname, pval.umr) %>%
    group_by(fid, fname) %>% slice(1) %>% ungroup() %>%
    select(fid, mid)
x = tk %>% filter(!is.na(re)) %>% select(lid, mid) %>%
    inner_join(mtf[,c('mid','fid')], by='mid') %>% rename(mid0=mid) %>%
    inner_join(x6, by='fid') %>%
    inner_join(tk0, by='mid') %>%
    inner_join(tl0[,c('lid','bat_mid')], by='lid') %>%
    distinct(bat_mid, fid, mid, kmers)
#
x2 = x %>%
    inner_join(tln2, by='bat_mid') %>%
    select(fid, mid, kmers, bat_mid, ng, ng_c, sids, sids_c)

x3 = x2 %>% unnest(kmers) %>% rename(kmer=kmers) %>%
    inner_join(kl, by='kmer') %>%
    mutate(tl = pmap(list(sids, sids_c, loc), filter_kmer_loc)) %>%
    select(fid, mid, kmer, bat_mid, ng, ng_c, tl) %>% unnest(tl) %>%
    group_by(fid, mid, bat_mid, ng, ng_c) %>%
    nest() %>% rename(tl=data) %>% ungroup() %>%
    mutate(tp = pmap(list(tl, ng, ng_c), loc_freq)) %>%
    mutate(tf = map(tl, loc_ftype))

x4a = x3 %>% select(bat_mid,mid,tl) %>% unnest(tl) %>%
    distinct(bat_mid,mid,sid,type,umr) %>%
    group_by(bat_mid,mid) %>%
    summarise(nh = length(unique(sid[type!='control'])),
              nh_c = length(unique(sid[type=='control'])),
              nhu = length(unique(sid[type!='control' & umr > 0])),
              nhu_c = length(unique(sid[type=='control' & umr >0]))) %>%
    ungroup()
x4 = x3 %>% inner_join(x4a, by=c('bat_mid','mid')) %>%
    mutate(fc = (nh/ng)/(nh_c/ng_c), fcu = (nhu/ng)/(nhu_c/ng_c)) %>%
    mutate(pval = phyper(nh-1, nh_c, ng_c-nh_c, ng, lower.tail=F)) %>%
    mutate(pval.umr = phyper(nhu-1, nhu_c, ng_c-nhu_c, ng, lower.tail=F)) %>%
    inner_join(grp[,c('fid','fname')], by='fid')
mer = x4

fo = file.path(dirw, '07.mtf.grp.enrich.rds')
saveRDS(mer, fo)
#}}}

fi = file.path(dirw, '07.mtf.grp.enrich.rds')
mer = readRDS(fi)

#{{{ heatmap
tp = mer %>% filter(nh_c >= 25) %>%
    mutate(score = -log10(pval.umr)) %>%# filter(!is.infinite(score)) %>%
    filter(score >= 10) %>%
    mutate(score2 = -log10(pval)) %>%# filter(!is.infinite(score2)) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} ({ng0}) {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    arrange(ylab, fid, fname, desc(score)) %>%
    group_by(ylab, fid, fname) %>% slice(1) %>% ungroup() %>%
    mutate(lab = sprintf("%.0f %.0f %.1f %.1f", score2, score, fc, fcu)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
#fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
fids = rev(tpm$gid)
tp = tp %>% mutate(fid=factor(fid, levels=fids))
tps = tp %>% distinct(fid, name_s) %>% arrange(fid)
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=45, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/12.mtf.grp.enrich.pdf", dirw)
p %>% ggexport(filename=fo, width=10, height=80)
#}}}
#}}}



#{{{ ##assess significance of all found motifs [obsolete]
#{{{ run st.23.mmd.1.R
fi = file.path(dirr, '23_kmer_count.rds')
kmer_cnt = readRDS(fi)
kmer_hits <- function(kmers, kmer_cnt)
    kmer_cnt %>% filter(kmer %in% kmers) %>% pull(sids) %>% unlist() %>% unique()

fi = file.path(dirw, '01.motif.grp.rds')
res = readRDS(fi)

rename_mid <- function(lid, mid) sprintf("%s_%03d", lid, as.numeric(str_replace(mid,'DREME-','')))
fi = sprintf("%s/dreme_kmer.tsv", dirr)
tk = read_tsv(fi) %>% mutate(mid = map2_chr(lid, mid, rename_mid))
tk0 = tk %>% filter(is.na(re)) %>% select(-re,-seq_rc) %>%
    group_by(mid) %>% summarise(kmers = list(seq)) %>% ungroup()

x = tk %>% filter(!is.na(re)) %>% select(lid, mid) %>%
    inner_join(tk0, by='mid') %>%
    mutate(hits = map(kmers,kmer_hits, kmer_cnt=kmer_cnt))
x1 = tl0 %>% filter(pick) %>% select(bat_mid, bin_epi, lid)

x2 = x %>% inner_join(x1, by='lid') %>% rename(lid0=lid,bin_epi0=bin_epi) %>%
    inner_join(x1, by='bat_mid') %>%
    select(bat_mid, mid, lid0, bin_epi0, lid, bin_epi, hits) %>%
    inner_join(tl0[,c('lid','ng0','ng','tg','ng_c','tg_c')], by='lid') %>%
    mutate(sids = map(tg, 'sid'), sids_c = map(tg_c, 'sid')) %>%
    select(-tg, -tg_c)

x3 = x2 %>%# slice(1:300) %>%
    mutate(sids_p = map2(sids, hits, intersect)) %>%
    mutate(sids_n = map2(sids_c, hits, intersect)) %>%
    mutate(pos = map_int(sids_p, length)) %>%
    mutate(neg = map_int(sids_n, length)) %>%
    mutate(pval = phyper(pos-1, neg, ng_c-neg, ng, lower.tail=F)) %>%
    select(-sids,-sids_c,-hits)
#}}}

fi = file.path(dirw, '11.phyper.rds')
t_hy = readRDS(fi)

fi = file.path(dirw, '01.motif.grp.rds')
tg = readRDS(fi)
mtf = tg$mtf; grp=tg$grp
#
te = t_hy %>% inner_join(mtf[,c('mid','fid')],by='mid') %>%
    inner_join(grp[,c('fid','fname')], by='fid')

#{{{ all enriched motif: +/-2k:umr
tls = tl0 %>% distinct(bat_mid, note)
tp0 = te %>% filter(bin_epi=='+/-2k:umr') %>%
    mutate(score = -log10(pval)) %>% filter(!is.infinite(score))
tp = tp0 %>%
    filter(score >= 10) %>%
    inner_join(tls, by='bat_mid') %>%
    mutate(ylab = glue("{bat_mid} {note}"),
           ylab = fct_reorder(ylab, as.numeric(bat_mid))) %>%
    group_by(ylab, fid, fname) %>%
    summarise(score = max(score)) %>% ungroup() %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(name_s = str_sub(str_c(fid, fname, sep=' '), 1, 50))
tpm = tp %>% select(fid, ylab, score) %>% complete(fid,ylab,fill=list(score=0)) %>%
    spread(ylab, score) %>% rename(gid=fid)
nrow(tpm)
fids = hc_order_row(tpm, cor.opt='euclidean', hc.opt='average')
tps = tp %>% distinct(fid, name_s) %>% arrange(desc(fid))
tp = tp %>% mutate(fid=factor(fid, levels=fids))
swit = (min(tp$score) + max(tp$score)) / 2
tpy = tibble(x=cumsum(c(5,3,3,4))+.5)
p = ggplot(tp, aes(x=ylab,y=fid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    geom_vline(xintercept=tpy$x, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tps$fid, labels=tps$name_s, expand=c(0,0)) +
    scale_fill_gradientn(name='-log10(phyper)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1.3,.3,.3), ygrid=T,
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=75, hjust=0, vjust=.5, size=7.5)) +
    theme(axis.text.y = element_text(size=7)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
fo = sprintf("%s/14.all.motif.grps.pdf", dirw)
p %>% ggexport(filename=fo, width=7, height=80)
#}}}

grp %>% filter(grp=='g0119')
grp %>% filter(str_detect(mid,'l0790_001')) %>% pull(mids)
grp %>% filter(str_detect(mid,'f137'))
mtf0 = rk$fam %>% filter(fid=='f137') %>% pull(pwm)
mtf1 = mtf %>% filter(mid=='l0790_001') %>% pull(pwm)
mtf0a = rk$mtf %>% filter(mid=='M06942_2.00') %>% pull(pwm)
mtf2 = mtf %>% filter(mid=='l1086_001') %>% pull(pwm)
mtf3 = mtf %>% filter(mid=='l0585_012') %>% pull(pwm)
compare_motifs(c(mtf0,mtf1,mtf2), method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)
compare_motifs(c(mtf3,mtfs_known), compare.to=1, method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)
compare_motifs(c(mtf2,mtfs_all), compare.to=1, method="PCC", min.mean.ic=0, score.strat="a.mean", min.overlap=6)

fo=file.path(dirw,'tmp.pdf')
pdf(file=fo, width=6,height=16)
midsd = grp %>% filter(str_detect(mid,'l0790_001')) %>% pull(mids)
mtfsd = mtf %>% filter(str_detect(midsd, mid)) %>% pull(pwm)
view_motifs(c(mtf0,mtf0a,mtfsd), method="PCC",min.mean.ic=0,min.overlap=6)
dev.off()
#}}}

#{{{ #single bat_mid motif plot
i = 1
tg0 = tl0 %>% filter(bat_mid == bat_mids[i]) %>% pull(tg)
tg0 = tg0[[1]] %>% select(-size)
tp0 = r$mtf %>% inner_join(ts, by=c('bat_mid','bin_epi')) %>%
    filter(bat_mid == bat_mids[i]) %>%
    select(grp,grpname,conseq,cnt)
tp1 = tp0 %>% mutate(nhits = map_int(cnt, nrow)) %>%
    arrange(grp,desc(nhits)) %>% group_by(grp) %>%
    summarise(conseq = conseq[1]) %>% ungroup()
tp = tp0 %>% select(-conseq) %>% inner_join(tp1, by='grp') %>%
    unnest(cnt) %>%
    group_by(grp, grpname, conseq, sid) %>%
    summarise(nseq=max(nseq)) %>% ungroup() %>%
    inner_join(tg0, by='sid') %>%
    mutate(score=nseq) %>%
    mutate(lab = number(score, accuracy=1)) %>%
    mutate(grpname = ifelse(grpname=='', str_c(grp,conseq,sep=' '), grpname)) %>%
    mutate(grpname_s = str_sub(grpname, 1, 15))
#
tps = tp %>% distinct(grpname, grpname_s)
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=grpname,y=gid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=1.5) +
    scale_x_discrete(breaks=tps$grpname, labels=tps$grpname_s, expand=c(0,0), position='top') +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,3.3,.3,.3),
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=25, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/09.%d.pdf", dirw, i), width = 8, height = 12)
#}}}


fi = file.path(dirw, '07.mtf.grp.enrich.rds')
mer = readRDS(fi)

# motif group meta plot
fids = c(
'g0437',
'g0098'
'g0054',
'g0645',
'g0024',
'g0130',
'g0014',
'g0050',
'g0009',
'g0290',
'g0722',
'g0191',
'g0154'
)
mtf_meta_plot <- function(fid, fname, mer, dirw=dirw) {
# {{{ plot
fname=grp %>% filter(grp == grp1) %>% pull(fname)
tp = mer %>% filter(fid==grp1) %>%
    mutate(pan = sprintf("%s n=%d fc=%.1f fcu=%.1f pval=%.2g pal.umr=%.2g", bat_mid, ng, fc, fcu, pval, pval.umr)) %>%
    select(pan, tp) %>% unnest(tp) #%>%
    #mutate(p = ifelse(srd == '-', -p, p))
#tp1 = tp %>% filter(srd=='+')
#tp2 = tp %>% filter(srd=='-')
types = c("cold/heat-responsive",'control')
tpx = tibble(x=c(.5,10.5,20.5),lab=c('-2kb','TSS','+2kb'))
p = ggplot(tp, aes(x=bin,y=p,col=epi)) +
    geom_line(aes(linetype=type), size=.5, na.rm = F) +
    geom_point(aes(shape=type), size=1, na.rm = F) +
    scale_x_continuous(expand=expansion(mult=c(.05,.05)),breaks=tpx$x,labels=tpx$lab) +
    scale_y_continuous(name="Proportion of TFBS", expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    #scale_shape(labels=types) +
    #scale_linetype(labels=types) +
    facet_wrap(~pan, scale='free', ncol=2) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           strip.style='white',margin = c(.3,.3,.3,.3),
           xgrid=T, xtick=T, ytick=T, ytitle=T,xtext=T, ytext=T) +
    guides(fill=F)
fp = sprintf("%s/%s.%s.pdf", dirw, grp1, fname)
p %>% ggexport(filename=fp, width = 10, height = 10)
#}}}
}
tibble(fid=fids) %>% inner_join(grp, by='fid') %>%
    mutate(p=pmap(list(fid,fname,mtf_meta_plot, mer=mer, dirw=dirw)))

