source('functions.R')
require(universalmotif)
#require(dynamicTreeCut)
dirw = file.path(dird, '23_mmd')
#
fl = file.path(dird, '21_seq/15.rds')
tl = readRDS(fl)

#{{{ read DREME motifs/kmers
fk = sprintf("%s/23.kmer.tsv", dirr)
fm = sprintf("%s/23.kmer.motif.tsv", dirr)
tkk = read_tsv(fk) %>% select(-seq_rc) %>% rename(kmer=seq)
tkm = read_tsv(fm) %>% select(-seq_rc) %>% rename(kmer=seq)
make_kmer_motif <- function(kmer) create_motif(kmer, name=kmer, type='PPM', alphabet='DNA')
kmers = tkk %>% distinct(kmer) %>% mutate(mtf = map(kmer, make_kmer_motif))
tk0 = tkk %>% group_by(mid) %>% summarise(kmers = list(kmer)) %>% ungroup()
tk = tkm %>% inner_join(tk0, by='mid')
#
fi = sprintf("%s/23.dreme.rds", dirr)
rd = readRDS(fi)
#
km = tl %>% filter(pick) %>% select(lid, bat, note) %>%
    inner_join(tk, by='lid') %>%
    inner_join(rd, by=c('lid','mid')) %>%
    select(bat, note, lid, mid, kmer, kmers, pval, mtf)

fo = file.path(dirw, "00.kmers.rds")
saveRDS(km, fo)
#}}}

#{{{ collapse DREME motifs w. known cisbp motifs
#{{{ read in & merge
fi = file.path(dirw, "00.kmers.rds")
km = readRDS(fi)
diri = '~/projects/cre/data/01_tfbs'
fi = file.path(diri, '10.fam.rds')
rk = readRDS(fi)
mtfs_known = rk$fam$pwm
mtfs_known = rk$mtf$pwm

mtfs_dreme = km$mtf
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
    mutate(fname=str_replace_all(fname, ".*\\|", "")) %>%
    mutate(fname=str_replace_all(fname, "\\[.*", "")) %>%
    mutate(fname=str_replace_all(fname, "-", "")) %>%
    mutate(pwm = map2(pwm, fid, rename_mtf))
#}}}

mtf1 = mtf %>% select(mid, fid)
grp1 = grp %>% mutate(cisbp=n_cisbp>0) %>% select(fid,fname,cisbp)
km2 = km %>% inner_join(mtf1, by='mid') %>%
    inner_join(grp1, by='fid')
fo = file.path(dirw, '01.kmer.grp.rds')
r = list(kmer=km2, mtf = mtf, grp = grp)
saveRDS(r, fo)
#}}}

#{{{ summerize kmers found in each module/searching parameter
fr = file.path(dirw, '01.kmer.grp.rds')
r = readRDS(fr)

#{{{ num. motifs found in each module
tps = r$kmer %>% group_by(bat,note) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup() %>%
    inner_join(tl %>% distinct(bat,note,bat_mid), by=c('bat','note'))
tp1 = r$kmer %>% group_by(lid) %>%
    summarise(n_grp=length(unique(fid))) %>% ungroup()
tp = tl %>% select(lid,bat,mid,bat_mid, bin_epi, ng0, note) %>%
    inner_join(tp1, by='lid') %>% rename(score = n_grp) %>%
    mutate(lab = score)
tps2 = tp %>% distinct(bat_mid) %>% arrange(desc(bat_mid)) %>% mutate(i=1:n())
tps = tps %>% inner_join(tps2, by='bat_mid')
tp = tp %>%
    inner_join(tps2, by='bat_mid') %>%
    mutate(bin_epi = factor(bin_epi, levels=unique(tp$bin_epi)))
tpy = tp %>% distinct(i, bat, mid, ng0, note) %>%
    mutate(ylab = sprintf("%s: %s (%d)", bat, note, ng0))
tp1 = tibble(o=cumsum(c(5,4,5))+.5)
tp2 = tibble(o=seq(4,to=20,by=4)+.5)

#{{{ all motifs
swit = (min(tp$score) + max(tp$score)) / 2
p_all = ggplot(tp, aes(x=bin_epi,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_hline(yintercept=tp1$o, color='purple') +
    geom_vline(xintercept=tp2$o, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_continuous(breaks=tpy$i, labels=tpy$ylab, expand=c(0,0),
                     sec.axis=sec_axis(~., breaks=as.numeric(tps$i), labels=tps$n_grp)) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2,.3,.3),
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
#}}}
#{{{ selected motifs
tpf = tp %>% filter(str_detect(bin_epi, "\\+/\\-2k"))
swit = (min(tpf$score) + max(tpf$score)) / 2
p_sel = ggplot(tpf, aes(x=bin_epi,y=i)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_hline(yintercept=tp1$o, color='purple') +
    geom_vline(xintercept=tp2$o, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_continuous(breaks=tpy$i, labels=tpy$ylab, expand=c(0,0),
                     sec.axis=sec_axis(~., breaks=as.numeric(tps$i), labels=tps$n_grp)) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,1,.3,.3),
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
#}}}
p_sel %>% ggexport(filename =sprintf("%s/10.n.mtf.2.pdf", dirw), width=4, height=5)
saveRDS(p_sel, file.path(dirf, 'f.3b.rds'))
p_all %>% ggexport(filename =file.path(dirf, "sf6.pdf"), width=8, height=5)
#}}}

#{{{ top 40 motifs found in each module
r1 = r$kmer %>% select(mid,fid,fname,cisbp, lid,pval) %>%
    inner_join(tl[,c('lid','bat_mid','bat','note','bin_epi','ng0')], by='lid')
tp = r1 %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid, bat, note, ng0, fid, fname, cisbp) %>%
    summarise(pval = min(pval)) %>% ungroup() %>%
    arrange(bat_mid, pval) %>%
    group_by(bat_mid, bat, note) %>%
    slice(1:40) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = ifelse(cisbp, fname, '')) %>%
    mutate(i = factor(i, levels=50:1))
tpx = tp %>% distinct(bat_mid, bat, ng0, note) %>%
    mutate(xlab = sprintf("%s: %s (%d)", bat, note, ng0))
tp1 = tibble(o=cumsum(c(6,5,4))+.5)

swit = (min(tp$score) + max(tp$score)) / 2
p4 = ggplot(tp, aes(x=bat_mid,y=i)) +
    geom_tile(aes(fill=score), na.rm = F, size=.5, color='white') +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_vline(xintercept=tp1$o, color='blue') +
    scale_x_discrete(breaks=tpx$bat_mid, labels=tpx$xlab, expand=c(0,0), position='top') +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='-log10(Pval)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=T,
           margin = c(.3,4.3,.3,.3), legend.vjust=-1.7,
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(legend.position = c(1,0), legend.justification = c(0,0)) +
    theme(axis.text.x = element_text(angle=25, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p4 %>% ggexport(filename=glue("{dirf}/sf7.pdf"), width=10, height=7)
#}}}
#}}}




