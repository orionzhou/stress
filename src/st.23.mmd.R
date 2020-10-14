source('functions.R')
require(universalmotif)
#require(dynamicTreeCut)
dirw = file.path(dird, '23_mmd')
#
fl = file.path(dird, '21_seq/15.rds')
tl = readRDS(fl)

#{{{ read DREME kmers
fk = sprintf("%s/23.kmer.tsv", dirr)
fm = sprintf("%s/23.kmer.motif.tsv", dirr)
tkk = read_tsv(fk) %>% select(-seq_rc) %>% rename(kmer=seq)
tkm = read_tsv(fm) %>% select(-seq_rc) %>% rename(kmer=seq)
make_kmer_motif <- function(kmer) create_motif(kmer, name=kmer, type='PPM', alphabet='DNA')
kmers = tkk %>% distinct(kmer) %>% mutate(mtf = map(kmer, make_kmer_motif))
tk0 = tkk %>% group_by(mid) %>% summarise(kmers = list(kmer)) %>% ungroup()
tk = tkm %>% inner_join(tk0, by='mid')

kn = tl %>% filter(pick) %>% select(lid, bat, note) %>%
    inner_join(tk, by='lid') %>%
    select(bat, note, kmers)

fo = file.path(dirw, "04.kmer.note.rds")
saveRDS(kn, fo)
#}}}

#{{{ collapse DREME motifs w. known cisbp motifs
#{{{ read in & merge
fi = sprintf("%s/23.dreme.rds", dirr)
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
    mutate(fname=str_replace_all(fname, ".*\\|", "")) %>%
    mutate(fname=str_replace_all(fname, "\\[.*", "")) %>%
    mutate(fname=str_replace_all(fname, "-", "")) %>%
    mutate(pwm = map2(pwm, fid, rename_mtf))
#}}}

fo = file.path(dirw, '01.motif.grp.rds')
res = list(mtf = mtf, grp = grp)
saveRDS(res, fo)
#}}}

#{{{ summerize kmers found in each module/searching parameter
fi = file.path(dirw, '01.motif.grp.rds')
r = readRDS(fi)

#{{{ num. motifs found in each module
tk1 = tk %>% select(lid,kmer,pval)
tp = tl %>% select(lid,bat,mid,bat_mid, bin_epi, ng0, note) %>%
    inner_join(tk1, by='lid') %>%
    count(lid,bat,mid,bat_mid, bin_epi,ng0,note) %>% rename(score=n) %>%
    mutate(lab = score) %>%
    mutate(bat_mid = factor(bat_mid, levels=rev(levels(tl$bat_mid)))) %>%
    mutate(bin_epi = factor(bin_epi, levels=levels(tl$bin_epi)))
tpy = tp %>% distinct(bat_mid, bat, mid, ng0, note) %>%
    mutate(ylab = sprintf("%s: %s (%d)", bat, note, ng0))
tp1 = tibble(o=cumsum(c(5,4,5))+.5)
tp2 = tibble(o=seq(4,to=20,by=4)+.5)

swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=bat_mid)) +
    geom_tile(aes(fill=score), na.rm = F) +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_hline(yintercept=tp1$o, color='purple') +
    geom_vline(xintercept=tp2$o, color='blue') +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(breaks=tpy$bat_mid, labels=tpy$ylab, expand=c(0,0)) +
    scale_fill_gradientn(name='# hits',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='none', legend.dir='v', legend.title=F,
           margin = c(.3,2,.3,.3),
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/10.n.mtf.pdf", dirw), width = 8, height = 6)
#}}}

#{{{ top 50 motifs found in each module
r1 = r$mtf %>% filter(ctag=='dreme') %>% select(mid,fid) %>%
    inner_join(r$grp[c('fid','fname','n_cisbp')], by='fid') %>%
    select(mid, fid, fname, n_cisbp) %>%
    inner_join(tk, by='mid') %>%
    select(mid,fid,fname,n_cisbp, lid,pval) %>%
    inner_join(tl[,c('lid','bat_mid','bat','note','bin_epi','ng0')], by='lid')
tp = r1 %>% arrange(bat_mid, pval) %>%
    group_by(bat_mid, bat, note, ng0, fid, fname, n_cisbp) %>%
    summarise(pval = min(pval)) %>% ungroup() %>%
    arrange(bat_mid, pval) %>%
    group_by(bat_mid, bat, note) %>% slice(1:50) %>%
    mutate(i = 1:n()) %>% ungroup() %>%
    mutate(score = -log10(pval)) %>%
    mutate(lab = ifelse(n_cisbp>0, fname, '')) %>%
    mutate(i = factor(i, levels=50:1))
tpx = tp %>% distinct(bat_mid, bat, ng0, note) %>%
    mutate(xlab = sprintf("%s: %s (%d)", bat, note, ng0))
tp1 = tibble(o=cumsum(c(6,5,4))+.5)

swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bat_mid,y=i)) +
    geom_tile(aes(fill=score), na.rm = F, size=.5, color='white') +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2.5) +
    geom_vline(xintercept=tp1$o, color='blue') +
    scale_x_discrete(breaks=tpx$bat_mid, labels=tpx$xlab, expand=c(0,0), position='top') +
    scale_y_discrete(expand=expansion(mult=c(0,0))) +
    scale_fill_gradientn(name='-log10(Pval)',colors=cols100v) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           margin = c(.5,4.3,.3,.3), legend.vjust=-1.7,
           xtick=T, ytick=F, xtitle=F, xtext=T, ytext=F) +
    theme(axis.text.x = element_text(angle=25, hjust=0, vjust=.5, size=7.5)) +
    #theme(axis.text.y = element_markdown(size=7.5)) +
    guides(color=F)
p %>% ggexport(filename =sprintf("%s/10.top50.mtf.pdf", dirw), width = 10, height = 10)
#}}}
#}}}




