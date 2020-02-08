source('functions.R')
dirw = file.path(dird, '11_brb')

#{{{ [obsolete] read1 UMI
#{{{ read umi_cnt and save
ti = crossing(x=1:12, y=1:8) %>%
    mutate(yl = LETTERS[y]) %>%
    mutate(fn = str_c(yl,x,sep='')) %>%
    rbind(tibble(x=13, y=9, yl='I',fn='undetermined')) %>%
    mutate(fi = sprintf("%s/04_umi_cnt/%s/%s.tsv", diri, sid, fn)) %>%
    mutate(res = map(fi, read_tsv, col_names=c('umi','cnt'))) %>%
    select(x, y, yl, fn, res) %>%
    unnest()

fo = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
saveRDS(ti, file=fo)
#}}}

fi = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
ti = readRDS(fi)
n_read = sum(ti$cnt); n_umi = nrow(ti); pct_dedup = percent(n_umi/n_read)
tit = sprintf("%s UMIs / %s total reads = %s", number(n_umi,big.mark=','), number(n_read,big.mark=','), pct_dedup)

#{{{ plot reads + umi
tp = ti %>% group_by(x,y,fn) %>%
    summarise(cnt = sum(cnt), umi=n()) %>% ungroup() %>%
    mutate(lab = str_c(number(cnt), number(umi), sep="\n"))
mid = (min(tp$cnt) + max(tp$cnt)) / 2
tp = tp %>% mutate(color = ifelse(cnt < mid, 'white','black'))
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=cnt), color='black') +
    geom_text(aes(x,y, label=lab, color=cnt<mid), hjust=1, size=2.5, nudge_x=.4) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/01.nseq.%s.pdf", dirw, sid)
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ fastqc
f_qc = sprintf("%s/02_trim/%s_2_fastqc.zip", diri, sid)
qc = qc_read(f_qc, modules=c("Basic Statistics", "Sequence Duplication Levels"))
n_read2 = as.integer(qc$basic_statistics$Value[[4]])
pct_dedup2 = qc$total_deduplicated_percentage
cat(str_c(number(n_read2), pct_dedup2, '\n', sep=', '))
#}}}
#}}}

#{{{ read DGE and store
read_dge <- function(ft)
    read_tsv(ft) %>% rename(gid=1) %>% gather(coord, cnt, -gid)

sids = str_c("batch", c('1','2a','2b','3'), sep="")
td = tibble(sid=sids) %>%
    mutate(fd1=sprintf("%s/cache/21_dge/%s/output.dge.reads.txt", dird, sid)) %>%
    mutate(fd2=sprintf("%s/cache/21_dge/%s/output.dge.umis.txt", dird, sid)) %>%
    mutate(td1 = map(fd1, read_dge)) %>%
    mutate(td2 = map(fd2, read_dge))
td1 = td %>% select(sid, td1) %>% unnest(td1) %>% rename(nr=cnt)
td2 = td %>% select(sid, td2) %>% unnest(td2) %>% rename(nu=cnt)
cnts = td1 %>% inner_join(td2, by=c('sid','coord','gid')) %>%
    select(sid, coord, gid, nr, nu)

fo = file.path(dirw, '01.cnts.rds')
saveRDS(cnts, fo)

#{{{ ngene umi
l2n = 1:8; names(l2n) = LETTERS[l2n]
min_nr = 2; min_nu = 2 #min_nr = 10; min_nu = 10
tp = cnts %>% dplyr::rename(SampleID = coord) %>%
    #filter(sid == !!sid) %>%
    rename(n_read = nr, n_umi = nu) %>%
    group_by(sid, SampleID) %>%
    summarise(nr = sum(n_read >= min_nr), nu = sum(n_umi >= min_nu)) %>%
    ungroup() %>%
    mutate(yl = str_sub(SampleID, 1, 1), x = str_sub(SampleID, 2)) %>%
    mutate(y = as.integer(l2n[yl]), x = as.integer(x))

mid = (min(tp$nu) + max(tp$nu)) / 2
tit=sprintf('# genes with >= %d UMIs', min_nu)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nu), color='black') +
    geom_text(aes(x,y, label=number(nu), color=nu<mid), hjust=1, size=2.5, nudge_x=.3) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    facet_wrap(~sid, nrow=2) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/05.ngene.umi.pdf", dirw)
ggsave(fo, p, width=10, height=7)
#
mid = (min(tp$nr) + max(tp$nr)) / 2
tit=sprintf('# genes with >= %d reads', min_nr)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nr), color='black') +
    geom_text(aes(x,y, label=number(nr), color=nr<mid), hjust=1, size=2.5, nudge_x=.3) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    facet_wrap(~sid, nrow=2) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/05.ngene.read.pdf", dirw)
ggsave(fo, p, width=10, height=7)
#}}}
#}}}

fh = file.path(dird, '01_exp_design/05.brb.local.tsv')
th = read_tsv(fh) %>% select(sid=Tissue,coord=file_prefix,SampleID=Treatment)
#{{{ normalize read counts
fi = file.path(dirw, '01.cnts.rds')
cnts = readRDS(fi)

t_rc = cnts %>%
    inner_join(th, by=c('sid','coord')) %>%
    mutate(SampleID = str_c(sid,coord,sep='_')) %>%
    select(SampleID, gid, ReadCount = nu)
res = readcount_norm(t_rc)
tm_u = res$tm

t_rc = cnts %>%
    inner_join(th, by=c('sid','coord')) %>%
    mutate(SampleID = str_c(sid,coord,sep='_')) %>%
    select(SampleID, gid, ReadCount = nr)
res = readcount_norm(t_rc)
tm_r = res$tm

res = list(cnts=cnts, tm_u=tm_u, tm_r = tm_r)
fo = file.path(dirw, '03.cpm.rds')
saveRDS(res, fo)
#}}}

fi = file.path(dirw, '03.cpm.rds')
res = readRDS(fi)

tid = 'set1'; sids = 'b01'
th = res$th %>% filter(batch %in% sids) %>%
    mutate(lab = str_c(Genotype, Treatment, Timepoint, sep='_')) %>%
    mutate(grp = str_c(Genotype, Treatment, sep='_'))
tm = res$tm_r %>% filter(SampleID %in% th$SampleID)


