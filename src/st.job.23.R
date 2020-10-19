source('functions.R')
dirw = file.path(dird, '23_mmd')
#
fl = file.path(dirw, '../21_seq/15.rds')
tl0 = readRDS(fl)
tl = tl0 %>% filter(pick) %>% select(lid,ng0,ng,tg,clid,ng_c,tg_c)

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

fo = file.path(dirw, '11.phyper.rds')
saveRDS(x3, fo)
