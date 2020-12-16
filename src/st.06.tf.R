source('functions.R')
dirw = file.path(dird, '06_tf_list')

#{{{ read cbf/hsf list
tac = read_gspread(book='stress_TF_list', sheet='cold')
tah = read_gspread(book='stress_TF_list', sheet='heat')

tac1 = tac %>% select(gid=Zm_gid, fam=family, name=Zm_alias, name2=At_alias, At_gid) %>%
    mutate(stress='cold')
tah1 = tah %>% select(gid=Zm_gid, fam=family, name=Zm_alias, name2=At_alias, At_gid) %>%
    mutate(stress='heat')
ta = tac1 %>% bind_rows(tah1) %>%
    filter(!is.na(gid), gid != 'N/A') %>%
    mutate(gid = map(gid, str_split, pattern='[,/]')) %>%
    mutate(gid = map(gid, 1)) %>%
    unnest(gid) %>%
    mutate(At_gid = str_to_upper(str_replace(At_gid, "\\.[1-9]+$", ''))) %>%
    #mutate(name = map_chr(name, get_first <- function(x) str_split(x,',')[[1]][1])) %>%
    select(stress, gid, fam, name, name2, At_gid)
ta %>% distinct(stress,gid,fam) %>% count(stress, gid) %>% filter(n>1)

ta2 = ta %>% arrange(stress, gid, fam, At_gid) %>%
    group_by(stress, gid, fam) %>%
    summarise(name=str_c(name, collapse=','),
              name_At=str_c(name2, collapse=','),
              At_gid=At_gid[1]) %>% ungroup()
to = ta2 %>% separate(gid, c('gid','tid'), sep='_', fill='right') %>% select(-tid)
to %>% count(stress, gid)
to %>% count(stress, fam) %>% arrange(desc(n)) %>%
    print(n=30)

fo = file.path(dirw, '10.stress.tf.tsv')
write_tsv(to, fo, na='')
#}}}

fi = '~/projects/grn/data/08_y1h/tf45.xlsx'
ti = read_xlsx(fi)
tf45 = ti %>% select(tf=1,gid=3) %>% distinct(tf,gid)

tag = 'degB'
fd = glue("{dird}/17_cluster//50_modules/{tag}.rds")
md = readRDS(fd)
tmd = md %>% select(cid, cond,note, gids) %>%
    unnest(gids) %>% rename(gid=gids) %>%
    separate(gid, c("gt",'gid'), sep='_') %>% select(-gt)

to = tmd %>% inner_join(tf45, by='gid') %>%
    select(cond,note,gid,tf)
