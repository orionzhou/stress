source('functions.R')
dirw = file.path(dird, '06_tf_list')

#{{{ read cbf/hsf list
tac = read_gspread(book='stress_TF_list', sheet='cold')
tah = read_gspread(book='stress_TF_list', sheet='heat')

tac1 = tac %>% select(gid=Zm_gid, fam=family, name=Zm_alias, name2=At_alias) %>%
    mutate(stress='cold')
tah1 = tah %>% select(gid=Zm_gid, fam=family, name=Zm_alias, name2=At_alias) %>%
    mutate(stress='heat')
ta = tac1 %>% bind_rows(tah1) %>%
    filter(!is.na(gid), gid != 'N/A') %>%
    mutate(gid = map(gid, str_split, pattern='[,/]')) %>%
    mutate(gid = map(gid, 1)) %>%
    unnest(gid) %>%
    #mutate(name = map_chr(name, get_first <- function(x) str_split(x,',')[[1]][1])) %>%
    select(stress, gid, fam, name, name2)
ta %>% distinct(stress,gid,fam) %>% count(stress, gid) %>% filter(n>1)

ta2 = ta %>% group_by(stress, gid, fam) %>%
    summarise(name=str_c(name, collapse=','),
              name_At=str_c(name2, collapse=',')) %>% ungroup()
to = ta2
to %>% count(stress, gid)
to %>% count(stress, fam) %>% arrange(desc(n)) %>%
    print(n=30)

fo = file.path(dirw, '10.stress.tf.tsv')
write_tsv(to, fo, na='')
#}}}

