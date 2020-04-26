source('functions.R')
dirw = file.path(dird, '22_mmm')

#{{{ get known cold/heat TF motifs
#{{{ read cbf/hsf list
tac = read_gspread(book='stress_TF_list', sheet='cold')
tah = read_gspread(book='stress_TF_list', sheet='heat')
tac1 = tac %>% select(At_gid, At_alias, Zm_gid, Zm_alias) %>% mutate(stress='cold')
tah1 = tah %>% select(At_gid, At_alias, Zm_gid, Zm_alias) %>% mutate(stress='heat')
ta = tac1 %>% bind_rows(tah1) %>%
    filter(!is.na(At_gid)) %>%
    mutate(name = map_chr(At_alias, get_first <- function(x) str_split(x,',')[[1]][1])) %>%
    mutate(name2 = map_chr(Zm_alias, get_first <- function(x) str_split(x,',')[[1]][1])) %>%
    distinct(stress,At_gid)
#}}}

fm = '/home/springer/zhoux379/projects/cre/data/01_tfbs/01.motifs.tsv'
tm0 = read_tsv(fm)
#
tm = tm0 %>% select(-ctag) %>% inner_join(ta, by=c('gid'='At_gid'))
#}}}
tm1 = tm %>% group_by(stress,motif,src_type,src_id) %>%
    summarise(name = str_c(unique(name), collapse=',')) %>% ungroup()

fi = file.path(dirw, '../21_seq/10.rds')
tl = readRDS(fi) %>%
    mutate(db_size=map_dbl(tg, get_db_size <- function(x) sum(x$size))) %>%
    select(-tg)

fi = '/home/springer/zhoux379/projects/nf/raw/mmm/fimo.tsv'
ti = read_tsv(fi, col_names=c("lid",'motif','bp','score')) %>%
    inner_join(tl, by='lid') %>%
    separate(bat, c('stress','drc'), sep="_", remove=F) %>%
    mutate(pscore=score/db_size) %>%
    mutate(pbp=bp/db_size) %>%
    inner_join(tm1, by=c('motif'))

names = c('DREB1A', 'DREB1B', 'DREB1C', 'DREB1D', 'DREB2A',
          'ERF5','ERF054','ERF011','CRF3','RAV1',
          'HSFA6A','HSFB2A','HSFB2B')
#{{{
for (name in names) {
tp = ti %>% mutate(score = pscore) %>%
    filter(name == !!name) %>%
    mutate(motif = sprintf("%s %s %s", name, src_type, src_id)) %>%
    mutate(bat_mid = factor(bat_mid, levels=rev(levels(tl$bat_mid))))
swit = (min(tp$score) + max(tp$score)) / 2
p = ggplot(tp, aes(x=bin_epi,y=bat_mid)) +
    geom_tile(aes(fill=score)) +
    #geom_text(aes(label=score, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expansion(mult=c(0,0)), position='top') +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name='number genes in cluster',colors=cols100) +
    #scale_fill_viridis(name='normalized eigengene value') +
    scale_color_manual(values=c('black','white')) +
    facet_wrap(.~motif, nrow=2, dir='v', scale='free_y') +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.3,1.9,.3,.3),
           xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=0, vjust=0, size=7.5)) +
    theme(axis.text.y = element_text(size=7.5)) +
    guides(color = F, fill=F)
#
fo = sprintf("%s/03.%s.pdf", dirw, name)
p %>% ggexport(filename = fo, width = 12, height = 9)
}
#}}}



