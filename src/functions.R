#{{{ load & read
require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
require(lubridate)
dirg = '~/data/genome'
dirp = "~/projects/stress"
dird = file.path(dirp, 'data')
dirr = file.path(dird, 'raw')
dirf = file.path(dird, '95_figures')
gcfg = read_genome_conf()
fh = file.path(dird, 'samples.xlsx')
th = read_xlsx(fh, sheet='merged') %>%
    mutate(Time = sprintf("%02d:%02d", hour(Time), minute(Time)))
th2 = th %>% mutate(Tissue = 'leaf') %>%
    mutate(Treatment=str_c(Experiment, Treatment, Timepoint, sep='_')) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate=Rep)
#}}}


