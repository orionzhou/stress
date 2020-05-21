args = commandArgs(trailingOnly=T)
source('functions.R')
dirw = file.path(dird, '16_ase')
#{{{ functions
require(stats4)
require(multidplyr)
require(bbmle)
calc_bic <- function(i, p1, p2, h1, h2, disp, method='L-BFGS-B') {
    #{{{
    cat(i, "\n", sep = '')
    #{{{ LLs
LL1 <- function(mu1, mu3) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu1, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu3, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL2 <- function(mu1, mu3, mu4) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu1, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL3 <- function(mu1, mu2, mu3) {
    #{{{
    mu4 = (mu2 / mu1) * mu3
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL4 <- function(mu1, mu2, mu3) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu3, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL5 <- function(mu1, mu2, mu3, mu4) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
    #}}}
    size = 1 / disp
    n_obs = length(p1)
    m1 = round(mean(p1)); m2 = round(mean(p2)); m3 = round(mean(h1)); m4 = round(mean(h2))
    p1 = round(p1); p2 = round(p2); h1 = round(h1); h2 = round(h2)
    min_mu = 1e-2; max_mu = 1e8
    fit1 = mle2(LL1, start = list(mu1=(m1+m2)/2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 2), upper=rep(max_mu, 2),
          method = method)#, nobs = n_obs)
    fit2 = mle2(LL2, start = list(mu1=(m1+m2)/2, mu3=m3, mu4=m4),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)#, control = list(trace=3, maxit=1000))
    fit3 = mle2(LL3, start = list(mu1=m1, mu2=m2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)
    fit4 = mle2(LL4, start = list(mu1=m1, mu2=m2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)
    fit5 = mle2(LL5, start = list(mu1=m1, mu2=m2, mu3=m3, mu4=m4),
          lower = rep(min_mu, 4), upper=rep(max_mu, 4),
          method = method)#, nobs = n_obs)
    #coef(fitc)
    bic = AICtab(fit1, fit2, fit3, fit4, fit5, k = log(n_obs), sort=F)
    tb = as_tibble(bic) %>%
        mutate(reg = c('conserved','unexpected','cis','trans','cis+trans')) %>%
        arrange(dAIC)
    tb$reg[1]
    #}}}
}
calc_bic_2 <- function(i, p1,p2,h1,h2,bp1,bp2,bh1,bh2,disp, method='L-BFGS-B') {
    #{{{
    cat(i, "\n", sep = '')
    #{{{ LLs
LL1 <- function(mu1,mu2,mu3,mu4,mu5,mu7) {
    #{{{
    mu6 = mu5 - mu1 + mu2
    mu8 = mu7 - mu3 + mu4
    mu6 = ifelse(mu6 <= 0, 1, mu6)
    mu8 = ifelse(mu8 <= 0, 1, mu8)
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l) | is.infinite(l), 1e6, l)
    #}}}
}
LL2 <- function(mu1,mu2,mu3,mu4,mu5,mu7,mu8) {
    #{{{
    mu6 = mu5 - mu1 + mu2
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l) | is.infinite(l), 1e6, l)
    #}}}
}
LL3 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7) {
    #{{{
    mu8 = ifelse(mu5==mu1, m8, mu4 - (mu3-mu7) * (mu2-mu6) / (mu1-mu5))
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l) | is.infinite(l), 1e6, l)
    #}}}
}
LL4 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7) {
    #{{{
    mu8 = mu7 - mu3 + mu4
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l) | is.infinite(l), 1e6, l)
    #}}}
}
LL5 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l) | is.infinite(l), 1e6, l)
    #}}}
}
    #}}}
    size = 1 / disp
    n_obs = length(p1)
    m1 = round(mean(p1)); m2 = round(mean(p2)); m3 = round(mean(h1)); m4 = round(mean(h2))
    m5 = round(mean(bp1)); m6 = round(mean(bp2)); m7 = round(mean(bh1)); m8 = round(mean(bh2))
    dm1=m1-m5; dm2=m2-m6; dm3=m3-m7; dm4=m4-m8
    p1 = round(p1); p2 = round(p2); h1 = round(h1); h2 = round(h2)
    bp1 = round(bp1); bp2 = round(bp2); bh1 = round(bh1); bh2 = round(bh2)
    min_mu = 1e-2; max_mu = 1e8
    dm12 = (dm1+dm2)/2; dm34 = (dm3+dm4)/2
    m1.f = m5 + dm12; m2.f = m6 + dm12
    m3.f = m7 + dm34; m4.f = m8 + dm34
    #mu1=m1.f;mu2=m2.f;mu3=m3.f;mu4=m4.f;mu5=m5;mu6=m6;mu7=m7;mu8=m8
    fit1 = mle2(LL1, start = list(mu1=m1.f,mu2=m2.f,mu3=m3.f,mu4=m4.f,
                                  mu5=m5,mu7=m7),
          lower = rep(min_mu, 6), upper=rep(max_mu, 6),
          method = method)#, nobs = n_obs)
    fit2 = mle2(LL2, start = list(mu1=m1.f,mu2=m2.f,mu3=m3,mu4=m4,
                                  mu5=m5,mu7=m7,mu8=m8),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    pr1.p = ifelse(dm1+dm2==0, .5, dm1/(dm1+dm2))
    pr1.h = ifelse(dm3+dm4==0, .5, dm3/(dm3+dm4))
    pr2.p=1-pr1.p; pr2.h=1-pr1.h
    mpr1 = (pr1.p+pr1.h)/2; mpr2 = (pr2.p+pr2.h)/2
    m1.f = (dm1+dm2) * mpr1 + m5; m2.f = (dm1+dm2) * mpr2 + m6
    m3.f = (dm3+dm4) * mpr1 + m7; m4.f = (dm3+dm4) * mpr2 + m8
    mu1=m1.f;mu2=m2.f;mu3=m3.f;mu4=m4.f;mu5=m5;mu6=m6;mu7=m7;mu8=m8
    fit3 = mle2(LL3, start = list(mu1=m1.f,mu2=m2.f,mu3=m3.f,mu4=m4.f,
                                  mu5=m5,mu6=m6,mu7=m7),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    m3.f = m7 + dm34; m4.f = m8 + dm34
    fit4 = mle2(LL4, start = list(mu1=m1,mu2=m2,mu3=m3.f,mu4=m4.f,
                                  mu5=m5,mu6=m6,mu7=m7),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    fit5 = mle2(LL5, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,
                                  mu5=m5,mu6=m6,mu7=m7,mu8=m8),
          lower = rep(min_mu, 8), upper=rep(max_mu, 8),
          method = method)#, nobs = n_obs)
    #coef(fitc)
    bic = AICtab(fit1, fit2, fit3, fit4, fit5, k = log(n_obs), sort=F, base=T)
    tb = as_tibble(bic) %>%
        mutate(reg = c('conserved','unexpected','cis','trans','cis+trans')) %>%
        arrange(dAIC)
    tb
    #}}}
}
#}}}

n_cpu = as.integer(args[1])
#{{{ prepare for parallel computing
n_thread = n_cpu
cluster = new_cluster(n_thread)
cluster_library(cluster, "tidyverse")
cluster_library(cluster, "stats4")
cluster_library(cluster, "bbmle")
cluster_copy(cluster, 'calc_bic_2')
#}}}

gid = 'Zm00001d034205'
#i = ti %>% filter(cond=='Cold25',condB=='Control25',cross=='B73xMo17') %>%
    #filter(gid == !!gid) %>% pull(i)
#i=7549
#disp=ti$disp[i]
#p1=ti$rc.p1[[i]]; p2=ti$rc.p2[[i]]; h1=ti$rc.h1[[i]]; h2=ti$rc.h2[[i]];
#bp1=ti$rc0.p1[[i]]; bp2=ti$rc0.p2[[i]]; bh1=ti$rc0.h1[[i]]; bh2=ti$rc0.h2[[i]];

#{{{ stress cis-trans
fi = file.path(dirw, '01.raw.rds')
ra = readRDS(fi)
tmt = ra$tmt

#{{{ prepare & filtering
min_rc = 10
ti = ra$th.cmp %>% inner_join(tmt, by=c('condB'='cond')) %>% select(-trc.h) %>%
    rename(trc0.p1=trc.p1,trc0.p2=trc.p2,trc0.h1=trc.h1,trc0.h2=trc.h2,
           rc0.p1=rc.p1,rc0.p2=rc.p2,rc0.h1=rc.h1,rc0.h2=rc.h2) %>%
    inner_join(tmt, by=c("cond",'cross','gid','disp')) %>%
    filter(trc.p1 + trc.p2 >= 2*min_rc, trc.h1+trc.h2 >= min_rc) %>%
    filter(trc0.p1 + trc0.p2 >= 2*min_rc, trc0.h1+trc0.h2 >= min_rc) %>%
    mutate(n.p1 = map_int(rc.p1, length), n.p2=map_int(rc.p2, length),
           n.h1 = map_int(rc.h1, length), n.h2=map_int(rc.h2, length),
           n0.p1 = map_int(rc0.p1, length), n0.p2=map_int(rc0.p2, length),
           n0.h1 = map_int(rc0.h1, length), n0.h2=map_int(rc0.h2, length)) %>%
    mutate(mrc.p1 = trc.p1/n.p1, mrc.p2 = trc.p2/n.p2,
           mrc.h1 = trc.h1/n.h1, mrc.h2 = trc.h2/n.h2,
           mrc0.p1 = trc0.p1/n0.p1, mrc0.p2 = trc0.p2/n0.p2,
           mrc0.h1 = trc0.h1/n0.h1, mrc0.h2 = trc0.h2/n0.h2) %>%
    mutate(dmrc.p1 = mrc.p1 - mrc0.p1, dmrc.p2 = mrc.p2 - mrc0.p2,
           dmrc.h1 = mrc.h1 - mrc0.h1, dmrc.h2 = mrc.h2 - mrc0.h2) %>%
    filter((dmrc.p1>=0 & dmrc.p2>=0 & dmrc.p1+dmrc.p2>0 &
            dmrc.h1>=0 & dmrc.h2>=0 & dmrc.h1+dmrc.h2>0) |
           (dmrc.p1<=0 & dmrc.p2<=0 & dmrc.p1+dmrc.p2<0 &
            dmrc.h1<=0 & dmrc.h2<=0 & dmrc.h1+dmrc.h2<0)) %>%
    mutate(prop.p=dmrc.p1/(dmrc.p1+dmrc.p2), prop.h=dmrc.h1/(dmrc.h1+dmrc.h2)) %>%
    select(cond,condB,cross,gid,mrc.p1,mrc.p2,mrc.h1,mrc.h2,
           mrc0.p1,mrc0.p2,mrc0.h1,mrc0.h2,
           dmrc.p1,dmrc.p2,dmrc.h1,dmrc.h2,prop.p,prop.h,
           rc.p1,rc.p2,rc.h1,rc.h2,
           rc0.p1,rc0.p2,rc0.h1,rc0.h2,disp) %>%
    mutate(i= 1:n())
ti %>% count(cond,condB,cross)
#}}}

tw = ti %>%# dplyr::slice(101:110) %>%
    partition(cluster) %>%
    mutate(reg = pmap(list(i,rc.p1,rc.p2,rc.h1,rc.h2,
                           rc0.p1,rc0.p2,rc0.h1,rc0.h2,disp), calc_bic_2)) %>%
    collect() %>%
    select(-rc.p1,-rc.p2,-rc.h1,-rc.h2,
        -rc0.p1,-rc0.p2,-rc0.h1,-rc0.h2) %>%
    print(n=40)

fo = file.path(dirw, '05.modes.x.rds')
saveRDS(tw, fo)
#}}}


