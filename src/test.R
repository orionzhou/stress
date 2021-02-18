source('functions.R')
dirw = dird
setwd(dirw)
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
    mu8 = ifelse(m5==m1, mu6-mu2+mu4, (mu6-mu2) / (mu5-mu1) * (mu7-mu3) + mu4)
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
    dm12 = mean(dm1,dm2); dm34 = mean(dm3,dm4)
    m1.f = m5 + dm12; m2.f = m6 + dm12
    m3.f = m7 + dm34; m4.f = m8 + dm34
    fit1 = mle2(LL1, start = list(mu1=m1.f,mu2=m2.f,mu3=m3.f,mu4=m4.f,
                                  mu5=m5,mu7=m7),
          lower = rep(min_mu, 6), upper=rep(max_mu, 6),
          method = method)#, nobs = n_obs)
    fit2 = mle2(LL2, start = list(mu1=m1.f,mu2=m2.f,mu3=m3,mu4=m4,
                                  mu5=m5,mu7=m7,mu8=m8),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    pr1.p = dm1/(dm1+dm2); pr1.h = dm3/(dm3+dm4); pr2.p=1-pr1.p; pr2.h=1-pr1.h
    mpr1 = mean(pr1.p, pr1.h); mpr2 = mean(pr2.p, pr2.h)
    m1.f = (dm1+dm2) * mpr1 + m5; m2.f = (dm1+dm2) * mpr2 + m6
    m3.f = (dm3+dm4) * mpr1 + m7; m4.f = (dm3+dm4) * mpr2 + m8
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
    bic = AICtab(fit1, fit2, fit3, fit4, fit5, k = log(n_obs), sort=F)
    tb = as_tibble(bic) %>%
        mutate(reg = c('conserved','unexpected','cis','trans','cis+trans')) %>%
        arrange(dAIC)
    tb$reg[[1]]
    #}}}
}
#}}}

args = commandArgs(trailingOnly=T)
n_cpu = as.integer(args[1])

cat(n_cpu, "\n")
cat(args, "\n")

diri = '~/projects/cre/data/01_tfbs'
fi = file.path(diri, '05.motifs.meme')
all_mtfs = read_meme(fi)
fi = file.path(diri, '10.fam.meme')
mtfs_known = read_meme(fi)

lid='l1069'
lid='l1070'
lid='l1086'
fi=sprintf("%s/22a_dreme_mtf/%s.dreme", dirr, lid)
mtfs0 = read_meme(fi)
fi=sprintf("%s/22a_dreme_mtf/%s.txt", dirr, lid)
fi=sprintf("%s/dreme_out/dreme.txt", dirr)
mtfs0 = read_meme(fi)

#mtfs1 = c(mtfs0[1], all_mtfs)
mtfs1 = c(mtfs_known[137], mtfs0)
cmp1 = compare_motifs(mtfs1, 1, method="PCC", min.mean.ic=.0,
                      min.overlap=6, score.strat="a.mean")


