source('functions.R')
dirw = file.path(dird, '10_qc')

yid = 'rn20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m


