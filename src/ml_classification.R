#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

ps <- ArgumentParser(description = 'Run machine learning classification on given dataset')
ps$add_argument("fi", nargs=1, help = "input dataset")
ps$add_argument("fo1", nargs=1, help = "output metrics file")
#ps$add_argument("fo2", nargs=1, help = "output best model file")
ps$add_argument("--perm", default=1, type='integer',
                help = "number permutations [default: %(default)s]")
ps$add_argument("--alg", default='rf',
                help = "ML algorithm [default: %(default)s]")
ps$add_argument("--holdout", default=.8, type='double',
                help = "proportion data to hold out for test [default: %(default)s]")
ps$add_argument("--fold", default=5, type='integer',
                help = "cv fold [default: %(default)s]")
ps$add_argument("--fold_repeat", default=1, type='integer',
                help = "repeat in each cv [default: %(default)s]")
ps$add_argument("--nlevel", default=3, type='integer',
                help = "levels of hyperparameters to tune [default: %(default)s]")
ps$add_argument("--downsample", action='store_true',
                help = "downsample to get balanced [default: %(default)s]")
ps$add_argument("--seed", default=26, type='integer',
                help = "random seed [default: %(default)s]")
ps$add_argument("--response", default='status',
                help = "response variable name [default: %(default)s]")
ps$add_argument("--cpu", default=1, type='integer',
                help = "num. processors to use [default: %(default)s]")
args <- ps$parse_args()

require(tidyverse)
require(tidymodels)
require(vip)

#{{{ functions
downsample <- function(ti, seed=1, colname='status') {
  #{{{
  if(colname != 'status')
    ti = ti %>% rename(status = get('colname'))
  levs = levels(ti$status)
  stopifnot(length(levs) == 2)
  tis = ti %>% count(status) %>% arrange(n)
  lev1 = tis %>% pluck('status', 1)
  lev2 = tis %>% pluck('status', 2)
  n1 = tis %>% pluck('n', 1)
  n2 = tis %>% pluck('n', 2)
  ti1 = ti %>% filter(status == lev1)
  set.seed(seed)
  ti2 = ti %>% filter(status == lev2) %>% slice(sample(n2, n1))
  ti1 %>% bind_rows(ti2)
  #}}}
}

fit_split <- function(object, model, split, ...) {
  #{{{
  if (inherits(object, "formula")) {
    object <- add_model(add_formula(workflow(), object, blueprint = hardhat::default_formula_blueprint(indicators = FALSE)), model)
  }
  tune::last_fit(object, split, ...)
  #}}}
}

get_tuner <- function(alg='rf', cpu=1, mode='classification') {
    #{{{
    if(alg == 'svm') {
      svm_poly(cost=tune(), degree=tune(),
               scale_factor=tune(), margin=tune()) %>%
      set_engine("kernlab") %>%
      set_mode(mode)
    } else if (alg == 'xgb') {
      boost_tree(trees=tune(), min_n=tune(),
                 tree_depth=tune(), learn_rate=tune(), loss_reduction=tune()) %>%
      set_engine("xgboost", importance="permutation", nthread=cpu) %>%
      set_mode(mode)
    } else if (alg == 'rf') {
      rand_forest(trees=tune(), min_n=tune()) %>%
      set_engine("ranger", importance="permutation", num.threads=cpu) %>%
      set_mode(mode)
    }
    #}}}
}

get_grid <- function(alg='rf', nlevel=3) {
    #{{{
    if(alg == 'svm') {
       grid_regular(cost(), degree(), scale_factor(), svm_margin(), levels=nlevel)
    } else if (alg == 'xgb') {
       grid_regular(
                    trees(),
                    min_n(),
                    tree_depth(),
                    learn_rate(),
                    loss_reduction(),
                    levels = nlevel)
    } else if (alg == 'rf') {
       grid_regular(
                    #finalize(mtry(trans=sqrt_trans()), data),
                    trees(),
                    min_n(),
                    levels = nlevel)
    }
    #}}}
}

get_model <- function(alg='rf', cpu=1, mode='classification') {
    #{{{
    if(alg == 'svm') {
      svm_poly() %>%
      set_engine("kernlab") %>%
      set_mode(mode)
    } else if (alg == 'xgb') {
      boost_tree() %>%
      set_engine("xgboost", importance="permutation", nthread=cpu) %>%
      set_mode(mode)
    } else if (alg == 'rf') {
      rand_forest() %>%
      set_engine("ranger", importance="permutation", num.threads=cpu) %>%
      set_mode(mode)
    }
    #}}}
}

metrics6 <- metric_set(sens,spec,precision,accuracy, f_meas, roc_auc, pr_auc)

wfl_fit <- function(wfl) wfl %>% pluck(".workflow", 1) %>% pull_workflow_fit()
wfl_metric <- function(wfl) wfl %>% collect_metrics() %>% select(metric=.metric,estimate=.estimate)
wfl_pred <- function(wfl) wfl %>% pluck('.predictions', 1) %>% rename(pred=.pred_class,truth=status)

recipe_ml <- function(train_data) {
  #{{{
  recipe(status ~ ., data = train_data) %>%
    #update_role(sid, new_role = "ID") %>%
    #step_medianimpute(all_numeric(), -all_outcomes()) %>%
    #step_modeimpute(all_nominal(), -all_outcomes()) %>%
    #step_novel(all_nominal()) %>%
    step_zv(all_predictors())
    #step_dummy(all_nominal()) %>%
    #themis::step_downsample(all_outcomes())
  #}}}
}

ml0 <- function(ti, alg='rf', split_prop=.8, cpu=1, seed=26) {
  #{{{
  require(vip)
  set.seed(seed)
  data_split <- initial_split(ti, prop = split_prop)
  train_data <- training(data_split)
  test_data  <- testing(data_split)
  #
  rec = recipe_ml(train_data)
  #
  wfl = workflow() %>%
    add_recipe(rec) %>%
    add_model(get_model(alg=alg, cpu=cpu)) %>%
    fit_split(split = data_split, metrics = metrics6)
  # metrics
  ft = wfl_fit(wfl)
  metric = wfl_metric(wfl)
  pred = wfl_pred(wfl)
  vis = vi(ft) %>% as_tibble()
  #metricB = metric_balanced(pred)
  #
  #tibble(fit=list(ft), vis=list(vis), pred=list(pred), metric=list(metric), metricB=list(metricB))
  tibble(fit=list(ft), vis=list(vis), pred=list(pred), metric=list(metric))
  #}}}
}

ml1 <- function(ti, alg='rf', split_prop=.8, fold=10, fold_repeat=1,
                nlevel=3, cpu=1, seed=26) {
  #{{{
  require(vip)
  set.seed(seed)
  #cat("seed =", seed, "\n")
  data_split <- initial_split(ti, prop = split_prop, strata = status)
  train_data <- training(data_split)
  test_data  <- testing(data_split)
  #
  rec = recipe_ml(train_data)
  #
  wfl = workflow() %>%
    add_recipe(rec) %>%
    add_model(get_tuner(alg=alg, cpu=cpu))
  folds = vfold_cv(train_data, v = fold, repeats=fold_repeat)
  mfit = wfl %>%
      tune_grid(resamples = folds,
                grid = get_grid(alg=alg, nlevel=nlevel),
                metrics = metric_set(f_meas))
  mres = mfit %>% collect_metrics(summarize = F)
  mfit %>% show_best(metric='f_meas') %>% print(width=Inf)
  param = mfit %>% select_best(metric='f_meas')
  wfl = wfl %>% finalize_workflow(param) %>%
      fit_split(split = data_split, metrics = metrics6)
  # metrics
  ft = wfl_fit(wfl)
  metric = wfl_metric(wfl)
  pred = wfl_pred(wfl)
  vis = tibble()
  if(alg %in% c("rf",'xgb')) vis = vi(ft) %>% as_tibble()
  #
  tibble(fit=list(ft), param=list(param), vis=list(vis), pred=list(pred), metric=list(metric), mres=list(mres))
  #}}}
}

#}}}

ti = read_tsv(args$fi)
if('gid' %in% colnames(ti)) ti = ti %>% select(-gid)
if(args$response != 'status') ti = ti %>% rename(status=args$response)
ti = ti %>% mutate(status = factor(status, levels=c(1,0)))
to0 = tibble(perm = 1:args$perm, ti0 = list(ti))
if(args$downsample) {
    to = to0 %>% mutate(ti = map2(ti0, perm, downsample)) %>% select(-ti0)
} else {
    to = to0 %>% rename(ti=ti0)
}

to = to %>%
    mutate(r = map(ti, ml1, alg=args$alg, split_prop=args$holdout,
        fold=args$fold, fold_repeat=args$fold_repeat,
        nlevel=args$nlevel,
        cpu=args$cpu, seed=args$seed)) %>%
    select(-ti) %>% unnest(r)

tb = to %>% select(perm, metric) %>% unnest(metric) %>%
    filter(metric=='f_meas') %>% arrange(desc(estimate))
cat("best model F1-score: ",tb %>% pluck('estimate', 1), "\n")
perm_best = tb %>% pluck('perm', 1)


res = to %>% mutate(fit = ifelse(perm==perm_best, fit, NA))
saveRDS(res, args$fo1)
#saveRDS(r2, args$fo2)

