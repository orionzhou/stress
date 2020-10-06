require(tidyverse)
require(tidymodels)
require(skimr)
require(themis)
require(vip)
fit_split <- function(object, model, split, ...) {
  #{{{
  if (inherits(object, "formula")) {
    object <- add_model(add_formula(workflow(), object, blueprint = hardhat::default_formula_blueprint(indicators = FALSE)), model)
  }
  tune::last_fit(object, split, ...)
  #}}}
}

dirp = "~/projects/ml"
dird = file.path(dirp, 'data')
dirw = file.path(dird, '03_workshop')

fi = file.path(dirw, 'data.txt')
ti = read_tsv(fi) %>% rename(id=X1)

features = c(
"Transferase",
"p450",
"UDPGT",
"tandemDupGenes",
"Nicotiana_tabacum.TN90_AYMY.SS.csv",
"Ppatens_318_v3.3.csv",
"Atrichopoda_291_v1.0.csv",
"Coffea_canephora.csv",
"Nicotiana_tomen.csv",
"BrapaFPsc_277_v1.3.csv",
"FamilySize_cat"
)
ti2 = ti %>%
    select(-id) %>%
    #select(any_of(c('Class',features))) %>%
    filter(Class != 'unknown') %>%
    mutate(Class=factor(Class))

set.seed(555)
data_split <- initial_split(ti2, prop = .9)
train_data <- training(data_split)
test_data  <- testing(data_split)

rec = recipe(Class ~ ., data = train_data) %>%
  #update_role(id, new_role = "ID") %>%
  #step_novel(all_nominal()) %>%
  #step_knnimpute(all_predictors(), neighbors = 3) %>%
  step_medianimpute(all_numeric(), -all_outcomes()) %>%
  step_modeimpute(all_nominal(), -all_outcomes()) %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_zv(all_predictors()) %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 10) %>%
  step_downsample(all_outcomes())
summary(rec)

wfl = workflow() %>%
  add_recipe(rec) %>%
  add_model(xgb_tuner)
  #add_model(rf_tuner)
  #add_model(svm_tuner)
folds = vfold_cv(train_data, v = 10)
mfit = wfl %>%
    tune_grid(resamples = folds,
        metrics = metric_set(accuracy, precision, f_meas, roc_auc, pr_auc))
mres = mfit %>% collect_metrics(summarize = F)

mfit %>% show_best(metric='f_meas') %>% print(width=Inf)
mfit1 = mfit %>% select_best(metric='f_meas')
mfit1
wfl1 = wfl %>% finalize_workflow(mfit1)
lfit = wfl1 %>% fit_split(split = data_split,
        metrics = metric_set(accuracy, precision, f_meas, roc_auc, pr_auc))
lfit %>% collect_metrics()

fo = file.path(dirw, 'vip.pdf')
pdf(fo, width=6, height=7)
lmod = lfit %>% pluck(".workflow", 1)
vip(lmod, num_features=20)
dev.off()


