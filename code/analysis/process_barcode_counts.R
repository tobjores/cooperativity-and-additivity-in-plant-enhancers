library(readr)
library(dplyr)
library(tidyr)
library(tibble)


### load and preprocess experiment data (barcode count files) ###
# input:        character or connection; data file with input read counts
# output:       named vector with conditon and output read count file; <condition 1> = <data file 1>, <condition 2> = <data file 2>, ...
# subassembly:  R object; subassembly data frame
# rcc:          numeric; read count cutoff to be applied to the input and output counts
load_experiment_bc <- function(input, output, subassembly, rcc = 1) {
  # load and merge files (only keep barcodes detected in both input and output)
  data_inp <- read_table(input, col_names = c('count', 'barcode'), col_types = 'ic') |>
    filter(count >= rcc)
  
  data_out <- tibble(condition = names(output), file = unlist(output)) |>
    nest_by(condition) |>
    reframe(read_table(data$file, col_names = c('count', 'barcode'), col_types = 'ic')) |>
    filter(count >= rcc) |>
    inner_join(data_inp, by = 'barcode', suffix = c('_out', '_inp'))
  
  # merge with subassembly
  data_out <- inner_join(data_out, subassembly, by = 'barcode', relationship = 'many-to-many')
  
  # return result
  return(data_out)
}


### PEfl library ###
## load subassembly
subassembly_PEfl <- read_tsv('data/subassembly/PEfl/subassembly_pPSup_PEfl.tsv.gz')

## load experiment data
data_PEfl <- expand_grid(rep = c(1, 2), conditions = c('dark', 'light')) |>
  group_by(rep) |>
  summarise(
    conditions = list(setNames(paste0('data/barcode_counts/PEfl/rep', rep,'/barcode_pPSup_PEfl_', as.roman(rep), '_', conditions, '.count.gz'), conditions))
  ) |>
  nest_by(rep) |>
  reframe(
    load_experiment_bc(
      input = paste0('data/barcode_counts/PEfl/rep', rep, '/barcode_pPSup_PEfl_', as.roman(rep), '_input.count.gz'),
      output = unlist(data$conditions),
      subassembly = subassembly_PEfl,
      rcc = 5
    )
  )

## remove barcodes detected in only one replicate or condition
data_filtered_PEfl <- data_PEfl |>
  group_by(barcode) |>
  filter(n() == 4) |>
  ungroup()

## calculate enrichment and normalize to library median
data_norm_PEfl <- data_filtered_PEfl |>
  group_by(condition, rep) |>
  mutate(
    enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup()

## save data
save(data_norm_PEfl, file = 'data/RData/PEfl_data_main.Rdata')


### PEV library ###
## load subassembly
subassembly_PEV <- read_tsv('data/subassembly/PEV/subassembly_pPSup_PEV.tsv.gz')

## load experiment data
# light and dark (ld)
data_PEV_ld <- expand_grid(rep = c(1, 2, 3), conditions = c('dark', 'light')) |>
  group_by(rep) |>
  summarise(
    conditions = list(setNames(paste0('data/barcode_counts/PEV/light_dark/rep', rep,'/barcode_pPSup_PEV_', as.roman(rep), '_', conditions, '.count.gz'), conditions))
  ) |>
  nest_by(rep) |>
  reframe(
    load_experiment_bc(
      input = paste0('data/barcode_counts/PEV/light_dark/rep', rep, '/barcode_pPSup_PEV_', as.roman(rep), '_input.count.gz'),
      output = unlist(data$conditions),
      subassembly = subassembly_PEV,
      rcc = 5
    )
  )

# circadian rhythm (cr)
data_PEV_cr <- expand_grid(rep = c(1, 2), conditions = c(0, 6, 12, 18, 24)) |>
  group_by(rep) |>
  summarise(
    conditions = list(setNames(paste0('data/barcode_counts/PEV/circadian_rhythm/rep', rep,'/barcode_pPSup_PEV_', as.roman(rep), '_CR-', conditions, '.count.gz'), conditions))
  ) |>
  nest_by(rep) |>
  reframe(
    load_experiment_bc(
      input = paste0('data/barcode_counts/PEV/circadian_rhythm/rep', rep, '/barcode_pPSup_PEV_', as.roman(rep), '_CR-input.count.gz'),
      output = unlist(data$conditions),
      subassembly = subassembly_PEV,
      rcc = 5
    )
  ) |>
  mutate(
    condition = as.numeric(condition)
  )

## remove barcodes detected in only one replicate
data_filtered_PEV_ld <- data_PEV_ld |>
  mutate(
    duplSeq = replace_na(duplSeq, FALSE)
  ) |>
  group_by(barcode, condition) |>
  filter(sum(! duplSeq) > 1) |>
  ungroup()

data_filtered_PEV_cr <- data_PEV_cr |>
  mutate(
    duplSeq = replace_na(duplSeq, FALSE)
  ) |>
  group_by(barcode, condition) |>
  filter(sum(! duplSeq) > 1) |>
  ungroup()

## aggregate barcodes (sum of individual counts)
data_ag_PEV_ld <- data_filtered_PEV_ld |>
  group_by(across(-c(count_inp, count_out, barcode))) |>
  summarise(
    n_bc = n(),
    count_inp = sum(count_inp),
    count_out = sum(count_out)
  ) |>
  ungroup()

data_ag_PEV_cr <- data_filtered_PEV_cr |>
  group_by(across(-c(count_inp, count_out, barcode))) |>
  summarise(
    n_bc = n(),
    count_inp = sum(count_inp),
    count_out = sum(count_out)
  ) |>
  ungroup()

## calculate enrichment and normalize to library median
data_norm_PEV_ld <- data_ag_PEV_ld |>
  group_by(condition, rep) |>
  mutate(
    enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup()

data_norm_PEV_cr <- data_ag_PEV_cr |>
  group_by(condition, rep) |>
  mutate(
    enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup()

## remove duplicated sequences (for replicate correlations)
data_reps_PEV_ld <- data_norm_PEV_ld |>
  filter(! duplSeq)

data_reps_PEV_cr <- data_norm_PEV_cr |>
  filter(! duplSeq)

## calculate mean enrichment across replicates
data_mean_PEV_ld <- data_norm_PEV_ld |>
  group_by(across(-c(rep, n_bc, count_inp, count_out, enrichment))) |>
  summarise(
    n_experiments = n(),
    min_bc = min(n_bc),
    min_ci = min(count_inp),
    min_co = min(count_out),
    enrichment = mean(enrichment)
  ) |>
  ungroup()

data_mean_PEV_cr <- data_norm_PEV_cr |>
  group_by(across(-c(rep, n_bc, count_inp, count_out, enrichment))) |>
  summarise(
    n_experiments = n(),
    min_bc = min(n_bc),
    min_ci = min(count_inp),
    min_co = min(count_out),
    enrichment = mean(enrichment)
  ) |>
  ungroup()

## save data
save(data_reps_PEV_ld, file = 'data/RData/PEV_ld_data_reps.Rdata')
save(data_mean_PEV_ld, file = 'data/RData/PEV_ld_data_main.Rdata')

save(data_reps_PEV_cr, file = 'data/RData/PEV_cr_data_reps.Rdata')
save(data_mean_PEV_cr, file = 'data/RData/PEV_cr_data_main.Rdata')


### PEF library ###
## load subassembly
subassembly_PEF <- read_tsv('data/subassembly/PEF/subassembly_pPSup_PEF.tsv.gz', col_types = 'ccccccci')

## load experiment data
data_PEF <- expand_grid(rep = c(1, 2, 3), conditions = c('dark', 'light')) |>
  group_by(rep) |>
  summarise(
    conditions = list(setNames(paste0('data/barcode_counts/PEF/rep', rep,'/barcode_pPSup_PEF_', as.roman(rep), '_', conditions, '.count.gz'), conditions))
  ) |>
  nest_by(rep) |>
  reframe(
    load_experiment_bc(
      input = paste0('data/barcode_counts/PEF/rep', rep, '/barcode_pPSup_PEF_', as.roman(rep), '_input.count.gz'),
      output = unlist(data$conditions),
      subassembly = subassembly_PEF,
      rcc = 5
    )
  )

## remove barcodes detected in only one replicate
data_filtered_PEF <- data_PEF |>
  group_by(barcode, condition) |>
  filter(n() > 1) |>
  ungroup()

## aggregate barcodes (sum of individual counts)
data_ag_PEF <- data_filtered_PEF |>
  group_by(across(-c(count_inp, count_out, barcode))) |>
  summarise(
    n_bc = n(),
    count_inp = sum(count_inp),
    count_out = sum(count_out)
  ) |>
  ungroup()

## calculate enrichment and normalize to library median
data_reps_PEF <- data_ag_PEF |>
  group_by(condition, rep) |>
  mutate(
    enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup()

## calculate mean enrichment across replicates and normalize to control construct without enhancer
data_mean_PEF <- data_reps_PEF |>
  group_by(across(-c(rep, n_bc, count_inp, count_out, enrichment))) |>
  summarise(
    n_experiments = n(),
    min_bc = min(n_bc),
    min_ci = min(count_inp),
    min_co = min(count_out),
    enrichment = mean(enrichment)
  ) |>
  group_by(condition) |>
  mutate(
    enrichment = enrichment - enrichment[type == 'noEnh-control']
  ) |>
  ungroup()

## save data
save(data_reps_PEF, file = 'data/RData/PEF_data_reps.Rdata')
save(data_mean_PEF, file = 'data/RData/PEF_data_main.Rdata')


### PEVdouble and PEFval library (aka PEval library) ###
# the two libraries were pooled for the Plant STARR-seq experiments

## load subassembly
subassembly_PEval <- bind_rows(
    'PEVdouble' = read_tsv('data/subassembly/PEVdouble/subassembly_pPSup_PEVdouble.tsv.gz'),
    'PEFval' = read_tsv('data/subassembly/PEFval/subassembly_pPSup_PEFval.tsv.gz', col_types = 'ccccci'),
    .id = 'library'
  ) |>
  filter(! (library == 'PEFval' & grepl('control', type, fixed = TRUE)))

# remove duplicated barcodes
subassembly_PEval <- subassembly_PEval |>
  filter(! barcode %in% barcode[duplicated(barcode)])

## load experiment data
data_PEval <- expand_grid(rep = c(1, 2), conditions = c('dark', 'light')) |>
  group_by(rep) |>
  summarise(
    conditions = list(setNames(paste0('data/barcode_counts/PEVdouble+PEFval/rep', rep,'/barcode_pPSup_PEVdouble+PEFval_', as.roman(rep), '_', conditions, '.count.gz'), conditions))
  ) |>
  nest_by(rep) |>
  reframe(
    load_experiment_bc(
      input = paste0('data/barcode_counts/PEVdouble+PEFval/rep', rep, '/barcode_pPSup_PEVdouble+PEFval_', as.roman(rep), '_input.count.gz'),
      output = unlist(data$conditions),
      subassembly = subassembly_PEval,
      rcc = 5
    )
  )

## remove barcodes detected in only one replicate
data_filtered_PEval <- data_PEval |>
  group_by(barcode, condition) |>
  filter(n() > 1) |>
  ungroup()

## aggregate barcodes (sum of individual counts)
data_ag_PEval <- data_filtered_PEval |>
  group_by(across(-c(count_inp, count_out, barcode))) |>
  summarise(
    n_bc = n(),
    count_inp = sum(count_inp),
    count_out = sum(count_out)
  ) |>
  ungroup()

## calculate enrichment and normalize to library median
data_norm_PEval <- data_ag_PEval |>
  group_by(condition, rep) |>
  mutate(
    variant = replace_na(variant, 'fragments'),
    enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup()

## calculate mean enrichment across replicates and normalize to control construct without enhancer
data_mean_PEval <- data_norm_PEval |>
  group_by(across(-c(rep, n_bc, count_inp, count_out, enrichment))) |>
  summarise(
    n_experiments = n(),
    min_bc = min(n_bc),
    min_ci = min(count_inp),
    min_co = min(count_out),
    enrichment = mean(enrichment)
  ) |>
  group_by(condition) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'noEnh-control']
  ) |>
  ungroup()

## select relevant columns and split by library
data_reps_PEval <- data_norm_PEval |>
  unite(
    col = 'fragments',
    starts_with('fragment'),
    na.rm = TRUE
  ) |>
  mutate(
    variant = if_else(library == 'PEFval', fragments, variant)
  ) |>
  select(library, rep, condition, enhancer, type, variant, orientation, enrichment)

data_mean_PEVdouble <- data_mean_PEval |>
  filter(library == 'PEVdouble') |>
  select(-library, -where(~ all(is.na(.x))))

data_mean_PEFval <- data_mean_PEval |>
  filter(library == 'PEFval' | type %in% c('control', 'WT')) |>
  mutate(
    type = variant
  ) |>
  select(-library, -where(~ all(is.na(.x))), -variant)

## save data
save(data_reps_PEval, file = 'data/RData/PEval_data_reps.Rdata')
save(data_mean_PEVdouble, file = 'data/RData/PEVdouble_data_main.Rdata')
save(data_mean_PEFval, file = 'data/RData/PEFval_data_main.Rdata')
