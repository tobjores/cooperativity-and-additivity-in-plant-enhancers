library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(RcppRoll)
library(hexbin)
library(universalmotif)
library(openxlsx)


### load export functions & define output folder ###
source('code/analysis/export_functions.R')

export_dir <- 'figures/rawData/'
supp_data_dir <- 'data/supplemental_data/'

paper_ref <- '' # update once we have a DOI (e.g. "Plant Cell (2015). 10.1105/tpc.18.00001.")
paper_id <- str_split_1(paper_ref, fixed('.')) |> tail(2) |> head(1)


### load data ###
sapply(dir(path = 'data/RData/', pattern = '*.Rdata', full.names = TRUE), load, environment())

enhancer_sequences <- read_tsv('data/sequence_files/plant-enhancer_sequences.tsv')
enhancer_variants <- read_tsv('data/sequence_files/plant-enhancer_variants.tsv')
enhancer_fragments <- read_tsv('data/sequence_files/enhancer-fragments.tsv')

TF_motifs <- read_meme('data/miscellaneous/TF-clusters.meme')


### PEfl library ###
## replicate correlations
# transform data
data_wide_PEfl <- data_norm_PEfl |>
  select(rep, condition, barcode, enhancer, enrichment) |>
  pivot_wider(
    names_from = rep,
    names_prefix = 'rep',
    values_from = enrichment
  ) |>
  drop_na(starts_with('rep'))

# export plot data
for (cond in unique(data_wide_PEfl$condition)) {
  data_export <- data_wide_PEfl |>
    filter(condition == cond) |>
    select(starts_with('rep'), enhancer) |>
    slice_sample(prop = 1)
  
  write_tsv(data_export, paste0(export_dir, 'PEfl_cor_reps_', cond, '_points.tsv'))
  
  LaTeX_cor_stats(
    x = data_export$rep1,
    y = data_export$rep2,
    file = paste0(export_dir, 'PEfl_cor_reps_', cond)
  )
}

# export axis limits
data_wide_PEfl |>
  summarise(
    xmin = min(rep1, rep2),
    xmax = max(rep1, rep2),
    ymin = xmin,
    ymax = xmax
  ) |>
  write_tsv(paste0(export_dir, 'PEfl_cor_reps_axes.tsv'))


## enhancer strength (fwd orientation only)
# normalize to control construct without enhancer
data_norm_noEnh_PEfl <- data_norm_PEfl |>
  group_by(condition) |>
  mutate(
    enrichment = enrichment - median(enrichment[enhancer == 'none'])
  ) |>
  ungroup()

# order samples
sample_order <- expand_grid(c('FL', 'A', 'B'), c('light', 'dark')) |>
  unite(
    col = 'sample',
    everything(),
    sep = '.'
  ) |>
  pull()

data_ordered_PEfl_fwd <- data_norm_noEnh_PEfl |>
  filter(orientation %in% c('fwd', NA)) |>
  unite(
    col = 'sample',
    c(part, condition),
    sep = '.',
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    sample = ordered(sample, levels = sample_order)
  ) |>
  drop_na(sample)

# export plot data
for (enh in unique(data_ordered_PEfl_fwd$enhancer)) {
  data_export <- data_ordered_PEfl_fwd |>
    filter(enhancer == enh)
  
  LaTeX_boxplot(
    data = data_export,
    samples_from = sample,
    values_from = enrichment,
    file = paste0(export_dir, 'PEfl_strength_fwd_', enh),
    part,
    outliers = TRUE,
    p_values = FALSE
  )
}

# export axis limits
data_ordered_PEfl_fwd |>
  rename(y = enrichment) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEfl_strength_fwd_axes.tsv'))


## enhancer strength by orientation
## order samples
sample_order <- expand_grid(c('FL', 'A', 'B'), c('fwd', 'rev')) |>
  unite(
    col = 'sample',
    everything(),
    sep = '.'
  ) |>
  pull()

data_ordered_PEfl_ori <- data_norm_noEnh_PEfl |>
  unite(
    col = 'sample',
    c(part, orientation),
    sep = '.',
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    sample = ordered(sample, levels = sample_order)
  ) |>
  drop_na(sample)

# export plot data
enh_cond <- data_ordered_PEfl_fwd |>
  distinct(enhancer, condition)

for (i in seq_len(nrow(enh_cond))) {
  enh <- enh_cond[i, 'enhancer', drop = TRUE]
  cond <- enh_cond[i, 'condition', drop = TRUE]
  
  data_export <- data_ordered_PEfl_ori |>
    filter(enhancer == enh & condition == cond)
  
  LaTeX_boxplot(
    data = data_export,
    samples_from = sample,
    values_from = enrichment,
    file = paste0(export_dir, 'PEfl_strength_', cond, '_', enh),
    part,
    outliers = TRUE,
    p_values = FALSE
  )
}

# export axis limits
for (cond in unique(data_ordered_PEfl_ori$condition)) {
  data_ordered_PEfl_ori |>
    filter(condition == cond) |>
    rename(y = enrichment) |>
    summarise(
      across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
      yhelper = 0.05 * (ymax - ymin)
    ) |>
    ungroup() |>
    mutate(
      ymin = ymin - 3 * yhelper,
      ymax = ymax + yhelper
    ) |>
    select(-yhelper) |>
    write_tsv(paste0(export_dir, 'PEfl_strength_', cond, '_axes.tsv')) 
}


## calculate light-responsiveness for each barcode
data_lightResp_PEfl <- data_norm_noEnh_PEfl |>
  select(rep, condition, barcode, enhancer, part, orientation, enrichment) |>
  group_by(rep, barcode, enhancer, part, orientation) |>
  summarise(
    lightResp = enrichment[condition == 'light'] - enrichment[condition == 'dark']
  ) |>
  ungroup()

## light-responsiveness (fwd orientation only)
# order samples
sample_order <- expand_grid(enhancer = c('35S', 'AB80', 'Cab-1', 'rbcS-E9'), part = c('FL', 'A', 'B')) |>
  filter(enhancer != '35S' | part == 'FL') |>
  unite(
    col = 'sample',
    everything(),
    sep = '.'
  ) |> 
  pull()

data_lightResp_PEfl_fwd <- data_lightResp_PEfl |>
  filter(orientation == 'fwd') |>
  unite(
    col = 'sample',
    c(enhancer, part),
    sep = '.',
    remove = FALSE
  ) |>
  mutate(
    sample = ordered(sample, levels = sample_order)
  ) |>
  drop_na(sample)

# export plot data
LaTeX_boxplot(
  data = data_lightResp_PEfl_fwd,
  samples_from = sample,
  values_from = lightResp,
  file = paste0(export_dir, 'PEfl_lightResp_fwd'),
  enhancer,
  part,
  outliers = TRUE,
  p_values = FALSE
)

# export axis limits
data_lightResp_PEfl_fwd |>
  rename(y = lightResp) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEfl_lightResp_fwd_axes.tsv'))


## light-responsiveness by orientation
# order samples
sample_order <- expand_grid(c('FL', 'A', 'B'), c('fwd', 'rev')) |>
  unite(
    col = 'sample',
    everything(),
    sep = '.'
  ) |>
  pull()

data_lightResp_PEfl_ori <- data_lightResp_PEfl |>
  unite(
    col = 'sample',
    c(part, orientation),
    sep = '.',
    remove = FALSE
  ) |>
  mutate(
    sample = ordered(sample, levels = sample_order)
  ) |>
  drop_na(sample)

# export plot data
for (enh in unique(data_lightResp_PEfl_ori$enhancer)) {
  data_export <- data_lightResp_PEfl_ori |>
    filter(enhancer == enh)
  
  LaTeX_boxplot(
    data = data_export,
    samples_from = sample,
    values_from = lightResp,
    file = paste0(export_dir, 'PEfl_lightResp_', enh),
    enhancer,
    part,
    outliers = TRUE,
    p_values = FALSE
  ) 
}

# export axis limits
data_lightResp_PEfl_ori |>
  rename(y = lightResp) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEfl_lightResp_axes.tsv'))


### PEV library - light/dark ###
## replicate correlations
# calculate axis limits
axis_limits <- data_reps_PEV_ld |>
  summarise(
    xmin = min(enrichment),
    xmax = max(enrichment),
    ymin = xmin,
    ymax = xmax
  )

# get all possible combinations
rep_combos <- expand_grid(rep1 = unique(data_reps_PEV_ld$rep), rep2 = unique(data_reps_PEV_ld$rep)) |>
  filter(rep2 > rep1)

max_count <- 0

for (i in seq_len(nrow(rep_combos))) {
  rep1 <- rep_combos[i, 'rep1', drop = TRUE]
  rep2 <- rep_combos[i, 'rep2', drop = TRUE]
  
  # transform data
  data_wide_PEV_ld <- data_reps_PEV_ld |>
    mutate(
      rep = case_when(
        rep == rep1 ~ 'rep1',
        rep == rep2 ~ 'rep2'
      )
    ) |>
    drop_na(rep) |>
    select(rep, condition, enhancer, part, orientation, variant, enrichment) |>
    pivot_wider(
      names_from = rep,
      values_from = enrichment
    ) |>
    drop_na(starts_with('rep'))
  
  # export plot data
  for (cond in unique(data_wide_PEV_ld$condition)) {
    data_export <- data_wide_PEV_ld |>
      filter(condition == cond) |>
      select(starts_with('rep'))
    
    this_count <- LaTeX_hexbin(
      data = data_export,
      x_values = rep1,
      y_values = rep2,
      file = paste0(export_dir, 'PEV_ld_cor_reps_', rep1, 'vs', rep2, '_', cond),
      xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
      axis_limits = FALSE
    )
    
    max_count <- max(max_count, this_count)
  }
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEV_ld_cor_reps_axes.tsv')) 


## heatmaps
# get WT nucleotides
WT_sequences <- enhancer_sequences |>
  group_by(enhancer) |>
  mutate(
    offset = str_locate(sequence[part == 'FL'], sequence)[, 1] - 1
  ) |>
  ungroup() |>
  filter(part != 'FL') |>
  nest_by(enhancer, part) |>
  reframe(
    varNuc = str_split(data$sequence, ''),
    position = list(data$offset + seq_along(varNuc))
  ) |>
  unnest(everything()) |>
  ungroup() |>
  mutate(
    varNuc = ordered(varNuc, levels = c('X', 'A', 'C', 'G', 'T')),
    type = 'WT',
    condition = 'light/dark',
    enrichment = 0
  ) |>
  separate_longer_delim(
    condition,
    delim = '/'
  )

# get all possible variants
all_variants <- enhancer_variants |>
  select(-sequence) |>
  mutate(
    condition = list(unique(data_mean_PEV_ld$condition))
  ) |>
  unnest_longer(condition)

# merge data with all possible variants and normalize enrichment to WT
data_heatmap_PEV_ld <- data_mean_PEV_ld |>
  filter(type != 'control') |>
  group_by(condition, enhancer, part, orientation) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'WT']
  ) |>
  group_by(condition, enhancer, part, type, variant, position, varNuc) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  ungroup() |>
  right_join(all_variants) |>
  filter(type != 'WT') |>
  bind_rows(WT_sequences) |>
  mutate(
    varNuc = ordered(varNuc, levels = c('X', 'A', 'C', 'G', 'T'))
  )

# export heatmap data
enh_part_cond <- data_heatmap_PEV_ld |>
  distinct(enhancer, part, condition)

for (i in seq_len(nrow(enh_part_cond))) {
  enh = enh_part_cond[i, 'enhancer', drop = TRUE]
  pt = enh_part_cond[i, 'part', drop = TRUE]
  cond = enh_part_cond[i, 'condition', drop = TRUE]
  
  # filter data
  data_export <- data_heatmap_PEV_ld |>
    filter(enhancer == enh & part == pt & condition == cond)
  
  # export insertion variants
  data_export |>
    filter(type == 'insertion') |>
    select(position, varNuc, enrichment) |>
    arrange(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_ld_heatmap_', enh, '_', pt, '_', cond, '_ins.tsv'), na = 'NaN')
  
  # export substitution and deletion variants
  data_export |>
    filter(type %in% c('WT', 'substitution', 'deletion')) |>
    select(position, varNuc, enrichment) |>
    arrange(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_ld_heatmap_', enh, '_', pt, '_', cond, '_sub+del.tsv'), na = 'NaN')
}

# export WT positions and axis limits
for (enh in unique(data_heatmap_PEV_ld$enhancer)) {
  # WT positions
  WT_sequences |>
    filter(enhancer == enh) |>
    distinct(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_ld_heatmap_', enh, '_WT.tsv'), na = 'X')
  
  # axis limits
  data_heatmap_PEV_ld |>
    filter(enhancer == enh) |>
    drop_na(enrichment) |>
    summarise(
      `point meta max` = max(abs(enrichment)),
      `point meta min` = -`point meta max`
    ) |>
    write_tsv(paste0(export_dir, 'PEV_ld_heatmap_', enh, '_axes.tsv'))
}

# export summary data
data_summary_PEV_ld <- data_heatmap_PEV_ld |>
  filter(variant != 'WT' & ! duplSeq) |>
  drop_na(enrichment) |>
  mutate(
    effect = case_when(
      enrichment > 1 ~ 'increase',
      enrichment < -1 ~ 'decrease',
      .default = 'neutral'
    )
  )

data_effects_PEV_ld <- data_summary_PEV_ld |>
  count(condition, enhancer, part, effect) |>
  pivot_wider(
    names_from = effect,
    values_from = n
  ) |>
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  )

data_range <- data_summary_PEV_ld |>
  pull(enrichment) |>
  range()

axis_limits <- tibble()

for (pt in unique(data_summary_PEV_ld$part)) {
  data_part <- data_summary_PEV_ld |>
    filter(part == pt)
  
  enh_cond <- data_part |>
    distinct(enhancer, condition)
  
  for (i in seq_len(nrow(enh_cond))) {
    enh <- enh_cond[i, 'enhancer', drop = TRUE]
    cond <- enh_cond[i, 'condition', drop = TRUE]
    
    data_export <- data_part |>
      filter(enhancer == enh & condition == cond)
    
    axis_limits <- axis_limits |>
      bind_rows(
        LaTeX_histogram(
          data = data_export,
          values_from = enrichment,
          file = paste0(export_dir, 'PEV_ld_summary_', enh, '_', pt, '_', cond),
          data_range = data_range,
          bins = 30,
          value_axis = 'y',
          save_axis_limits = FALSE
        ) |>
          bind_cols(list(condition = cond, part = pt))
      )
    
    data_effects_PEV_ld |>
      filter(enhancer == enh & part == pt & condition == cond) |>
      write_tsv(paste0(export_dir, 'PEV_ld_summary_', enh, '_', pt, '_', cond, '_stats.tsv'))
  }
}

axis_limits <- axis_limits |>
  group_by(condition, part) |>
  summarise(
    across(everything(), max)
  ) |>
  ungroup()

pt_cond <- axis_limits |>
  distinct(part, condition)

for (i in seq_len(nrow(pt_cond))) {
  pt <- pt_cond[i, 'part', drop = TRUE]
  cond <- pt_cond[i, 'condition', drop = TRUE]
  
  axis_limits |>
    filter(part == pt, condition == cond) |>
    select(-c(part, condition)) |>
    write_tsv(paste0(export_dir, 'PEV_ld_summary_', pt, '_', cond, '_axes.tsv')) 
}


## create supplemental data set with single-nucleotide variant data in light and dark
# prepare data
WT_nucs <- WT_sequences |>
  unite(
    col = 'enh_pos',
    enhancer,
    position
  ) |>
  distinct(enh_pos, varNuc) |>
  mutate(
    varNuc = as.character(varNuc)
  ) |>
  deframe()

variant_sequences <- enhancer_variants |>
  filter(variant != 'WT') |>
  group_by(sequence) |>
  arrange(position) |>
  mutate(
    first_var = first(variant),
    unique = case_when(
      n() == 1 ~ 'yes',
      variant == first_var & n() == 2 ~ 'no; duplicated once',
      variant == first_var ~ paste0('no; duplicated ', n() - 1, ' times'),
      .default = paste0('no; same as ', first_var)
    )
  ) |>
  ungroup() |>
  select(enhancer, part, variant, unique, sequence)

data_table <- data_heatmap_PEV_ld |>
  filter(type != 'WT') |>
  mutate(
    segment = if_else(part == 'A', "5'", "3'"),
    wtNuc = WT_nucs[paste0(enhancer, '_', position)],
    wtNuc = replace_na(wtNuc, 'X'),
    across(ends_with('Nuc'), ~ if_else(.x == 'X', '-', .x))
  ) |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment
  ) |>
  full_join(
    variant_sequences
  ) |>
  arrange(enhancer, part, position, varNuc) |>
  select(enhancer, segment, type, variant, position, 'wild-type base' = wtNuc, 'variant base' = varNuc, unique, dark, light, sequence)

# export data to an excel file
DS_id <- 1

wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'single-nucleotide variants (ld)')

writeData(
  wb,
  sheet = 1,
  paste('Supplemental Data. Jores et al.', paper_ref),
  startCol = 1,
  startRow = 1
)

writeData(
  wb,
  sheet = 1,
  paste0('Supplemental Data Set ', DS_id,' Enhancer strength of single-nucleotide variants of the AB80, Cab-1, and rbcS-E9 enhancers in the light or dark.'),
  startCol = 1,
  startRow = 3
)

makeNameBold(wb, DS_id)
mergeCells(wb, sheet = 1, rows = 3, cols = seq_len(ncol(data_table)))

writeData(
  wb,
  sheet = 1,
  paste(
    "All possible single-nucleotide substitution, deletion, and insertion variants of the 5' and 3' segments of the AB80, Cab-1, and rbcS-E9",
    "enhancers were subjected to Plant STARR-seq in tobacco plants grown in normal light/dark cycles (light) or completely in the dark (dark)",
    "for two days prior to RNA extraction. The enhancer strength of the individual insertion, substitution, and deletion variants was measured",
    "and normalized to the wild-type variant (log2 set to 0)."
  ),
  startCol = 1,
  startRow = 4
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 4, cols = 1)
mergeCells(wb, sheet = 1, rows = 4, cols = seq_len(ncol(data_table)))
setRowHeights(wb, sheet = 1, rows = 4, heights = 12.5 * 3)

xlsx_header <- names(data_table)
enrichment_ids <- which(xlsx_header %in% c('light', 'dark'))
xlsx_header[enrichment_ids] <- 'log2(enhancer strength)'

writeData(
  wb,
  sheet = 1,
  matrix(xlsx_header, nrow = 1),
  startCol = 1,
  startRow = 6,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_boldwrap, cols = seq_along(xlsx_header), rows = 6)

writeData(
  wb,
  sheet = 1,
  data_table,
  startCol = 1,
  startRow = 7,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

mergeCells(wb, sheet = 1, rows = 6, cols = enrichment_ids)
addStyle(wb, sheet = 1, style = xlsx_center, rows = 6, cols = enrichment_ids, stack = TRUE)

for (i in seq_along(xlsx_header)[-enrichment_ids]) {
  mergeCells(wb, sheet = 1, rows = 6:7, cols = i)
}

addStyle(wb, sheet = 1, style = xlsx_2digit, cols = enrichment_ids, rows = 8:(nrow(data_table) + 8), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_seq_font, cols = which(xlsx_header == 'sequence'), rows = 8:(nrow(data_table) + 8))

setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'unique'), widths = 19)

freezePane(wb, sheet = 1, firstActiveRow = 8)

saveWorkbook(wb, paste0(supp_data_dir, if_else(paper_id == '', '', paste0('tpc', paper_id, '-')), 'SupplementalDS', DS_id, '.xlsx'), overwrite = TRUE)


## scan enhancer sequences for transcription factor binding sites
FL_sequences <- enhancer_sequences |>
  filter(part == 'FL') |>
  select(enhancer, sequence) |>
  deframe() |>
  Biostrings::DNAStringSet()

TF_scan <- scan_sequences(TF_motifs, FL_sequences, RC = TRUE) |>
  as_tibble() |>
  select(sequence, 'id' = motif.i, start, stop, strand) |>
  mutate(
    id = paste0('TF-', id),
    orientation = if_else(strand == '+', 'fwd', 'rev'),
    tmp = if_else(strand == '-', start, stop),
    start = if_else(strand == '-', stop, start),
    stop = tmp,
    start = start - 0.5,
    stop = stop + 0.5,
    center = 0.5 * (start + stop),
  ) |>
  select(-strand, -tmp) |>
  distinct(sequence, id, start, stop, center, .keep_all = TRUE) |>
  arrange(center)

# number of mutation-sensitive regions with/without a contained TFBS
TF_content <- function(enhancer, reg_start, reg_end) {
  TF_scan |>
    filter(sequence == enhancer) |>
    mutate(
      start = start + 0.5,
      stop = stop - 0.5,
      in_region = between(start, reg_start, reg_end) & between(stop, reg_start, reg_end)
    ) |>
    filter(in_region) |>
    nrow()
}

enhancer_fragments |>
  filter(nchar(fragment) == 1) |>
  rowwise() |>
  mutate(
    TFBSs = TF_content(enhancer, start, end)
  ) |>
  ungroup() |>
  count(TFBSs > 0)

# export data
for (enhancer in unique(TF_scan$sequence)) {
  TF_scan |>
    filter(sequence == enhancer) |>
    select(-sequence) |>
    write_tsv(paste0(export_dir, 'PEV_ld_TF-scan_', enhancer, '.tsv'))
}

# function to lowercase non-ACGT letters
lower_ambiguous <- function(seq) {
  ACGT <- str_split(seq, '[^ACGT]')
  ambiguous <- str_split(seq, '[ACGT]') |>
    lapply(str_to_lower)
  
  sapply(
    seq_along(ACGT),
    function(i) {
      if (ACGT[[i]][1] == '') 
        interleave_vectors(ACGT[[i]], ambiguous[[i]])
      else
        interleave_vectors(ambiguous[[i]], ACGT[[i]])
    }
  ) |>
    sapply(paste0, collapse = '')
}

interleave_vectors <- function(vec_1, vec_2) {
  c(vec_1, vec_2)[order(c(seq_along(vec_1)*2 - 1, seq_along(vec_2)*2))]
}

# export transcription factor families
TF_scan |>
  distinct(id) |>
  rowwise() |>
  mutate(
    numID = as.integer(str_replace(id, 'TF-', '')),
    family = TF_motifs[[numID]]['altname'],
    consensus = TF_motifs[[numID]]['consensus']
  ) |>
  ungroup() |>
  arrange(numID) |>
  select(-numID) |>
  mutate(
    family = str_replace_all(family, '_', ' '),
    consensus = lower_ambiguous(consensus)
  ) |>
  rename('ID' = id) |>
  write_tsv(paste0(export_dir, 'PEV_ld_TF-scan_legend.tsv'))


## generate sequence logos for mutation-sensitive regions
# process data
mutSens_regions <- enhancer_fragments |>
  filter(nchar(fragment) == 1)

data_regions_PEV_ld <- data_heatmap_PEV_ld |>
  filter(part == 'B' & (condition == 'light' | enhancer == 'rbcS-E9') & type %in% c('substitution', 'WT')) |>
  inner_join(
    mutSens_regions,
    by = join_by(enhancer, between(x$position, y$start, y$end))
  ) |>
  filter(condition == 'light' | fragment %in% letters[1:3]) |>
  unite(
    col = 'region',
    enhancer,
    fragment
  )

# generate sequence logos for each region
cond_reg <- data_regions_PEV_ld |>
  distinct(condition, region)

all_motifs <- vector('list', length = nrow(cond_reg))

for (i in seq_along(all_motifs)) {
  cond <- cond_reg[i, 'condition', drop = TRUE]
  reg <- cond_reg[i, 'region', drop = TRUE]
  
  # filter data
  data_filtered <- data_regions_PEV_ld |>
    filter(condition == cond & region == reg)
  
  # get enrichment matrix
  motif_data <- data_filtered |>
    select(position, varNuc, enrichment) |>
    arrange(varNuc) |>
    pivot_wider(
      names_from = varNuc,
      values_from = enrichment
    ) |>
    arrange(position) |>
    select(-position) |>
    as.matrix() |>
    t()
  
  dimnames(motif_data) <- NULL
  
  # replace NA values with -1
  motif_data[is.na(motif_data)] <- -1
  
  # find scaling factor to yield an average information content of 1 per base
  beta <- 1
  ic <- create_motif(motif_data * beta, alphabet = 'DNA')[]$icscore
  
  while (ic < ncol(motif_data)) {
    beta <- beta + 0.1
    ic <- create_motif(motif_data * beta, alphabet = 'DNA')[]$icscore
  }
  
  beta <- max(beta - 0.1, 1)
  
  # create and export motifs and store in list
  motif <- create_motif(motif_data * beta, alphabet = 'DNA', type = 'ICM', name = paste0(reg, '_', cond))
  
  PWM_to_LaTeX(
    motif = motif,
    file = paste0(export_dir, 'PEV_ld_seqLogo_', reg, '_', cond, '.tsv'),
    offset = data_filtered |> distinct(start) |> pull() - 1,
    WT_seq = data_filtered |> distinct(sequence) |> pull()
  )
  
  all_motifs[[i]] <- motif
}

# region d of AB80 seems to contain two transcription factor binding sites, so we manually add its first half
AB80_d_trunc <- filter_motifs(all_motifs, name = 'AB80_d_light')[[1]][]$motif[, 1:8]
AB80_d_trunc <- create_motif(AB80_d_trunc, name = 'AB80_d_trunc')

all_motifs <- c(all_motifs, AB80_d_trunc)

# find matching transcription factor motifs
TF_hits <- tibble()

for (motif in all_motifs) {
  TF_hits <- TF_hits |>
    bind_rows(
      as_tibble(compare_motifs(c(TF_motifs, motif), compare.to = 73))
    )
}

# function to get alignment parameters for two motifs
motif_alignment <- function(motif_1, motif_2) {
  alignment <- view_motifs(c(motif_1, motif_2), return.raw = TRUE)
  
  start_pos <- sapply(alignment, function(x) min(which(colSums(x) > 0)))
  
  length <- unname(length(which(colSums(alignment[[2]]) > 0)))
  
  offset <- unname(diff(start_pos))
  
  orientation <- if_else(grepl('RC', names(alignment)[2], fixed = TRUE), 'rev', 'fwd')
  
  return(list(length = length, offset = offset, orientation = orientation))
}

# process transcription factor matches
TF_hits_filtered <- TF_hits |>
  filter(score > 0.85) |>
  rowwise() |>
  mutate(
    TFfamily = TF_motifs[[target.i]]['altname'],
    TFfamily = case_match(
      target.i,
      6 ~ 'G-box',
      13 ~ 'MYB',
      .default = TFfamily
    ),
    alignment = list(motif_alignment(filter_motifs(all_motifs, name = subject), TF_motifs[target.i]))
  ) |>
  ungroup() |>
  unnest_wider(alignment) |>
  separate_wider_delim(
    cols = subject,
    delim = '_',
    names = c('enhancer', 'region', 'condition')
  ) |>
  mutate(
    condition = if_else(condition == 'trunc', 'light', condition)
  )

# export data
enh_cond <- TF_hits_filtered |>
  distinct(enhancer, condition)

for (i in seq_len(nrow(enh_cond))) {
  enh <- enh_cond[i, 'enhancer', drop = TRUE]
  cond <- enh_cond[i, 'condition', drop = TRUE]
  
  TF_hits_filtered |> 
    filter(enhancer == enh & condition == cond) |>
    mutate(
      TF = paste0('TF-', target.i)
    ) |>
    select(region, TF, orientation, length, offset, orientation, TFfamily) |>
    write_tsv(paste0(export_dir, 'PEV_ld_TF-matches_', enh, '_', cond, '.tsv'))
}

# export sequence logos
TF_hits_PWMs <- TF_hits_filtered |>
  distinct(target.i, orientation) |>
  nest_by(across(everything())) |>
  reframe(
    motif = if_else(orientation == 'rev', motif_rc(TF_motifs[target.i]), TF_motifs[target.i]),
    LaTeX_PWM = list(PWM_to_LaTeX(motif, paste0(export_dir, 'seqLogo_TF-', target.i, '_', orientation, '.tsv')))
  )


## correlation between variant enhancer strength and TF motif match
# combine TF hits from the scanning and mutagenesis data approach
TF_hits_pos <- TF_hits_filtered |>
  filter(condition == 'light') |>
  inner_join(
    enhancer_fragments |> select(enhancer, 'region' = fragment, 'region_start' = start),
    by = c('enhancer', 'region')
  ) |>
  mutate(
    start = region_start + offset,
    stop = start + length - 1,
    start = start - 0.5,
    stop = stop + 0.5,
    center = 0.5 * (start + stop),
    method = 'mutagenesis'
  ) |>
  select(enhancer, 'id' = target, start, stop, center, orientation, method) |>
  bind_rows(
    TF_scan |>
      rename('enhancer' = sequence) |>
      mutate(
        id = str_replace(id, '-', '-cluster_'),
        method = 'scanning'
      )
  )

# calculate correlation
cor_TF_match <- data_heatmap_PEV_ld |>
  filter(type == 'substitution' & condition == 'light' & part == 'B') |>
  inner_join(
    TF_hits_pos,
    join_by(enhancer, between(position, start, stop))
  ) |>
  group_by(method, enhancer, offset, id, start, stop, center, orientation) |>
  group_modify(
    ~ add_row(.x, type = 'WT', variant = 'WT', duplSeq = FALSE, enrichment = 0)
  ) |>
  ungroup() |>
  inner_join(
    enhancer_variants,
    by = join_by(enhancer, type, variant, position, varNuc, offset, duplSeq)
  ) |>
  mutate(
    start = start + 0.5,
    stop = stop - 0.5,
    across(c(start, stop), ~ .x - offset),
    matched_seq = str_sub(sequence, start, stop),
    matched_seq = if_else(orientation == 'rev', as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(matched_seq))), matched_seq)
  ) |>
  filter(start > 0) |>
  group_by(id) |>
  mutate(
    match_score = score_match(filter_motifs(TF_motifs, name = id), matched_seq)
  ) |>
  group_by(method, enhancer, id, start, stop, center, orientation) |>
  summarise(
    pearson = cor(enrichment, match_score, use = 'complete.obs')
  ) |>
  ungroup() |>
  mutate(
    enhancer = ordered(enhancer, levels = c('AB80', 'Cab-1', 'rbcS-E9'))
  )

# export data
for (meth in unique(cor_TF_match$method)) {
  cor_TF_match |>
    filter(method == meth) |>
    arrange(enhancer) |>
    mutate(
      sample = paste0('points.', as.numeric(enhancer))
    ) |>
    group_by(enhancer) |>
    mutate(
      id = seq_len(n())
    ) |>
    ungroup() |>
    select(id, sample, pearson) |>
    pivot_wider(
      names_from = sample,
      values_from = pearson
    ) |>
    write_tsv(paste0(export_dir, 'PEV_cor_TF_match_', meth, '_points.tsv'))
  
  cor_TF_match |>
    filter(method == meth) |>
    group_by(enhancer) |>
    summarise(
      mean = mean(pearson)
    ) |>
    ungroup() |>
    select(enhancer, mean) |>
    mutate(
      id = as.numeric(enhancer)
    ) |>
    arrange(id) |>
    write_tsv(paste0(export_dir, 'PEV_cor_TF_match_', meth, '_mean.tsv'))
}

# export axis limits
cor_TF_match |>
  select('y' = pearson) |>
  summarise(
    across(y, list('min' = min, 'max' = max), .names = '{.col}{.fn}')
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cor_TF_match_axes.tsv'))


## consistency of mutation effect directions
# prepare data
data_TF_mutations <- TF_hits_pos |>
  mutate(
    start = start + 0.5,
    stop = stop - 0.5
  ) |>
  inner_join(
    data_heatmap_PEV_ld |>
      filter(part == 'B' & condition == 'light' & ! duplSeq),
    by = join_by(enhancer, between(y$position, x$start, x$stop))
  ) |>
  drop_na(enrichment) |>
  unite(
    'TF',
    id,
    center
  ) |>
  mutate(
    TF = str_replace_all(TF, 'cluster_', ''),
    TF = str_replace_all(TF, '_', '@')
  )

data_TF_mutations_count <- data_TF_mutations |>
  mutate(
    activity = if_else(enrichment > 0, 'increased', 'decreased')
  ) |>
  count(method, TF, activity) |>
  group_by(method, TF) |>
  mutate(
    percent = n / sum(n) * 100,
    order = percent[activity == 'decreased']
  ) |>
  ungroup() |>
  arrange(order) |>
  mutate(
    TF = ordered(TF, levels = unique(TF))
  )

# export data
for (met in unique(data_TF_mutations_count$method)) {
  data_TF_mutations_count |>
    filter(method == met) |>
    select(TF, activity, n, 'p' = percent) |>
    pivot_wider(
      names_from = activity,
      values_from = c(n, p)
    ) |>
    mutate(
      across(-TF, ~ replace_na(.x, 0)),
      id = as.numeric(droplevels(TF))
    ) |>
    write_tsv(paste0(export_dir, 'PEV_effect_directions_', met, '_bars.tsv'))
}

# export summary data
data_summary <- data_TF_mutations_count |>
  filter(activity == 'decreased') |>
  mutate(
    method = ordered(method, levels = c('mutagenesis', 'scanning'))
  )

LaTeX_boxplot(
  data = data_summary,
  samples_from = method,
  values_from = percent,
  file = paste0(export_dir, 'PEV_effect_directions_summary'),
  p_values = FALSE
)

data_summary |>
  rename(y = percent) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEV_effect_directions_summary_axes.tsv'))

# sample random regions
random <- vector('numeric', length = 100)

for (i in seq_along(random)) {
  set.seed(i)
  
  random_pos <- data_TF_mutations |>
    group_by(enhancer, offset, method, TF) |>
    summarise(
      length = unique(stop - start + 1),
      .groups = 'drop'
    ) |>
    group_by(enhancer) |>
    mutate(
      start = sample(seq(unique(offset), unique(offset + 168 - max(length))), n(), replace = TRUE),
      stop = start + length - 1,
      method = 'random'
    )
  
  random[i] <- random_pos |>
    inner_join(
      data_heatmap_PEV_ld |>
        filter(part == 'B' & condition == 'light' & ! duplSeq),
      by = join_by(enhancer, between(y$position, x$start, x$stop))
    ) |>
    drop_na(enrichment) |>
    mutate(
      activity = if_else(enrichment > 0, 'increased', 'decreased')
    ) |>
    count(method, TF, activity) |>
    group_by(method, TF) |>
    filter(n() == 2) |>
    mutate(
      percent = n[activity == 'decreased'] / sum(n) * 100
    ) |>
    ungroup() |>
    summarise(
      percent = median(percent)
    ) |>
    pull(percent)
}

tibble(
  sample = c('mean-sd', 'mean', 'mean+sd'),
  percent = c(mean(random) - sd(random), mean(random), mean(random) + sd(random)),
  linetype = c('dashed', 'solid', 'dashed')
) |>
  write_tsv(paste0(export_dir, 'PEV_effect_directions_summary_lines.tsv'))


## mean positional sensitivity to mutations
# process data
data_mutSens_PEV_ld <- data_mean_PEV_ld |>
  filter(type != 'control') |>
  group_by(enhancer, part, condition, orientation) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'WT']
  ) |>
  filter(variant != 'WT') |>
  group_by(enhancer, part, condition, variant, position) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  group_by(enhancer, part, position, condition) |>
  arrange(position) |>
  summarise(
    enrichment = mean(enrichment, na.rm = TRUE)
  ) |>
  group_by(enhancer, part, condition) |>
  mutate(
    enrichment = roll_mean(enrichment, n = 12, fill = NA, align = 'center')
  ) |> 
  ungroup() |>
  drop_na()

# export data
enh_part <- data_mutSens_PEV_ld |>
  distinct(enhancer, part)

for (i in seq_len(nrow(enh_part))) {
  enh = enh_part[i, 'enhancer', drop = TRUE]
  pt = enh_part[i, 'part', drop = TRUE]
  
  data_export <- data_mutSens_PEV_ld |>
    filter(enhancer == enh & part == pt)
  
  # export data
  data_export |>
    select(condition, position, enrichment) |>
    drop_na(enrichment) |>
    pivot_wider(
      names_from = condition,
      values_from = enrichment
    ) |>
    arrange(position) |>
    write_tsv(paste0(export_dir, 'PEV_ld_mutSens_', enh, '_', pt, '_lines.tsv'))
}

# export axis limits
axis_limits <- data_mutSens_PEV_ld |>
  select(enhancer, 'y' = enrichment) |>
  group_by(enhancer) |>
  summarise(
    across(everything(), list(min = min, max = max), .names = '{.col}{.fn}')
  ) |>
  ungroup()

axis_limits_recap <- axis_limits |>
  mutate(
    ratio = abs(ymin / ymax),
    fix_max = max(ratio) / ratio * ymin,
    fix_min = ratio / max(ratio) * ymax,
    ymin = if_else(ratio == max(ratio) | ymin <= fix_max, ymin, fix_max),
    ymax = if_else(ratio == max(ratio) | ymax >= fix_min, ymax, fix_min)
  ) 

axis_limits_enlarged <- axis_limits |>
  group_by(enhancer) |>
  mutate(
    yhelper = 0.05 * (ymax - ymin),
    ymin = ymin - yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  ungroup()

for (enh in unique(axis_limits_enlarged$enhancer)) {
  axis_limits_enlarged |>
    filter(enhancer == enh) |>
    select(-enhancer) |>
    write_tsv(paste0(export_dir, 'PEV_ld_mutSens_', enh, '_axes.tsv'))
  
  axis_limits_recap |>
    filter(enhancer == enh) |>
    select(ymin, ymax) |>
    write_tsv(paste0(export_dir, 'PEV_ld_mutSens_', enh, '_recap_axes.tsv'))
}


## correlation between variants in overlapping part of enhancer segments
data_overlap_PEV_ld <- data_heatmap_PEV_ld |>
  filter(! duplSeq & type != 'WT') |>
  select(condition, enhancer, part, variant, enrichment) |>
  pivot_wider(
    names_from = part,
    values_from = enrichment,
    names_prefix = 'part'
  ) |>
  drop_na(starts_with('part'))

# export plot data
enh_cond <- data_overlap_PEV_ld |>
  distinct(enhancer, condition)

for (i in seq_len(nrow(enh_cond))) {
  enh <- enh_cond[i, 'enhancer', drop = TRUE]
  cond <- enh_cond[i, 'condition', drop = TRUE]
  
  data_export <- data_overlap_PEV_ld |>
    filter(enhancer == enh & condition == cond) |>
    select(starts_with('part')) |>
    slice_sample(prop = 1)
  
  write_tsv(data_export, paste0(export_dir, 'PEV_cor_overlap_', enh, '_', cond, '_points.tsv'))
  
  LaTeX_cor_stats(
    x = data_export$partA,
    y = data_export$partB,
    file = paste0(export_dir, 'PEV_cor_overlap_', enh, '_', cond)
  )
}

# export axis limits
data_overlap_PEV_ld |>
  summarise(
    xmin = min(partA, partB),
    xmax = max(partA, partB),
    ymin = xmin,
    ymax = xmax
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cor_overlap_axes.tsv'))


### PEV library - circadian rhythm ###
## replicate correlations
# calculate axis limits
axis_limits <- data_reps_PEV_cr |>
  summarise(
    xmin = min(enrichment),
    xmax = max(enrichment),
    ymin = xmin,
    ymax = xmax
  )

max_count <- 0

# transform data
data_wide_PEV_cr <- data_reps_PEV_cr |>
  select(rep, condition, enhancer, part, orientation, variant, enrichment) |>
  pivot_wider(
    names_from = rep,
    names_prefix = 'rep',
    values_from = enrichment
  ) |>
  drop_na(starts_with('rep'))

# export plot data
for (cond in unique(data_wide_PEV_cr$condition)) {
  data_export <- data_wide_PEV_cr |>
    filter(condition == cond) |>
    select(starts_with('rep'))
  
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = rep1,
    y_values = rep2,
    file = paste0(export_dir, 'PEV_cr_cor_reps_', cond),
    xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
    axis_limits = FALSE
  )
  
  max_count <- max(max_count, this_count)
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cr_cor_reps_axes.tsv'))


## correlation between time points
# calculate axis limits
axis_limits <- data_reps_PEV_cr |>
  summarise(
    xmin = min(enrichment),
    xmax = max(enrichment),
    ymin = xmin,
    ymax = xmax
  )

max_count <- 0

# transform data
data_wide_PEV_cr_time <- data_mean_PEV_cr |>
  filter(! duplSeq) |>
  select(condition, enhancer, part, orientation, variant, enrichment) |>
  pivot_wider(
    names_from = condition,
    names_prefix = 't_',
    values_from = enrichment
  )

# export plot data
for (timepoint in c(6, 12, 18, 24)) {
  data_export <- data_wide_PEV_cr_time |>
    select(all_of(paste0('t_', c(0, timepoint)))) |>
    drop_na()
  
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = paste0('t_', timepoint),
    y_values = t_0,
    file = paste0(export_dir, 'PEV_cr_cor_0vs', timepoint),
    xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
    axis_limits = FALSE
  )
  
  max_count <- max(max_count, this_count)
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cr_cor_timepoints_axes.tsv')) 


## Correlation with PEV light
# transform data
data_wide_PEV_cr_ld <- data_mean_PEV_cr |>
  filter(! duplSeq) |>
  select(condition, enhancer, part, orientation, variant, enrichment) |>
  inner_join(
    data_mean_PEV_ld |>
      filter(condition == 'light') |>
      select(enhancer, orientation, part, variant, enrichment),
    by = c('enhancer', 'orientation', 'part', 'variant'),
    suffix = c('.cr', '.ld')
  )

# calculate axis limits
axis_limits <- data_wide_PEV_cr_ld |>
  summarise(
    xmin = min(enrichment.cr, enrichment.ld),
    xmax = max(enrichment.cr, enrichment.ld),
    ymin = xmin,
    ymax = xmax
  )

max_count <- 0

# export plot data
for (cond in unique(data_wide_PEV_cr_ld$condition)) {
  data_export <- data_wide_PEV_cr_ld |>
    filter(condition == cond) |>
    select(starts_with('enrichment'))
  
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = enrichment.cr,
    y_values = enrichment.ld,
    file = paste0(export_dir, 'PEV_cr_cor_PEV_ld_', cond),
    xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
    axis_limits = FALSE
  )
  
  max_count <- max(max_count, this_count)
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cr_cor_PEV_ld_axes.tsv')) 


## circadian rhythm of enhancer variants
# filter data
data_filtered_PEV_cr <- data_mean_PEV_cr |>
  filter(part == 'B' & ! duplSeq) |>
  group_by(condition, enhancer, variant) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  group_by(enhancer, variant) |>
  filter(n() == 5) |>
  ungroup()

# model circadian rhythm
models_PEV_cr <- data_filtered_PEV_cr |>
  nest_by(enhancer, variant) |>
  summarise(
    model = list(lm(enrichment ~ sin(2*pi*condition/24) + cos(2*pi*condition/24), data = data))
  ) |>
  ungroup()

# generate predictions for WT enhancers
data_predWT_PEV_cr <- models_PEV_cr |>
  filter(variant == 'WT') |>
  nest_by(enhancer, model) |>
  reframe(
    condition = seq(-2, 26, 0.1),
    enrichment = predict(model, tibble(condition)),
    rsquare = summary(model)$r.squared,
    amplitude = .5 * diff(range(enrichment)),
    equilibrium = max(enrichment) - amplitude,
    peak = condition[which.max(enrichment)],
    type = 'lines'
  ) |>
  ungroup() |>
  select(-model)

equilibrium_WT <- data_predWT_PEV_cr |>
  distinct(enhancer, equilibrium) |>
  deframe()

# merge WT predictions with measurements and process data 
data_WT_PEV_cr <- data_filtered_PEV_cr |>
  filter(variant == 'WT') |>
  select(-variant) |>
  mutate(
    type = 'points'
  ) |>
  bind_rows(
    data_predWT_PEV_cr |>
      select(enhancer, condition, enrichment, type)
  ) |>
  mutate(
    enrichment = enrichment - equilibrium_WT[enhancer]
  ) |>
  pivot_wider(
    names_from = enhancer,
    values_from = enrichment
  ) |>
  mutate(
    ZT = condition + 8
  ) |>
  arrange(ZT)

# export WT data
for (t in unique(data_WT_PEV_cr$type)) {
  data_WT_PEV_cr |>
    filter(type == t) |>
    select(-type) |>
    write_tsv(paste0(export_dir, 'PEV_cr_WT_', t, '.tsv'))
}

data_predWT_PEV_cr |>
  distinct(enhancer, rsquare, amplitude, peak) |>
  mutate(
    peak = peak + 8,
    peak = if_else(peak > 24, peak - 24, peak)
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cr_WT_stats.tsv'))

data_WT_PEV_cr |>
  summarise(
    xmin = min(ZT),
    xmax = max(ZT),
    ymin = min(AB80, `Cab-1`, `rbcS-E9`),
    ymax = max(AB80, `Cab-1`, `rbcS-E9`),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEV_cr_WT_axes.tsv'))

# calculate variant effects
data_variants_PEV_cr <- models_PEV_cr |>
  nest_by(across(everything())) |>
  summarise(
    condition = list(seq(0, 24, 0.1)),
    enrichment = list(predict(model, tibble(condition))),
    rsquare = summary(model)$r.squared,
    amplitude = .5 * diff(range(enrichment)),
    peak = condition[which.max(enrichment)]
  ) |>
  ungroup() |>
  select(enhancer, variant, rsquare, amplitude, peak)

# add enrichment data from timepoint 0 and normalize to WT
data_variants_norm_PEV_cr <- data_variants_PEV_cr |>
  inner_join(
    data_filtered_PEV_cr |>
      filter(condition == 0) |>
      select(-condition)
  ) |>
  group_by(enhancer) |>
  mutate(
    across(
      c(enrichment, amplitude, peak),
      ~ .x - .x[variant == 'WT']
    ),
    peak = case_when(
      peak < -12 ~ peak + 24,
      peak > 12 ~ peak - 24,
      .default = peak
    ),
  ) |>
  ungroup()

# export histogram data for peak shifts
axis_limits <- tibble()

for (enh in unique(data_variants_norm_PEV_cr$enhancer)) {
  data_export <- data_variants_norm_PEV_cr |>
    filter(enhancer == enh) |>
    mutate(
      fit = if_else(rsquare >= mean(rsquare), 'goodFit', 'badFit')
    )
  
  axis_limits <- axis_limits |>
    bind_rows(
      LaTeX_histogram(
        data = data_export,
        values_from = peak,
        group = fit,
        file = paste0(export_dir, 'PEV_cr_peak_', enh),
        data_range = c(-12, 12),
        bins = 30,
        value_axis = 'y',
        save_axis_limits = FALSE
      )
    )
}

axis_limits |>
  summarise(
    across(everything(), max)
  ) |>
  write_tsv(paste0(export_dir, 'PEV_cr_peak_axes.tsv'))

# export histogram data for amplitude and enhancer strength changes
data_range <- data_variants_norm_PEV_cr |>
  pull(enrichment) |>
  range()
data_range[2] <- data_range[2] + 1e-10 # the highest enrichment of Cab-1 ends up as `NA` without this (maybe due to a rounding error?)

axis_limits <- tibble()

for (enh in unique(data_variants_norm_PEV_cr$enhancer)) {
  data_export <- data_variants_norm_PEV_cr |>
    filter(enhancer == enh)
  
  axis_limits <- axis_limits |>
    bind_rows(
      LaTeX_histogram(
        data = data_export,
        values_from = amplitude,
        file = paste0(export_dir, 'PEV_cr_amplitude_', enh),
        data_range = data_range,
        bins = 30,
        value_axis = 'y',
        save_axis_limits = FALSE
      ) |>
        bind_cols(metric = 'amplitude')
    )
  
  axis_limits <- axis_limits |>
    bind_rows(
      LaTeX_histogram(
        data = data_export,
        values_from = enrichment,
        file = paste0(export_dir, 'PEV_cr_enrichment_', enh),
        data_range = data_range,
        bins = 30,
        value_axis = 'y',
        save_axis_limits = FALSE
      ) |>
        bind_cols(metric = 'enrichment')
    )
}

axis_limits <- axis_limits |>
  group_by(metric) |>
  summarise(
    across(everything(), max)
  )

for (m in unique(axis_limits$metric)) {
  axis_limits |>
    filter(metric == m) |>
    select(-metric) |>
    write_tsv(paste0(export_dir, 'PEV_cr_', m, '_axes.tsv')) 
}


## heatmaps
# get all possible variants
all_variants <- enhancer_variants |>
  select(-sequence) |>
  mutate(
    condition = list(unique(data_mean_PEV_cr$condition))
  ) |>
  unnest_longer(condition)

# merge data with all possible variants and normalize enrichment to WT
data_heatmap_PEV_cr <- data_mean_PEV_cr |>
  filter(type != 'control') |>
  group_by(condition, enhancer, part, orientation) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'WT']
  ) |>
  group_by(condition, enhancer, part, type, variant, position, varNuc) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  ungroup() |>
  right_join(all_variants) |>
  mutate(
    varNuc = ordered(varNuc, levels = c('X', 'A', 'C', 'G', 'T'))
  )

# export data
enh_part_cond <- data_heatmap_PEV_cr |>
  distinct(enhancer, part, condition)

for (i in seq_len(nrow(enh_part_cond))) {
  enh = enh_part_cond[i, 'enhancer', drop = TRUE]
  pt = enh_part_cond[i, 'part', drop = TRUE]
  cond = enh_part_cond[i, 'condition', drop = TRUE]
  
  # filter data
  data_export <- data_heatmap_PEV_cr |>
    filter(enhancer == enh & part == pt & condition == cond)
  
  # export insertion variants
  data_export |>
    filter(type == 'insertion') |>
    select(position, varNuc, enrichment) |>
    arrange(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_cr_heatmap_', enh, '_', pt, '_', cond, '_ins.tsv'), na = 'NaN')
  
  # export substitution and deletion variants
  data_export |>
    filter(type %in% c('substitution', 'deletion')) |>
    bind_rows(
      WT_sequences |>
        filter(enhancer == enh & part == pt) |>
        distinct(across(-condition))
    ) |>
    select(position, varNuc, enrichment) |>
    arrange(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_cr_heatmap_', enh, '_', pt, '_', cond, '_sub+del.tsv'), na = 'NaN')
}

# export WT positions and axis limits
for (enh in unique(data_heatmap_PEV_cr$enhancer)) {
  # WT positions
  WT_sequences |>
    filter(enhancer == enh) |>
    distinct(position, varNuc) |>
    write_tsv(paste0(export_dir, 'PEV_cr_heatmap_', enh, '_WT.tsv'), na = 'X')
  
  # axis limits
  data_heatmap_PEV_cr |>
    filter(enhancer == enh) |>
    drop_na(enrichment) |>
    summarise(
      `point meta max` = max(abs(enrichment)),
      `point meta min` = -`point meta max`
    ) |>
    write_tsv(paste0(export_dir, 'PEV_cr_heatmap_', enh, '_axes.tsv'))
}

## create supplemental data set with single-nucleotide variant data
# prepare data
data_table <- data_heatmap_PEV_cr |>
  filter(type != 'WT') |>
  mutate(
    segment = if_else(part == 'A', "5'", "3'"),
    wtNuc = WT_nucs[paste0(enhancer, '_', position)],
    wtNuc = replace_na(wtNuc, 'X'),
    across(ends_with('Nuc'), ~ if_else(.x == 'X', '-', .x)),
    time = paste0('ZT ', condition + 8)
  ) |>
  select(enhancer, part, segment, type, variant, position, wtNuc, varNuc, time, enrichment) |>
  pivot_wider(
    names_from = time,
    values_from = enrichment
  ) |>
  full_join(
    variant_sequences,
    by = join_by(enhancer, part, variant)
  ) |>
  full_join(
    data_variants_norm_PEV_cr |>
      filter(variant != 'WT') |>
      mutate(part = 'B') |>
      select(-enrichment),
    by = join_by(enhancer, part, variant)
  ) |>
  group_by(
    sequence
  ) |>
  arrange(position) |>
  mutate(
    across(
      c(amplitude, peak, rsquare),
      first
    ),
    peak = round(peak, 1)
  ) |>
  ungroup() |>
  arrange(enhancer, part, position, varNuc) |>
  select(
    enhancer, segment, type, variant, position,
    'wild-type base' = wtNuc, 'variant base' = varNuc, unique,
    starts_with('ZT'), 'Δ amplitude' = amplitude,
    'Δ peak time (h)' = peak, 'R²' = rsquare, sequence
  )

# export data to an excel file
DS_id <- 2

wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'single-nucleotide variants (cr)')

writeData(
  wb,
  sheet = 1,
  paste('Supplemental Data. Jores et al.', paper_ref),
  startCol = 1,
  startRow = 1
)

writeData(
  wb,
  sheet = 1,
  paste0('Supplemental Data Set ', DS_id,' Enhancer strength of single-nucleotide variants of the AB80, Cab-1, and rbcS-E9 enhancers in time course experiment.'),
  startCol = 1,
  startRow = 3
)

makeNameBold(wb, DS_id)
mergeCells(wb, sheet = 1, rows = 3, cols = seq_len(ncol(data_table)))

writeData(
  wb,
  sheet = 1,
  paste(
    "All possible singlenucleotide variants of the AB80, Cab-1, and rbcS-E9 enhancers were subjected to Plant STARR-seq in tobacco leaves.",
    "On the morning of the third day after transformation (ZT 0), the plants were shifted to constant light. Leaves were harvested for RNA",
    "extraction starting at mid-day (ZT 8) and in 6 hour intervals (ZT 14, 20, 26, and 32) afterwards. The enhancer strength of the individual",
    "insertion, substitution, and deletion variants was measured and normalized to the wild-type variant (log2 set to 0). For the 3' enhancer",
    "segments, a sine wave with a period of 24 h was fitted to the enhancer strength of a given variant across all sampled time points. The",
    "amplitude and time of highest strength (peak time) was extracted from the fit and the difference (Δ) between the variant and the wild-type",
    "enhancer is listed. The goodness-of-fit (R²) is also indicated."
  ),
  startCol = 1,
  startRow = 4
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 4, cols = 1)
mergeCells(wb, sheet = 1, rows = 4, cols = seq_len(ncol(data_table)))
setRowHeights(wb, sheet = 1, rows = 4, heights = 12.5 * 4)

xlsx_header <- names(data_table)
enrichment_ids <- grep('ZT', xlsx_header, fixed = TRUE)
xlsx_header[enrichment_ids] <- 'log2(enhancer strength)'

writeData(
  wb,
  sheet = 1,
  matrix(xlsx_header, nrow = 1),
  startCol = 1,
  startRow = 6,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_boldwrap, cols = seq_along(xlsx_header), rows = 6)

writeData(
  wb,
  sheet = 1,
  data_table,
  startCol = 1,
  startRow = 7,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

mergeCells(wb, sheet = 1, rows = 6, cols = enrichment_ids)
addStyle(wb, sheet = 1, style = xlsx_center, rows = 6, cols = enrichment_ids, stack = TRUE)

for (i in seq_along(xlsx_header)[-enrichment_ids]) {
  mergeCells(wb, sheet = 1, rows = 6:7, cols = i)
}

addStyle(wb, sheet = 1, style = xlsx_2digit, cols = c(enrichment_ids, which(xlsx_header %in% c('Δ amplitude', 'R²'))), rows = 8:(nrow(data_table) + 8), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_seq_font, cols = which(xlsx_header == 'sequence'), rows = 8:(nrow(data_table) + 8))

setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'unique'), widths = 19)
setColWidths(wb, sheet = 1, cols = grep('peak', xlsx_header, fixed = TRUE), widths = 13.5)

freezePane(wb, sheet = 1, firstActiveRow = 8)

saveWorkbook(wb, paste0(supp_data_dir, if_else(paper_id == '', '', paste0('tpc', paper_id, '-')), 'SupplementalDS', DS_id, '.xlsx'), overwrite = TRUE)


## mean positional sensitivity to mutations
# process data (normalize to WT variant at ZT 8)
data_mutSens_PEV_cr <- data_mean_PEV_cr |>
  filter(type != 'control') |>
  group_by(enhancer, part, orientation) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'WT' & condition == 0]
  ) |>
  filter(variant != 'WT') |>
  group_by(enhancer, part, condition, variant, position) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  group_by(enhancer, part, position, condition) |>
  arrange(position) |>
  summarise(
    enrichment = mean(enrichment, na.rm = TRUE)
  ) |>
  group_by(enhancer, part, condition) |>
  mutate(
    enrichment = roll_mean(enrichment, n = 12, fill = NA, align = 'center')
  ) |> 
  ungroup() |>
  drop_na()

# export data
enh_part <- data_mutSens_PEV_cr |>
  distinct(enhancer, part)

for (i in seq_len(nrow(enh_part))) {
  enh = enh_part[i, 'enhancer', drop = TRUE]
  pt = enh_part[i, 'part', drop = TRUE]
  
  data_export <- data_mutSens_PEV_cr |>
    filter(enhancer == enh & part == pt)
  
  # export data
  data_export |>
    select(condition, position, enrichment) |>
    drop_na(enrichment) |>
    pivot_wider(
      names_from = condition,
      values_from = enrichment
    ) |>
    arrange(position) |>
    write_tsv(paste0(export_dir, 'PEV_cr_mutSens_', enh, '_', pt, '_lines.tsv'))
}

# export axis limits
for (enh in unique(data_mutSens_PEV_cr$enhancer)) {
  data_mutSens_PEV_cr |>
    filter(enhancer == enh) |>
    select('y' = enrichment) |>
    summarise(
      across(everything(), list(min = min, max = max), .names = '{.col}{.fn}')
    ) |>
    write_tsv(paste0(export_dir, 'PEV_cr_mutSens_', enh, '_axes.tsv'))
}


### PEF library ###
## replicate correlations
# calculate axis limits
axis_limits <- data_reps_PEF |>
  summarise(
    xmin = min(enrichment),
    xmax = max(enrichment),
    ymin = xmin,
    ymax = xmax
  )

# get all possible combinations
rep_combos <- expand_grid(rep1 = unique(data_reps_PEF$rep), rep2 = unique(data_reps_PEF$rep)) |>
  filter(rep2 > rep1)

max_count <- 0

for (i in seq_len(nrow(rep_combos))) {
  rep1 <- rep_combos[i, 'rep1', drop = TRUE]
  rep2 <- rep_combos[i, 'rep2', drop = TRUE]
  
  # transform data
  data_wide_PEF <- data_reps_PEF |>
    mutate(
      rep = case_when(
        rep == rep1 ~ 'rep1',
        rep == rep2 ~ 'rep2'
      )
    ) |>
    drop_na(rep) |>
    select(rep, condition, type, contains('frag'), enrichment) |>
    pivot_wider(
      names_from = rep,
      values_from = enrichment
    ) |>
    drop_na(starts_with('rep'))
  
  # export plot data
  for (cond in unique(data_wide_PEF$condition)) {
    data_export <- data_wide_PEF |>
      filter(condition == cond) |>
      select(starts_with('rep'))
    
    this_count <- LaTeX_hexbin(
      data = data_export,
      x_values = rep1,
      y_values = rep2,
      file = paste0(export_dir, 'PEF_cor_reps_', rep1, 'vs', rep2, '_', cond),
      xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
      axis_limits = FALSE
    )
    
    max_count <- max(max_count, this_count)
  }
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEF_cor_reps_axes.tsv'))


## correlation between single fragment activity and DMS results
# area under the curve from the DMS experiment
DMS_AUCs <- data_mutSens_PEV_ld |>
  filter(part == 'B') |>
  inner_join(
    enhancer_fragments,
    by = join_by(
      enhancer,
      between(x$position, y$start, y$end)
    )
  ) |>
  group_by(condition, enhancer, fragment) |>
  arrange(position) |>
  summarise(
    AUC = -sum(diff(position) * roll_mean(enrichment, 2))
  ) |>
  ungroup()

# merge with fragment data
data_frags_AUCs <- data_mean_PEF |>
  filter(n_frags == 1) |>
  separate_wider_delim(
    cols = fragment_1,
    delim = '_',
    names = c('enhancer', 'fragment'),
    too_many = 'merge'
  ) |>
  select(condition, enhancer, fragment, enrichment) |>
  inner_join(DMS_AUCs) |>
  filter(fragment != 'd_rand')

# export plot data
for (cond in unique(data_frags_AUCs$condition)) {
  data_export <- data_frags_AUCs |>
    filter(condition == cond) |>
    select(enrichment, AUC, enhancer) |>
    slice_sample(prop = 1)
  
  write_tsv(data_export, paste0(export_dir, 'PEF_cor_AUC_', cond, '_points.tsv'))
  
  LaTeX_cor_stats(
    x = data_export$enrichment,
    y = data_export$AUC,
    file = paste0(export_dir, 'PEF_cor_AUC_', cond)
  )
  
  # export correlation for each enhancer individually (for rebuttal letter)
  data_export |>
    group_by(enhancer) |>
    summarise(
      intercept = coef(lm(enrichment ~ AUC))[1],
      slope = coef(lm(enrichment ~ AUC))[2],
      rsquare = cor(enrichment, AUC)^2
    ) |>
    ungroup() |>
    pivot_wider(
      names_from = enhancer,
      values_from = c(intercept, slope, rsquare)
    ) |>
    write_tsv(paste0(export_dir, 'PEF_cor_AUC_', cond, '_stats_per_enhancer.tsv'))
  
  # export axis limits
  data_frags_AUCs |>
    summarise(
      xmin = min(AUC[condition == cond]),
      xmax = max(AUC[condition == cond]),
      ymin = min(enrichment),
      ymax = max(enrichment)
    ) |>
    write_tsv(paste0(export_dir, 'PEF_cor_AUC_', cond, '_axes.tsv'))
}


## enhancer strength by number of fragments
data_frags_number <- data_mean_PEF |>
  filter(between(n_frags, 1, 3)) |>
  mutate(
    n_frags = ordered(n_frags)
  )

# export axis limits
axis_limits <- data_frags_number |>
  rename(y = enrichment) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'PEF_strength_by_frags_axes.tsv'))

# export data
for (cond in unique(data_frags_number$condition)) {
  data_export <- data_frags_number |>
    filter(condition == cond)
  
  LaTeX_violinplot(
    data = data_export,
    samples_from = n_frags,
    values_from = enrichment,
    file = paste0(export_dir, 'PEF_strength_by_frags_', cond),
    p_values = FALSE
  )
}


## enhancer strength in dark and light
# transform data
data_mean_PEF_wide <- data_mean_PEF |>
  filter(type == 'fragments') |>
  select(condition, enrichment, contains('frag')) |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment
  ) |>
  mutate(
    lightResp = light - dark,
    type = case_when(
      light <= 1 & dark <= 1 ~ 'inactive',
      lightResp >= 1 ~ 'light-activated',
      lightResp <= -1 ~ 'dark-activated',
      .default = 'constitutive'
    )
  ) |>
  mutate(
    type = ordered(type, levels = c('inactive', 'constitutive', 'light-activated', 'dark-activated')),
    type_id = as.numeric(type)
  )

# export scatter plot data
data_mean_PEF_wide |>
  drop_na(lightResp) |>
  select(dark, light, type, type_id) |>
  slice_sample(prop = 1) |>
  write_tsv(paste0(export_dir, 'PEF_strength_points.tsv'))

# export axis limits
data_mean_PEF_wide |>
  drop_na(lightResp) |>
  summarise(
    xmin = min(dark, light),
    xmax = max(dark, light),
    ymin = xmin,
    ymax = xmax
  ) |>
  write_tsv(paste0(export_dir, 'PEF_strength_axes.tsv'))

# export legend entries
data_mean_PEF_wide |>
  drop_na(lightResp) |>
  count(type) |>
  pivot_wider(
    names_from = type,
    values_from = n
  ) |>
  write_tsv(paste0(export_dir, 'PEF_strength_legend.tsv'))


## create supplemental data set with enhancer fragment combinations data in light and dark
# prepare data
fragment_sequences <- enhancer_fragments |>
  unite(
    col = 'fragment',
    enhancer,
    fragment
  ) |>
  select(fragment, sequence) |>
  deframe()

data_table <- data_mean_PEF_wide |>
  filter(n_frags <= 3) |>
  select(all_of(paste0('fragment_', seq(3, 1))), n_frags, type, dark, light) |>
  unite(
    col = 'name',
    starts_with('fragment'),
    sep = '+',
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    across(
      starts_with('fragment'),
      ~ fragment_sequences[.x],
      .names = '{.col}_seq'
    )
  ) |>
  unite(
    col = 'sequence',
    ends_with('_seq'),
    sep = 'gtgatg',
    na.rm = TRUE
  ) |>
  mutate(
    across(
      starts_with('fragment'),
      ~ replace_na(.x, '')
    )
  ) |>
  arrange(n_frags) |>
  select(name, 'number of fragments' = n_frags, all_of(paste0('fragment_', seq(3, 1))), 'category' = type, dark, light, sequence) |>
  rename_with(~ str_replace(.x, '_', ' '))

# export data to an excel file
DS_id <- 3

wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'enhancer fragment combinations')

writeData(
  wb,
  sheet = 1,
  paste('Supplemental Data. Jores et al.', paper_ref),
  startCol = 1,
  startRow = 1
)

writeData(
  wb,
  sheet = 1,
  paste0('Supplemental Data Set ', DS_id,' Enhancer strength of combinations of fragments of the AB80, Cab-1, and rbcS-E9 enhancers in the light or dark.'),
  startCol = 1,
  startRow = 3
)

makeNameBold(wb, DS_id)
mergeCells(wb, sheet = 1, rows = 3, cols = seq_len(ncol(data_table)))

writeData(
  wb,
  sheet = 1,
  paste(
    "Fragments of the AB80, Cab-1, and rbcS-E9 enhancers spanning 1–3 mutation-sensitive regions (a–e, ab, abc, bc, de; see Figure 6A) as",
    "well as a control fragment from a mutation-insensitive region in Cab-1 (ctrl) and a shuffled version of the AB80 fragment d (d_rand)",
    "were ordered as oligonucleotides. These fragments were randomly combined (separated by a 6-bp spacer: gtgatg) to create synthetic enhancer",
    "with up to three fragments which were then subjected to Plant STARR-seq in tobacco plants grown in normal light/dark cycles (light) or",
    "completely in the dark (dark) for two days prior to RNA extraction. Enhancer strength was normalized to a control construct without an",
    "enhancer (log2 set to 0). The synthetic enhancers were grouped into four categories: inactive, log2(enhancer strength) ≤ 1 in both conditions;",
    "constitutive, similar strength in both conditions; light-activated, at least two-fold more active in the light; dark-activated, at least",
    "two-fold more active in the dark."
  ),
  startCol = 1,
  startRow = 4
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 4, cols = 1)
mergeCells(wb, sheet = 1, rows = 4, cols = seq_len(ncol(data_table)))
setRowHeights(wb, sheet = 1, rows = 4, heights = 12.5 * 5)

xlsx_header <- names(data_table)
enrichment_ids <- which(xlsx_header %in% c('light', 'dark'))
xlsx_header[enrichment_ids] <- 'log2(enhancer strength)'

writeData(
  wb,
  sheet = 1,
  matrix(xlsx_header, nrow = 1),
  startCol = 1,
  startRow = 6,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_boldwrap, cols = seq_along(xlsx_header), rows = 6)

writeData(
  wb,
  sheet = 1,
  data_table,
  startCol = 1,
  startRow = 7,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

mergeCells(wb, sheet = 1, rows = 6, cols = enrichment_ids)
addStyle(wb, sheet = 1, style = xlsx_center, rows = 6, cols = enrichment_ids, stack = TRUE)

for (i in seq_along(xlsx_header)[-enrichment_ids]) {
  mergeCells(wb, sheet = 1, rows = 6:7, cols = i)
}

addStyle(wb, sheet = 1, style = xlsx_2digit, cols = enrichment_ids, rows = 8:(nrow(data_table) + 8), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_seq_font, cols = which(xlsx_header == 'sequence'), rows = 8:(nrow(data_table) + 8))

setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'name'), widths = 36)
setColWidths(wb, sheet = 1, cols = grep('fragment |category', xlsx_header), widths = 12)
setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'sequence'), widths = 53)

freezePane(wb, sheet = 1, firstActiveRow = 8, firstActiveCol = 2)

saveWorkbook(wb, paste0(supp_data_dir, if_else(paper_id == '', '', paste0('tpc', paper_id, '-')), 'SupplementalDS', DS_id, '.xlsx'), overwrite = TRUE)


## predict enhancer strength from single fragments
# add enhancer strength/light-responsiveness of single fragments
PEF_single_frags <- data_mean_PEF |>
  filter(n_frags == 1) |>
  unite(
    col = 'fragment',
    fragment_1,
    condition
  ) |>
  select(fragment, enrichment) |>
  deframe()

data_single_PEF <- data_mean_PEF |>
  filter(between(n_frags, 2, 3)) |>
  select(-where(~ all(is.na(.x)))) |>
  mutate(
    across(
      starts_with('fragment'),
      ~ PEF_single_frags[paste(.x, condition, sep = '_')]
    ),
    across(
      starts_with('fragment'),
      ~ replace_na(.x, 0)
    )
  )

# train linear models and use for prediction
data_predict_strength_PEF <- data_single_PEF |>
  nest_by(condition) |>
  mutate(
    model = list(lm('enrichment ~ fragment_1 + fragment_2 + fragment_3', data = data)),
    prediction = list(predict(model, data))
  ) |>
  ungroup() |>
  select(-model) |>
  unnest(everything())

# calculate axis limits
axis_limits <- data_predict_strength_PEF |>
  summarise(
    xmin = min(enrichment, prediction),
    xmax = max(enrichment, prediction),
    ymin = xmin,
    ymax = xmax
  )

max_count <- 0

# export plot data
for (cond in unique(data_predict_strength_PEF$condition)) {
  data_export <- data_predict_strength_PEF |>
    filter(condition == cond) |>
    select(enrichment, prediction)
  
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = prediction,
    y_values = enrichment,
    file = paste0(export_dir, 'PEF_cor_prediction_strength_', cond),
    xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
    axis_limits = FALSE
  )
  
  max_count <- max(max_count, this_count)
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEF_cor_prediction_strength_axes.tsv')) 


## predict light-responsiveness from single fragments
# calculate light-responsiveness
data_lr_PEF <- data_mean_PEF |>
  filter(between(n_frags, 2, 3)) |>
  select(-where(~ all(is.na(.x)))) |>
  group_by(across(contains('frag'))) |>
  filter(n() == 2) |>
  summarise(
    lightResp = enrichment[condition == 'light'] - enrichment[condition == 'dark']
  ) |>
  ungroup() |>
  mutate(
    across(
      starts_with('fragment'),
      ~ PEF_single_frags[paste0(.x, '_light')] - PEF_single_frags[paste0(.x, '_dark')]
    ),
    across(
      starts_with('fragment'),
      ~ replace_na(.x, 0)
    )
  )
  
# train linear models and use for prediction
data_predict_lr_PEF <- data_lr_PEF |>
  nest(data = everything()) |>
  rowwise() |>
  mutate(
    model = list(lm('lightResp ~ fragment_1 + fragment_2 + fragment_3', data = data)),
    prediction = list(predict(model, data))
  ) |>
  ungroup() |>
  select(-model) |>
  unnest(everything())

# export plot data
this_count <- LaTeX_hexbin(
  data = data_predict_lr_PEF,
  x_values = prediction,
  y_values = lightResp,
  file = paste0(export_dir, 'PEF_cor_prediction_lightResp'),
  xy_range = data_predict_lr_PEF |> select(prediction, lightResp) |> range()
)


## strength of individual and combined fragments
# set up axis limits
axis_limits <- tibble()

# select combination fragments
combi_frags <- enhancer_fragments |>
  filter(nchar(fragment) == 2)

for (i in seq_len(nrow(combi_frags))) {
  frag <- combi_frags[i, 'fragment', drop = TRUE]
  enh <- combi_frags[i, 'enhancer', drop = TRUE]
  
  # select and order all relevant fragments 
  all_fragments <- paste(enh, c(unlist(str_split(frag, '')), frag), sep = '_')
  
  fragment_order <- expand_grid(f1 = c(NA, all_fragments), f2 = all_fragments) |>
    filter(is.na(f1) | f1 != f2) |>
    unite(
      col = 'construct',
      everything(),
      sep = '+',
      na.rm = TRUE
    ) |>
    mutate(
      construct = str_replace_all(construct, paste0(enh, '_'), '')
    ) |>
    filter(nchar(construct) <= 3) |>
    pull(construct)
  
  # get enrichment and normalize to control construct without enhancer
  data_fragments <- data_reps_PEF |>
    group_by(condition) |>
    mutate(
      enrichment = enrichment - median(enrichment[type == 'noEnh-control'])
    ) |>
    ungroup() |>
    filter(n_frags <= 2 & fragment_1 %in% all_fragments & fragment_2 %in% c(NA, all_fragments)) |>
    unite(
      col = 'construct',
      fragment_2,
      fragment_1,
      sep = '+',
      na.rm = TRUE
    ) |>
    mutate(
      construct = str_replace_all(construct, paste0(enh, '_'), ''),
      construct = ordered(construct, levels = fragment_order)
    ) |>
    select(rep, condition, construct, enrichment) |>
    drop_na()
  
  # update axis limits
  axis_limits <- data_fragments |>
    bind_rows(axis_limits) |>
    group_by(condition) |>
    reframe(
      enrichment = range(enrichment),
      type = c('ymin', 'ymax')
    ) |>
    ungroup()
    
  
  # calculate mean enrichment across replicates
  data_summary <- data_fragments |>
    group_by(condition, construct) |>
    summarise(
      enrichment = mean(enrichment)
    ) |>
    ungroup()
  
  # export data
  for (cond in unique(data_fragments$condition)) {
    data_fragments |>
      filter(condition == cond) |>
      arrange(construct) |>
      mutate(
        sample = paste0('points.', as.numeric(construct))
      ) |>
      select(rep, sample, enrichment) |>
      pivot_wider(
        names_from = sample,
        values_from = enrichment
      ) |>
      write_tsv(paste0(export_dir, 'PEF_fragments_', enh, '_', frag, '_', cond, '_points.tsv'))
    
    data_summary |>
      filter(condition == cond) |>
      select(construct, 'mean' = enrichment) |>
      mutate(
        id = as.numeric(construct)
      ) |>
      arrange(id) |>
      write_tsv(paste0(export_dir, 'PEF_fragments_', enh, '_', frag, '_', cond, '_mean.tsv'))
  }
}

# export axis limits
for (cond in unique(axis_limits$condition)) {
  axis_limits |>
    filter(condition == cond) |>
    select(type, enrichment) |>
    pivot_wider(
      names_from = type,
      values_from = enrichment
    ) |>
    write_tsv(paste0(export_dir, 'PEF_fragments_', cond, '_axes.tsv'))
}

axis_limits |>
  pivot_wider(
    names_from = type,
    values_from = enrichment
  ) |>
  summarise(
    ymin = min(ymin),
    ymax = max(ymax)
  ) |>
  write_tsv(paste0(export_dir, 'PEF_fragments_axes.tsv'))


## effect of fragment order on enhancer strength
data_export <- data_mean_PEF |>
  filter(n_frags == 2 & fragment_1 != fragment_2) |>
  rowwise() |>
  mutate(
    fragments = paste0(sort(c(fragment_1, fragment_2)), collapse = '+')
  ) |>
  group_by(condition, fragments) |>
  filter(n() == 2) |>
  summarise(
    diff = abs(diff(enrichment)),
  ) |>
  ungroup() |>
  mutate(
    condition = ordered(condition, levels = c('light', 'dark'))
  )

LaTeX_violinplot(
  data = data_export,
  samples_from = condition,
  values_from = diff,
  file = paste0(export_dir, 'PEF_order_diff')
)


### PEval library (pooled PEFval and PEVdouble libraries) ###
## replicate correlations
# calculate axis limits
axis_limits <- data_reps_PEval |>
  summarise(
    xmin = min(enrichment),
    xmax = max(enrichment),
    ymin = xmin,
    ymax = xmax
  )

max_count <- 0

# transform data
data_wide_PEval <- data_reps_PEval |>
  pivot_wider(
    names_from = rep,
    names_prefix = 'rep',
    values_from = enrichment
  ) |>
  drop_na(starts_with('rep'))

# export plot data
for (cond in unique(data_wide_PEval$condition)) {
  data_export <- data_wide_PEval |>
    filter(condition == cond) |>
    select(starts_with('rep'))
    
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = rep1,
    y_values = rep2,
    file = paste0(export_dir, 'PEval_cor_reps')
  ) 
  
  this_count <- LaTeX_hexbin(
    data = data_export,
    x_values = rep1,
    y_values = rep2,
    file = paste0(export_dir, 'PEval_cor_reps_', cond),
    xy_range = axis_limits |> select(starts_with('x')) |> as.numeric(),
    axis_limits = FALSE
  )
  
  max_count <- max(max_count, this_count)
}

# export axis limits
axis_limits |>
  bind_cols(
    'point meta max' = max_count
  ) |>
  write_tsv(paste0(export_dir, 'PEval_cor_reps_axes.tsv'))


## correlation with main library results - PEVdouble
# combine data
data_wide_PEVdouble <- data_mean_PEV_ld |>
  filter(type %in% c('WT', 'deletion') & part == 'B') |>
  select(condition, type, enhancer, variant, orientation, enrichment) |>
  inner_join(
    data_mean_PEVdouble |>
      filter(type %in% c('WT', 'deletion') & is.na(mut2)) |>
      select(condition, type, enhancer, variant, orientation, enrichment),
    by = c('condition', 'type', 'enhancer', 'variant', 'orientation'),
    suffix = c('_PEV_ld', '_PEVdouble')
  ) |>
  group_by(condition) |>
  mutate(
    across(
      starts_with('enrichment'),
      ~ .x - median(.x)
    )
  ) |>
  ungroup()

# export plot data
for (cond in unique(data_wide_PEVdouble$condition)) {
  data_export <- data_wide_PEVdouble |>
    filter(condition == cond) |>
    select(starts_with('enrichment'), enhancer) |>
    slice_sample(prop = 1)
  
  write_tsv(data_export, paste0(export_dir, 'PEVdouble_cor_PEV_', cond, '_points.tsv'))
  
  LaTeX_cor_stats(
    x = data_export$enrichment_PEV_ld,
    y = data_export$enrichment_PEVdouble,
    file = paste0(export_dir, 'PEVdouble_cor_PEV_', cond)
  )
}

# export axis limits
data_wide_PEVdouble |>
  summarise(
    xmin = min(enrichment_PEV_ld, enrichment_PEVdouble),
    xmax = max(enrichment_PEV_ld, enrichment_PEVdouble),
    ymin = xmin,
    ymax = xmax
  ) |>
  write_tsv(paste0(export_dir, 'PEVdouble_cor_PEV_axes.tsv'))


## correlation with main library results - PEFval
# combine data
data_wide_PEFval <- data_mean_PEF |>
  filter(n_frags <= 3 & type == 'fragments') |>
  select(condition, fragment_1, fragment_2, fragment_3, n_frags, enrichment) |>
  inner_join(
    data_mean_PEFval |>
      filter(n_frags <= 3 & type == 'fragments') |>
      select(condition, fragment_1, fragment_2, fragment_3, n_frags, enrichment),
    by = c('condition', 'fragment_1', 'fragment_2', 'fragment_3', 'n_frags'),
    suffix = c('_PEF', '_PEFval')
  )

# export plot data
for (cond in unique(data_wide_PEFval$condition)) {
  data_export <- data_wide_PEFval |>
    filter(condition == cond) |>
    select(starts_with('enrichment'), n_frags) |>
    slice_sample(prop = 1)
  
  write_tsv(data_export, paste0(export_dir, 'PEFval_cor_PEF_', cond, '_points.tsv'))
  
  LaTeX_cor_stats(
    x = data_export$enrichment_PEF,
    y = data_export$enrichment_PEFval,
    file = paste0(export_dir, 'PEFval_cor_PEF_', cond)
  )
}

# export axis limits
data_wide_PEFval |>
  summarise(
    xmin = min(enrichment_PEF, enrichment_PEFval),
    xmax = max(enrichment_PEF, enrichment_PEFval),
    ymin = xmin,
    ymax = xmax
  ) |>
  write_tsv(paste0(export_dir, 'PEFval_cor_PEF_axes.tsv'))


## effect of double mutation relative to single ones
# single mutations
data_del_single <- data_mean_PEVdouble |>
  filter(type %in% c('deletion', 'WT') & is.na(mut2)) |>
  group_by(condition, enhancer, orientation) |>
  mutate(
    enrichment = enrichment - enrichment[variant == 'WT']
  ) |>
  ungroup() |>
  unite(
    col = 'variant',
    enhancer,
    variant,
    orientation,
    condition
  ) |>
  select(variant, enrichment) |>
  deframe()

# use sum of effects of single deletions to predict effect of double deletion
data_del_doubles_pred <- data_mean_PEVdouble |>
  filter(type == 'WT' | (type == 'deletion' & ! is.na(mut2))) |>
  unite(
    col = 'variant1',
    enhancer,
    mut1,
    orientation,
    condition,
    remove = FALSE
  ) |>
  unite(
    col = 'variant2',
    enhancer,
    mut2,
    orientation,
    condition,
    remove = FALSE
  ) |>
  mutate(
    eff1 = data_del_single[variant1],
    eff2 = data_del_single[variant2],
    add_eff = eff1 + eff2,
    distance = pos2 - pos1
  ) |>
  group_by(enhancer, orientation, condition) |>
  mutate(
    prediction = add_eff + enrichment[variant == 'WT']
  ) |>
  ungroup() |>
  drop_na(enrichment, add_eff, pos2)

# only keep deletion pairs at least 8 bp apart
data_del_doubles_filtered <- data_del_doubles_pred |>
  filter(distance > 7)

# annotate deletions within core TF binding sites (for rebuttal letter)
dels_in_TF <- data_mean_PEVdouble |>
  filter(type %in% c('deletion', 'WT')) |>
  distinct(enhancer, mut1, pos1) |>
  drop_na(mut1) |>
  inner_join(
    TF_hits_pos |>
      filter(method == 'mutagenesis') |>
      bind_rows( # add mutation-sensitive CCAAT and TGTGG sequences
        tibble(
          enhancer = c('AB80', 'Cab-1', 'rbcS-E9'),
          start = c(204.5, 217.5, 170.5),
          stop = c(209.5, 222.5, 175.5)
        )
      ),
    by = join_by(enhancer, between(pos1, start, stop))
  ) |>
  distinct(enhancer, mut1) |>
  unite(
    'mut',
    enhancer,
    mut1
  ) |>
  pull()

data_del_doubles_filtered_TF <- data_del_doubles_filtered |>
  rowwise() |>
  mutate(
    in_TF = sum(paste(enhancer, mut1, sep = '_') %in% dels_in_TF, paste(enhancer, mut2, sep = '_') %in% dels_in_TF)
  ) |>
  ungroup()

# export data
enh_cond <- data_del_doubles_filtered_TF |>
  distinct(enhancer, condition) |>
  filter(condition == 'light' | enhancer == 'rbcS-E9')

for (i in seq_len(nrow(enh_cond))) {
  enh <- enh_cond[i, 'enhancer', drop = TRUE]
  cond <- enh_cond[i, 'condition', drop = TRUE]

  data_export <- data_del_doubles_filtered_TF |>
    filter(enhancer == enh & condition == cond) |>
    select(enrichment, prediction, distance, in_TF) |>
    slice_sample(prop = 1)

  write_tsv(data_export, paste0(export_dir, 'PEVdouble_prediction_', enh, '_', cond, '_points.tsv'))

  LaTeX_cor_stats(
    x = data_export$prediction,
    y = data_export$enrichment,
    file = paste0(export_dir, 'PEVdouble_prediction_', enh, '_', cond)
  )
}

# export stats split by number of deletions in TFBSs
enh_cond_TF <- data_del_doubles_filtered_TF |>
  distinct(enhancer, condition, in_TF) |>
  filter(condition == 'light' | enhancer == 'rbcS-E9')

for (i in seq_len(nrow(enh_cond_TF))) {
  enh <- enh_cond_TF[i, 'enhancer', drop = TRUE]
  cond <- enh_cond_TF[i, 'condition', drop = TRUE]
  TF <- enh_cond_TF[i, 'in_TF', drop = TRUE]
  
  data_export <- data_del_doubles_filtered_TF |>
    filter(enhancer == enh & condition == cond & in_TF == TF)
  
  LaTeX_cor_stats(
    x = data_export$prediction,
    y = data_export$enrichment,
    file = paste0(export_dir, 'PEVdouble_prediction_', enh, '_', cond, '_', TF)
  )
}

# axis limits
axis_limits <- data_del_doubles_filtered_TF |>
  inner_join(enh_cond) |>
  group_by(condition) |>
  mutate(
    xmin = min(enrichment, prediction),
    xmax = max(enrichment, prediction),
    ymin = xmin,
    ymax = xmax
  ) |>
  group_by(pick(enhancer, condition, ends_with(c('min', 'max')))) |>
  summarise(
    `point meta min` = min(distance),
    `point meta max` = max(distance)
  ) |>
  ungroup()

enh_cond <- axis_limits |>
  distinct(enhancer, condition)

for (i in seq_len(nrow(enh_cond))) {
  enh <- enh_cond[i, 'enhancer', drop = TRUE]
  cond <- enh_cond[i, 'condition', drop = TRUE]
  
  axis_limits |>
    filter(enhancer == enh & condition == cond) |>
    select(-enhancer, -condition) |>
    write_tsv(paste0(export_dir, 'PEVdouble_prediction_', enh, '_', cond, '_axes.tsv'))
}


### dual-luciferase experiments ###
## 3' enhancer segments (light only)
# order samples
sample_order <- c('none', '35S', 'AB80', 'Cab-1', 'rbcS-E9')

data_dl_enhancers <- data_dl_norm |>
  filter(set == 1 & condition == 'light') |>
  mutate(
    enhancer = ordered(enhancer, levels = sample_order)
  ) |>
  drop_na(enhancer)

# export plot data
LaTeX_boxplot(
  data = data_dl_enhancers,
  samples_from = enhancer,
  values_from = l2ratio,
  file = paste0(export_dir, 'DL_enhancers'),
  outliers = TRUE,
  p_values = FALSE
)

# export axis limits
data_dl_enhancers |>
  rename(y = l2ratio) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'DL_enhancers_axes.tsv'))

# correlation with Plant STARR-seq results
data_dl_enh_summary <- data_dl_enhancers |>
  group_by(enhancer) |>
  summarise(
    CI = 0.5 * diff(t.test(l2ratio)$conf.int),
    l2ratio = mean(l2ratio)
  ) |>
  ungroup()

data_dl_enh_cor <- data_norm_PEfl |>
  mutate(
    part = if_else(enhancer %in% c('none', '35S'), 'B', part),
    orientation = replace_na(orientation, 'fwd')
  ) |>
  filter(condition == 'light' & part == 'B' & orientation == 'fwd') |>
  mutate(
    enrichment = enrichment - median(enrichment[enhancer == 'none'])
  ) |>
  group_by(enhancer) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  inner_join(data_dl_enh_summary)

# export correlation data
data_dl_enh_cor |> write_tsv(paste0(export_dir, 'DL_enhancers_cor_points.tsv'))

LaTeX_cor_stats(
  x = data_dl_enh_cor$enrichment,
  y = data_dl_enh_cor$l2ratio,
  file = paste0(export_dir, 'DL_enhancers_cor')
)

data_dl_enh_cor |>
  rename('x' = enrichment) |>
  summarise(
    across(x, list('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    ymin = min(l2ratio - CI),
    ymax = max(l2ratio + CI)
  ) |>
  write_tsv(paste0(export_dir, 'DL_enhancers_cor_axes.tsv'))


## 3' enhancer segments (light vs. dark)
# order samples
sample_order <- expand_grid(
  enhancer = c('none', '35S', 'AB80', 'Cab-1', 'rbcS-E9'),
  condition = c('light', 'dark')
) |>
  unite(
    'sample',
    enhancer,
    condition,
    sep = '.',
  ) |>
  pull(sample)

data_dl_enhancers_ld <- data_dl_norm |>
  filter(set == 1) |>
  unite(
    'sample',
    enhancer,
    condition,
    sep = '.',
    remove = FALSE
  ) |>
  mutate(
    sample = ordered(sample, levels = sample_order)
  ) |>
  drop_na(sample)

# export plot data
LaTeX_boxplot(
  data = data_dl_enhancers_ld,
  samples_from = sample,
  values_from = l2ratio,
  file = paste0(export_dir, 'DL_enhancers_ld'),
  enhancer,
  outliers = TRUE,
  p_values = FALSE
)

# export axis limits
data_dl_enhancers_ld |>
  rename(y = l2ratio) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'DL_enhancers_ld_axes.tsv'))


## synthetic enhancers
# order samples
synEnh  <- c(
  'Ce+Cb+Ce', 'Cb+Aa+Ad', 'Ce+Aa+Ac', 'Ca+Abc+Cb', 'Cb+Cb+Ca',
  'Cd+Cb+Cb', 'Aa+Cb+Cb', 'Cbc+Cde', 'Cbc+Abc', 'Cde+Ra', 'Rb+Abc+Rd'
)

sample_order <- data_dl_norm |>
  filter(enhancer %in% synEnh & condition == 'light') |>
  group_by(enhancer) |>
  summarise(
    l2ratio = mean(l2ratio)
  ) |>
  arrange(l2ratio) |>
  pull(enhancer)

sample_order <- c('none', '35S', sample_order)

synEnh_id <- setNames(c('none', '35S', paste0('syn', seq_len(length(sample_order) - 2))), sample_order)

data_dl_synEnh <- data_dl_norm |>
  filter((! enhancer %in% c('none', '35S') | set == 2) & condition == 'light') |>
  mutate(
    enhancer = ordered(enhancer, levels = sample_order, labels = synEnh_id)
  ) |>
  drop_na(enhancer)

# export plot data
LaTeX_boxplot(
  data = data_dl_synEnh,
  samples_from = enhancer,
  values_from = l2ratio,
  file = paste0(export_dir, 'DL_synEnh'),
  outliers = TRUE,
  p_values = FALSE
)

# export axis limits
data_dl_synEnh |>
  rename(y = l2ratio) |>
  summarise(
    across(y, c('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    yhelper = 0.05 * (ymax - ymin)
  ) |>
  ungroup() |>
  mutate(
    ymin = ymin - 3 * yhelper,
    ymax = ymax + yhelper
  ) |>
  select(-yhelper) |>
  write_tsv(paste0(export_dir, 'DL_synEnh_axes.tsv'))

# correlation with Plant STARR-seq results
data_dl_synEnh_summary <- data_dl_synEnh |>
  group_by(enhancer) |>
  summarise(
    CI = 0.5 * diff(t.test(l2ratio)$conf.int),
    l2ratio = mean(l2ratio)
  ) |>
  ungroup()

data_dl_synEnh_cor <- data_mean_PEF |>
  filter(condition == 'light' & (grepl('control', type, fixed = TRUE) | between(n_frags, 2, 3))) |>
  select(-fragment_5, -fragment_4) |>
  unite(
    'enhancer',
    fragment_3,
    fragment_2,
    fragment_1,
    sep = '+',
    na.rm = TRUE,
    remove = FALSE
  ) |>
  mutate(
    enhancer = case_when(
      type == 'noEnh-control' ~ 'none',
      type == '35Senh-control' ~ '35S',
      .default = enhancer
    ),
    enhancer = str_replace_all(enhancer, 'B80_|ab-1_|bcS-E9_', ''),
    enhancer = str_replace_all(enhancer, 'r', 'R'),
    synEnh = enhancer,
    enhancer = ordered(synEnh_id[enhancer], levels = synEnh_id)
  ) |>
  drop_na(enhancer) |>
  inner_join(data_dl_synEnh_summary) |>
  arrange(enhancer)

# export correlation data
data_dl_synEnh_cor |> 
  select(enhancer, synEnh, enrichment, l2ratio, CI) |>
  write_tsv(paste0(export_dir, 'DL_synEnh_cor_points.tsv'))

LaTeX_cor_stats(
  x = data_dl_synEnh_cor$enrichment,
  y = data_dl_synEnh_cor$l2ratio,
  file = paste0(export_dir, 'DL_synEnh_cor')
)

data_dl_synEnh_cor |>
  rename('x' = enrichment) |>
  summarise(
    across(x, list('min' = min, 'max' = max), .names = '{.col}{.fn}'),
    ymin = min(l2ratio - CI),
    ymax = max(l2ratio + CI)
  ) |>
  write_tsv(paste0(export_dir, 'DL_synEnh_cor_axes.tsv'))
