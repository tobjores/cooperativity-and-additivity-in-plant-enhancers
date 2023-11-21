library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)


### function to generate all possible single base variants of a given sequence
# WTseq:  character; wild-type sequence for which to generate variants
# offset: numeric; offset to add to the position in the variant specification
generate_variants <- function(WTseq, offset = 0) {
  if (WTseq == '') {
    stop('A wild-type sequence must be provided')
  }
  
  WTseq <- toupper(WTseq)
  
  if (grepl('[^ACGT]', WTseq)) {
    stop('The sequence can only contain the letters "A", "C", "G", and "T".')
  }
  
  WTlen <- nchar(WTseq)
  
  variants <- bind_rows(
    tibble(
      type = 'WT',
      sequence = WTseq,
      variant = 'WT'
    ),
    tibble(
      type = 'deletion',
      position = seq_len(WTlen),
      sequence = paste0(str_sub(WTseq, 1, position - 1), str_sub(WTseq, position + 1, WTlen)),
      variant = paste0(position + offset, 'del', str_sub(WTseq, position, position)),
      varNuc = 'X'
    ),
    tibble(
      type = 'insertion',
      position = rep(seq_len(WTlen - 1) + 0.5, each = 4),
      varNuc = rep(c('A', 'C', 'G', 'T'), times = WTlen - 1),
      sequence = paste0(str_sub(WTseq, 1, position), varNuc, str_sub(WTseq, position + 1, WTlen)),
      variant = paste0(position + offset, 'ins', varNuc)
    ),
    tibble(
      type = 'substitution',
      position = rep(seq_len(WTlen), each = 4),
      varNuc = rep(c('A', 'C', 'G', 'T'), times = WTlen),
      wtNuc = str_sub(WTseq, position, position),
      sequence = paste0(str_sub(WTseq, 1, position - 1), varNuc, str_sub(WTseq, position + 1, WTlen)),
      variant = paste0(position + offset, wtNuc, '>', varNuc)
    ) |>
      filter(wtNuc != varNuc) |>
      select(-wtNuc)
  ) |>
    mutate(
      position = position + offset,
      duplSeq = duplicated(sequence)
    )
  
  return(variants)
}


### PEfl library ###
## load raw subassembly
PEfl_subassembly_raw <- read_tsv('data/subassembly/PEfl/barcodes_pPSup_PEfl.tsv.gz') |>
  bind_rows(
    read_tsv('data/sequence_files/control_BCs_set1.tsv'),
    read_tsv('data/sequence_files/BCs_full-length_plant-enhancers.tsv')
  ) |>
  mutate(
    part = replace(part, is.na(part) & enhancer != 'none', 'FL')
  )

## remove duplicated barcodes
PEfl_subassembly_noDup <- PEfl_subassembly_raw |>
  filter(! barcode %in% barcode[duplicated(barcode)])

## save final subassembly
write_tsv(PEfl_subassembly_noDup, 'data/subassembly/PEfl/subassembly_pPSup_PEfl.tsv.gz')


### PEV library ###
## load enhancer sequences and generate variants
WT_sequences <- read_tsv('data/sequence_files/plant-enhancer_sequences.tsv') |>
  group_by(enhancer) |>
  mutate(
    offset = str_locate(sequence[part == 'FL'], sequence)[, 1] - 1
  ) |>
  ungroup() |>
  filter(part != 'FL')

variant_sequences <- WT_sequences |> 
  group_by(across(everything())) |>
  reframe(
    generate_variants(sequence, offset)
  )

write_tsv(variant_sequences, 'data/sequence_files/plant-enhancer_variants.tsv')

## load raw subassemblies
PEV_subassembly_raw <- expand_grid(enhancer = c('AB80', 'Cab-1', 'rbcS-E9'), orientation = c('fwd', 'rev')) |>
  group_by(across(everything())) |>
  reframe(
    read_tsv(paste0('data/subassembly/PEV/raw_subassembly_pPSup_', enhancer, '-var_', orientation,'.tsv.gz'))
  ) |>
  bind_rows(read_tsv('data/sequence_files/control_BCs_set2.tsv')) |>
  mutate(
    assembly_count = if_else(enhancer %in% c('none', '35S'), Inf, assembly_count)
  )

## remove duplicated barcodes and barcodes seen less than 5 times
PEV_subassembly_noDup <- PEV_subassembly_raw |>
  filter(! barcode %in% barcode[duplicated(barcode)] & assembly_count >= 5)

## add variant information
PEV_subassembly_final <- PEV_subassembly_noDup |>
  left_join(
    variant_sequences,
    by = c('enhancer', 'sequence'),
    relationship = 'many-to-many'
  ) |>
  mutate(
    type = if_else(enhancer %in% c('none', '35S'), 'control', type),
    variant = case_when(
      enhancer == 'none' ~ 'noEnh-control',
      enhancer == '35S' ~ '35Senh-control',
      .default = variant
    )
  ) |>
  drop_na(type) |>
  select(-sequence , -assembly_count)

## save final subassembly
write_tsv(PEV_subassembly_final, 'data/subassembly/PEV/subassembly_pPSup_PEV.tsv.gz')


### PEVdouble library ###
## load library sequences
PEVdouble_seqs <- read_tsv('data/sequence_files/PEVdouble_sequences.tsv') 

## load raw subassemblies
PEVdouble_subassembly_raw <- read_tsv('data/subassembly/PEVdouble/raw_subassembly_pPSup_PEVdouble.tsv.gz') |>
  bind_rows(read_tsv('data/sequence_files/control_BCs_set1.tsv')) |>
  mutate(
    assembly_count = if_else(enhancer %in% c('none', '35S'), Inf, assembly_count)
  ) |>
  filter(is.na(orientation) | orientation == 'fwd') |>
  select(-orientation)

## remove duplicated barcodes and barcodes seen less than 5 times
PEVdouble_subassembly_noDup <- PEVdouble_subassembly_raw |>
  filter(! barcode %in% barcode[duplicated(barcode)] & assembly_count >= 5)

## add variant information
PEVdouble_subassembly_final <- PEVdouble_subassembly_noDup |>
  left_join(PEVdouble_seqs, by = 'sequence') |>
  filter(! if_all(starts_with('enhancer'), is.na)) |>
  mutate(
    enhancer = coalesce(enhancer.x, enhancer.y),
    type = replace_na(type, 'control'),
    variant = case_when(
      enhancer == 'none' ~ 'noEnh-control',
      enhancer == '35S' ~ '35Senh-control',
      .default = variant
    )
  ) |>
  select(-c(ends_with(c('.x', '.y')), sequence, assembly_count))

## save final subassembly
write_tsv(PEVdouble_subassembly_final, 'data/subassembly/PEVdouble/subassembly_pPSup_PEVdouble.tsv.gz')


### PEF library ###
## load raw subassembly
PEF_subassembly_raw <- read_tsv('data/subassembly/PEF/raw_subassembly_pPSup_PEF.tsv.gz') |>
  bind_rows(read_tsv('data/sequence_files/control_BCs_set2.tsv')) |>
  mutate(
    type = case_when(
      enhancer == 'none' ~ 'noEnh-control',
      enhancer == '35S' ~ '35Senh-control',
      .default = 'fragments'
    ),
    assembly_count = if_else(is.na(enhancer), assembly_count, Inf)
  ) |>
  select(-enhancer)

## remove duplicated barcodes and barcodes seen less than 5 times
PEF_subassembly_noDup <- PEF_subassembly_raw |>
  filter(! barcode %in% barcode[duplicated(barcode)] & assembly_count >= 5)

## identify individual fragments
fragment_names <- read_tsv('data/sequence_files/enhancer-fragments.tsv') |>
  unite(
    col = 'fragment',
    c(enhancer, fragment)
  ) |>
  select(sequence, fragment) |>
  deframe()

PEF_subassembly_frags <- PEF_subassembly_noDup |>
  mutate(
    sequence = str_sub(sequence, 7, -7)
  ) |>
  separate_wider_delim(
    cols = sequence,
    delim = 'GTGATG',
    names_sep = '_',
    too_few = 'align_end'
  ) |>
  mutate(
    across(
      starts_with('sequence'),
      ~ if_else(.x %in% c(names(fragment_names), NA), unname(fragment_names[.x]), 'unknown fragment')
    )
  ) |>
  filter(
    if_all(
      starts_with('sequence'),
      ~ is.na(.x) | .x != 'unknown fragment'
    )
  ) |>
  select(
    where(~ any(! is.na(.x)))
  ) |>
  rename_with(
    ~ paste0('fragment_', seq(length(.x), 1)),
    starts_with('sequence')
  )

## count number of fragments
frags_helper <- PEF_subassembly_frags |>
  select(starts_with('fragment')) |>
  is.na() |>
  rowSums()

PEF_subassembly_final <- PEF_subassembly_frags |>
  mutate(
    n_frags = max(frags_helper) - frags_helper
  ) |>
  select(-assembly_count)

## save final subassembly
write_tsv(PEF_subassembly_final, 'data/subassembly/PEF/subassembly_pPSup_PEF.tsv.gz')


### PEFval library ###
## load raw subassembly
PEFval_subassembly_raw <- read_tsv('data/subassembly/PEFval/raw_subassembly_pPSup_PEFval.tsv.gz') |>
  bind_rows(read_tsv('data/sequence_files/control_BCs_set1.tsv')) |>
  mutate(
    type = case_when(
      enhancer == 'none' ~ 'noEnh-control',
      enhancer == '35S' ~ '35Senh-control',
      .default = 'fragments'
    ),
    assembly_count = if_else(is.na(enhancer), assembly_count, Inf)
  ) |>
  filter(is.na(orientation) | orientation == 'fwd') |>
  select(-enhancer, -orientation)

## remove duplicated barcodes and barcodes seen less than 5 times
PEFval_subassembly_noDup <- PEFval_subassembly_raw |>
  filter(! barcode %in% barcode[duplicated(barcode)] & assembly_count >= 5)

## identify individual fragments
fragment_names <- read_tsv('data/sequence_files/enhancer-fragments.tsv') |>
  unite(
    col = 'fragment',
    c(enhancer, fragment)
  ) |>
  select(sequence, fragment) |>
  deframe()

PEFval_subassembly_frags <- PEFval_subassembly_noDup |>
  mutate(
    sequence = str_sub(sequence, 7, -7)
  ) |>
  separate_wider_delim(
    cols = sequence,
    delim = 'GTGATG',
    names_sep = '_',
    too_few = 'align_end'
  ) |>
  mutate(
    across(
      starts_with('sequence'),
      ~ if_else(.x %in% c(names(fragment_names), NA), unname(fragment_names[.x]), 'unknown fragment')
    )
  ) |>
  filter(
    if_all(
      starts_with('sequence'),
      ~ is.na(.x) | .x != 'unknown fragment'
    )
  ) |>
  select(
    where(~ any(! is.na(.x)))
  ) |>
  rename_with(
    ~ paste0('fragment_', seq(length(.x), 1)),
    starts_with('sequence')
  )

## count number of fragments and remove construct with more than 3 fragments
frags_helper <- PEFval_subassembly_frags |>
  select(starts_with('fragment')) |>
  is.na() |>
  rowSums()

PEFval_subassembly_final <- PEFval_subassembly_frags |>
  mutate(
    n_frags = max(frags_helper) - frags_helper
  ) |>
  filter(n_frags <= 3) |>
  select(-assembly_count, -where(~ all(is.na(.x))))

## save final subassembly
write_tsv(PEFval_subassembly_final, 'data/subassembly/PEFval/subassembly_pPSup_PEFval.tsv.gz')
