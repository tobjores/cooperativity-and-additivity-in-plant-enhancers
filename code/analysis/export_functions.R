### define functions to export data for plotting in LaTeX with pgfplots

## function to export data for hexbin plots
# data:             R object; data frame to export
# x_values:         unquoted column specification; column containing the x values
# y_values:         unquoted column specification; column containing the y values
# file:             character; base name of the export file without extension (full name is: "<file>_hexbin.tsv")
# stats:            logical; whether to export correlation statistics (file name: "<file>_stats.tsv")
# x_bins:           numeric; number of bins along the x range
# x_range:          numeric vector, length 2; force x range to x_range[1] - x_range[2] (default is to use data range)
# y_range:          numeric vector, length 2; force y range to y_range[1] - y_range[2] (default is to use data range)
# xy_range:         numeric vector, length 2; set x_range and y_range simultaneously
# axis_limits:      logical; whether to export a file with axis limits (file name: "<file>_axes.tsv")
# count_fun:        function; function to apply to the counts per hexbin (default is `log2`; set to `NULL` to get raw counts)
LaTeX_hexbin <- function(data, x_values, y_values, file, stats = TRUE, x_bins = 50, x_range = NULL, y_range = NULL, xy_range = NULL, axis_limits = TRUE, count_fun = log2) {
  x_col <- enquo(x_values)
  y_col <- enquo(y_values)
  
  x_data <- data |> pull(!! x_col)
  y_data <- data |> pull(!! y_col)
  
  # get axis limits
  if (! is.null(xy_range)) {
    x_range <- xy_range
    y_range <- xy_range
  }
  
  if (is.null(x_range)) {
    x_range <- range(x_data) 
  }
  if (is.null(y_range)) {
    y_range <- range(y_data) 
  }
  
  if (axis_limits) {
    axis_limits <- tibble(type = c('min', 'max'), x_range, y_range) |>
      pivot_wider(
        names_from = type,
        values_from = -type,
        names_glue = '{substr(.value, 1, 1)}{type}'
      )
    
    write_tsv(axis_limits, paste0(file, '_axes.tsv'))
  }
  
  # generate hexbin coordinates
  hbin <- hexbin(
    x = x_data,
    y = y_data,
    xbins = x_bins,
    xbnds = x_range,
    ybnds = y_range
  )
  
  data_hbin <- as_tibble(hcell2xy(hbin)) |>
    mutate(
      count = slot(hbin, 'count')
    )
  
  if (! is.null(count_fun)) {
    data_hbin <- data_hbin |>
      mutate(
        count = count_fun(count)
      )
  }
  
  write_tsv(data_hbin, paste0(file, '_hexbin.tsv'))
  
  # caluclate correlation
  if (stats) {
    LaTeX_cor_stats(
      x = x_data,
      y = y_data,
      file = file
    )
  }
  
  return(max(data_hbin$count))
}


## function to export correlation statistics
# x:    numeric; x values for the correlation
# y:    numeric; y values for the correlation
# file: character; base name of the export file without extension (full name is: "<file>_stats.tsv")
LaTeX_cor_stats <- function(x, y, file) {
  tibble(x = x, y = y) |>
    drop_na() |>
    summarise(
      n = n(),
      spearman = cor(x, y, method = 'spearman'),
      rsquare = cor(x, y, method = 'pearson')^2
    ) |>
    write_tsv(paste0(file,  '_stats.tsv')) 
}


## function to export data for boxplots
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_boxplot.tsv")
# ...:              unquoted column specifications; additional columns to include in the export; should contain only a single value per sample!
# outliers:         logical; whether to create a file for the outlier points (will be saved as: "<file>_outliers.tsv")
# p_values:         logical; whether to calculate all possible pairwise p values using a Mann-Whitney-Wilcox test (will be saved as: "<file>_pvalues.tsv")
# use_sample_name:  logical; whether to use the sample name instead of a numeric index as column name for the outliers
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)
LaTeX_boxplot <- function(data, samples_from, values_from, file, ..., outliers = TRUE, p_values = TRUE, use_sample_name = FALSE, data_range = NULL) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  additional_cols <- enquos(...)
  
  # save boxplot data (median, quantiles, whiskers, sample size, ...)
  box <- data |>
    select(sample = !! sample_col, value =  !! value_col, !!! additional_cols) |>
    group_by(sample, !!! additional_cols) |>
    summarise(
      id = cur_group_id(),
      med = median(value, na.rm = TRUE),
      lq = quantile(value, 0.25, na.rm = TRUE),
      uq = quantile(value, 0.75, na.rm = TRUE),
      iqr = 1.5 * IQR(value, na.rm = TRUE),
      lw = max(lq - iqr, min(value)),
      uw = min(uq + iqr, max(value)),
      outliers = sum(outliers(value) != 0),
      n = n()
    ) |>
    select(-iqr) |>
    ungroup() |>
    arrange(sample)
  
  if (length(box$sample) != length(unique(box$sample))) {
    warning('The additional columns contain more than one value per sample. The results might be unexpected!')
  }
  
  write_tsv(box, paste0(file, '_boxplot.tsv'), na = 'NaN')
  
  # save outlier coordinates (optional)
  if (outliers) {
    outlier_data <- data |>
      select(value =  !! value_col, sample = !! sample_col) |>
      mutate(
        sample = droplevels(sample)
      ) |>
      group_by(sample) |>
      summarise(
        outlier = list(value[outliers(value) != 0]),
        id = list(seq_along(unlist(outlier)))
      ) |>
      unnest(c(outlier, id), keep_empty = FALSE) |>
      ungroup() |>
      arrange(sample)
    
    if (nrow(outlier_data) > 0) {
      if (use_sample_name) {
        outlier_data  <- outlier_data |>
          mutate(
            sample = paste0('outlier.', sample)
          )
      } else {
        outlier_data  <- outlier_data |>
          mutate(
            sample = paste0('outlier.', as.numeric(sample))
          )
      }
      
      outlier_data <- outlier_data |>
        pivot_wider(
          names_from = sample,
          values_from = outlier,
          id_cols = id
        ) |>
        select(-id)
      
      write_tsv(outlier_data, paste0(file, '_outliers.tsv'), na = 'NaN')
    }
  }
  
  # generate pairwise p-values (optional)
  if (p_values) {
    wilcox_all(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      data_range = data_range,
      use_sample_name = use_sample_name
    )
  }
}


## function to export data for violin plots
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_violin.tsv")
# boxplot:          logical; whether to create a file for a boxplot (will be saved as: "<file>_boxplot.tsv")
# ...:              unquoted column specifications; additional columns to include in the boxplot export; should contain only a single value per sample!
# outliers:         logical; whether to create a file for the outlier points (will be saved as: "<file>_outliers.tsv")
# p_values:         logical; whether to calculate all possible pairwise p values using a Mann-Whitney-Wilcox test (will be saved as: "<file>_pvalues.tsv")
# half:             logical; whether to export only one half of the violinplot for split violins
# use_sample_name:  logical; whether to use the sample name instead of a numeric index as column name
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)
LaTeX_violinplot <- function(
  data, samples_from, values_from, file, ..., boxplot = TRUE, outliers = FALSE, p_values = TRUE, half = FALSE, use_sample_name = FALSE, data_range = NULL
) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  
  # create boxplot (optional)
  if (boxplot) {
    LaTeX_boxplot(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      ...,
      p_values = FALSE,
      outliers = outliers,
      use_sample_name = use_sample_name,
      data_range = data_range
    )
  }
  
  # save violin plot data
  violin <- data |>
    select(value =  !! value_col, sample = !! sample_col) |>
    group_by(sample) |>
    summarise(
      x = list(density(value, from = min(value), to = max(value), n = 256)$y),
      y = list(density(value, from = min(value), to = max(value), n = 256)$x)
    ) |>
    arrange(sample) |>
    ungroup()
  
  if (! use_sample_name) {
    violin <- violin |>
      mutate(
        sample = ordered(seq_len(n()))
      )
  }
  
  violin <- violin |>
    pivot_wider(names_from = sample, values_from = c(x, y), names_sep = '.') |>
    unnest(everything()) |>
    ungroup() |>
    mutate(
      id = row_number()
    )
  
  # normalize width (width of widest violin is set to 1)
  max_x <- max(select_at(violin, vars(starts_with('x.'))))
  
  violin <- violin |>
    mutate_at(
      vars(starts_with('x.')),
      function(x) {x / max_x * 0.5}
    )
  
  if (half) {
    max_y <- sapply(colnames(violin), function(x) max(pull(violin, x)))
    max_y[grepl('x.', names(max_y), fixed = TRUE)] <- 0
    min_y <- sapply(colnames(violin), function(x) min(pull(violin, x)))
    min_y[grepl('x.', names(min_y), fixed = TRUE)] <- 0
    
    violin <- bind_rows(min_y, violin, max_y) |>
      select(-id)
  } else {
    violin <- bind_rows(violin, arrange(mutate_at(violin, vars(starts_with('x.')), function(x) {-x}), desc(id))) |>
      select(-id)
  }
  
  write_tsv(violin, paste0(file, '_violin.tsv'), na = 'NaN')
  
  # generate pairwise p-values (optional)
  if (p_values) {
    wilcox_all(
      data = data,
      samples_from = !! sample_col,
      values_from = !! value_col,
      file = file,
      data_range = data_range,
      use_sample_name = use_sample_name,
      half = half
    )
  }
}


## calculate and export all pairwise p values using a Mann-Whitney-Wilcox test
# data:             R object; data frame to export
# samples_from:     unquoted column specification; column containing the sample names by which to summarise
# values_from:      unquoted column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_pvalues.tsv")
# use_sample_name:  logical; whether to use the sample name instead of a numeric index as column name
# data_range:       numeric vector, length 2; min and max coordinates to be considered for pvalue coordinates (the range of <data> is used if NULL)
# half:             logical; whether to export a simplified version of the pvalues for half violin plots
wilcox_all <- function(data, samples_from, values_from, file, data_range = NULL, use_sample_name = FALSE, half = FALSE) {
  sample_col <- enquo(samples_from)
  value_col <- enquo(values_from)
  
  if (is.null(data_range)) {
    data_range <- range(pull(data, !! value_col))
  } else if (length(data_range) != 2 | ! is.numeric(data_range)) {
    stop('"data_range" must be a two component numerical vector')
  }
  
  samples <- data |> distinct(!! sample_col) |> arrange(!! sample_col) |> pull(!! sample_col)
  
  if (! is.ordered(samples)) {
    samples <- ordered(samples) 
  }
  
  
  pvalues <- tibble(
    exp1 = combn(samples, 2)[1, ],
    exp2 = combn(samples, 2)[2, ]
  ) |>
    rowwise() |>
    mutate(
      p.value = wilcox.test(
        pull(filter(data, !! sample_col == exp1), !! value_col),
        pull(filter(data, !! sample_col == exp2), !! value_col)
      )$p.value,
      max1 = max(pull(filter(data, !! sample_col == exp1), !! value_col)),
      max2 = max(pull(filter(data, !! sample_col == exp2), !! value_col)),
      x1 = as.numeric(exp1),
      x2 = as.numeric(exp2)
    )
  
  if (half) {
    pvalues <- pvalues |>
      filter(as.numeric(exp1) + 1 == as.numeric(exp2) & as.numeric(exp1) %% 2 == 1) |>
      mutate(
        max1 = max(max1, max2) - 0.1 * diff(data_range),
        max2 = max1
      )
  }
  
  if (! use_sample_name) {
    pvalues <- pvalues |>
      mutate(
        exp1 = x1,
        exp2 = x2
      )
  }
  
  pvalues <- pvalues |>
    mutate(
      x = list(c(rep(x1, 2), (x1 + x2) / 2, rep(x2, 2))),
      y = list(c(max1, rep(data_range[2] + 0.1 * diff(data_range), 3), max2)),
      p.value = list(c(rep(NA_real_, 2), p.value, rep(NA_real_, 2))),
      comparison = paste(exp1, exp2, sep = '_')
    ) |>
    select(comparison, x, y, p.value) |>
    pivot_wider(
      names_from = comparison,
      values_from = c(x, y, p.value),
      names_sep = '.'
    ) |>
    unnest(everything())
  
  write_tsv(pvalues, paste0(file, '_pvalues.tsv'), na = 'NaN')
}


## function to export a PWM logo
# PWM:      R object; the PWM as a universalmotif class object
# file:     character; name of the export file
# offset:   numeric; adjust the positions (motif position + offset)
# WT_seq:   character; sequence of the wild-type motif; if given wild-type bases are stored as capital letters; alternative bases are lowercase letters
# remove_WT:  logical; whether to remove the column with the WT nucleotides in the output file
PWM_to_LaTeX <- function(motif, file, offset = 0, WT_seq = NULL, remove_WT = FALSE) {
  if (class(motif) != 'universalmotif') {
    stop('The supplied `motif` is not a valid universalmotif class object')
  }
  
  if (motif['type'] != 'ICM') {
    motif <- convert_type(motif, 'ICM')
  }
  
  icm <- as_tibble(t(motif['motif'])) |>
    mutate(
      pos = seq_len(n()) + offset
    ) |>
    pivot_longer(
      -pos,
      names_to = 'base',
      values_to = 'IC'
    ) |>
    group_by(pos) |>
    arrange(IC, .by_group = TRUE) |>
    mutate(
      plot = seq_len(n())
    ) |>
    pivot_wider(
      names_from = plot,
      values_from = c(IC, base)
    )
  
  if (! is.null(WT_seq)) {
    icm <- icm |>
      bind_cols(
        WT = unlist(str_split(WT_seq, ''))
      ) |>
      mutate(
        across(starts_with('base'), ~ if_else(.x != WT, tolower(.x), .x))
      )
    
    if(remove_WT) {
      icm <- icm |>
        select(-WT)
    }
  }
  
  write_tsv(icm, file, na = 'NaN')
}


## annotate outliers (high outliers -> 1, low outliers -> -1, no outlier -> 0)
# data: numeric vector; data for which to annotate outliers
outliers <- function(data) {
  case_when(
    data > quantile(data, 0.75, na.rm = TRUE) + 1.5 * IQR(data, na.rm = TRUE) ~ 1,
    data < quantile(data, 0.25, na.rm = TRUE) - 1.5 * IQR(data, na.rm = TRUE) ~ -1,
    TRUE ~ 0
  )
}


## function to export data for histograms
# data:             R object; data frame to export
# values_from:      column specification; column containing the values to be summarised and exported
# file:             character; base name of the export file without extension (full name is: "<file>_hist.tsv")
# group:            column specification; column to group the data by (count values will be saved in colums named after the "group" contents)
# bins:             numeric; number of bins
# data_range:       numeric vector, length 2; min and max values used to determine bins and axis limits (the range of <data>$<values_from> is used if NULL)
# save_axis_limits: logical; whether to save axis limits to a file (file name: "<file>_axes.tsv") or to use them as the return value of the function
# value_axis:       character, "x" or "y"; axis along which the values are shown
# stacked:          logical; whether to stack the histograms of every group or plot them overlapping each other
LaTeX_histogram <- function(data, values_from, file, group = NULL, bins = 50, data_range = NULL, save_axis_limits = TRUE, value_axis = 'x', stacked = TRUE) {
  # get data range and calculate bin width and breaks
  if (is.null(data_range)) {
    data_range <- range(pull(data, {{values_from}}))
  } else if (length(data_range) != 2 | ! is.numeric(data_range)) {
    stop('"data_range" must be a two component numerical vector')
  } else if(! all(between(range(pull(data, {{values_from}})), data_range[1], data_range[2]))) {
    warning('The supplied "data_range" is narrower that the actual data range. Data outside of "data_range" will be ignored.')
  }
  
  bin_width <- diff(data_range) / bins
  breaks <- seq(data_range[1], data_range[2], bin_width)
  
  # assign axes
  if (value_axis == 'x') {
    count_axis  <-  'y'
  } else if (value_axis == 'y') {
    count_axis  <-  'x'
  } else {
    stop('"value_axis" must be either "x" or "y"')
  }
  
  # calculate histogram data
  histogram <- data |>
    select(value =  {{values_from}}, group = {{group}}) |>
    mutate(
      {{value_axis}} := cut(value, breaks, labels = FALSE, include.lowest = TRUE),
      across({{value_axis}}, ~ breaks[.x] + 0.5 * bin_width)
    ) |>
    count(across(-value), name = 'count')
  
  # calculate axis limits
  if (stacked) {
    hist_type <- paste0(count_axis, 'bar stacked')
    axis_limits <- histogram |>
      group_by(across(all_of(value_axis))) |>
      summarise(
        {{count_axis}} := sum(count)
      ) |>
      ungroup()
  } else {
    hist_type <- paste0(count_axis, 'bar')
    axis_limits <- histogram |>
      rename({{count_axis}} := count)
  }
  
  axis_limits <- axis_limits |>
    drop_na() |>
    reframe(
      across(all_of(count_axis), range),
      {{value_axis}} := data_range
    ) |>
    summarise(
      across(c(x, y), list('min' = min, 'max' = max), .names = '{.col}{.fn}'),
      "{count_axis}min" := 0,
      `bar width` = bin_width,
      `bar shift` = 0,
      `histogram type` = hist_type
    )
  
  # create columns for each group
  if (length(histogram) == 3) {
    histogram <- histogram |>
      pivot_wider(
        names_from = group,
        values_from = count,
        values_fill = 0
      ) |>
      drop_na({{value_axis}})
  } 
  
  # save histogram data to file
  write_tsv(histogram, paste0(file, '_hist.tsv'))
  
  # save or return axis limits
  if (save_axis_limits) {
    write_tsv(axis_limits, paste0(file, '_axes.tsv'))
  } else {
    return(axis_limits)
  }
}


### define styles for exported excel tables
xlsx_bold <- createStyle(textDecoration = 'bold')
xlsx_wrap <- createStyle(wrapText = TRUE)
xlsx_1digit <- createStyle(numFmt = '0.0')
xlsx_2digit <- createStyle(numFmt = '0.00')
xlsx_seq_font <- createStyle(fontName = 'Courier New')
xlsx_boldwrap <- createStyle(textDecoration = 'bold', wrapText = TRUE)
xlsx_center <- createStyle(halign = 'center')


### function to make the the name of a supplemental data set bold
# wb:   Workbook object; the workbook to format
# id:   numeric or character; the number of the current supplemental data set
makeNameBold <- function(wb, id) {
  # make the data set name bold
  for (i in grep(paste0('Supplemental Data Set S', id), wb$sharedStrings, fixed = TRUE)) {
    # insert additional formatting in shared string
    wb$sharedStrings[[i]] <- gsub('<si>', '<si><r>', gsub('</si>', '</r></si>', wb$sharedStrings[[i]]))
    wb$sharedStrings[[i]] <- gsub(
      '(Supplemental Data Set S[0-9]+)',
      '</t></r><r><rPr><b val=\"true\"/></rPr><t xml:space=\"preserve\">\\1</t></r><r><t xml:space=\"preserve\">',
      wb$sharedStrings[[i]]
    )
  }
}
