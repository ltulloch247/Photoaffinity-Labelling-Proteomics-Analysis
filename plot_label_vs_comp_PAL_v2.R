# Load packages
library(tools)
library(readxl)
library(dplyr)
library(plotly)
library(writexl)

# Define analysis parameters
data_path <- "C:/Users/path-to-data"                              # Replace path with your owb path. Use / not \
pulldown_file <- "DNDI-6690-AxAmPAL-rep1.csv"                     # Include extension (.xlsx/.csv)
pulldown_sheet <- ""                                              # Leave as "" if only one sheet
highlight_proteins_red <- "A0A3Q8IRJ0"                            # Protein IDs or "kinases" (pre-processed data only)
highlight_proteins_blue <- ""                                     # Protein IDs or "host kinases" (pre-processed data only)
highlight_interactors <- FALSE                                    # TRUE = yes, FALSE = no (pre-processed data only)
remove_contaminants <- TRUE                                       # TRUE = yes, FALSE = no
remove_low_peptides <- TRUE                                       # TRUE = yes, FALSE = no
add_imputation <- FALSE                                           # TRUE = yes, FALSE = no
add_imputation_DMSO_only <- TRUE                                  # TRUE = yes, FALSE = no
add_hit_threshold <- TRUE                                         # TRUE = yes, FALSE = no
hit_threshold_method_Probe_DMSO <- "median+4MAD"                  # Choose format: "2", "mean+2SD", "median+3MAD"
hit_threshold_method_Probe_ProbeComp <- "median+4MAD"             # Choose format: "2", "mean+2SD", "median+3MAD"
export_results <- TRUE                                            # TRUE = yes, FALSE = no
results_filename <- "DNDI-6690-AxAmPAL-rep1_result.xlsx"          # Include .xlsx extension

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          PART 1: FORMAT DATASET                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load data (supports .xlsx and .csv)
setwd(data_path)
ext <- file_ext(pulldown_file)
pulldown <- switch(
  tolower(ext),
  "xlsx" = if (pulldown_sheet == "") {read_excel(pulldown_file)} else 
  {read_excel(pulldown_file, sheet = pulldown_sheet)},
  "csv"  = read.csv(pulldown_file, stringsAsFactors = FALSE))

# Rename columns
pulldown <- pulldown %>%
  rename(any_of(c(
    Name = "Name",
    Protein_ID = "Protein_ID",
    Gene_ID = "Gene_ID",          # Use "" if no gene information
    Peptides = "Peptides",
    Intensity = "Intensity",      # Use "" if no intensity column - it will be calculated automatically
    DMSO = "DMSO",
    Probe = "Probe",
    ProbeComp = "ProbeComp"
  )))

# Define graph limits:
xmin <- -6    # add limits or leave "" for autofit
xmax <- 6     # add limits or leave "" for autofit
ymin <- -6    # add limits or leave "" for autofit
ymax <- 12    # add limits or leave "" for autofit

# Calculate graph limits
xmin_plot <- if (xmin == "") min(pulldown$log2FC_Probe_ProbeComp, na.rm = TRUE) else xmin
xmax_plot <- if (xmax == "") max(pulldown$log2FC_Probe_ProbeComp, na.rm = TRUE) else xmax
ymin_plot <- if (ymin == "") min(pulldown$log2FC_Probe_DMSO, na.rm = TRUE) else ymin
ymax_plot <- if (ymax == "") max(pulldown$log2FC_Probe_DMSO, na.rm = TRUE) else ymax

# Clean and prepare data
pulldown <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO", "Probe", "ProbeComp", "Intensity", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
pulldown <- pulldown %>% dplyr::mutate(across(any_of(c("DMSO", "Probe", "ProbeComp")), as.numeric))

# Determine data type
cols <- intersect(c("DMSO", "Probe", "ProbeComp"), colnames(pulldown))
if (length(cols) > 0) {
  vals <- unlist(pulldown[cols])
  max_val <- max(vals, na.rm = TRUE)
  data_type <- if (max_val < 1e5) "normalised" else "raw"
}

# Create function to handle NA/NaN values
is.any.na <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  is.na(x) | (is.character(x) & tolower(x) %in% c("na", "nan", "NA", "NaN"))
}

# Remove any previous imputation
remove_imputation <- TRUE  # TRUE = yes, FALSE = no
if (remove_imputation && length(cols) > 0 && all(!is.any.na(pulldown[cols]))) {
  val_remove <- case_when(
    any(vals == 11111, na.rm = TRUE) ~ 11111,
    data_type == "normalised" & max_val > 2000 ~ max_val,
    TRUE ~ min(vals, na.rm = TRUE)
  )
  pulldown <- pulldown %>%
    dplyr::mutate(across(all_of(cols), ~ na_if(., val_remove)))
}

# Calculate max intensity
if (!"Intensity" %in% colnames(pulldown) && length(cols) > 0) {
  pulldown$Intensity <- do.call(pmax, c(pulldown[cols], na.rm = TRUE))
}

# Freeze raw data for hit threshold calculations
raw_pulldown <- pulldown

# Remove contaminants
if (remove_contaminants) {
  pulldown <- dplyr::filter(pulldown, !grepl("cRAP", Protein_ID))
  raw_pulldown <- dplyr::filter(raw_pulldown, !grepl("cRAP", Protein_ID))
}

# Remove proteins with low peptide counts (≤ 2)
if (remove_low_peptides) {
  pulldown <- pulldown %>% dplyr::filter(Peptides > 2)
  raw_pulldown <- raw_pulldown %>% dplyr::filter(Peptides > 2)
}

# Impute missing values
impute_val <- FALSE
if (add_imputation) {
  if (length(cols) > 0) {
    vals <- unlist(pulldown[cols])
    impute_val <- min(vals, na.rm = TRUE) * 0.999999
    pulldown <- pulldown %>%
      dplyr::mutate(across(all_of(cols), ~ replace(., is.any.na(.), impute_val)))
  }
} else if (add_imputation_DMSO_only) {
  if (length(cols) > 0) {
    vals <- unlist(pulldown$DMSO)
    impute_val_DMSO_only <- min(vals, na.rm = TRUE) * 0.999999
    pulldown <- pulldown %>%
      dplyr::mutate(across(DMSO, ~ replace(., is.any.na(.), impute_val_DMSO_only)))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                        PART 2: PLOT ENRICHMENT                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Calculate log2 FC for Probe/DMSO and Probe/Probe+Comp and filter zero values
pulldown$log2FC_Probe_DMSO <- log2(pulldown$Probe / pulldown$DMSO)
pulldown$log2FC_Probe_ProbeComp <- log2(pulldown$Probe / pulldown$ProbeComp)
pulldown <- pulldown %>% dplyr::filter(log2FC_Probe_DMSO != 0 | log2FC_Probe_ProbeComp != 0)

# Calculate log2 FC values from observed values only
raw_pulldown$log2FC_Probe_DMSO <- log2(raw_pulldown$Probe / raw_pulldown$DMSO)
raw_pulldown$log2FC_Probe_ProbeComp <- log2(raw_pulldown$Probe / raw_pulldown$ProbeComp)
raw_pulldown <- raw_pulldown %>% dplyr::filter(log2FC_Probe_DMSO != 0 | log2FC_Probe_ProbeComp != 0)

# Reformat highlighted protein list
highlight_kinases_red <- highlight_proteins_red %in% c("kinases", "kinase")
highlight_host_kinases_blue <- highlight_proteins_blue %in% c("host kinases", "host kinase")
highlight_IDs_red <- if (highlight_kinases_red) { character(0) } else { unlist(strsplit(highlight_proteins_red, "\\s*,\\s*")) }
highlight_IDs_blue <- if (highlight_host_kinases_blue) { character(0) } else { unlist(strsplit(highlight_proteins_blue, "\\s*,\\s*")) }

# Define visual styling
pulldown <- pulldown %>%
  plotly::mutate(Imputed = if (add_imputation) { 
    pulldown$DMSO == impute_val | pulldown$Probe == impute_val | pulldown$ProbeComp == impute_val
  } else if (add_imputation_DMSO_only) { 
    pulldown$DMSO == impute_val_DMSO_only 
  } else { FALSE },
  Contaminant = grepl("^cRAP", Protein_ID),
  LowPeptide = Peptides <= 2,
  Shape = case_when(
    Imputed & LowPeptide ~ "diamond-open",
    Imputed & !LowPeptide ~ "diamond",
    !Imputed & LowPeptide ~ "circle-open",
    TRUE ~ "circle"
  ),
  Colour = case_when(
    if ("Kinase" %in% colnames(pulldown)) { highlight_kinases_red & Kinase ~ "red" },
    if ("Host_kinase" %in% colnames(pulldown)) { highlight_host_kinases_blue & Host_kinase ~ "blue" },
    if ("Kinase_interactor" %in% colnames(pulldown) & highlight_interactors) { highlight_kinases_red & Kinase_interactor ~ "goldenrod" },
    if ("Host_kinase_interactor" %in% colnames(pulldown) & highlight_interactors) { highlight_host_kinases_blue & Host_kinase_interactor ~ "green" },
    Protein_ID %in% highlight_IDs_red ~ "red",
    Protein_ID %in% highlight_IDs_blue ~ "blue",
    Contaminant ~ "grey",
    TRUE ~ "black"
  ),
  StyleGroup = paste(Shape, Colour, sep = "_")
  )

# Create a function to define the hit threshold
parse_hit_threshold <- function(threshold_input, values) {
  if (grepl("^median([+-])\\d*\\.?\\d*MAD$", threshold_input)) {
    matches <- regmatches(threshold_input, regexec("^median([+-])(\\d*\\.?\\d*)MAD$", threshold_input))[[1]]
    operator <- matches[2]
    multiplier <- as.numeric(matches[3])
    med_val <- median(values, na.rm = TRUE)
    mad_val <- mad(values, constant = 1, na.rm = TRUE)
    if (operator == "+") {
      return(med_val + multiplier * mad_val)
    } else {
      return(med_val - multiplier * mad_val)
    }
  } else if (grepl("^mean([+-])\\d*\\.?\\d*SD$", threshold_input)) {
    matches <- regmatches(threshold_input, regexec("^mean([+-])(\\d*\\.?\\d*)SD$", threshold_input))[[1]]
    operator <- matches[2]
    multiplier <- as.numeric(matches[3])
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    if (operator == "+") {
      return(mean_val + multiplier * sd_val)
    } else {
      return(mean_val - multiplier * sd_val)
    }
  } else if (!is.na(suppressWarnings(as.numeric(threshold_input)))) {
    return(as.numeric(threshold_input))
  }
}

# Calculate hit threshold
log2FC_threshold_Probe_DMSO <- FALSE
log2FC_Probe_ProbeComp_threshold <- FALSE
if (add_hit_threshold) {
  log2FC_threshold_Probe_DMSO <- parse_hit_threshold(hit_threshold_method_Probe_DMSO, raw_pulldown$log2FC_Probe_DMSO)
  log2FC_threshold_Probe_ProbeComp <- parse_hit_threshold(hit_threshold_method_Probe_ProbeComp, raw_pulldown$log2FC_Probe_ProbeComp)
}

# Initialise plot
fig <- plot_ly()

# Define legend labels
legend_labels <- c(
  "diamond_black" = "imputed",
  "circle_black" = "> 2 peptides",
  "circle-open_black" = "≤ 2 peptides",
  "circle_grey" = "contaminant",
  "circle_blue" = "highlighted",
  "circle_red" = "highlighted",
  "circle_goldenrod" = "kinase interactor",
  "circle_green" = "host kinase interactor"
)

# Define display order (top = foreground, bottom = background)
style_order <- c(
  "diamond-open_red",        # highlighted_red, ≤ 2 peptides, imputed
  "diamond_red",             # highlighted_red, > 2 peptides, imputed
  "circle-open_red",         # highlighted_red, ≤ 2 peptides
  "circle_red",              # highlighted_red, ≤ 2 peptides
  "diamond-open_blue",       # highlighted_blue, ≤ 2 peptides, imputed
  "diamond_blue",            # highlighted_blue, > 2 peptides, imputed
  "circle-open_blue",        # highlighted_blue, ≤ 2 peptides
  "circle_blue",             # highlighted_blue, ≤ 2 peptides
  "diamond-open_goldenrod",  # highlighted_goldenrod, ≤ 2 peptides, imputed
  "diamond_goldenrod",       # highlighted_goldenrod, > 2 peptides, imputed
  "circle-open_goldenrod",   # highlighted_goldenrod, ≤ 2 peptides
  "circle_goldenrod",        # highlighted_goldenrod, ≤ 2 peptides
  "diamond-open_green",      # highlighted_green, ≤ 2 peptides, imputed
  "diamond_green",           # highlighted_green, > 2 peptides, imputed
  "circle-open_green",       # highlighted_green, ≤ 2 peptides
  "circle_green",            # highlighted_green, ≤ 2 peptides
  "diamond-open_grey",       # contaminant, ≤ 2 peptides, imputed
  "diamond_grey",            # contaminant, > 2 peptides, imputed
  "circle-open_grey",        # contaminant, ≤ 2 peptides
  "circle_grey",             # contaminant, > 2 peptides
  "diamond-open_black",      # ≤ 2 peptides, imputed
  "diamond_black",           # > 2 peptides, imputed
  "circle-open_black",       # ≤ 2 peptides
  "circle_black"             # > 2 peptides
)

# Add hit threshold lines
if (add_hit_threshold) {
  fig <- fig %>%
    add_trace(
      x = c(log2FC_threshold_Probe_ProbeComp, xmax),
      y = c(log2FC_threshold_Probe_DMSO, log2FC_threshold_Probe_DMSO),
      type = "scatter",
      mode = "lines",
      line = list(color = "black", dash = "dash", width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    ) %>%
    add_trace(
      x = c(log2FC_threshold_Probe_ProbeComp, log2FC_threshold_Probe_ProbeComp),
      y = c(log2FC_threshold_Probe_DMSO, ymax),
      type = "scatter",
      mode = "lines",
      line = list(color = "black", dash = "dash", width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    )
}

# Add data points
for(group in rev(style_order)) {
  group_data <- pulldown %>% filter(StyleGroup == group)
  if (nrow(group_data) == 0) next
  
  parts <- strsplit(group, "_")[[1]]
  Shape <- parts[1]
  Color <- parts[2]
  
  marker_size <- ifelse(
    (Shape == "circle" && Color == "black") || 
      (Shape == "diamond" && Color == "black") || 
      (Shape == "circle" && Color == "red") ||
      (Shape == "diamond" && Color == "red") ||
      (Shape == "circle" && Color == "blue") ||
      (Shape == "diamond" && Color == "blue") ||
      (Shape == "circle" && Color == "goldenrod") ||
      (Shape == "diamond" && Color == "goldenrod") ||
      (Shape == "circle" && Color == "green") ||
      (Shape == "diamond" && Color == "green"),
    10, 8
  )
  
  legend_name <- legend_labels[group]
  show_legend <- !is.na(legend_name)
  
  if ("Gene_ID" %in% colnames(group_data)) {
    hover_text <- paste0(
      "Protein ID: ", group_data$Protein_ID,
      "<br>Gene ID: ", group_data$Gene_ID,
      "<br>Name: ", group_data$Name,
      "<br>Peptides: ", group_data$Peptides,
      "<br>log2FC label: ", signif(group_data$log2FC_Probe_DMSO, 3),
      "<br>log2FC comp: ", signif(group_data$log2FC_Probe_ProbeComp, 3)
    )
  } else {
    hover_text <- paste0(
      "Protein ID: ", group_data$Protein_ID,
      "<br>Name: ", group_data$Name,
      "<br>Peptides: ", group_data$Peptides,
      "<br>log2FC label: ", signif(group_data$log2FC_Probe_DMSO, 3),
      "<br>log2FC comp: ", signif(group_data$log2FC_Probe_ProbeComp, 3)
    )
  }
  
  fig <- fig %>%
    add_trace(
      data = group_data,
      x = ~log2FC_Probe_ProbeComp,
      y = ~log2FC_Probe_DMSO,
      type = 'scatter',
      mode = 'markers',
      name = ifelse(show_legend, legend_name, ""),
      showlegend = show_legend,
      marker = list(
        symbol = Shape,
        size = marker_size,
        color = Color,
        line = list(color = Color, width = 1)
      ),
      text = hover_text,
      hoverinfo = 'text'
    )
}

# Finalise layout
fig <- fig %>%
  layout(
    xaxis = list(
      title = list(
        text = "Competition (log<sub>2</sub> FC probe / probe + competitor)",
        font = list(size = 20)
      ),
      range = c(xmin, xmax), # <- Remove '#' sign before 'range' to set x-axis limits
      showline = TRUE,
      linecolor = "black",
      ticks = "outside",
      tick0 = 0,
      dtick = 2,
      ticklen = 5,
      tickwidth = 1,
      tickcolor = "black",
      tickfont = list(size = 18),
      showgrid = FALSE
    ),
    yaxis = list(
      title = list(
        text = "Labelling (log<sub>2</sub> FC probe / DMSO)",
        font = list(size = 20)
      ),
      range = c(ymin, ymax), # <- Remove '#' sign before 'range' to set x-axis limits
      showline = TRUE,
      linecolor = "black",
      ticks = "outside",
      tick0 = 0,
      dtick = 2,
      ticklen = 5,
      tickwidth = 1,
      tickcolor = "black",
      tickfont = list(size = 18),
      showgrid = FALSE
    ),
    legend = list(
      title = list(text = "Legend"),
      x = 1.05,
      y = 1,
      bordercolor = "grey",
      borderwidth = 1
    ),
    margin = list(l = 50, r = 50, t = 10, b = 60)
  )

# Display plot
print(fig)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                         PART 3: EXPORT RESULTS                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Format final dataset and create parameter dataframe
final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO", "Probe", "ProbeComp", "Intensity", "log2FC_Probe_DMSO", "log2FC_Probe_ProbeComp", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
parameters <- data.frame(
  Parameter = c("data path", "data file", "data sheet", "script", "remove contaminants", "remove low peptides", "add imputation (all samples)", "add imputation (DMSO only)", "imputed value (LOD)", "add hit threshold"),
  Value = if (add_imputation_DMSO_only) {
    c(data_path, pulldown_file, pulldown_sheet, "plot_label_vs_comp_PAL.R", remove_contaminants, remove_low_peptides, add_imputation, add_imputation_DMSO_only, impute_val_DMSO_only, add_hit_threshold)
  } else {
    c(data_path, pulldown_file, pulldown_sheet, "plot_label_vs_comp_PAL.R", remove_contaminants, remove_low_peptides, add_imputation, add_imputation_DMSO_only, impute_val, add_hit_threshold)
  },
  stringsAsFactors = FALSE
)

# Apply hit thresholds
if (add_hit_threshold) {
  pulldown$Hit <- FALSE
  pulldown$Hit <- pulldown$log2FC_Probe_DMSO > log2FC_threshold_Probe_DMSO & pulldown$log2FC_Probe_ProbeComp > log2FC_threshold_Probe_ProbeComp
  final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO", "Probe", "ProbeComp", "Intensity", "log2FC_Probe_DMSO", "log2FC_Probe_ProbeComp", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor", "Hit")))
  parameters <- data.frame(
    Parameter = c("data path", "data file", "data sheet", "script", "remove contaminants", "remove low peptides", "add imputation (all samples)", "add imputation (DMSO only)", "imputed value (LOD)", "add hit threshold", "hit threshold method Probe/DMSO", "hit threshold Probe/DMSO", "hit threshold method Probe/Probe+Comp", "hit threshold Probe/Probe+Comp"),
    Value = if (add_imputation_DMSO_only) {
      c(data_path, pulldown_file, pulldown_sheet, "plot_label_vs_comp_PAL.R", remove_contaminants, remove_low_peptides, add_imputation, add_imputation_DMSO_only, impute_val_DMSO_only, add_hit_threshold, hit_threshold_method_Probe_DMSO, log2FC_threshold_Probe_DMSO, hit_threshold_method_Probe_ProbeComp, log2FC_threshold_Probe_ProbeComp)
    } else {
      c(data_path, pulldown_file, pulldown_sheet, "plot_label_vs_comp_PAL.R", remove_contaminants, remove_low_peptides, add_imputation, add_imputation_DMSO_only, impute_val, add_hit_threshold, hit_threshold_method_Probe_DMSO, log2FC_threshold_Probe_DMSO, hit_threshold_method_Probe_ProbeComp, log2FC_threshold_Probe_ProbeComp)
    },
    stringsAsFactors = FALSE
  )
}

# Export to file
if (export_results) {
  write_xlsx(list("Results" = final_dataset, "Parameters" = parameters), results_filename)
}