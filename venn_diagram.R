# Load packages
library(readxl)
library(dplyr)
library(VennDiagram)
library(writexl)

# Define analysis parameters
data_path <- "C:/Users/path-to-data"                        # Replace path with your owb path. Use / not \
pulldown_file_1 <- "DNDI-6690-AxAmPAL-rep1_result.xlsx"     # Include .xlsx extension
pulldown_file_2 <- "DNDI-6690-AxAmPAL-rep2_result.xlsx"     # Include .xlsx extension
pulldown_file_3 <- "DNDI-6690-AxAmPAL-rep3_result.xlsx"     # Include .xlsx extension
pulldown_file_4 <- ""                                       # Include .xlsx extension
pulldown_file_5 <- ""                                       # Include .xlsx extension
mode <- "hits"                                              # "all" or "hits"
export_results <- TRUE                                      # TRUE = yes, FALSE = no
results_filename <- "DNDI-6690-AxAmPAL-common_hits.xlsx"    # Include .xlsx extension

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          PART 1: FORMAT DATASETS                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Remove empty entries
setwd(data_path)
files <- c(pulldown_file_1, pulldown_file_2, pulldown_file_3, pulldown_file_4, pulldown_file_5)
files <- files[files != ""]

# Load data (supports .xlsx only)
pulldowns <- lapply(seq_along(files), function(i) { read_excel(files[i]) })

# Select data according to analysis mode
if (mode == "hits") {
  analysis_data <- lapply(pulldowns, function(df) { filter(df, Hit == TRUE) })
} else if (mode == "all") {
  analysis_data <- pulldowns
}

# Create sets for each condition
condition_sets <- setNames(
  lapply(analysis_data, function(df) unique(df$Protein_ID)),
  paste0("dataset ", seq_along(analysis_data))
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          PART 2: GENERATE VENN DIAGRAM                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_sets <- length(condition_sets)
venn.plot <- venn.diagram(
  x = condition_sets,
  category.names = names(condition_sets),
  filename = NULL,
  output = TRUE,
  fill = grDevices::rainbow(n_sets),
  alpha = 0.5,
  lwd = 2,
  margin = 0.1,
  cat.cex = 1,
  cex = 1,
  fontface = "bold",
  fontfamily = "Arial",      # numbers
  cat.fontfamily = "Arial"   # labels
)

grid.newpage()
grid.draw(venn.plot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           PART 3: EXPORT RESULTS                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Columns to merge across datasets
merge_fields <- c("Name", "Gene_ID", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")

# Function to combine text fields
collapse_values <- function(x) {
  x <- as.character(x)
  vals <- unique(na.omit(x))
  if (length(vals) == 0) { return(NA_character_) }
  vals <- unique(unlist(strsplit(vals, "[,;]\\s*")))
  paste(vals, collapse = ";")
}

# Rename columns by dataset
renamed_datasets <- lapply(seq_along(analysis_data), function(i) {
  analysis_data[[i]] %>%
    rename_with(~ paste0(.x, "_dataset", i), -Protein_ID)
})

# Merge datasets
merged <- Reduce(
  function(x, y) full_join(x, y, by = "Protein_ID"),
  renamed_datasets
)

# Combine selected columns
for (field in merge_fields) {
  cols <- grep(paste0("^", field, "_dataset"), names(merged), value = TRUE)
  if (length(cols) > 0) {
    merged[[field]] <- apply(merged[, cols, drop = FALSE], 1, collapse_values)
    merged <- merged %>% dplyr::select(-all_of(cols))
  }
}

# Collapse duplicate Protein_ID entries
merged <- merged %>%
  group_by(Protein_ID) %>%
  summarise(across(where(is.character), collapse_values),
            across(where(is.numeric), ~ first(na.omit(.))),
            across(where(is.logical), ~ any(.x, na.rm = TRUE)),
            .groups = "drop"
            )

# Move key columns to front
priority_cols <- intersect(c("Name", "Protein_ID", "Gene_ID", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor"), names(merged))
merged <- merged %>% dplyr::select(all_of(priority_cols), everything())

# Export
if (export_results) {
  write_xlsx(merged, results_filename)
}
