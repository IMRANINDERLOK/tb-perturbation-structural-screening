############################################################
# 09_prepare_CMap_CLUE_input_and_result_template.R
# Prepare CMap/CLUE query gene lists and result-entry template
############################################################

# -----------------------------
# 1. Load packages
# -----------------------------

packages <- c(
  "dplyr", "readr", "stringr", "openxlsx", "tidyr"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 2. Project directory
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

out_dir <- file.path(project_dir, "results", "CMap_CLUE_validation")
gene_dir <- file.path(out_dir, "gene_query_files")
table_dir <- file.path(out_dir, "tables")
template_dir <- file.path(out_dir, "templates")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:\n")
cat(project_dir, "\n\n")

# -----------------------------
# 3. Helper functions
# -----------------------------

find_file_recursive <- function(project_dir, file_pattern, required = TRUE) {
  
  matched <- list.files(
    project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matched <- matched[!grepl("CMap_CLUE_validation", matched)]
  
  if (length(matched) == 0 && required) {
    stop(paste("File not found:", file_pattern))
  }
  
  if (length(matched) == 0 && !required) {
    return(NA_character_)
  }
  
  return(matched[1])
}

clean_gene_list <- function(file_path) {
  
  genes <- readLines(file_path, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  genes <- toupper(genes)
  genes <- unique(genes)
  
  # keep only likely gene symbols
  genes <- genes[!grepl("\\s", genes)]
  genes <- genes[!is.na(genes)]
  
  return(genes)
}

normalize_name <- function(x) {
  
  x <- as.character(x)
  x <- tolower(x)
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "")
  x <- stringr::str_trim(x)
  
  return(x)
}

# -----------------------------
# 4. Find gene-list files
# -----------------------------

up_file <- find_file_recursive(
  project_dir,
  "CMap_UP_active_TB_top150\\.txt$",
  required = FALSE
)

if (is.na(up_file)) {
  up_file <- find_file_recursive(
    project_dir,
    "active_TB_up_top150\\.txt$",
    required = TRUE
  )
}

down_file <- find_file_recursive(
  project_dir,
  "CMap_DOWN_posttherapy_top150\\.txt$",
  required = FALSE
)

if (is.na(down_file)) {
  down_file <- find_file_recursive(
    project_dir,
    "posttherapy_up_top150\\.txt$",
    required = TRUE
  )
}

cat("Using CMap UP gene file:\n")
cat(up_file, "\n\n")

cat("Using CMap DOWN gene file:\n")
cat(down_file, "\n\n")

# -----------------------------
# 5. Read and clean gene lists
# -----------------------------

up_genes <- clean_gene_list(up_file)
down_genes <- clean_gene_list(down_file)

# Remove overlap if any gene appears in both lists
overlap_genes <- intersect(up_genes, down_genes)

if (length(overlap_genes) > 0) {
  cat("Overlapping genes found between UP and DOWN lists:\n")
  print(overlap_genes)
  
  up_genes <- setdiff(up_genes, overlap_genes)
  down_genes <- setdiff(down_genes, overlap_genes)
}

cat("Final UP genes:", length(up_genes), "\n")
cat("Final DOWN genes:", length(down_genes), "\n\n")

# -----------------------------
# 6. Save CMap/CLUE query files
# -----------------------------

write.table(
  up_genes,
  file.path(gene_dir, "CMap_query_UP_genes_active_TB_signature.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  down_genes,
  file.path(gene_dir, "CMap_query_DOWN_genes_posttherapy_signature.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Also save CSV versions
write_csv(
  data.frame(UP_Genes = up_genes),
  file.path(gene_dir, "CMap_query_UP_genes_active_TB_signature.csv")
)

write_csv(
  data.frame(DOWN_Genes = down_genes),
  file.path(gene_dir, "CMap_query_DOWN_genes_posttherapy_signature.csv")
)

# Save combined long-format query
query_long <- bind_rows(
  data.frame(Direction = "UP_in_active_TB", Gene = up_genes),
  data.frame(Direction = "DOWN_in_active_TB_or_UP_posttherapy", Gene = down_genes)
)

write_csv(
  query_long,
  file.path(gene_dir, "CMap_query_gene_signature_long_format.csv")
)

# -----------------------------
# 7. Gene query summary
# -----------------------------

gene_query_summary <- data.frame(
  Query_Component = c(
    "UP genes",
    "DOWN genes",
    "Removed overlapping genes",
    "UP source file",
    "DOWN source file"
  ),
  Description = c(
    "Genes elevated in active TB and used as the UP query.",
    "Genes relatively higher after therapy and used as the DOWN query.",
    "Genes appearing in both UP and DOWN query files were removed to avoid directional ambiguity.",
    up_file,
    down_file
  ),
  Count = c(
    length(up_genes),
    length(down_genes),
    length(overlap_genes),
    NA,
    NA
  )
)

write_csv(
  gene_query_summary,
  file.path(table_dir, "01_CMap_query_gene_summary.csv")
)

# -----------------------------
# 8. Find candidate files
# -----------------------------

all_candidate_file <- find_file_recursive(
  project_dir,
  "01_all_61_candidates_final_clean_labels\\.csv$",
  required = FALSE
)

if (is.na(all_candidate_file)) {
  all_candidate_file <- find_file_recursive(
    project_dir,
    "01_all_61_refined_candidate_perturbagens_final\\.csv$",
    required = TRUE
  )
}

validation_candidate_file <- find_file_recursive(
  project_dir,
  "02_candidates_for_validation_final_clean_labels\\.csv$",
  required = FALSE
)

if (is.na(validation_candidate_file)) {
  validation_candidate_file <- find_file_recursive(
    project_dir,
    "02_candidates_for_CMap_CTD_target_mapping\\.csv$",
    required = FALSE
  )
}

cat("Using all-candidate file:\n")
cat(all_candidate_file, "\n\n")

if (!is.na(validation_candidate_file)) {
  cat("Using validation-candidate file:\n")
  cat(validation_candidate_file, "\n\n")
}

# -----------------------------
# 9. Read candidate files
# -----------------------------

all_candidates <- readr::read_csv(all_candidate_file, show_col_types = FALSE)

if (!"Refined_Candidate_Name" %in% names(all_candidates)) {
  stop("Refined_Candidate_Name column missing in all-candidate file.")
}

if ("Display_Name" %in% names(all_candidates)) {
  candidate_name_col <- "Display_Name"
} else {
  candidate_name_col <- "Refined_Candidate_Name"
}

all_candidates <- all_candidates %>%
  mutate(
    Candidate_Display_Name = .data[[candidate_name_col]],
    Candidate_Match_Key = normalize_name(Candidate_Display_Name)
  )

if (!is.na(validation_candidate_file)) {
  
  validation_candidates <- readr::read_csv(validation_candidate_file, show_col_types = FALSE)
  
  if ("Display_Name" %in% names(validation_candidates)) {
    validation_name_col <- "Display_Name"
  } else {
    validation_name_col <- "Refined_Candidate_Name"
  }
  
  validation_candidates <- validation_candidates %>%
    mutate(
      Candidate_Display_Name = .data[[validation_name_col]],
      Candidate_Match_Key = normalize_name(Candidate_Display_Name)
    )
  
} else {
  
  validation_candidates <- all_candidates %>%
    filter(Candidate_Group %in% c(
      "Candidate compound",
      "Mapped LINCS compound",
      "Biologic / reference"
    ))
}

write_csv(
  all_candidates,
  file.path(table_dir, "02_all_candidate_universe_for_CMap_tracking.csv")
)

write_csv(
  validation_candidates,
  file.path(table_dir, "03_candidate_validation_set_for_CMap_tracking.csv")
)

# -----------------------------
# 10. Candidate name mapping template
# -----------------------------

candidate_mapping_template <- validation_candidates %>%
  transmute(
    Candidate_Display_Name,
    Candidate_Match_Key,
    Candidate_Group,
    Broad_Target_Axis = if ("Broad_Target_Axis_Clean" %in% names(validation_candidates)) {
      Broad_Target_Axis_Clean
    } else {
      Broad_Target_Axis
    },
    Suitability_Adjusted_Score,
    Supporting_Signatures,
    Max_Anchor_Hit_Count,
    Max_Core8_Hit_Count,
    Anchor_Hits,
    Core8_Hits,
    Possible_CMap_Name = "",
    CMap_Perturbagen_ID = "",
    Manual_Name_Check = "",
    Notes = ""
  ) %>%
  arrange(desc(Suitability_Adjusted_Score))

write_csv(
  candidate_mapping_template,
  file.path(template_dir, "CMap_candidate_name_mapping_template.csv")
)

# -----------------------------
# 11. CMap result manual-entry template
# -----------------------------

CMap_result_template <- data.frame(
  CMap_Result_Rank = integer(),
  CMap_Perturbagen_Name = character(),
  CMap_Perturbagen_ID = character(),
  CMap_Cell_ID = character(),
  CMap_Dose = character(),
  CMap_Time = character(),
  Tau = numeric(),
  Connectivity_Score = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  Matched_Candidate_Name = character(),
  Match_Type = character(),
  Reversal_Category = character(),
  Notes = character()
)

write_csv(
  CMap_result_template,
  file.path(template_dir, "CMap_result_manual_entry_template.csv")
)

# -----------------------------
# 12. Tau classification table
# -----------------------------

tau_classification <- data.frame(
  Tau_Category = c(
    "Strong reversal",
    "Moderate reversal",
    "Weak or no reversal",
    "Disease-mimicking"
  ),
  Criterion = c(
    "tau <= -90",
    "-90 < tau <= -70",
    "-70 < tau <= 0",
    "tau > 0"
  ),
  Interpretation = c(
    "Strong negative connectivity against the active-TB signature.",
    "Moderate negative connectivity and worth retaining for secondary validation.",
    "Weak reversal evidence; keep only if supported by other evidence.",
    "Positive connectivity suggests similarity to the active-TB signature and should not be prioritized."
  )
)

write_csv(
  tau_classification,
  file.path(table_dir, "04_CMap_tau_classification_criteria.csv")
)

# -----------------------------
# 13. Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "CMap_CLUE_query_and_result_templates.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "UP_genes")
writeData(wb, "UP_genes", data.frame(Gene = up_genes))

addWorksheet(wb, "DOWN_genes")
writeData(wb, "DOWN_genes", data.frame(Gene = down_genes))

addWorksheet(wb, "Query_summary")
writeData(wb, "Query_summary", gene_query_summary)

addWorksheet(wb, "Candidate_universe")
writeData(wb, "Candidate_universe", all_candidates)

addWorksheet(wb, "Validation_candidates")
writeData(wb, "Validation_candidates", validation_candidates)

addWorksheet(wb, "Name_mapping_template")
writeData(wb, "Name_mapping_template", candidate_mapping_template)

addWorksheet(wb, "CMap_result_entry")
writeData(wb, "CMap_result_entry", CMap_result_template)

addWorksheet(wb, "Tau_criteria")
writeData(wb, "Tau_criteria", tau_classification)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(wb)) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:100, gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:100, widths = "auto")
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# -----------------------------
# 14. Write CMap instructions
# -----------------------------

instructions <- c(
  "CMap/CLUE query preparation",
  "",
  "Use the following files for CLUE/CMap query:",
  "1. CMap_query_UP_genes_active_TB_signature.txt",
  "2. CMap_query_DOWN_genes_posttherapy_signature.txt",
  "",
  "Interpretation:",
  "UP genes are genes elevated in active TB.",
  "DOWN genes are genes relatively higher after therapy or lower in active TB.",
  "",
  "Candidate interpretation after CMap:",
  "tau <= -90: strong reversal",
  "-90 < tau <= -70: moderate reversal",
  "-70 < tau <= 0: weak or no reversal",
  "tau > 0: disease-mimicking or non-reversal direction",
  "",
  "After downloading CMap results, paste them into the CMap_result_entry sheet or save as a CSV.",
  "Then run the next integration script to merge tau scores with the 61-candidate table."
)

writeLines(
  instructions,
  file.path(log_dir, "CMap_CLUE_query_instructions.txt")
)

# -----------------------------
# 15. Method note
# -----------------------------

method_note <- data.frame(
  Item = c(
    "CMap query purpose",
    "UP query",
    "DOWN query",
    "Candidate universe",
    "Tau interpretation",
    "Next script"
  ),
  Description = c(
    "The CMap/CLUE query is used to test whether candidate perturbagens show negative connectivity against the active-TB host-response signature.",
    "The UP query contains genes elevated in active TB.",
    "The DOWN query contains genes relatively higher after therapy, representing the opposite direction of the active-TB state.",
    "The 61 refined perturbagens are retained as the candidate universe and tracked against CMap perturbagen names.",
    "Negative tau indicates potential reversal of the active-TB signature, whereas positive tau indicates similarity to the active-TB signature.",
    "After obtaining CMap results, run 10_integrate_CMap_tau_with_candidate_scores.R."
  )
)

write_csv(
  method_note,
  file.path(log_dir, "CMap_CLUE_preparation_method_note.csv")
)

sink(file.path(log_dir, "session_info_CMap_CLUE_preparation.txt"))
sessionInfo()
sink()

# -----------------------------
# 16. Console output
# -----------------------------

cat("\nCMap/CLUE preparation completed successfully.\n\n")

cat("UP gene query file:\n")
cat(file.path(gene_dir, "CMap_query_UP_genes_active_TB_signature.txt"), "\n\n")

cat("DOWN gene query file:\n")
cat(file.path(gene_dir, "CMap_query_DOWN_genes_posttherapy_signature.txt"), "\n\n")

cat("Candidate mapping template:\n")
cat(file.path(template_dir, "CMap_candidate_name_mapping_template.csv"), "\n\n")

cat("CMap result entry template:\n")
cat(file.path(template_dir, "CMap_result_manual_entry_template.csv"), "\n\n")

cat("Excel workbook:\n")
cat(xlsx_file, "\n\n")

cat("Gene query summary:\n")
print(gene_query_summary)