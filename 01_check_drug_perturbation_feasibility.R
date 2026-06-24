############################################################
# 01_check_drug_perturbation_feasibility.R
# Check drug perturbation feasibility using Enrichr/LINCS libraries
############################################################

# -----------------------------
# 1. Load packages
# -----------------------------

packages <- c(
  "httr", "jsonlite", "dplyr", "stringr",
  "readr", "purrr", "ggplot2", "openxlsx"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 2. Set project directory
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

out_dir <- file.path(project_dir, "results", "drug_perturbation")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:\n")
cat(project_dir, "\n\n")

# -----------------------------
# 3. Read and clean gene lists
# -----------------------------

gene_list_dir <- file.path(project_dir, "gene_lists")

active_file <- file.path(gene_list_dir, "active_TB_up_top150.txt")
post_file   <- file.path(gene_list_dir, "posttherapy_up_top150.txt")

if (!file.exists(active_file)) stop("Missing file: active_TB_up_top150.txt")
if (!file.exists(post_file)) stop("Missing file: posttherapy_up_top150.txt")

clean_gene_list <- function(file_path) {
  
  genes <- readLines(file_path, warn = FALSE)
  genes <- genes[!is.na(genes)]
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  genes <- toupper(genes)
  
  # Keep standard gene-symbol characters
  genes <- gsub("[^A-Z0-9._-]", "", genes)
  genes <- genes[genes != ""]
  genes <- unique(genes)
  
  return(genes)
}

active_genes <- clean_gene_list(active_file)
post_genes <- clean_gene_list(post_file)

cat("Active-TB-up genes loaded:", length(active_genes), "\n")
cat("Post-therapy-up genes loaded:", length(post_genes), "\n\n")

cat("First active-TB-up genes:\n")
print(head(active_genes, 10))

cat("\nFirst post-therapy-up genes:\n")
print(head(post_genes, 10))

# Save cleaned gene lists
write.table(
  active_genes,
  file.path(gene_list_dir, "active_TB_up_top150_cleaned.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  post_genes,
  file.path(gene_list_dir, "posttherapy_up_top150_cleaned.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

input_summary <- data.frame(
  Gene_Set = c("Active_TB_up_top150", "Posttherapy_up_top150"),
  Gene_Count = c(length(active_genes), length(post_genes))
)

write.csv(
  input_summary,
  file.path(table_dir, "input_gene_list_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# 4. Define Enrichr helper functions
# -----------------------------

enrichr_add_url <- "https://maayanlab.cloud/Enrichr/addList"
enrichr_enrich_url <- "https://maayanlab.cloud/Enrichr/enrich"
enrichr_stats_url <- "https://maayanlab.cloud/Enrichr/datasetStatistics"

submit_gene_list <- function(genes, description, max_tries = 5) {
  
  genes <- genes[!is.na(genes)]
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  genes <- unique(genes)
  
  genes_string <- paste(genes, collapse = "\n")
  
  for (attempt in 1:max_tries) {
    
    cat("Submitting:", description, "| attempt", attempt, "\n")
    
    response <- httr::POST(
      enrichr_add_url,
      body = list(
        list = genes_string,
        description = description
      ),
      encode = "multipart",
      httr::timeout(60)
    )
    
    status <- httr::status_code(response)
    response_text <- httr::content(response, as = "text", encoding = "UTF-8")
    
    if (status == 200) {
      result <- jsonlite::fromJSON(response_text)
      cat("Submission successful:", description, "\n")
      return(result$userListId)
    }
    
    cat("Submission failed. Status code:", status, "\n")
    cat("Response text:\n")
    cat(response_text, "\n")
    
    Sys.sleep(5)
  }
  
  stop(paste("Failed to submit gene list after retries:", description))
}

run_enrichr <- function(user_list_id, library_name) {
  
  response <- httr::GET(
    enrichr_enrich_url,
    query = list(
      userListId = user_list_id,
      backgroundType = library_name
    ),
    httr::timeout(60)
  )
  
  if (httr::status_code(response) != 200) {
    warning(paste("Failed library:", library_name))
    return(NULL)
  }
  
  response_text <- httr::content(response, as = "text", encoding = "UTF-8")
  result <- jsonlite::fromJSON(response_text, simplifyVector = FALSE)
  
  if (!library_name %in% names(result)) {
    warning(paste("No result for:", library_name))
    return(NULL)
  }
  
  raw_results <- result[[library_name]]
  
  if (length(raw_results) == 0) {
    return(NULL)
  }
  
  df_list <- lapply(raw_results, function(x) {
    
    overlapping_genes <- x[[6]]
    
    if (length(overlapping_genes) > 1) {
      overlapping_genes <- paste(unlist(overlapping_genes), collapse = ";")
    } else {
      overlapping_genes <- as.character(overlapping_genes)
    }
    
    data.frame(
      Rank = as.numeric(x[[1]]),
      Term = as.character(x[[2]]),
      P_value = as.numeric(x[[3]]),
      Z_score = as.numeric(x[[4]]),
      Combined_score = as.numeric(x[[5]]),
      Overlapping_genes = overlapping_genes,
      Adjusted_P_value = as.numeric(x[[7]]),
      Old_P_value = as.numeric(x[[8]]),
      Old_Adjusted_P_value = as.numeric(x[[9]]),
      stringsAsFactors = FALSE
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  return(df)
}

# -----------------------------
# 5. Check available Enrichr libraries
# -----------------------------

cat("\nChecking available Enrichr libraries...\n")

stats_response <- httr::GET(enrichr_stats_url, httr::timeout(60))

if (httr::status_code(stats_response) != 200) {
  stop("Could not access Enrichr dataset statistics.")
}

stats_json <- jsonlite::fromJSON(
  httr::content(stats_response, as = "text", encoding = "UTF-8")
)

library_table <- as.data.frame(stats_json$statistics)

write.csv(
  library_table,
  file.path(table_dir, "Enrichr_available_libraries.csv"),
  row.names = FALSE
)

drug_related_libraries <- library_table %>%
  filter(str_detect(
    libraryName,
    regex("LINCS|L1000|Drug|Perturb|CMap|DSigDB", ignore_case = TRUE)
  )) %>%
  select(libraryName, numTerms)

write.csv(
  drug_related_libraries,
  file.path(table_dir, "Enrichr_drug_related_libraries.csv"),
  row.names = FALSE
)

cat("\nDrug-related libraries found:\n")
print(drug_related_libraries)

# -----------------------------
# 6. Select libraries for analysis
# -----------------------------

candidate_libraries <- c(
  "LINCS_L1000_Chem_Pert_down",
  "LINCS_L1000_Chem_Pert_up",
  "Drug_Perturbations_from_GEO_down",
  "Drug_Perturbations_from_GEO_up",
  "DSigDB"
)

available_libraries <- library_table$libraryName
candidate_libraries <- candidate_libraries[candidate_libraries %in% available_libraries]

cat("\nLibraries selected for analysis:\n")
print(candidate_libraries)

if (length(candidate_libraries) == 0) {
  stop("None of the selected drug perturbation libraries are available.")
}

write.csv(
  data.frame(Selected_Library = candidate_libraries),
  file.path(table_dir, "selected_drug_perturbation_libraries.csv"),
  row.names = FALSE
)

# -----------------------------
# 7. Submit gene lists
# -----------------------------

cat("\nSubmitting gene lists to Enrichr...\n")

active_id <- submit_gene_list(
  active_genes,
  "Active_TB_up_genes"
)

Sys.sleep(5)

post_id <- submit_gene_list(
  post_genes,
  "Posttherapy_up_genes"
)

list_ids <- data.frame(
  Gene_Set = c("Active_TB_up", "Posttherapy_up"),
  Enrichr_UserListID = c(active_id, post_id)
)

write.csv(
  list_ids,
  file.path(log_dir, "Enrichr_user_list_ids.csv"),
  row.names = FALSE
)

cat("\nActive-TB list ID:", active_id, "\n")
cat("Post-therapy list ID:", post_id, "\n\n")

# -----------------------------
# 8. Run enrichment against selected libraries
# -----------------------------

all_results <- list()

for (lib in candidate_libraries) {
  
  cat("Running library:", lib, "\n")
  
  active_result <- run_enrichr(active_id, lib)
  
  if (!is.null(active_result)) {
    active_result$Query <- "Active_TB_up"
    active_result$Library <- lib
    all_results[[paste0("Active_", lib)]] <- active_result
  }
  
  Sys.sleep(2)
  
  post_result <- run_enrichr(post_id, lib)
  
  if (!is.null(post_result)) {
    post_result$Query <- "Posttherapy_up"
    post_result$Library <- lib
    all_results[[paste0("Post_", lib)]] <- post_result
  }
  
  Sys.sleep(2)
}

if (length(all_results) == 0) {
  stop("No enrichment results returned.")
}

all_enrichr_results <- bind_rows(all_results)

write.csv(
  all_enrichr_results,
  file.path(table_dir, "all_drug_perturbation_results.csv"),
  row.names = FALSE
)

saveRDS(
  all_enrichr_results,
  file.path(table_dir, "all_drug_perturbation_results.rds")
)

cat("\nAll Enrichr results saved.\n")
cat("Total result rows:", nrow(all_enrichr_results), "\n")

# -----------------------------
# 9. Filter significant terms
# -----------------------------

significant_results <- all_enrichr_results %>%
  filter(Adjusted_P_value < 0.05) %>%
  arrange(Adjusted_P_value, desc(Combined_score))

write.csv(
  significant_results,
  file.path(table_dir, "significant_drug_perturbation_results.csv"),
  row.names = FALSE
)

saveRDS(
  significant_results,
  file.path(table_dir, "significant_drug_perturbation_results.rds")
)

cat("\nSignificant drug perturbation terms:", nrow(significant_results), "\n")

# -----------------------------
# 10. Feasibility summary
# -----------------------------

if (nrow(significant_results) > 0) {
  
  feasibility_summary <- significant_results %>%
    group_by(Library, Query) %>%
    summarise(
      Significant_Terms = n(),
      Best_Adjusted_P = min(Adjusted_P_value, na.rm = TRUE),
      Max_Combined_Score = max(Combined_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Library, Query)
  
} else {
  
  feasibility_summary <- data.frame(
    Library = character(),
    Query = character(),
    Significant_Terms = integer(),
    Best_Adjusted_P = numeric(),
    Max_Combined_Score = numeric()
  )
}

write.csv(
  feasibility_summary,
  file.path(table_dir, "drug_perturbation_feasibility_summary.csv"),
  row.names = FALSE
)

cat("\nFeasibility summary:\n")
print(feasibility_summary)

# -----------------------------
# 11. Extract reversal-related hits
# -----------------------------

# Reversal logic:
# Active-TB-up genes matching drug DOWN signatures = active program reduction
# Post-therapy-up genes matching drug UP signatures = post-therapy mimicry

active_down_hits <- significant_results %>%
  filter(
    Query == "Active_TB_up",
    str_detect(Library, regex("down", ignore_case = TRUE))
  ) %>%
  arrange(Adjusted_P_value, desc(Combined_score))

posttherapy_up_hits <- significant_results %>%
  filter(
    Query == "Posttherapy_up",
    str_detect(Library, regex("up", ignore_case = TRUE))
  ) %>%
  arrange(Adjusted_P_value, desc(Combined_score))

write.csv(
  active_down_hits,
  file.path(table_dir, "active_TB_up_matching_drug_DOWN_signatures.csv"),
  row.names = FALSE
)

write.csv(
  posttherapy_up_hits,
  file.path(table_dir, "posttherapy_up_matching_drug_UP_signatures.csv"),
  row.names = FALSE
)

cat("\nActive-TB-up genes matching drug DOWN signatures:", nrow(active_down_hits), "\n")
cat("Post-therapy-up genes matching drug UP signatures:", nrow(posttherapy_up_hits), "\n")

# -----------------------------
# 12. Find bidirectional reversal terms
# -----------------------------

clean_term <- function(x) {
  x <- tolower(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  return(x)
}

active_down_hits$Term_Clean <- clean_term(active_down_hits$Term)
posttherapy_up_hits$Term_Clean <- clean_term(posttherapy_up_hits$Term)

common_clean_terms <- intersect(
  active_down_hits$Term_Clean,
  posttherapy_up_hits$Term_Clean
)

bidirectional_reversal_terms <- data.frame(
  Term_Clean = common_clean_terms
)

write.csv(
  bidirectional_reversal_terms,
  file.path(table_dir, "bidirectional_reversal_terms_clean_match.csv"),
  row.names = FALSE
)

cat("\nBidirectional clean-match reversal terms:", length(common_clean_terms), "\n")

if (length(common_clean_terms) > 0) {
  print(head(common_clean_terms, 20))
}

# -----------------------------
# 13. Save top 30 review tables
# -----------------------------

top30_active_down <- active_down_hits %>%
  select(Query, Library, Term, Adjusted_P_value, Combined_score, Overlapping_genes) %>%
  head(30)

top30_post_up <- posttherapy_up_hits %>%
  select(Query, Library, Term, Adjusted_P_value, Combined_score, Overlapping_genes) %>%
  head(30)

top30_all_significant <- significant_results %>%
  select(Query, Library, Term, Adjusted_P_value, Combined_score, Overlapping_genes) %>%
  head(30)

write.csv(
  top30_active_down,
  file.path(table_dir, "TOP30_active_TB_reversal_hits.csv"),
  row.names = FALSE
)

write.csv(
  top30_post_up,
  file.path(table_dir, "TOP30_posttherapy_mimic_hits.csv"),
  row.names = FALSE
)

write.csv(
  top30_all_significant,
  file.path(table_dir, "TOP30_all_significant_drug_perturbation_hits.csv"),
  row.names = FALSE
)

cat("\nTop 30 review tables saved.\n")

# -----------------------------
# 14. Decision summary
# -----------------------------

decision_summary <- data.frame(
  Metric = c(
    "Total significant drug perturbation terms",
    "Active-TB-up genes matching drug DOWN signatures",
    "Post-therapy-up genes matching drug UP signatures",
    "Bidirectional clean-match reversal terms"
  ),
  Value = c(
    nrow(significant_results),
    nrow(active_down_hits),
    nrow(posttherapy_up_hits),
    length(common_clean_terms)
  )
)

write.csv(
  decision_summary,
  file.path(log_dir, "decision_summary.csv"),
  row.names = FALSE
)

cat("\nDecision summary:\n")
print(decision_summary)

# -----------------------------
# 15. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "TB_drug_perturbation_feasibility_results.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "Input_summary")
writeData(wb, "Input_summary", input_summary)

addWorksheet(wb, "Selected_libraries")
writeData(wb, "Selected_libraries", data.frame(Selected_Library = candidate_libraries))

addWorksheet(wb, "Feasibility_summary")
writeData(wb, "Feasibility_summary", feasibility_summary)

addWorksheet(wb, "Decision_summary")
writeData(wb, "Decision_summary", decision_summary)

addWorksheet(wb, "Top30_active_reversal")
writeData(wb, "Top30_active_reversal", top30_active_down)

addWorksheet(wb, "Top30_posttherapy_mimic")
writeData(wb, "Top30_posttherapy_mimic", top30_post_up)

addWorksheet(wb, "Top30_all")
writeData(wb, "Top30_all", top30_all_significant)

addWorksheet(wb, "Bidirectional_terms")
writeData(wb, "Bidirectional_terms", bidirectional_reversal_terms)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

cat("\nExcel workbook saved:\n")
cat(xlsx_file, "\n")

# -----------------------------
# 16. Make summary figures
# -----------------------------

if (nrow(feasibility_summary) > 0) {
  
  p1 <- ggplot(
    feasibility_summary,
    aes(x = Library, y = Significant_Terms, fill = Query)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw(base_size = 12) +
    labs(
      title = "Drug perturbation feasibility summary",
      x = "Enrichr library",
      y = "Number of significant terms"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    )
  
  ggsave(
    filename = file.path(figure_dir, "drug_perturbation_feasibility_summary.png"),
    plot = p1,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(figure_dir, "drug_perturbation_feasibility_summary.pdf"),
    plot = p1,
    width = 9,
    height = 5
  )
}

if (nrow(significant_results) > 0) {
  
  plot_top <- significant_results %>%
    mutate(
      Label = paste0(Query, " | ", Library),
      MinusLog10AdjP = -log10(Adjusted_P_value)
    ) %>%
    group_by(Query, Library) %>%
    slice_min(Adjusted_P_value, n = 5, with_ties = FALSE) %>%
    ungroup()
  
  p2 <- ggplot(
    plot_top,
    aes(x = reorder(Term, MinusLog10AdjP), y = MinusLog10AdjP, fill = Query)
  ) +
    geom_col() +
    coord_flip() +
    facet_wrap(~Library, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      title = "Top drug perturbation terms",
      x = "",
      y = "-log10 adjusted p-value"
    ) +
    theme(
      legend.title = element_blank(),
      strip.text = element_text(size = 8)
    )
  
  ggsave(
    filename = file.path(figure_dir, "top_drug_perturbation_terms.png"),
    plot = p2,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(figure_dir, "top_drug_perturbation_terms.pdf"),
    plot = p2,
    width = 12,
    height = 8
  )
}

# -----------------------------
# 17. Save session information
# -----------------------------

sink(file.path(log_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\nInterpretation guide:\n")
cat(">100 significant terms = strong feasibility\n")
cat("50-100 significant terms = moderate feasibility\n")
cat("<50 significant terms = weak feasibility for signature-reversal QSAR\n")

cat("\nStep completed successfully.\n")
cat("All tables saved in:\n")
cat(table_dir, "\n\n")
cat("All figures saved in:\n")
cat(figure_dir, "\n\n")
cat("Excel workbook saved in:\n")
cat(xlsx_file, "\n")
############################################
# DRUG PERTURBATION VISUALIZATION
# Clean summary + publication-quality figures
############################################

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(tidyr)

# ==========================================
# 1. Project directory
# ==========================================

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

# Input files
active_file <- file.path(project_dir, "active_TB_up_matching_drug_DOWN_signatures.csv")
post_file   <- file.path(project_dir, "posttherapy_up_matching_drug_UP_signatures.csv")

# Optional files (not mandatory for plotting below)
allsig_file <- file.path(project_dir, "significant_drug_perturbation_results.csv")
libs_file   <- file.path(project_dir, "selected_drug_perturbation_libraries.csv")

# Output folders
out_table_dir <- file.path(project_dir, "results", "drug_perturbation", "tables")
out_fig_dir   <- file.path(project_dir, "results", "drug_perturbation", "figures")
out_log_dir   <- file.path(project_dir, "results", "drug_perturbation", "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 2. Check required files
# ==========================================

required_files <- c(active_file, post_file)

for (f in required_files) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}

# ==========================================
# 3. Helper functions
# ==========================================

find_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_library_name <- function(x) {
  recode(
    x,
    "Drug_Perturbations_from_GEO_down" = "GEO drug perturbations (down)",
    "Drug_Perturbations_from_GEO_up"   = "GEO drug perturbations (up)",
    "DSigDB"                           = "DSigDB",
    "LINCS_L1000_Chem_Pert_down"      = "LINCS chemical perturbations (down)",
    "LINCS_L1000_Chem_Pert_up"        = "LINCS chemical perturbations (up)",
    .default = x
  )
}

clean_term_name <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("&circ;", "^") %>%
    str_replace_all("&sup2;", "2") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

standardize_perturbation_df <- function(df, group_name) {
  
  term_col <- find_first_col(df, c(
    "Term", "term", "Description", "Geneset", "Gene_set"
  ))
  
  adj_col <- find_first_col(df, c(
    "Adjusted.P.value", "Adjusted P-value", "Adjusted_P_value",
    "adjusted_p_value", "adj.P.Val", "P.adjust", "padj"
  ))
  
  lib_col <- find_first_col(df, c(
    "Library", "library", "Gene_set_library", "GeneSetLibrary"
  ))
  
  overlap_col <- find_first_col(df, c(
    "Overlap", "overlap"
  ))
  
  combined_col <- find_first_col(df, c(
    "Combined.Score", "Combined Score", "combined_score"
  ))
  
  pval_col <- find_first_col(df, c(
    "P.value", "P-value", "P_Value", "p_value"
  ))
  
  if (is.na(term_col) || is.na(adj_col)) {
    stop("Could not detect required columns (Term and adjusted p-value). Please check column names.")
  }
  
  out <- data.frame(
    Term = as.character(df[[term_col]]),
    Adjusted_P = suppressWarnings(as.numeric(df[[adj_col]])),
    Library = if (!is.na(lib_col)) as.character(df[[lib_col]]) else "Unknown",
    Overlap = if (!is.na(overlap_col)) as.character(df[[overlap_col]]) else NA_character_,
    Combined_Score = if (!is.na(combined_col)) suppressWarnings(as.numeric(df[[combined_col]])) else NA_real_,
    P_Value = if (!is.na(pval_col)) suppressWarnings(as.numeric(df[[pval_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  out <- out %>%
    filter(!is.na(Term), !is.na(Adjusted_P)) %>%
    mutate(
      Signature_Group = group_name,
      Library_clean = clean_library_name(Library),
      Term_clean = clean_term_name(Term),
      MinusLog10AdjP = -log10(Adjusted_P)
    )
  
  out
}

save_plot_versions <- function(plot_obj, file_base, width = 12, height = 8) {
  
  # PNG: high DPI
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  # PDF: vector quality
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
}

# ==========================================
# 4. Read input files
# ==========================================

active_raw <- read_csv(active_file, show_col_types = FALSE)
post_raw   <- read_csv(post_file, show_col_types = FALSE)

active_df <- standardize_perturbation_df(active_raw, "Active_TB_up")
post_df   <- standardize_perturbation_df(post_raw, "Posttherapy_up")

all_df <- bind_rows(active_df, post_df)

# Keep only significant rows
all_sig <- all_df %>%
  filter(Adjusted_P <= 0.05)

# Save cleaned tables
write_csv(active_df, file.path(out_table_dir, "active_TB_up_cleaned_drug_perturbation.csv"))
write_csv(post_df,   file.path(out_table_dir, "posttherapy_up_cleaned_drug_perturbation.csv"))
write_csv(all_df,    file.path(out_table_dir, "all_cleaned_drug_perturbation_results.csv"))
write_csv(all_sig,   file.path(out_table_dir, "all_significant_drug_perturbation_results_cleaned.csv"))

# ==========================================
# 5. Summary counts by library
# ==========================================

library_order <- c(
  "GEO drug perturbations (down)",
  "GEO drug perturbations (up)",
  "DSigDB",
  "LINCS chemical perturbations (down)",
  "LINCS chemical perturbations (up)"
)

summary_counts <- all_sig %>%
  count(Signature_Group, Library_clean, name = "Significant_Terms") %>%
  complete(
    Signature_Group = c("Active_TB_up", "Posttherapy_up"),
    Library_clean = library_order,
    fill = list(Significant_Terms = 0)
  ) %>%
  mutate(
    Signature_Group = factor(Signature_Group,
                             levels = c("Active_TB_up", "Posttherapy_up")),
    Library_clean = factor(Library_clean, levels = library_order)
  )

write_csv(summary_counts, file.path(out_table_dir, "drug_perturbation_summary_counts_by_library.csv"))

# ==========================================
# 6. Figure 1: Feasibility summary barplot
# ==========================================

p_summary <- ggplot(summary_counts,
                    aes(x = Library_clean,
                        y = Significant_Terms,
                        fill = Signature_Group)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.68,
           color = "black",
           linewidth = 0.25) +
  geom_text(aes(label = Significant_Terms),
            position = position_dodge(width = 0.75),
            hjust = -0.15,
            size = 5.2,
            fontface = "bold") +
  coord_flip() +
  scale_fill_manual(
    values = c("Active_TB_up" = "#F8766D",
               "Posttherapy_up" = "#00BFC4"),
    labels = c("Active-TB-up signature", "Post-therapy-up signature")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Drug perturbation feasibility summary",
    x = NULL,
    y = "Number of significant terms",
    fill = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 15),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 18, 12, 12)
  )

save_plot_versions(
  p_summary,
  file.path(out_fig_dir, "Figure_1_drug_perturbation_feasibility_summary"),
  width = 12,
  height = 7.5
)

# ==========================================
# 7. Function for top-term plots
# ==========================================

make_top_terms_plot <- function(df, group_label, top_n = 5, fill_color = "#F8766D") {
  
  plot_df <- df %>%
    filter(Signature_Group == group_label, Adjusted_P <= 0.05) %>%
    group_by(Library_clean) %>%
    arrange(Adjusted_P, desc(MinusLog10AdjP)) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    filter(Library_clean %in% library_order) %>%
    mutate(
      Library_clean = factor(Library_clean, levels = library_order)
    ) %>%
    group_by(Library_clean) %>%
    arrange(MinusLog10AdjP, .by_group = TRUE) %>%
    mutate(
      Term_wrapped = str_wrap(Term_clean, width = 38),
      Plot_ID = paste0(Term_wrapped, "__", row_number())
    ) %>%
    ungroup()
  
  if (nrow(plot_df) == 0) return(NULL)
  
  write_csv(
    plot_df,
    file.path(out_table_dir,
              paste0("top_terms_", group_label, ".csv"))
  )
  
  p <- ggplot(plot_df,
              aes(x = fct_reorder(Plot_ID, MinusLog10AdjP),
                  y = MinusLog10AdjP)) +
    geom_col(fill = fill_color,
             color = "black",
             linewidth = 0.25,
             width = 0.72) +
    geom_text(aes(label = sprintf("%.2f", MinusLog10AdjP)),
              hjust = -0.08,
              size = 4.5,
              fontface = "bold") +
    coord_flip() +
    facet_wrap(~ Library_clean, scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = function(x) str_replace(x, "__\\d+$", "")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = ifelse(group_label == "Active_TB_up",
                     "Top active-TB-up drug perturbation terms",
                     "Top post-therapy-up drug perturbation terms"),
      x = NULL,
      y = expression(-log[10]("adjusted p-value"))
    ) +
    theme_bw(base_size = 17) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12.5),
      strip.text = element_text(size = 16, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  p
}

# ==========================================
# 8. Figure 2: Top terms for active-TB-up
# ==========================================

p_active <- make_top_terms_plot(
  df = all_df,
  group_label = "Active_TB_up",
  top_n = 5,
  fill_color = "#F8766D"
)

if (!is.null(p_active)) {
  save_plot_versions(
    p_active,
    file.path(out_fig_dir, "Figure_2_top_active_TB_up_drug_perturbation_terms"),
    width = 16,
    height = 11
  )
}

# ==========================================
# 9. Figure 3: Top terms for post-therapy-up
# ==========================================

p_post <- make_top_terms_plot(
  df = all_df,
  group_label = "Posttherapy_up",
  top_n = 5,
  fill_color = "#00BFC4"
)

if (!is.null(p_post)) {
  save_plot_versions(
    p_post,
    file.path(out_fig_dir, "Figure_3_top_posttherapy_up_drug_perturbation_terms"),
    width = 16,
    height = 10
  )
}

# ==========================================
# 10. Optional: Combined top-term table
# ==========================================

top_combined <- all_df %>%
  filter(Adjusted_P <= 0.05) %>%
  group_by(Signature_Group, Library_clean) %>%
  arrange(Adjusted_P, desc(MinusLog10AdjP)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(Signature_Group, Library_clean, Term_clean, Adjusted_P, MinusLog10AdjP, Overlap, Combined_Score)

write_csv(
  top_combined,
  file.path(out_table_dir, "top5_drug_perturbation_terms_by_group_and_library.csv")
)

# ==========================================
# 11. Save session info
# ==========================================

sink(file.path(out_log_dir, "session_info_drug_perturbation_visualization.txt"))
sessionInfo()
sink()

cat("\nDrug perturbation visualization completed successfully.\n")
cat("All cleaned tables saved in:\n", out_table_dir, "\n")
cat("All figures saved in:\n", out_fig_dir, "\n")
project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

all_csv <- list.files(
  project_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

matched_files <- all_csv[grepl("active_TB_up_matching_drug_DOWN_signatures", basename(all_csv))]

cat("Matched files:\n")
print(matched_files)