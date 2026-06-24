############################################################
# 10B_integrate_CLUE_query_result_with_candidates.R
# Use CLUE query_result.gct for final NCS/FDR candidate integration
############################################################

packages <- c(
  "dplyr", "readr", "stringr", "tidyr",
  "openxlsx", "ggplot2", "forcats", "scales"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Paths
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

clue_query_result <- file.path(
  project_dir,
  "results",
  "CLUE",
  "arfs",
  "TAG",
  "query_result.gct"
)

candidate_file <- file.path(
  project_dir,
  "results",
  "final_publication_figures",
  "tables",
  "01_all_61_candidates_final_clean_labels.csv"
)

out_dir <- file.path(project_dir, "results", "CLUE_query_result_candidate_integration")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(clue_query_result)) {
  stop("query_result.gct not found. Please check the CLUE folder path.")
}

if (!file.exists(candidate_file)) {
  stop("Candidate table not found. Please check candidate_file path.")
}

cat("Using CLUE query result:\n")
cat(clue_query_result, "\n\n")

cat("Using candidate table:\n")
cat(candidate_file, "\n\n")

# -----------------------------
# 2. Helper functions
# -----------------------------

normalize_name <- function(x) {
  x <- tolower(as.character(x))
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "")
  x <- stringr::str_trim(x)
  return(x)
}

classify_ncs <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "No CMap match",
    x <= -1.50 ~ "Strong opposing NCS",
    x > -1.50 & x <= -1.00 ~ "Moderate opposing NCS",
    x > -1.00 & x <= 0 ~ "Weak opposing NCS",
    x > 0 ~ "Similar / disease-mimicking NCS",
    TRUE ~ "Unclassified"
  )
}

save_plot_versions <- function(plot_obj, file_base, width = 12, height = 8) {
  
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

# -----------------------------
# 3. Read GCT 1.3 file
# -----------------------------

read_gct_13 <- function(gct_file) {
  
  lines <- readLines(gct_file, warn = FALSE)
  version <- trimws(lines[1])
  
  if (version != "#1.3") {
    stop("This script expects GCT version #1.3.")
  }
  
  dims <- as.integer(strsplit(lines[2], "\t")[[1]])
  
  n_data_cols <- dims[2]
  n_col_meta  <- dims[4]
  
  header <- strsplit(lines[3], "\t")[[1]]
  data_start <- 4 + n_col_meta
  
  gct <- read.delim(
    gct_file,
    skip = data_start - 1,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  
  colnames(gct) <- header
  
  return(gct)
}

query_result <- read_gct_13(clue_query_result)

write_csv(
  query_result,
  file.path(table_dir, "01_CLUE_query_result_raw_parsed.csv")
)

cat("Rows in query_result:", nrow(query_result), "\n")
cat("Columns:\n")
print(colnames(query_result))

# -----------------------------
# 4. Check required CLUE columns
# -----------------------------

required_clue_cols <- c(
  "pert_id", "pert_iname", "cell_iname", "pert_type",
  "pert_idose", "pert_itime", "raw_cs",
  "fdr_q_nlog10", "norm_cs"
)

missing_cols <- setdiff(required_clue_cols, colnames(query_result))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

clue_clean <- query_result %>%
  mutate(
    raw_cs = suppressWarnings(as.numeric(raw_cs)),
    fdr_q_nlog10 = suppressWarnings(as.numeric(fdr_q_nlog10)),
    norm_cs = suppressWarnings(as.numeric(norm_cs)),
    CMap_Perturbagen_Name = pert_iname,
    CMap_Match_Key = normalize_name(CMap_Perturbagen_Name)
  ) %>%
  filter(!is.na(norm_cs))

write_csv(
  clue_clean,
  file.path(table_dir, "02_CLUE_query_result_clean_norm_cs.csv")
)

# -----------------------------
# 5. Save top opposing and similar CLUE perturbagens
# -----------------------------

top_opposing <- clue_clean %>%
  arrange(norm_cs, desc(fdr_q_nlog10)) %>%
  head(200)

top_similar <- clue_clean %>%
  arrange(desc(norm_cs), desc(fdr_q_nlog10)) %>%
  head(200)

write_csv(
  top_opposing,
  file.path(table_dir, "03_top200_opposing_CLUE_perturbagens_norm_cs.csv")
)

write_csv(
  top_similar,
  file.path(table_dir, "04_top200_similar_CLUE_perturbagens_norm_cs.csv")
)

# -----------------------------
# 6. Read 61-candidate table
# -----------------------------

cand <- read_csv(candidate_file, show_col_types = FALSE)

name_col <- if ("Display_Name" %in% colnames(cand)) {
  "Display_Name"
} else {
  "Refined_Candidate_Name"
}

cand <- cand %>%
  mutate(
    Candidate_Display_Name = .data[[name_col]],
    Candidate_Match_Key = normalize_name(Candidate_Display_Name)
  )

candidate_keys <- cand %>%
  select(Candidate_Display_Name, Candidate_Match_Key) %>%
  distinct()

# Add synonyms for better matching
extra_keys <- data.frame(
  Candidate_Display_Name = c(
    "1,25-dihydroxyvitamin D",
    "Vitamin D3",
    "Ursodeoxycholic acid",
    "PLX4720",
    "DMNQ",
    "PLX4032",
    "Cyclosporine"
  ),
  Candidate_Match_Key = normalize_name(c(
    "calcitriol",
    "cholecalciferol",
    "ursodiol",
    "plx-4720",
    "2,3-dimethoxy-1,4-naphthoquinone",
    "vemurafenib",
    "cyclosporin a"
  ))
)

candidate_keys <- bind_rows(candidate_keys, extra_keys) %>%
  distinct()

write_csv(
  candidate_keys,
  file.path(table_dir, "05_candidate_name_keys_used_for_matching.csv")
)

# -----------------------------
# 7. Match CLUE perturbagens with candidate list
# -----------------------------

matched_raw <- clue_clean %>%
  inner_join(
    candidate_keys,
    by = c("CMap_Match_Key" = "Candidate_Match_Key")
  )

write_csv(
  matched_raw,
  file.path(table_dir, "06_raw_candidate_matches_from_CLUE_query_result.csv")
)

candidate_best <- matched_raw %>%
  group_by(Candidate_Display_Name) %>%
  summarise(
    Best_Negative_NCS = min(norm_cs, na.rm = TRUE),
    Best_Positive_NCS = max(norm_cs, na.rm = TRUE),
    Mean_NCS = mean(norm_cs, na.rm = TRUE),
    Best_FDR_q_nlog10 = max(fdr_q_nlog10, na.rm = TRUE),
    Number_of_CLUE_Matches = n(),
    Best_Opposing_Perturbagen = CMap_Perturbagen_Name[which.min(norm_cs)][1],
    Best_Opposing_Cell = cell_iname[which.min(norm_cs)][1],
    Best_Opposing_Dose = pert_idose[which.min(norm_cs)][1],
    Best_Opposing_Time = pert_itime[which.min(norm_cs)][1],
    Best_Opposing_Pert_ID = pert_id[which.min(norm_cs)][1],
    .groups = "drop"
  ) %>%
  mutate(
    NCS_Category = classify_ncs(Best_Negative_NCS)
  )

write_csv(
  candidate_best,
  file.path(table_dir, "07_candidate_best_CLUE_NCS_summary.csv")
)

# -----------------------------
# 8. Integrate with full 61-candidate table
# -----------------------------

cand_integrated <- cand %>%
  left_join(candidate_best, by = "Candidate_Display_Name") %>%
  mutate(
    NCS_Category = ifelse(
      is.na(NCS_Category),
      "No CMap match",
      NCS_Category
    )
  ) %>%
  arrange(Best_Negative_NCS, desc(Suitability_Adjusted_Score))

write_csv(
  cand_integrated,
  file.path(table_dir, "08_all_61_candidates_integrated_with_CLUE_query_result.csv")
)

ncs_summary <- cand_integrated %>%
  count(NCS_Category, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  ncs_summary,
  file.path(table_dir, "09_NCS_category_summary.csv")
)

# -----------------------------
# 9. Create final CMap-supported candidate ranking
# -----------------------------

final_supported <- cand_integrated %>%
  filter(!is.na(Best_Negative_NCS)) %>%
  mutate(
    NCS_Reversal_Strength = pmax(0, -Best_Negative_NCS),
    NCS_Reversal_Scaled = scales::rescale(
      NCS_Reversal_Strength,
      to = c(0, 100),
      from = range(NCS_Reversal_Strength, na.rm = TRUE)
    ),
    Integrated_CLUE_Score = round(
      0.60 * NCS_Reversal_Scaled +
        0.25 * Suitability_Adjusted_Score +
        0.15 * scales::rescale(Best_FDR_q_nlog10, to = c(0, 100), na.rm = TRUE),
      2
    )
  ) %>%
  arrange(desc(Integrated_CLUE_Score), Best_Negative_NCS)

write_csv(
  final_supported,
  file.path(table_dir, "10_final_CMap_supported_candidate_ranking.csv")
)

# -----------------------------
# 10. Figures
# -----------------------------

top30_final <- final_supported %>%
  head(30) %>%
  mutate(
    Plot_Label = stringr::str_wrap(Candidate_Display_Name, width = 32),
    Plot_Label = forcats::fct_reorder(Plot_Label, Best_Negative_NCS)
  )

p_ncs <- ggplot(
  top30_final,
  aes(x = Plot_Label, y = Best_Negative_NCS, fill = NCS_Category)
) +
  geom_col(color = "black", linewidth = 0.3, width = 0.72) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Strong opposing NCS" = "#1B9E77",
      "Moderate opposing NCS" = "#66A61E",
      "Weak opposing NCS" = "#7570B3",
      "Similar / disease-mimicking NCS" = "#D95F02",
      "No CMap match" = "grey70"
    )
  ) +
  labs(
    title = "CMap-supported opposing perturbagens",
    x = NULL,
    y = "Best negative normalized connectivity score",
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, face = "plain"),
    axis.text.y = element_text(size = 11, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 35, 12)
  )

save_plot_versions(
  p_ncs,
  file.path(fig_dir, "Figure_1_CMap_supported_candidate_NCS"),
  width = 11,
  height = 9
)

p_summary <- ggplot(
  ncs_summary,
  aes(x = forcats::fct_reorder(NCS_Category, Candidate_Count), y = Candidate_Count)
) +
  geom_col(fill = "#4B8BBE", color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = Candidate_Count), hjust = -0.1, size = 5, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "CMap/NCS support summary",
    x = NULL,
    y = "Number of candidates"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )

save_plot_versions(
  p_summary,
  file.path(fig_dir, "Figure_2_CMap_NCS_category_summary"),
  width = 10,
  height = 6
)

# -----------------------------
# 11. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "CLUE_query_result_candidate_integration.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "CLUE_clean")
writeData(wb, "CLUE_clean", clue_clean)

addWorksheet(wb, "Top_opposing")
writeData(wb, "Top_opposing", top_opposing)

addWorksheet(wb, "Top_similar")
writeData(wb, "Top_similar", top_similar)

addWorksheet(wb, "Candidate_keys")
writeData(wb, "Candidate_keys", candidate_keys)

addWorksheet(wb, "Raw_matches")
writeData(wb, "Raw_matches", matched_raw)

addWorksheet(wb, "Candidate_best_NCS")
writeData(wb, "Candidate_best_NCS", candidate_best)

addWorksheet(wb, "All_61_integrated")
writeData(wb, "All_61_integrated", cand_integrated)

addWorksheet(wb, "NCS_summary")
writeData(wb, "NCS_summary", ncs_summary)

addWorksheet(wb, "Final_supported_rank")
writeData(wb, "Final_supported_rank", final_supported)

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
# 12. Save method note
# -----------------------------

method_note <- data.frame(
  Item = c(
    "Selected CLUE file",
    "Score used",
    "Why no tau threshold",
    "Candidate matching",
    "Next step"
  ),
  Description = c(
    "arfs/TAG/query_result.gct was used because it contains norm_cs, raw_cs and fdr_q_nlog10.",
    "norm_cs was used directly as the normalized connectivity score.",
    "The downloaded CLUE result provides normalized connectivity scores rather than tau; therefore NCS-specific thresholds were used.",
    "Candidate names were matched to perturbagen names using normalized names and selected synonyms.",
    "Use the final CMap-supported candidate ranking for CTD/DGIdb/DrugShot evidence mapping."
  )
)

write_csv(
  method_note,
  file.path(log_dir, "CLUE_query_result_integration_method_note.csv")
)

sink(file.path(log_dir, "session_info_CLUE_query_result_integration.txt"))
sessionInfo()
sink()

# -----------------------------
# 13. Console output
# -----------------------------

cat("\nCLUE query-result integration completed successfully.\n\n")

cat("Send me these files:\n")
cat(file.path(table_dir, "07_candidate_best_CLUE_NCS_summary.csv"), "\n")
cat(file.path(table_dir, "08_all_61_candidates_integrated_with_CLUE_query_result.csv"), "\n")
cat(file.path(table_dir, "09_NCS_category_summary.csv"), "\n")
cat(file.path(table_dir, "10_final_CMap_supported_candidate_ranking.csv"), "\n\n")

cat("NCS summary:\n")
print(ncs_summary)

cat("\nExcel workbook:\n")
cat(xlsx_file, "\n")
......................
/////////
  # =========================================================
# CLEAN CLUE / CMap FIGURE SCRIPT
# For: Imran_MTB_DRUG-MODULATION
# =========================================================

# -----------------------------
# 1. Load packages
# -----------------------------
required_pkgs <- c("ggplot2", "dplyr", "readr", "stringr", "forcats")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# -----------------------------
# 2. Set directories
# -----------------------------
project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"
clue_dir    <- file.path(project_dir, "results", "CLUE")
fig_dir     <- file.path(clue_dir, "final_clean_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Helper functions
# -----------------------------
find_file_by_name <- function(folder, filename_pattern) {
  files <- list.files(folder, pattern = filename_pattern, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop(paste("File not found:", filename_pattern))
  return(files[1])
}

clean_names_simple <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

wrap_axis_labels <- function(x, width = 20) {
  stringr::str_wrap(x, width = width)
}

# column detector
detect_col <- function(df, candidates) {
  nm_clean <- tolower(gsub("[^a-z0-9]", "", names(df)))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  idx <- which(nm_clean %in% cand_clean)
  if (length(idx) == 0) return(NA)
  names(df)[idx[1]]
}

# -----------------------------
# 4. Read main CLUE/CMap files
# -----------------------------
candidate_file <- find_file_by_name(clue_dir, "^10_final_CMap_supported_candidate_ranking\\.csv$")
ncs_summary_file <- find_file_by_name(clue_dir, "^09_NCS_category_summary\\.csv$")
group_summary_file <- find_file_by_name(clue_dir, "^05_final_candidate_group_summary\\.csv$|^05_candidate_group_summary\\.csv$")

cand_df <- read_csv(candidate_file, show_col_types = FALSE)
ncs_sum <- read_csv(ncs_summary_file, show_col_types = FALSE)
group_sum <- read_csv(group_summary_file, show_col_types = FALSE)

cat("Loaded files:\n")
cat(candidate_file, "\n")
cat(ncs_summary_file, "\n")
cat(group_summary_file, "\n\n")

# -----------------------------
# 5. Standardize candidate file
# -----------------------------
name_col <- detect_col(cand_df, c("candidate", "candidate_name", "compound", "compound_name",
                                  "perturbagen", "perturbagen_name", "drug", "drug_name",
                                  "matched_candidate_name", "candidate_clean"))
ncs_col  <- detect_col(cand_df, c("best_negative_norm_cs", "best_negative_ncs", "norm_cs",
                                  "ncs", "best_ncs", "best_negative_normalized_connectivity_score"))
cat_col  <- detect_col(cand_df, c("ncs_category", "category", "clue_ncs_category", "reversal_category"))
class_col <- detect_col(cand_df, c("candidate_class", "group", "candidate_group", "compound_group"))

if (is.na(name_col) | is.na(ncs_col) | is.na(cat_col)) {
  stop("Could not detect one or more required columns in 10_final_CMap_supported_candidate_ranking.csv")
}

cand_plot <- cand_df %>%
  dplyr::select(
    Candidate = all_of(name_col),
    Best_NCS = all_of(ncs_col),
    NCS_Category = all_of(cat_col),
    Candidate_Class = if (!is.na(class_col)) all_of(class_col) else NULL
  ) %>%
  mutate(
    Candidate = clean_names_simple(as.character(Candidate)),
    Best_NCS = as.numeric(Best_NCS),
    NCS_Category = as.character(NCS_Category)
  )

# keep only opposing candidates
cand_plot <- cand_plot %>%
  filter(!is.na(Best_NCS), Best_NCS < 0) %>%
  filter(str_detect(tolower(NCS_Category), "moderate|strong"))

# standardize category labels
cand_plot <- cand_plot %>%
  mutate(
    NCS_Category = case_when(
      str_detect(tolower(NCS_Category), "strong")   ~ "Strong opposing NCS",
      str_detect(tolower(NCS_Category), "moderate") ~ "Moderate opposing NCS",
      TRUE ~ NCS_Category
    )
  ) %>%
  arrange(Best_NCS)

# optional: keep top 30 only
top_n_show <- min(30, nrow(cand_plot))
cand_plot_top <- cand_plot %>%
  slice(1:top_n_show) %>%
  mutate(
    Candidate_wrapped = wrap_axis_labels(Candidate, width = 22),
    Candidate_wrapped = factor(Candidate_wrapped, levels = Candidate_wrapped)
  )

# save cleaned table used in plot
write_csv(cand_plot_top, file.path(fig_dir, "plot_input_top_opposing_candidates.csv"))

# -----------------------------
# 6. Color palettes
# -----------------------------
# Figure 1 colors
col_strong   <- "#D55E00"   # orange-red
col_moderate <- "#0072B2"   # blue

# Figure 2 colors
ncs_palette <- c(
  "Strong opposing NCS"   = "#D55E00",
  "Moderate opposing NCS" = "#0072B2",
  "Weak opposing NCS"     = "#E6AB02",
  "No CMap match"         = "#999999"
)

# Figure 3 colors
group_palette <- c(
  "Candidate compound"    = "#1B9E77",
  "Mapped LINCS compound" = "#7570B3",
  "Biologic / reference"  = "#E7298A",
  "Exposure / toxicant"   = "#D95F02",
  "Endogenous / metabolic"= "#66A61E"
)

# -----------------------------
# 7. Figure 1:
# CMap-supported opposing perturbagens
# -----------------------------
p1 <- ggplot(cand_plot_top,
             aes(x = Best_NCS, y = Candidate_wrapped, fill = NCS_Category)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = c(
    "Strong opposing NCS"   = col_strong,
    "Moderate opposing NCS" = col_moderate
  )) +
  labs(
    title = "CLUE/CMap-supported opposing perturbagens",
    x = "Best negative normalized connectivity score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12.5, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    plot.margin = margin(12, 25, 12, 12)
  )

ggsave(
  filename = file.path(fig_dir, "Figure_CLUE_opposing_perturbagens.png"),
  plot = p1, width = 11.5, height = 8.5, dpi = 600, bg = "white"
)
ggsave(
  filename = file.path(fig_dir, "Figure_CLUE_opposing_perturbagens.pdf"),
  plot = p1, width = 11.5, height = 8.5, bg = "white"
)

# -----------------------------
# 8. Standardize NCS summary file
# -----------------------------
sum_cat_col <- detect_col(ncs_sum, c("ncs_category", "category", "reversal_category"))
sum_n_col   <- detect_col(ncs_sum, c("count", "n", "frequency", "num_candidates"))

if (is.na(sum_cat_col) | is.na(sum_n_col)) {
  stop("Could not detect columns in 09_NCS_category_summary.csv")
}

ncs_plot <- ncs_sum %>%
  dplyr::select(Category = all_of(sum_cat_col), Count = all_of(sum_n_col)) %>%
  mutate(
    Category = case_when(
      str_detect(tolower(Category), "strong")   ~ "Strong opposing NCS",
      str_detect(tolower(Category), "moderate") ~ "Moderate opposing NCS",
      str_detect(tolower(Category), "weak")     ~ "Weak opposing NCS",
      str_detect(tolower(Category), "no")       ~ "No CMap match",
      TRUE ~ as.character(Category)
    ),
    Count = as.numeric(Count)
  ) %>%
  filter(!is.na(Category), !is.na(Count)) %>%
  mutate(Category = factor(Category, levels = c(
    "Strong opposing NCS",
    "Moderate opposing NCS",
    "Weak opposing NCS",
    "No CMap match"
  )))

write_csv(ncs_plot, file.path(fig_dir, "plot_input_ncs_summary.csv"))

# -----------------------------
# 9. Figure 2:
# CMap/NCS support summary
# -----------------------------
p2 <- ggplot(ncs_plot, aes(x = Count, y = fct_rev(Category), fill = Category)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = Count), hjust = -0.15, size = 5, fontface = "bold") +
  scale_fill_manual(values = ncs_palette, drop = FALSE) +
  expand_limits(x = max(ncs_plot$Count, na.rm = TRUE) * 1.15) +
  labs(
    title = "CLUE/CMap NCS support summary",
    x = "Number of candidates",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    legend.position = "none",
    plot.margin = margin(12, 20, 12, 12)
  )

ggsave(
  filename = file.path(fig_dir, "Figure_CLUE_NCS_summary.png"),
  plot = p2, width = 9.5, height = 6.8, dpi = 600, bg = "white"
)
ggsave(
  filename = file.path(fig_dir, "Figure_CLUE_NCS_summary.pdf"),
  plot = p2, width = 9.5, height = 6.8, bg = "white"
)

# -----------------------------
# 10. Standardize group summary file
# -----------------------------
grp_col <- detect_col(group_sum, c("candidate_group", "group", "class", "candidate_class"))
grp_ncol <- detect_col(group_sum, c("count", "n", "frequency", "num_candidates"))

if (is.na(grp_col) | is.na(grp_ncol)) {
  stop("Could not detect columns in candidate group summary file")
}

group_plot <- group_sum %>%
  dplyr::select(Group = all_of(grp_col), Count = all_of(grp_ncol)) %>%
  mutate(
    Group = as.character(Group),
    Count = as.numeric(Count)
  ) %>%
  filter(!is.na(Group), !is.na(Count))

write_csv(group_plot, file.path(fig_dir, "plot_input_candidate_group_summary.csv"))

# -----------------------------
# 11. Figure 3:
# Candidate group distribution
# -----------------------------
p3 <- ggplot(group_plot,
             aes(x = Count, y = fct_reorder(Group, Count), fill = Group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = Count), hjust = -0.15, size = 5, fontface = "bold") +
  scale_fill_manual(values = group_palette, na.value = "#4E79A7") +
  expand_limits(x = max(group_plot$Count, na.rm = TRUE) * 1.18) +
  labs(
    title = "Candidate group distribution",
    x = "Number of candidates",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    legend.position = "none",
    plot.margin = margin(12, 20, 12, 12)
  )

ggsave(
  filename = file.path(fig_dir, "Figure_candidate_group_distribution.png"),
  plot = p3, width = 9.5, height = 6.5, dpi = 600, bg = "white"
)
ggsave(
  filename = file.path(fig_dir, "Figure_candidate_group_distribution.pdf"),
  plot = p3, width = 9.5, height = 6.5, bg = "white"
)

# -----------------------------
# 12. Print output info
# -----------------------------
cat("\nAll clean figures saved in:\n")
cat(fig_dir, "\n")

cat("\nFiles created:\n")
print(list.files(fig_dir, full.names = TRUE))
.........................
############################################################
# 11_clean_CLUE_CMap_figures_FIXED.R
# Clean publication figures from CLUE/CMap integrated outputs
############################################################

packages <- c("ggplot2", "dplyr", "readr", "stringr", "forcats", "scales")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Project directory
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

out_dir <- file.path(project_dir, "results", "CLUE_CMap_clean_publication_figures")
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Robust file finder
# -----------------------------

find_project_file <- function(pattern, project_dir, required = TRUE) {
  
  files <- list.files(
    project_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  files <- files[!grepl("CLUE_CMap_clean_publication_figures", files)]
  
  if (length(files) == 0 && required) {
    stop(paste("File not found:", pattern))
  }
  
  if (length(files) == 0 && !required) {
    return(NA_character_)
  }
  
  return(files[1])
}

ranking_file <- find_project_file(
  "^10_final_CMap_supported_candidate_ranking\\.csv$",
  project_dir
)

integrated_file <- find_project_file(
  "^08_all_61_candidates_integrated_with_CLUE_query_result\\.csv$",
  project_dir
)

summary_file <- find_project_file(
  "^09_NCS_category_summary\\.csv$",
  project_dir
)

cat("Using ranking file:\n", ranking_file, "\n\n")
cat("Using integrated file:\n", integrated_file, "\n\n")
cat("Using summary file:\n", summary_file, "\n\n")

# -----------------------------
# 3. Read files
# -----------------------------

rank_df <- read_csv(ranking_file, show_col_types = FALSE)
all61_df <- read_csv(integrated_file, show_col_types = FALSE)
summary_df <- read_csv(summary_file, show_col_types = FALSE)

# -----------------------------
# 4. Check required columns
# -----------------------------

required_rank_cols <- c(
  "Candidate_Display_Name",
  "Best_Negative_NCS",
  "NCS_Category",
  "Suitability_Adjusted_Score",
  "Integrated_CLUE_Score"
)

missing_rank <- setdiff(required_rank_cols, colnames(rank_df))

if (length(missing_rank) > 0) {
  stop(paste("Missing columns in ranking file:", paste(missing_rank, collapse = ", ")))
}

required_all61_cols <- c(
  "Candidate_Display_Name",
  "Candidate_Group",
  "NCS_Category"
)

missing_all61 <- setdiff(required_all61_cols, colnames(all61_df))

if (length(missing_all61) > 0) {
  stop(paste("Missing columns in integrated file:", paste(missing_all61, collapse = ", ")))
}

# -----------------------------
# 5. Clean labels
# -----------------------------

clean_candidate_name <- function(x) {
  
  x <- as.character(x)
  
  x <- case_when(
    x == "1,25-dihydroxyvitamin D" ~ "1,25-dihydroxyvitamin D",
    x == "1,25 Dihydroxyvitamin D" ~ "1,25-dihydroxyvitamin D",
    x == "Ursodeoxycholic Acid" ~ "Ursodeoxycholic acid",
    x == "2,3-Dimethoxy-1,4-Naphtoquinone (Dmnq)" ~ "DMNQ",
    x == "2,3-Dimethoxy-1,4-Naphthoquinone (Dmnq)" ~ "DMNQ",
    TRUE ~ x
  )
  
  return(x)
}

rank_df <- rank_df %>%
  mutate(
    Candidate_Display_Name = clean_candidate_name(Candidate_Display_Name),
    Best_Negative_NCS = as.numeric(Best_Negative_NCS),
    Integrated_CLUE_Score = as.numeric(Integrated_CLUE_Score),
    NCS_Category = case_when(
      NCS_Category == "Strong opposing NCS" ~ "Strong opposing NCS",
      NCS_Category == "Moderate opposing NCS" ~ "Moderate opposing NCS",
      NCS_Category == "Weak opposing NCS" ~ "Weak opposing NCS",
      TRUE ~ NCS_Category
    )
  )

all61_df <- all61_df %>%
  mutate(
    Candidate_Display_Name = clean_candidate_name(Candidate_Display_Name),
    NCS_Category = ifelse(is.na(NCS_Category), "No CLUE/CMap match", NCS_Category),
    NCS_Category = ifelse(NCS_Category == "No CMap match", "No CLUE/CMap match", NCS_Category)
  )

# -----------------------------
# 6. Save plotting input tables
# -----------------------------

write_csv(
  rank_df,
  file.path(table_dir, "plot_input_final_CMap_supported_candidate_ranking.csv")
)

write_csv(
  all61_df,
  file.path(table_dir, "plot_input_all61_integrated_candidates.csv")
)

# -----------------------------
# 7. Colors
# -----------------------------

ncs_colors <- c(
  "Strong opposing NCS" = "#D95F02",
  "Moderate opposing NCS" = "#1F78B4",
  "Weak opposing NCS" = "#E6AB02",
  "No CLUE/CMap match" = "#999999",
  "Similar / disease-mimicking NCS" = "#7570B3"
)

group_colors <- c(
  "Candidate compound" = "#1B9E77",
  "Mapped LINCS compound" = "#7570B3",
  "Biologic / reference" = "#E7298A",
  "Exposure / toxicant" = "#D95F02",
  "Endogenous / metabolic" = "#66A61E"
)

# -----------------------------
# 8. Save function
# -----------------------------

save_plot <- function(plot, filename, width, height) {
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

# -----------------------------
# 9. Figure 1: Opposing candidates
# -----------------------------

top30 <- rank_df %>%
  filter(
    !is.na(Best_Negative_NCS),
    NCS_Category %in% c("Strong opposing NCS", "Moderate opposing NCS")
  ) %>%
  arrange(Best_Negative_NCS) %>%
  head(30) %>%
  mutate(
    Candidate_Label = stringr::str_wrap(Candidate_Display_Name, width = 28),
    Candidate_Label = factor(Candidate_Label, levels = rev(Candidate_Label))
  )

p1 <- ggplot(
  top30,
  aes(x = Best_Negative_NCS, y = Candidate_Label, fill = NCS_Category)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "grey35", linewidth = 0.65) +
  scale_fill_manual(values = ncs_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.06, 0.04))
  ) +
  labs(
    title = "CLUE/CMap opposing-connectivity candidates",
    x = "Best negative normalized connectivity score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12.5, color = "black"),
    axis.text.y = element_text(size = 11.5, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(0.45, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 30, 35, 12)
  )

save_plot(p1, "Figure_1_CLUE_CMap_opposing_connectivity_candidates", 11.5, 8.8)

# -----------------------------
# 10. Figure 2: NCS category summary
# -----------------------------

ncs_summary <- all61_df %>%
  count(NCS_Category, name = "Candidate_Count") %>%
  mutate(
    NCS_Category = factor(
      NCS_Category,
      levels = c(
        "Strong opposing NCS",
        "Moderate opposing NCS",
        "Weak opposing NCS",
        "Similar / disease-mimicking NCS",
        "No CLUE/CMap match"
      )
    )
  ) %>%
  arrange(NCS_Category)

write_csv(
  ncs_summary,
  file.path(table_dir, "plot_input_NCS_category_summary.csv")
)

p2 <- ggplot(
  ncs_summary,
  aes(x = Candidate_Count, y = fct_rev(NCS_Category), fill = NCS_Category)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.68) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.12,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = ncs_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "CLUE/CMap NCS category summary",
    x = "Number of candidates",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12.5, color = "black"),
    axis.text.y = element_text(size = 12.5, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )

save_plot(p2, "Figure_2_CLUE_CMap_NCS_category_summary", 9.8, 6.2)

# -----------------------------
# 11. Figure 3: Suitability score vs NCS
# -----------------------------

scatter_df <- all61_df %>%
  filter(!is.na(Best_Negative_NCS)) %>%
  mutate(
    Best_Negative_NCS = as.numeric(Best_Negative_NCS),
    Suitability_Adjusted_Score = as.numeric(Suitability_Adjusted_Score),
    Supporting_Signatures = as.numeric(Supporting_Signatures),
    NCS_Category = ifelse(NCS_Category == "No CMap match", "No CLUE/CMap match", NCS_Category)
  )

p3 <- ggplot(
  scatter_df,
  aes(
    x = Suitability_Adjusted_Score,
    y = Best_Negative_NCS,
    color = NCS_Category,
    size = Supporting_Signatures
  )
) +
  geom_hline(yintercept = -1.5, linetype = "dashed", color = "grey35", linewidth = 0.65) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = ncs_colors, drop = FALSE) +
  scale_size_continuous(range = c(3, 8)) +
  labs(
    title = "Relationship between prioritization score and CLUE/CMap connectivity",
    x = "Suitability-adjusted perturbation score",
    y = "Best negative normalized connectivity score",
    color = NULL,
    size = "Supporting\nsignatures"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 15.5, face = "bold"),
    axis.title.y = element_text(size = 15.5, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11.5, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 20, 12, 12)
  )

save_plot(p3, "Figure_3_suitability_score_vs_CLUE_CMap_NCS", 11, 7.2)

# -----------------------------
# 12. Figure 4: Integrated CLUE score ranking
# -----------------------------

top25_integrated <- rank_df %>%
  arrange(desc(Integrated_CLUE_Score)) %>%
  head(25) %>%
  mutate(
    Candidate_Label = stringr::str_wrap(Candidate_Display_Name, width = 28),
    Candidate_Label = factor(Candidate_Label, levels = rev(Candidate_Label))
  )

p4 <- ggplot(
  top25_integrated,
  aes(x = Integrated_CLUE_Score, y = Candidate_Label, fill = NCS_Category)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  scale_fill_manual(values = ncs_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  geom_text(
    aes(label = round(Integrated_CLUE_Score, 1)),
    hjust = -0.10,
    size = 4.1,
    fontface = "bold"
  ) +
  labs(
    title = "Integrated CLUE/CMap-supported candidate ranking",
    x = "Integrated CLUE/CMap score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 15.5, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11.5, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 35, 12)
  )

save_plot(p4, "Figure_4_integrated_CLUE_CMap_supported_candidate_ranking", 11, 8.2)

# -----------------------------
# 13. Figure 5: Candidate group distribution
# -----------------------------

group_summary <- all61_df %>%
  count(Candidate_Group, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  group_summary,
  file.path(table_dir, "plot_input_candidate_group_distribution.csv")
)

p5 <- ggplot(
  group_summary,
  aes(x = Candidate_Count, y = fct_reorder(Candidate_Group, Candidate_Count), fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.68) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.12,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "Candidate group distribution",
    x = "Number of candidates",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12.5, color = "black"),
    axis.text.y = element_text(size = 12.5, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )

save_plot(p5, "Figure_5_candidate_group_distribution_clean", 9.8, 6.2)

# -----------------------------
# 14. Finish
# -----------------------------

cat("\nClean CLUE/CMap figures completed successfully.\n\n")
cat("Figures saved in:\n")
cat(fig_dir, "\n\n")

cat("Tables saved in:\n")
cat(table_dir, "\n\n")

cat("Created figures:\n")
print(list.files(fig_dir, full.names = FALSE))