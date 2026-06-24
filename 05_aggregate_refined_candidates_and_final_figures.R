############################################################
# 05_aggregate_refined_candidates_and_final_figures.R
# Correct LINCS names, aggregate duplicate candidates,
# and generate clean final figures without background grid lines
############################################################

# -----------------------------
# 1. Load packages
# -----------------------------

packages <- c(
  "dplyr", "readr", "stringr", "ggplot2",
  "forcats", "scales", "openxlsx", "tidyr"
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

out_dir <- file.path(project_dir, "results", "candidate_final_aggregation")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:\n")
cat(project_dir, "\n\n")

# -----------------------------
# 3. Find input file
# -----------------------------

find_file_recursive <- function(project_dir, file_pattern, required = TRUE) {
  
  matched <- list.files(
    project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matched <- matched[!grepl("candidate_final_aggregation", matched)]
  
  if (length(matched) == 0 && required) {
    stop(paste("File not found:", file_pattern))
  }
  
  if (length(matched) == 0 && !required) {
    return(NA_character_)
  }
  
  return(matched[1])
}

input_file <- find_file_recursive(
  project_dir,
  "02_ranked_all_active_TB_signature_suppressor_candidates\\.csv$",
  required = TRUE
)

cat("Using input file:\n")
cat(input_file, "\n\n")

ranked <- readr::read_csv(input_file, show_col_types = FALSE)

if (!"Candidate_Name" %in% colnames(ranked)) {
  stop("Candidate_Name column not found.")
}

# -----------------------------
# 4. Required columns
# -----------------------------

needed_cols <- c(
  "Candidate_Class", "Priority_Score_100", "Supporting_Signatures",
  "Library_Count", "Best_Adjusted_P", "Max_Combined_Score",
  "Max_Overlap_Count", "Max_Anchor_Hit_Count", "Max_Core8_Hit_Count",
  "Anchor_Hits", "Core8_Hits", "Libraries", "Top_Term"
)

for (col in needed_cols) {
  if (!col %in% colnames(ranked)) {
    ranked[[col]] <- NA
  }
}

# -----------------------------
# 5. Helper functions
# -----------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  return(x)
}

extract_lincs_compound <- function(candidate_name, top_term = "") {
  
  x <- paste(candidate_name, top_term)
  x <- clean_text(x)
  
  hit <- stringr::str_match(
    x,
    stringr::regex(
      "(?:3H|6H|24H|48H)-(.+?)-[0-9]+(?:\\.[0-9]+)?(?:\\s|$)",
      ignore_case = TRUE
    )
  )[, 2]
  
  if (!is.na(hit) && hit != "") {
    
    hit <- stringr::str_replace_all(hit, "\\+", " ")
    hit <- stringr::str_replace_all(hit, "_", " ")
    hit <- stringr::str_trim(hit)
    hit_upper <- stringr::str_to_upper(hit)
    
    hit_final <- dplyr::case_when(
      stringr::str_detect(hit_upper, "DASATINIB") ~ "Dasatinib",
      stringr::str_detect(hit_upper, "WITHAFERIN-A") ~ "Withaferin-A",
      stringr::str_detect(hit_upper, "CRIZOTINIB") ~ "Crizotinib",
      stringr::str_detect(hit_upper, "PELITINIB") ~ "Pelitinib",
      stringr::str_detect(hit_upper, "GEFITINIB") ~ "Gefitinib",
      stringr::str_detect(hit_upper, "RADICICOL") ~ "Radicicol",
      stringr::str_detect(hit_upper, "CHELERYTHRINE") ~ "Chelerythrine Chloride",
      TRUE ~ hit_upper
    )
    
    return(hit_final)
  }
  
  return(NA_character_)
}

refine_candidate_name <- function(candidate_name, top_term = "") {
  
  lincs_hit <- extract_lincs_compound(candidate_name, top_term)
  
  if (!is.na(lincs_hit) && lincs_hit != "") {
    return(lincs_hit)
  }
  
  x <- clean_text(candidate_name)
  
  x <- stringr::str_replace(x, stringr::regex("\\s+CID$", ignore_case = TRUE), "")
  x <- stringr::str_replace(x, stringr::regex("\\s+[0-9]+\\s+(human|mouse|rat)$", ignore_case = TRUE), "")
  x <- stringr::str_replace(x, stringr::regex("\\s+(human|mouse|rat)$", ignore_case = TRUE), "")
  x <- stringr::str_replace(x, "\\s+[0-9]{3,}$", "")
  
  x <- stringr::str_replace(x, stringr::regex("^Curcumin\\s+Cid$", ignore_case = TRUE), "Curcumin")
  x <- stringr::str_replace(x, stringr::regex("^Palmitic Acid\\s+985\\s+Mouse$", ignore_case = TRUE), "Palmitic Acid")
  x <- stringr::str_replace(x, stringr::regex("^Formaldehyde\\s+712\\s+Rat$", ignore_case = TRUE), "Formaldehyde")
  x <- stringr::str_replace(x, stringr::regex("^Plx4720\\s+Cid$", ignore_case = TRUE), "PLX4720")
  
  x <- stringr::str_to_title(x)
  
  x <- stringr::str_replace_all(x, "Plx4720", "PLX4720")
  x <- stringr::str_replace_all(x, "Plx4032", "PLX4032")
  x <- stringr::str_replace_all(x, "Dmng", "DMNQ")
  x <- stringr::str_replace_all(x, "Dmno", "DMNQ")
  x <- stringr::str_replace_all(x, "Tpca-1", "TPCA-1")
  x <- stringr::str_replace_all(x, "Xmd-1150", "XMD-1150")
  x <- stringr::str_replace_all(x, "Bi-2536", "BI-2536")
  x <- stringr::str_replace_all(x, "Tg-101348", "TG-101348")
  x <- stringr::str_replace_all(x, "Zm-447439", "ZM-447439")
  x <- stringr::str_replace_all(x, "Hg-6-64-01", "HG-6-64-01")
  
  return(x)
}

assign_candidate_group <- function(refined_name, candidate_class = "") {
  
  x <- tolower(refined_name)
  cls <- tolower(candidate_class)
  
  toxic_terms <- c(
    "cigarette", "smoke", "nickel", "cadmium", "arsenic", "lead",
    "mercury", "hypochlorous acid", "nitric oxide", "formaldehyde",
    "chlorpyrifos", "ethanol", "vanadium", "cobalt dichloride",
    "heroin", "ozone", "diesel", "radiation", "ultraviolet"
  )
  
  endogenous_terms <- c(
    "adenosine triphosphate", "atp", "arachidonic acid",
    "palmitic acid", "cholesterol", "glucose", "oleic acid",
    "stearic acid", "pyruvate", "lactic acid", "hydrogen peroxide"
  )
  
  biologic_terms <- c(
    "etanercept", "rituximab", "adalimumab", "infliximab",
    "interferon", "antibody"
  )
  
  if (any(stringr::str_detect(x, toxic_terms)) || stringr::str_detect(cls, "exposure|toxicant")) {
    return("Exposure / toxicant")
  }
  
  if (any(stringr::str_detect(x, endogenous_terms)) || stringr::str_detect(cls, "endogenous|metabolic")) {
    return("Endogenous / metabolic")
  }
  
  if (any(stringr::str_detect(x, biologic_terms)) || stringr::str_detect(cls, "biologic")) {
    return("Biologic / reference")
  }
  
  if (stringr::str_detect(cls, "unmapped|genetic")) {
    return("Mapped LINCS compound")
  }
  
  return("Candidate compound")
}

assign_broad_axis <- function(refined_name) {
  
  x <- tolower(refined_name)
  
  if (stringr::str_detect(x, "alfacalcidol|calcitriol|vitamin d|cholecalciferol|dihydroxyvitamin")) {
    return("Vitamin D / VDR axis")
  }
  
  if (stringr::str_detect(x, "imatinib|dasatinib|sunitinib|plx4720|plx4032|tpca-1|xmd-1150|bi-2536|crizotinib|pazopanib|sb590885|zm-447439|nvp-auy922|tg-101348|as-601245|hg-6-64-01|gsk-461364|azd-7762|wz-4-145")) {
    return("Kinase signaling axis")
  }
  
  if (stringr::str_detect(x, "curcumin|celecoxib|dexamethasone|dmnq|withaferin")) {
    return("Inflammatory / oxidative signaling")
  }
  
  if (stringr::str_detect(x, "resveratrol|metformin|quercetin")) {
    return("Immunometabolic / antioxidant axis")
  }
  
  if (stringr::str_detect(x, "etanercept|tnf")) {
    return("TNF / cytokine axis")
  }
  
  if (stringr::str_detect(x, "methotrexate|cyclosporine")) {
    return("Broad immunomodulatory axis")
  }
  
  if (stringr::str_detect(x, "ursodeoxycholic")) {
    return("Bile acid / inflammatory regulation")
  }
  
  if (stringr::str_detect(x, "estradiol|anastrozole|medroxyprogesterone")) {
    return("Hormone-related axis")
  }
  
  if (stringr::str_detect(x, "citalopram|carbidopa")) {
    return("Neuroactive / off-axis pharmacology")
  }
  
  return("Needs biological mapping")
}

suitability_factor <- function(group) {
  
  dplyr::case_when(
    group == "Candidate compound" ~ 1.00,
    group == "Mapped LINCS compound" ~ 0.90,
    group == "Biologic / reference" ~ 0.75,
    group == "Endogenous / metabolic" ~ 0.45,
    group == "Exposure / toxicant" ~ 0.15,
    TRUE ~ 0.50
  )
}

choose_candidate_group <- function(groups) {
  
  groups <- unique(groups)
  
  priority <- c(
    "Candidate compound",
    "Mapped LINCS compound",
    "Biologic / reference",
    "Endogenous / metabolic",
    "Exposure / toxicant"
  )
  
  hit <- priority[priority %in% groups]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  return(groups[1])
}

choose_axis <- function(axes) {
  
  axes <- unique(axes)
  axes <- axes[!is.na(axes) & axes != ""]
  
  if (length(axes) == 0) {
    return("Needs biological mapping")
  }
  
  priority <- c(
    "Vitamin D / VDR axis",
    "Kinase signaling axis",
    "Inflammatory / oxidative signaling",
    "Immunometabolic / antioxidant axis",
    "TNF / cytokine axis",
    "Broad immunomodulatory axis",
    "Bile acid / inflammatory regulation",
    "Hormone-related axis",
    "Neuroactive / off-axis pharmacology",
    "Needs biological mapping"
  )
  
  hit <- priority[priority %in% axes]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  return(axes[1])
}

save_plot_versions <- function(plot_obj, file_base, width = 12, height = 8) {
  
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

# -----------------------------
# 6. Refine candidates before aggregation
# -----------------------------

refined_terms <- ranked %>%
  mutate(
    Refined_Candidate_Name = mapply(refine_candidate_name, Candidate_Name, Top_Term),
    Candidate_Group = mapply(assign_candidate_group, Refined_Candidate_Name, Candidate_Class),
    Broad_Target_Axis = sapply(Refined_Candidate_Name, assign_broad_axis),
    Suitability_Factor = suitability_factor(Candidate_Group),
    Suitability_Adjusted_Score = round(Priority_Score_100 * Suitability_Factor, 2)
  )

write_csv(
  refined_terms,
  file.path(table_dir, "01_refined_terms_before_aggregation.csv")
)

# -----------------------------
# 7. Aggregate duplicated candidates
# -----------------------------

aggregated <- refined_terms %>%
  group_by(Refined_Candidate_Name) %>%
  summarise(
    Candidate_Group = choose_candidate_group(Candidate_Group),
    Broad_Target_Axis = choose_axis(Broad_Target_Axis),
    Priority_Score_100 = max(Priority_Score_100, na.rm = TRUE),
    Suitability_Factor = suitability_factor(Candidate_Group),
    Suitability_Adjusted_Score = round(Priority_Score_100 * Suitability_Factor, 2),
    Supporting_Signatures = sum(Supporting_Signatures, na.rm = TRUE),
    Library_Count = max(Library_Count, na.rm = TRUE),
    Best_Adjusted_P = min(Best_Adjusted_P, na.rm = TRUE),
    Max_Combined_Score = max(Max_Combined_Score, na.rm = TRUE),
    Max_Overlap_Count = max(Max_Overlap_Count, na.rm = TRUE),
    Max_Anchor_Hit_Count = max(Max_Anchor_Hit_Count, na.rm = TRUE),
    Max_Core8_Hit_Count = max(Max_Core8_Hit_Count, na.rm = TRUE),
    Anchor_Hits = paste(unique(Anchor_Hits[!is.na(Anchor_Hits) & Anchor_Hits != ""]), collapse = ";"),
    Core8_Hits = paste(unique(Core8_Hits[!is.na(Core8_Hits) & Core8_Hits != ""]), collapse = ";"),
    Libraries = paste(unique(Libraries[!is.na(Libraries) & Libraries != ""]), collapse = ";"),
    Source_Terms = paste(unique(Top_Term), collapse = " || "),
    Source_Candidates = paste(unique(Candidate_Name), collapse = " || "),
    .groups = "drop"
  ) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

write_csv(
  aggregated,
  file.path(table_dir, "02_aggregated_refined_candidate_table.csv")
)

candidate_for_manual <- aggregated %>%
  filter(Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound", "Biologic / reference")) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

flagged_low_suitability <- aggregated %>%
  filter(Candidate_Group %in% c("Endogenous / metabolic", "Exposure / toxicant")) %>%
  arrange(desc(Priority_Score_100), Best_Adjusted_P)

write_csv(
  candidate_for_manual,
  file.path(table_dir, "03_final_candidate_perturbagens_for_manual_check.csv")
)

write_csv(
  flagged_low_suitability,
  file.path(table_dir, "04_flagged_low_suitability_perturbagens.csv")
)

group_summary <- aggregated %>%
  count(Candidate_Group, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

axis_summary <- aggregated %>%
  count(Broad_Target_Axis, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  group_summary,
  file.path(table_dir, "05_final_candidate_group_summary.csv")
)

write_csv(
  axis_summary,
  file.path(table_dir, "06_final_target_axis_summary.csv")
)

# -----------------------------
# 8. Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "final_aggregated_candidate_curation.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "Aggregated_candidates")
writeData(wb, "Aggregated_candidates", aggregated)

addWorksheet(wb, "Manual_check")
writeData(wb, "Manual_check", candidate_for_manual)

addWorksheet(wb, "Flagged_low_suitability")
writeData(wb, "Flagged_low_suitability", flagged_low_suitability)

addWorksheet(wb, "Refined_terms")
writeData(wb, "Refined_terms", refined_terms)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Target_axis_summary")
writeData(wb, "Target_axis_summary", axis_summary)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(wb)) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:100, gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:60, widths = "auto")
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# -----------------------------
# 9. Figures
# -----------------------------

group_colors <- c(
  "Candidate compound" = "#00B050",
  "Mapped LINCS compound" = "#7BC67B",
  "Biologic / reference" = "#00AEB3",
  "Endogenous / metabolic" = "#5B9BD5",
  "Exposure / toxicant" = "#F8766D"
)

top30_all <- aggregated %>%
  arrange(desc(Priority_Score_100), Best_Adjusted_P) %>%
  head(30) %>%
  mutate(
    Label = stringr::str_wrap(Refined_Candidate_Name, width = 35),
    Label = forcats::fct_reorder(Label, Priority_Score_100)
  )

p_before <- ggplot(
  top30_all,
  aes(x = Label, y = Priority_Score_100, fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(
    aes(label = round(Priority_Score_100, 2)),
    hjust = -0.08,
    size = 4.8,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Top perturbation signals before pharmacological filtering",
    x = NULL,
    y = "Priority score",
    fill = NULL
  ) +
  theme_classic(base_size = 19) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 21, face = "bold"),
    axis.text.x = element_text(size = 16, face = "plain"),
    axis.text.y = element_text(size = 15, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 13, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 60, 15, 15)
  )

save_plot_versions(
  p_before,
  file.path(fig_dir, "Figure_1_aggregated_top_perturbation_signals_before_filtering"),
  width = 13,
  height = 12
)

top25_candidates <- candidate_for_manual %>%
  head(25) %>%
  mutate(
    Label = stringr::str_wrap(Refined_Candidate_Name, width = 35),
    Label = forcats::fct_reorder(Label, Suitability_Adjusted_Score)
  )

p_after <- ggplot(
  top25_candidates,
  aes(x = Label, y = Suitability_Adjusted_Score, fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(
    aes(label = round(Suitability_Adjusted_Score, 2)),
    hjust = -0.08,
    size = 4.8,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Candidate perturbagens after suitability adjustment",
    x = NULL,
    y = "Suitability-adjusted score",
    fill = NULL
  ) +
  theme_classic(base_size = 19) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 21, face = "bold"),
    axis.text.x = element_text(size = 16, face = "plain"),
    axis.text.y = element_text(size = 15, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 13, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 60, 15, 15)
  )

save_plot_versions(
  p_after,
  file.path(fig_dir, "Figure_2_aggregated_candidate_perturbagens_after_suitability_adjustment"),
  width = 13,
  height = 11
)

p_group <- ggplot(
  group_summary,
  aes(x = forcats::fct_reorder(Candidate_Group, Candidate_Count), y = Candidate_Count)
) +
  geom_col(fill = "#4B8BBE", color = "black", linewidth = 0.35, width = 0.7) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.1,
    size = 5.8,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "Candidate group distribution",
    x = NULL,
    y = "Number of candidates"
  ) +
  theme_classic(base_size = 19) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 21, face = "bold"),
    axis.text.x = element_text(size = 16, face = "plain"),
    axis.text.y = element_text(size = 16, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 45, 15, 15)
  )

save_plot_versions(
  p_group,
  file.path(fig_dir, "Figure_3_aggregated_candidate_group_distribution"),
  width = 11,
  height = 7
)

# -----------------------------
# 10. Save notes
# -----------------------------

note <- data.frame(
  Item = c(
    "Duplicate aggregation",
    "LINCS correction",
    "Before filtering figure",
    "After suitability adjustment figure",
    "Next step"
  ),
  Description = c(
    "Candidates with the same refined name were merged before plotting.",
    "LINCS perturbagens such as AS-601245, XMD-1150, TPCA-1, BI-2536, TG-101348, HG-6-64-01 and ZM-447439 were extracted more accurately.",
    "This figure shows the raw perturbation-priority ranking and is best suited for transparency or supplementary material.",
    "This figure is more suitable for the main manuscript because toxicants and endogenous terms are down-weighted or removed from the displayed candidate list.",
    "Manually check the final candidate table before CMap/CLUE, CTD, DGIdb/DrugShot and target-axis mapping."
  )
)

write_csv(
  note,
  file.path(log_dir, "aggregation_and_plotting_notes.csv")
)

sink(file.path(log_dir, "session_info_final_candidate_aggregation.txt"))
sessionInfo()
sink()

# -----------------------------
# 11. Console output
# -----------------------------

cat("\nFinal aggregation completed successfully.\n\n")

cat("Main table for manual checking:\n")
cat(file.path(table_dir, "03_final_candidate_perturbagens_for_manual_check.csv"), "\n\n")

cat("Excel workbook:\n")
cat(xlsx_file, "\n\n")

cat("Figures saved in:\n")
cat(fig_dir, "\n\n")

cat("Group summary:\n")
print(group_summary)

cat("\nTop candidates after aggregation and suitability adjustment:\n")
print(head(candidate_for_manual, 25))
..............
############################################################
# 06_final_candidate_space_figures_and_tables.R
# Final candidate-space tables and clean publication figures
############################################################

packages <- c(
  "dplyr", "readr", "stringr", "ggplot2",
  "forcats", "scales", "openxlsx", "tidyr"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

out_dir <- file.path(project_dir, "results", "final_candidate_space")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

find_file_recursive <- function(project_dir, file_pattern, required = TRUE) {
  matched <- list.files(
    project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matched <- matched[!grepl("final_candidate_space", matched)]
  
  if (length(matched) == 0 && required) {
    stop(paste("File not found:", file_pattern))
  }
  
  if (length(matched) == 0 && !required) {
    return(NA_character_)
  }
  
  matched[1]
}

candidate_file <- find_file_recursive(
  project_dir,
  "03_final_candidate_perturbagens_for_manual_check\\.csv$",
  required = TRUE
)

cat("Using candidate file:\n")
cat(candidate_file, "\n\n")

cand <- readr::read_csv(candidate_file, show_col_types = FALSE)

needed_cols <- c(
  "Refined_Candidate_Name", "Candidate_Group", "Broad_Target_Axis",
  "Priority_Score_100", "Suitability_Factor", "Suitability_Adjusted_Score",
  "Supporting_Signatures", "Library_Count", "Best_Adjusted_P",
  "Max_Combined_Score", "Max_Overlap_Count", "Max_Anchor_Hit_Count",
  "Max_Core8_Hit_Count", "Anchor_Hits", "Core8_Hits", "Libraries",
  "Source_Terms", "Source_Candidates"
)

for (col in needed_cols) {
  if (!col %in% names(cand)) {
    cand[[col]] <- NA
  }
}

# -----------------------------
# 1. Correct VX classification
# -----------------------------

cand <- cand %>%
  mutate(
    Refined_Candidate_Name = ifelse(
      Refined_Candidate_Name == "Vx",
      "VX",
      Refined_Candidate_Name
    ),
    Candidate_Group = ifelse(
      Refined_Candidate_Name == "VX",
      "Exposure / toxicant",
      Candidate_Group
    ),
    Broad_Target_Axis = ifelse(
      Refined_Candidate_Name == "VX",
      "Excluded toxicant / nerve agent",
      Broad_Target_Axis
    ),
    Suitability_Factor = case_when(
      Candidate_Group == "Candidate compound" ~ 1.00,
      Candidate_Group == "Mapped LINCS compound" ~ 0.90,
      Candidate_Group == "Biologic / reference" ~ 0.75,
      Candidate_Group == "Endogenous / metabolic" ~ 0.45,
      Candidate_Group == "Exposure / toxicant" ~ 0.15,
      TRUE ~ 0.50
    ),
    Suitability_Adjusted_Score = round(Priority_Score_100 * Suitability_Factor, 2)
  ) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

# -----------------------------
# 2. Add manual-status columns
# -----------------------------

cand_final <- cand %>%
  mutate(
    Suggested_Use = case_when(
      Candidate_Group == "Exposure / toxicant" ~ "Exclude from therapeutic prioritization",
      Candidate_Group == "Endogenous / metabolic" ~ "Keep only as biological reference",
      Candidate_Group == "Biologic / reference" ~ "Reference / caution",
      Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound") ~ "Candidate for external validation",
      TRUE ~ "Manual review"
    ),
    Manual_Status = "",
    Manual_Reason = "",
    Keep_for_CMap = ifelse(
      Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound", "Biologic / reference"),
      "Yes",
      "No"
    ),
    Keep_for_CTD = ifelse(
      Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound", "Biologic / reference", "Endogenous / metabolic"),
      "Yes",
      "No"
    ),
    Keep_for_Target_Mapping = ifelse(
      Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound", "Biologic / reference"),
      "Yes",
      "No"
    )
  )

write_csv(
  cand_final,
  file.path(table_dir, "01_all_61_refined_candidate_perturbagens_final.csv")
)

candidate_validation <- cand_final %>%
  filter(Candidate_Group %in% c("Candidate compound", "Mapped LINCS compound", "Biologic / reference")) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

write_csv(
  candidate_validation,
  file.path(table_dir, "02_candidates_for_CMap_CTD_target_mapping.csv")
)

flagged <- cand_final %>%
  filter(Candidate_Group %in% c("Exposure / toxicant", "Endogenous / metabolic"))

write_csv(
  flagged,
  file.path(table_dir, "03_flagged_reference_or_excluded_perturbagens.csv")
)

group_summary <- cand_final %>%
  count(Candidate_Group, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

axis_summary <- cand_final %>%
  count(Broad_Target_Axis, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(group_summary, file.path(table_dir, "04_candidate_group_summary.csv"))
write_csv(axis_summary, file.path(table_dir, "05_target_axis_summary.csv"))

# -----------------------------
# 3. Optional gene input summary
# -----------------------------

active_gene_file <- find_file_recursive(
  project_dir,
  "CMap_UP_active_TB_top150\\.txt$",
  required = FALSE
)

post_gene_file <- find_file_recursive(
  project_dir,
  "CMap_DOWN_posttherapy_top150\\.txt$",
  required = FALSE
)

core8_file <- find_file_recursive(
  project_dir,
  "core8_host_response_module\\.txt$",
  required = FALSE
)

anchor_file <- find_file_recursive(
  project_dir,
  "anchor_STAT1_GBP1\\.txt$",
  required = FALSE
)

read_gene_list <- function(file_path) {
  if (is.na(file_path)) return("")
  genes <- readLines(file_path, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  paste(genes, collapse = ";")
}

gene_input_summary <- data.frame(
  Gene_Set = c(
    "Active-TB-up genes used as UP query",
    "Post-therapy-up genes used as DOWN query",
    "STAT1/GBP1 anchor genes",
    "Core 8-gene host-response module"
  ),
  File_Found = c(
    !is.na(active_gene_file),
    !is.na(post_gene_file),
    !is.na(anchor_file),
    !is.na(core8_file)
  ),
  Genes = c(
    read_gene_list(active_gene_file),
    read_gene_list(post_gene_file),
    read_gene_list(anchor_file),
    read_gene_list(core8_file)
  )
)

write_csv(
  gene_input_summary,
  file.path(table_dir, "06_gene_input_summary_for_perturbation_analysis.csv")
)

# -----------------------------
# 4. Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "final_candidate_space_tables.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "All_61_candidates")
writeData(wb, "All_61_candidates", cand_final)

addWorksheet(wb, "For_CMap_CTD_mapping")
writeData(wb, "For_CMap_CTD_mapping", candidate_validation)

addWorksheet(wb, "Flagged_reference_exclude")
writeData(wb, "Flagged_reference_exclude", flagged)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Target_axis_summary")
writeData(wb, "Target_axis_summary", axis_summary)

addWorksheet(wb, "Gene_input_summary")
writeData(wb, "Gene_input_summary", gene_input_summary)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(wb)) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:100, gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:80, widths = "auto")
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# -----------------------------
# 5. Clean figure settings
# -----------------------------

group_colors <- c(
  "Candidate compound" = "#00B050",
  "Mapped LINCS compound" = "#7BC67B",
  "Biologic / reference" = "#00AEB3",
  "Endogenous / metabolic" = "#5B9BD5",
  "Exposure / toxicant" = "#F8766D"
)

save_plot_versions <- function(plot_obj, file_base, width = 12, height = 8) {
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

# -----------------------------
# 6. Main top-25 candidate figure
# -----------------------------

top25 <- candidate_validation %>%
  head(25) %>%
  mutate(
    Label = stringr::str_wrap(Refined_Candidate_Name, width = 34),
    Label = forcats::fct_reorder(Label, Suitability_Adjusted_Score)
  )

p_top25 <- ggplot(
  top25,
  aes(x = Label, y = Suitability_Adjusted_Score, fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(
    aes(label = round(Suitability_Adjusted_Score, 2)),
    hjust = -0.08,
    size = 4.6,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "Candidate perturbagens after suitability adjustment",
    x = NULL,
    y = "Suitability-adjusted score",
    fill = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 15, face = "plain"),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.position = "right",
    legend.text = element_text(size = 13, face = "plain"),
    legend.key.size = unit(0.55, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 95, 15, 15)
  )

save_plot_versions(
  p_top25,
  file.path(fig_dir, "Figure_1_top25_candidate_perturbagens_no_cut_labels"),
  width = 15,
  height = 11
)

# -----------------------------
# 7. All-61 bubble plot
# -----------------------------

bubble_data <- cand_final %>%
  mutate(
    Candidate_Group = factor(
      Candidate_Group,
      levels = c(
        "Candidate compound",
        "Mapped LINCS compound",
        "Biologic / reference",
        "Endogenous / metabolic",
        "Exposure / toxicant"
      )
    ),
    Axis_Label = stringr::str_wrap(Broad_Target_Axis, width = 28)
  )

p_bubble <- ggplot(
  bubble_data,
  aes(
    x = Axis_Label,
    y = Candidate_Group,
    size = Supporting_Signatures,
    color = Suitability_Adjusted_Score
  )
) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "#D9EAF7", high = "#E95A4F") +
  labs(
    title = "Candidate perturbagen space derived from active-TB signature suppression",
    x = NULL,
    y = NULL,
    size = "Supporting\nsignatures",
    color = "Adjusted\nscore"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 13, face = "plain", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "plain"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(15, 35, 35, 15)
  )

save_plot_versions(
  p_bubble,
  file.path(fig_dir, "Figure_2_all61_candidate_perturbagen_space_bubble_plot"),
  width = 15,
  height = 8.5
)

# -----------------------------
# 8. Candidate group distribution
# -----------------------------

p_group <- ggplot(
  group_summary,
  aes(x = forcats::fct_reorder(Candidate_Group, Candidate_Count), y = Candidate_Count)
) +
  geom_col(fill = "#4B8BBE", color = "black", linewidth = 0.35, width = 0.7) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.1,
    size = 5.5,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "Candidate group distribution",
    x = NULL,
    y = "Number of candidates"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 15, face = "plain"),
    axis.text.y = element_text(size = 15, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 45, 15, 15)
  )

save_plot_versions(
  p_group,
  file.path(fig_dir, "Figure_3_candidate_group_distribution_final"),
  width = 11,
  height = 7
)

# -----------------------------
# 9. Save method note
# -----------------------------

method_note <- data.frame(
  Item = c(
    "Candidate universe",
    "VX correction",
    "Score rationale",
    "Gene-input record",
    "Next validation"
  ),
  Description = c(
    "All 61 refined perturbagens are retained as the working candidate universe.",
    "VX was reclassified as exposure/toxicant because CID 39793 corresponds to VX nerve agent.",
    "The base score combines enrichment strength, statistical significance, gene overlap, repeated support, core-module overlap, anchor-gene overlap and library support; the adjusted score applies a pharmacological suitability factor.",
    "Gene input files are saved to show the active-TB-up, post-therapy-up, anchor and core-module gene sets used in perturbation analysis.",
    "Use the candidate validation table for CMap/CLUE, CTD, DGIdb, DrugShot and target-axis mapping."
  )
)

write_csv(
  method_note,
  file.path(log_dir, "final_candidate_space_method_note.csv")
)

sink(file.path(log_dir, "session_info_final_candidate_space.txt"))
sessionInfo()
sink()

cat("\nFinal candidate-space outputs completed successfully.\n\n")
cat("All 61 candidate table:\n")
cat(file.path(table_dir, "01_all_61_refined_candidate_perturbagens_final.csv"), "\n\n")
cat("Candidate validation table:\n")
cat(file.path(table_dir, "02_candidates_for_CMap_CTD_target_mapping.csv"), "\n\n")
cat("Gene input summary:\n")
cat(file.path(table_dir, "06_gene_input_summary_for_perturbation_analysis.csv"), "\n\n")
cat("Figures saved in:\n")
cat(fig_dir, "\n")
...........
..............
############################################################
# 07_final_clean_candidate_figures_with_gene_support.R
# Clean final figures for 61 candidate perturbagens
# Includes top-25 bar plot, all-61 bubble plot, and gene-support heatmap
############################################################

# -----------------------------
# 1. Load packages
# -----------------------------

packages <- c(
  "dplyr", "readr", "stringr", "ggplot2",
  "forcats", "scales", "tidyr", "openxlsx"
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

out_dir <- file.path(project_dir, "results", "final_candidate_figures_gene_support")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Find input file
# -----------------------------

find_file_recursive <- function(project_dir, file_pattern, required = TRUE) {
  
  matched <- list.files(
    project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matched <- matched[!grepl("final_candidate_figures_gene_support", matched)]
  
  if (length(matched) == 0 && required) {
    stop(paste("File not found:", file_pattern))
  }
  
  if (length(matched) == 0 && !required) {
    return(NA_character_)
  }
  
  return(matched[1])
}

candidate_file <- find_file_recursive(
  project_dir,
  "01_all_61_refined_candidate_perturbagens_final\\.csv$",
  required = TRUE
)

cat("Using candidate file:\n")
cat(candidate_file, "\n\n")

cand <- readr::read_csv(candidate_file, show_col_types = FALSE)

# -----------------------------
# 4. Check required columns
# -----------------------------

required_cols <- c(
  "Refined_Candidate_Name",
  "Candidate_Group",
  "Broad_Target_Axis",
  "Priority_Score_100",
  "Suitability_Adjusted_Score",
  "Supporting_Signatures",
  "Max_Anchor_Hit_Count",
  "Max_Core8_Hit_Count",
  "Anchor_Hits",
  "Core8_Hits"
)

for (col in required_cols) {
  if (!col %in% names(cand)) {
    stop(paste("Required column missing:", col))
  }
}

# -----------------------------
# 5. Clean labels
# -----------------------------

clean_candidate_label <- function(x) {
  
  x <- as.character(x)
  
  x <- dplyr::case_when(
    x == "1,25 Dihydroxyvitamin D" ~ "1,25-dihydroxyvitamin D",
    x == "2,3-Dimethoxy-1,4-Naphtoquinone (Dmnq)" ~ "DMNQ",
    x == "Ursodeoxycholic Acid" ~ "Ursodeoxycholic acid",
    x == "Angiotensin Ii" ~ "Angiotensin II",
    x == "FTORAFUR" ~ "Ftorafur",
    x == "NORETHINDRONE" ~ "Norethindrone",
    x == "MEDROXYPROGESTERONE" ~ "Medroxyprogesterone",
    TRUE ~ x
  )
  
  return(x)
}

clean_axis_label <- function(x) {
  
  x <- as.character(x)
  
  x <- dplyr::case_when(
    is.na(x) ~ "Other / pending annotation",
    x == "Needs biological mapping" ~ "Other / pending annotation",
    x == "Excluded toxicant / nerve agent" ~ "Excluded toxicant",
    TRUE ~ x
  )
  
  return(x)
}

cand_clean <- cand %>%
  mutate(
    Display_Name = sapply(Refined_Candidate_Name, clean_candidate_label),
    Broad_Target_Axis_Clean = sapply(Broad_Target_Axis, clean_axis_label),
    Candidate_Group = case_when(
      Refined_Candidate_Name == "VX" ~ "Exposure / toxicant",
      TRUE ~ Candidate_Group
    ),
    Broad_Target_Axis_Clean = case_when(
      Refined_Candidate_Name == "VX" ~ "Excluded toxicant",
      TRUE ~ Broad_Target_Axis_Clean
    )
  ) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

# -----------------------------
# 6. Save cleaned all-61 table
# -----------------------------

write_csv(
  cand_clean,
  file.path(table_dir, "01_all_61_candidates_clean_labels.csv")
)

candidate_for_main <- cand_clean %>%
  filter(Candidate_Group %in% c(
    "Candidate compound",
    "Mapped LINCS compound",
    "Biologic / reference"
  )) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

write_csv(
  candidate_for_main,
  file.path(table_dir, "02_candidates_for_main_validation_clean_labels.csv")
)

# -----------------------------
# 7. Group and target-axis summaries
# -----------------------------

group_summary <- cand_clean %>%
  count(Candidate_Group, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

axis_summary <- cand_clean %>%
  count(Broad_Target_Axis_Clean, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  group_summary,
  file.path(table_dir, "03_candidate_group_summary_clean.csv")
)

write_csv(
  axis_summary,
  file.path(table_dir, "04_target_axis_summary_clean.csv")
)

# -----------------------------
# 8. Plot saving function
# -----------------------------

save_plot_versions <- function(plot_obj, file_base, width = 12, height = 8) {
  
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

# -----------------------------
# 9. Figure colors
# -----------------------------

group_colors <- c(
  "Candidate compound" = "#00B050",
  "Mapped LINCS compound" = "#7BC67B",
  "Biologic / reference" = "#00AEB3",
  "Endogenous / metabolic" = "#5B9BD5",
  "Exposure / toxicant" = "#F8766D"
)

# -----------------------------
# 10. Figure 1: Top-25 candidate perturbagens
# -----------------------------

top25 <- candidate_for_main %>%
  head(25) %>%
  mutate(
    Label = stringr::str_wrap(Display_Name, width = 32),
    Label = forcats::fct_reorder(Label, Suitability_Adjusted_Score)
  )

p_top25 <- ggplot(
  top25,
  aes(x = Label, y = Suitability_Adjusted_Score, fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(
    aes(label = round(Suitability_Adjusted_Score, 2)),
    hjust = -0.08,
    size = 4.4,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Candidate perturbagens after suitability adjustment",
    x = NULL,
    y = "Suitability-adjusted score",
    fill = NULL
  ) +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 14, face = "plain"),
    axis.text.y = element_text(size = 13.5, face = "plain"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12.5, face = "plain"),
    legend.key.size = unit(0.45, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 45, 35, 15)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

save_plot_versions(
  p_top25,
  file.path(fig_dir, "Figure_1_top25_candidate_perturbagens_clean_final"),
  width = 11.5,
  height = 10.5
)

# -----------------------------
# 11. Figure 2: All-61 candidate perturbagen space
# -----------------------------

bubble_data <- cand_clean %>%
  mutate(
    Candidate_Group = factor(
      Candidate_Group,
      levels = c(
        "Candidate compound",
        "Mapped LINCS compound",
        "Biologic / reference",
        "Endogenous / metabolic",
        "Exposure / toxicant"
      )
    ),
    Axis_Label = stringr::str_wrap(Broad_Target_Axis_Clean, width = 26)
  )

p_bubble <- ggplot(
  bubble_data,
  aes(
    x = Axis_Label,
    y = Candidate_Group,
    size = Supporting_Signatures,
    color = Suitability_Adjusted_Score
  )
) +
  geom_point(
    alpha = 0.80,
    position = position_jitter(width = 0.15, height = 0.10, seed = 10)
  ) +
  scale_size_continuous(range = c(3, 9)) +
  scale_color_gradient(low = "#D9EAF7", high = "#E95A4F") +
  labs(
    title = "Candidate perturbagen space from active-TB signature suppression",
    x = NULL,
    y = NULL,
    size = "Supporting\nsignatures",
    color = "Adjusted\nscore"
  ) +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12.5, face = "plain", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13.5, face = "plain"),
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 11.5),
    panel.grid = element_blank(),
    plot.margin = margin(15, 35, 45, 15)
  )

save_plot_versions(
  p_bubble,
  file.path(fig_dir, "Figure_2_all61_candidate_space_bubble_clean_final"),
  width = 14,
  height = 8
)

# -----------------------------
# 12. Figure 3: Candidate group distribution
# -----------------------------

p_group <- ggplot(
  group_summary,
  aes(x = forcats::fct_reorder(Candidate_Group, Candidate_Count), y = Candidate_Count)
) +
  geom_col(fill = "#4B8BBE", color = "black", linewidth = 0.35, width = 0.7) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.1,
    size = 5.2,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "Candidate group distribution",
    x = NULL,
    y = "Number of candidates"
  ) +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 14, face = "plain"),
    axis.text.y = element_text(size = 14, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(15, 45, 15, 15)
  )

save_plot_versions(
  p_group,
  file.path(fig_dir, "Figure_3_candidate_group_distribution_clean_final"),
  width = 10,
  height = 6.5
)

# -----------------------------
# 13. Gene-support heatmap
# -----------------------------

core_genes <- c(
  "STAT1", "GBP1", "IFIT3", "IFI35",
  "IRF7", "PARP9", "IFIT2", "IFI6"
)

extract_gene_hits <- function(x) {
  
  if (is.na(x) || x == "") {
    return(character(0))
  }
  
  genes <- unlist(stringr::str_split(x, ";|,|/|\\s+"))
  genes <- toupper(trimws(genes))
  genes <- genes[genes != ""]
  genes <- unique(genes)
  
  return(genes)
}

top30_gene <- cand_clean %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P) %>%
  head(30) %>%
  mutate(
    Display_Name = stringr::str_wrap(Display_Name, width = 28)
  )

heatmap_rows <- list()

for (i in seq_len(nrow(top30_gene))) {
  
  hits <- extract_gene_hits(top30_gene$Core8_Hits[i])
  
  for (g in core_genes) {
    heatmap_rows[[length(heatmap_rows) + 1]] <- data.frame(
      Candidate = top30_gene$Display_Name[i],
      Gene = g,
      Hit = ifelse(g %in% hits, 1, 0),
      Score = top30_gene$Suitability_Adjusted_Score[i]
    )
  }
}

heatmap_df <- bind_rows(heatmap_rows)

heatmap_df <- heatmap_df %>%
  mutate(
    Candidate = factor(Candidate, levels = rev(unique(top30_gene$Display_Name))),
    Gene = factor(Gene, levels = core_genes)
  )

write_csv(
  heatmap_df,
  file.path(table_dir, "05_top30_candidate_core_gene_support_matrix.csv")
)

p_heatmap <- ggplot(
  heatmap_df,
  aes(x = Gene, y = Candidate, fill = factor(Hit))
) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_manual(
    values = c("0" = "#F2F2F2", "1" = "#2C7FB8"),
    labels = c("No overlap", "Overlap")
  ) +
  labs(
    title = "Core host-response gene support among prioritized perturbagens",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12.5, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11.5, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(15, 25, 35, 15)
  )

save_plot_versions(
  p_heatmap,
  file.path(fig_dir, "Figure_4_top30_candidate_core_gene_support_heatmap"),
  width = 10,
  height = 10.5
)

# -----------------------------
# 14. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "final_candidate_figures_and_gene_support_tables.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "All_61_clean")
writeData(wb, "All_61_clean", cand_clean)

addWorksheet(wb, "Candidate_validation")
writeData(wb, "Candidate_validation", candidate_for_main)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Axis_summary")
writeData(wb, "Axis_summary", axis_summary)

addWorksheet(wb, "Gene_support_top30")
writeData(wb, "Gene_support_top30", heatmap_df)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(wb)) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:100, gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:80, widths = "auto")
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# -----------------------------
# 15. Save notes
# -----------------------------

method_note <- data.frame(
  Item = c(
    "All candidate universe",
    "Top-25 bar plot",
    "All-61 bubble plot",
    "Gene-support heatmap",
    "Next step"
  ),
  Description = c(
    "All 61 refined candidate perturbagens were retained in the final candidate table.",
    "The top-25 bar plot is used only for readability and does not replace the full 61-candidate table.",
    "The all-61 bubble plot summarizes the full candidate space by candidate group, target-axis annotation, supporting signatures and adjusted score.",
    "Core gene support is shown separately using STAT1, GBP1, IFIT3, IFI35, IRF7, PARP9, IFIT2 and IFI6.",
    "Use the candidate validation table for CMap/CLUE tau validation and CTD/DGIdb/DrugShot evidence mapping."
  )
)

write_csv(
  method_note,
  file.path(log_dir, "final_figure_and_gene_support_notes.csv")
)

sink(file.path(log_dir, "session_info_final_candidate_figures_gene_support.txt"))
sessionInfo()
sink()

cat("\nFinal clean figures and gene-support outputs completed successfully.\n\n")
cat("All 61 clean candidate table:\n")
cat(file.path(table_dir, "01_all_61_candidates_clean_labels.csv"), "\n\n")
cat("Candidate validation table:\n")
cat(file.path(table_dir, "02_candidates_for_main_validation_clean_labels.csv"), "\n\n")
cat("Figures saved in:\n")
cat(fig_dir, "\n\n")
cat("Excel workbook:\n")
cat(xlsx_file, "\n")