############################################################
# 08_final_publication_figures_no_cut_labels.R
# Final clean figures with corrected labels and no cut text
############################################################

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
# 1. Project directory
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"

out_dir <- file.path(project_dir, "results", "final_publication_figures")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Find input file
# -----------------------------

find_file_recursive <- function(project_dir, file_pattern, required = TRUE) {
  
  matched <- list.files(
    project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matched <- matched[!grepl("final_publication_figures", matched)]
  
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
# 3. Required columns
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
# 4. Clean labels
# -----------------------------

clean_candidate_label <- function(x) {
  
  x <- as.character(x)
  
  x <- dplyr::case_when(
    x == "1,25 Dihydroxyvitamin D" ~ "1,25-dihydroxyvitamin D",
    x == "1,25-dihydroxyvitamin D" ~ "1,25-dihydroxyvitamin D",
    x == "2,3-Dimethoxy-1,4-Naphtoquinone (Dmnq)" ~ "DMNQ",
    x == "2,3-Dimethoxy-1,4-Naphthoquinone (Dmnq)" ~ "DMNQ",
    x == "Ursodeoxycholic Acid" ~ "Ursodeoxycholic acid",
    x == "Angiotensin Ii" ~ "Angiotensin II",
    x == "FTORAFUR" ~ "Ftorafur",
    x == "NORETHINDRONE" ~ "Norethindrone",
    x == "MEDROXYPROGESTERONE" ~ "Medroxyprogesterone",
    x == "Vx" ~ "VX",
    TRUE ~ x
  )
  
  return(x)
}

clean_axis_label <- function(x) {
  
  x <- as.character(x)
  
  x <- dplyr::case_when(
    is.na(x) ~ "Other / unassigned axis",
    x == "Needs biological mapping" ~ "Other / unassigned axis",
    x == "Other / pending annotation" ~ "Other / unassigned axis",
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
      Display_Name == "VX" ~ "Exposure / toxicant",
      TRUE ~ Candidate_Group
    ),
    Broad_Target_Axis_Clean = case_when(
      Display_Name == "VX" ~ "Excluded toxicant",
      TRUE ~ Broad_Target_Axis_Clean
    )
  ) %>%
  arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P)

# -----------------------------
# 5. Save cleaned tables
# -----------------------------

write_csv(
  cand_clean,
  file.path(table_dir, "01_all_61_candidates_final_clean_labels.csv")
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
  file.path(table_dir, "02_candidates_for_validation_final_clean_labels.csv")
)

group_summary <- cand_clean %>%
  count(Candidate_Group, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

axis_summary <- cand_clean %>%
  count(Broad_Target_Axis_Clean, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  group_summary,
  file.path(table_dir, "03_candidate_group_summary_final.csv")
)

write_csv(
  axis_summary,
  file.path(table_dir, "04_target_axis_summary_final.csv")
)

# -----------------------------
# 6. Figure saving function
# -----------------------------

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
# 7. Colors
# -----------------------------

group_colors <- c(
  "Candidate compound" = "#00B050",
  "Mapped LINCS compound" = "#7BC67B",
  "Biologic / reference" = "#00AEB3",
  "Endogenous / metabolic" = "#5B9BD5",
  "Exposure / toxicant" = "#F8766D"
)

# -----------------------------
# 8. Figure 1: Top-25 candidate bar plot
# -----------------------------

top25 <- candidate_for_main %>%
  head(25) %>%
  mutate(
    Label = stringr::str_wrap(Display_Name, width = 30),
    Label = forcats::fct_reorder(Label, Suitability_Adjusted_Score)
  )

p_top25 <- ggplot(
  top25,
  aes(x = Label, y = Suitability_Adjusted_Score, fill = Candidate_Group)
) +
  geom_col(color = "black", linewidth = 0.32, width = 0.70) +
  geom_text(
    aes(label = round(Suitability_Adjusted_Score, 2)),
    hjust = -0.08,
    size = 4.2,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.24)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "Candidate perturbagens after suitability adjustment",
    x = NULL,
    y = "Suitability-adjusted score",
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.text.x = element_text(size = 13, face = "plain"),
    axis.text.y = element_text(size = 12.5, face = "plain"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 11.5, face = "plain"),
    legend.key.size = unit(0.40, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 45, 35, 12)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

save_plot_versions(
  p_top25,
  file.path(fig_dir, "Figure_1_top25_candidate_perturbagens_final_no_cut"),
  width = 12,
  height = 10
)

# -----------------------------
# 9. Figure 2: All-61 candidate space bubble plot
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
    Axis_Label = stringr::str_wrap(Broad_Target_Axis_Clean, width = 24)
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
    alpha = 0.82,
    position = position_jitter(width = 0.12, height = 0.08, seed = 10)
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
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 11.5, face = "plain", angle = 40, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12.5, face = "plain"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    panel.grid = element_blank(),
    plot.margin = margin(12, 35, 55, 12)
  )

save_plot_versions(
  p_bubble,
  file.path(fig_dir, "Figure_2_all61_candidate_space_bubble_final_no_cut"),
  width = 15.5,
  height = 8.5
)

# -----------------------------
# 10. Figure 3: Candidate group distribution
# -----------------------------

p_group <- ggplot(
  group_summary,
  aes(x = forcats::fct_reorder(Candidate_Group, Candidate_Count), y = Candidate_Count)
) +
  geom_col(fill = "#4B8BBE", color = "black", linewidth = 0.32, width = 0.68) +
  geom_text(
    aes(label = Candidate_Count),
    hjust = -0.10,
    size = 5.0,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.24)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "Candidate group distribution",
    x = NULL,
    y = "Number of candidates"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.text.x = element_text(size = 13, face = "plain"),
    axis.text.y = element_text(size = 13, face = "plain"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 45, 12, 12)
  )

save_plot_versions(
  p_group,
  file.path(fig_dir, "Figure_3_candidate_group_distribution_final_no_cut"),
  width = 10.5,
  height = 6.2
)

# -----------------------------
# 11. Gene-support heatmaps
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

make_gene_heatmap_data <- function(input_df, n_top = 30) {
  
  top_gene <- input_df %>%
    arrange(desc(Suitability_Adjusted_Score), Best_Adjusted_P) %>%
    head(n_top) %>%
    mutate(Display_Name = stringr::str_wrap(Display_Name, width = 28))
  
  heatmap_rows <- list()
  
  for (i in seq_len(nrow(top_gene))) {
    
    hits <- extract_gene_hits(top_gene$Core8_Hits[i])
    
    for (g in core_genes) {
      heatmap_rows[[length(heatmap_rows) + 1]] <- data.frame(
        Candidate = top_gene$Display_Name[i],
        Gene = g,
        Hit = ifelse(g %in% hits, 1, 0),
        Score = top_gene$Suitability_Adjusted_Score[i]
      )
    }
  }
  
  heatmap_df <- bind_rows(heatmap_rows)
  
  heatmap_df <- heatmap_df %>%
    mutate(
      Candidate = factor(Candidate, levels = rev(unique(top_gene$Display_Name))),
      Gene = factor(Gene, levels = core_genes)
    )
  
  return(heatmap_df)
}

heatmap_top30 <- make_gene_heatmap_data(cand_clean, n_top = 30)

write_csv(
  heatmap_top30,
  file.path(table_dir, "05_top30_candidate_core_gene_support_matrix.csv")
)

p_heatmap30 <- ggplot(
  heatmap_top30,
  aes(x = Gene, y = Candidate, fill = factor(Hit))
) +
  geom_tile(color = "white", linewidth = 0.40) +
  scale_fill_manual(
    values = c("0" = "#F2F2F2", "1" = "#2C7FB8"),
    labels = c("No overlap", "Overlap")
  ) +
  labs(
    title = "Core gene support of prioritized perturbagens",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 11.5, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10.5, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.margin = margin(12, 25, 35, 12)
  )

save_plot_versions(
  p_heatmap30,
  file.path(fig_dir, "Figure_4_top30_core_gene_support_heatmap_final_no_cut"),
  width = 10.5,
  height = 10.5
)

# -----------------------------
# 12. Supplementary all-61 heatmap
# -----------------------------

heatmap_all61 <- make_gene_heatmap_data(cand_clean, n_top = 61)

write_csv(
  heatmap_all61,
  file.path(table_dir, "06_all61_candidate_core_gene_support_matrix.csv")
)

p_heatmap61 <- ggplot(
  heatmap_all61,
  aes(x = Gene, y = Candidate, fill = factor(Hit))
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(
    values = c("0" = "#F2F2F2", "1" = "#2C7FB8"),
    labels = c("No overlap", "Overlap")
  ) +
  labs(
    title = "Core gene support across all refined perturbagens",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10.5, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8.5, face = "plain"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(12, 25, 35, 12)
  )

save_plot_versions(
  p_heatmap61,
  file.path(fig_dir, "Figure_S_all61_core_gene_support_heatmap_final_no_cut"),
  width = 10.5,
  height = 17
)

# -----------------------------
# 13. Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "final_publication_figure_tables.xlsx")

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
writeData(wb, "Gene_support_top30", heatmap_top30)

addWorksheet(wb, "Gene_support_all61")
writeData(wb, "Gene_support_all61", heatmap_all61)

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
# 14. Method note
# -----------------------------

method_note <- data.frame(
  Item = c(
    "All candidate universe",
    "Bar plot",
    "Bubble plot",
    "Gene support",
    "Axis label correction"
  ),
  Description = c(
    "All 61 refined candidate perturbagens were retained in the final candidate table.",
    "The top-25 bar plot is used for readability and does not replace the full candidate table.",
    "The bubble plot summarizes all 61 candidates by group, target-axis annotation, supporting signatures and adjusted score.",
    "Core gene support is shown separately using STAT1, GBP1, IFIT3, IFI35, IRF7, PARP9, IFIT2 and IFI6.",
    "Unmapped target-axis terms were labelled as Other / unassigned axis."
  )
)

write_csv(
  method_note,
  file.path(log_dir, "final_publication_figure_notes.csv")
)

sink(file.path(log_dir, "session_info_final_publication_figures.txt"))
sessionInfo()
sink()

# -----------------------------
# 15. Console output
# -----------------------------

cat("\nFinal publication figures completed successfully.\n\n")
cat("Clean tables saved in:\n")
cat(table_dir, "\n\n")
cat("Figures saved in:\n")
cat(fig_dir, "\n\n")
cat("Excel workbook:\n")
cat(xlsx_file, "\n")