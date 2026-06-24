############################################################
# 10_read_CLUE_ncs_and_match_candidates.R
# Read CLUE ncs.gct and match scores with 61 candidate perturbagens
############################################################

packages <- c("dplyr", "readr", "stringr", "tidyr", "openxlsx")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Folders
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"
clue_dir    <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION/results/CLUE"

out_dir   <- file.path(project_dir, "results", "CLUE_NCS_integration")
table_dir <- file.path(out_dir, "tables")
log_dir   <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Helper functions
# -----------------------------

find_file <- function(folder, pattern) {
  f <- list.files(folder, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) stop(paste("File not found:", pattern))
  f[1]
}

normalize_name <- function(x) {
  x <- tolower(as.character(x))
  x <- str_replace_all(x, "[^a-z0-9]+", "")
  x
}

classify_reversal <- function(x) {
  case_when(
    is.na(x) ~ "No CMap match",
    x <= -90 ~ "Strong reversal",
    x > -90 & x <= -70 ~ "Moderate reversal",
    x > -70 & x <= 0 ~ "Weak or no reversal",
    x > 0 ~ "Disease-mimicking",
    TRUE ~ "Unclassified"
  )
}

# -----------------------------
# 3. Find files
# -----------------------------

ncs_file <- find_file(clue_dir, "^ncs\\.gct$")

candidate_file <- find_file(
  project_dir,
  "01_all_61_candidates_final_clean_labels\\.csv$"
)

cat("Using ncs.gct:\n", ncs_file, "\n\n")
cat("Using candidate table:\n", candidate_file, "\n\n")

# -----------------------------
# 4. Read GCT file
# -----------------------------

read_gct <- function(gct_file) {
  
  lines <- readLines(gct_file, warn = FALSE)
  version <- trimws(lines[1])
  
  if (version == "#1.3") {
    
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
    score_cols <- tail(colnames(gct), n_data_cols)
    
    gct_long <- gct %>%
      pivot_longer(
        cols = all_of(score_cols),
        names_to = "Query_ID",
        values_to = "NCS_raw"
      ) %>%
      mutate(NCS_raw = suppressWarnings(as.numeric(NCS_raw)))
    
    return(gct_long)
  }
  
  if (version == "#1.2") {
    
    gct <- read.delim(
      gct_file,
      skip = 2,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
    
    score_cols <- setdiff(colnames(gct), c("Name", "Description"))
    
    gct_long <- gct %>%
      pivot_longer(
        cols = all_of(score_cols),
        names_to = "Query_ID",
        values_to = "NCS_raw"
      ) %>%
      mutate(NCS_raw = suppressWarnings(as.numeric(NCS_raw)))
    
    return(gct_long)
  }
  
  stop("Unsupported GCT version.")
}

ncs_long <- read_gct(ncs_file)

# If score is -1 to +1, convert to -100 to +100
max_abs <- max(abs(ncs_long$NCS_raw), na.rm = TRUE)

ncs_long <- ncs_long %>%
  mutate(
    NCS_100 = ifelse(max_abs <= 1.5, NCS_raw * 100, NCS_raw)
  )

write_csv(
  ncs_long,
  file.path(table_dir, "01_CLUE_NCS_long_format.csv")
)

cat("NCS rows extracted:", nrow(ncs_long), "\n")
cat("Columns found:\n")
print(colnames(ncs_long))

# -----------------------------
# 5. Detect perturbagen name column
# -----------------------------

possible_name_cols <- c(
  "pert_iname", "pert_desc", "pert_name", "cmap_name",
  "Name", "name", "Description", "description", "id"
)

name_col <- possible_name_cols[possible_name_cols %in% colnames(ncs_long)][1]

if (is.na(name_col)) {
  stop("Could not detect perturbagen name column. Check 01_CLUE_NCS_long_format.csv")
}

cat("\nUsing perturbagen name column:", name_col, "\n\n")

ncs_clean <- ncs_long %>%
  mutate(
    CMap_Perturbagen_Name = .data[[name_col]],
    CMap_Match_Key = normalize_name(CMap_Perturbagen_Name)
  ) %>%
  filter(!is.na(NCS_100))

write_csv(
  ncs_clean,
  file.path(table_dir, "02_CLUE_NCS_clean_with_match_keys.csv")
)

# -----------------------------
# 6. Top opposing and similar perturbagens
# -----------------------------

top_opposing <- ncs_clean %>%
  arrange(NCS_100) %>%
  head(100)

top_similar <- ncs_clean %>%
  arrange(desc(NCS_100)) %>%
  head(100)

write_csv(
  top_opposing,
  file.path(table_dir, "03_top100_opposing_negative_connectivity_perturbagens.csv")
)

write_csv(
  top_similar,
  file.path(table_dir, "04_top100_similar_positive_connectivity_perturbagens.csv")
)

# -----------------------------
# 7. Read candidate table
# -----------------------------

cand <- read_csv(candidate_file, show_col_types = FALSE)

name_column <- if ("Display_Name" %in% colnames(cand)) {
  "Display_Name"
} else {
  "Refined_Candidate_Name"
}

cand <- cand %>%
  mutate(
    Candidate_Display_Name = .data[[name_column]],
    Candidate_Match_Key = normalize_name(Candidate_Display_Name)
  )

candidate_keys <- cand %>%
  select(Candidate_Display_Name, Candidate_Match_Key) %>%
  distinct()

# Extra synonyms for better matching
extra_keys <- data.frame(
  Candidate_Display_Name = c(
    "1,25-dihydroxyvitamin D",
    "Vitamin D3",
    "Ursodeoxycholic acid",
    "PLX4720",
    "DMNQ"
  ),
  Candidate_Match_Key = normalize_name(c(
    "calcitriol",
    "cholecalciferol",
    "ursodiol",
    "plx-4720",
    "2,3-dimethoxy-1,4-naphthoquinone"
  ))
)

candidate_keys <- bind_rows(candidate_keys, extra_keys) %>%
  distinct()

# -----------------------------
# 8. Match CMap perturbagens with 61 candidates
# -----------------------------

matched_raw <- ncs_clean %>%
  inner_join(
    candidate_keys,
    by = c("CMap_Match_Key" = "Candidate_Match_Key")
  )

write_csv(
  matched_raw,
  file.path(table_dir, "05_raw_CMap_matches_to_candidate_names.csv")
)

candidate_best <- matched_raw %>%
  group_by(Candidate_Display_Name) %>%
  summarise(
    Best_Negative_NCS = min(NCS_100, na.rm = TRUE),
    Best_Positive_NCS = max(NCS_100, na.rm = TRUE),
    Mean_NCS = mean(NCS_100, na.rm = TRUE),
    Number_of_CMap_Matches = n(),
    Best_Opposing_Perturbagen = CMap_Perturbagen_Name[which.min(NCS_100)][1],
    Best_Opposing_Query_ID = Query_ID[which.min(NCS_100)][1],
    .groups = "drop"
  ) %>%
  mutate(
    CMap_Reversal_Category = classify_reversal(Best_Negative_NCS)
  )

cand_integrated <- cand %>%
  left_join(candidate_best, by = "Candidate_Display_Name") %>%
  mutate(
    CMap_Reversal_Category = ifelse(
      is.na(CMap_Reversal_Category),
      "No CMap match",
      CMap_Reversal_Category
    )
  ) %>%
  arrange(Best_Negative_NCS, desc(Suitability_Adjusted_Score))

write_csv(
  cand_integrated,
  file.path(table_dir, "06_candidate_table_integrated_with_CLUE_NCS.csv")
)

reversal_summary <- cand_integrated %>%
  count(CMap_Reversal_Category, name = "Candidate_Count") %>%
  arrange(desc(Candidate_Count))

write_csv(
  reversal_summary,
  file.path(table_dir, "07_CMap_reversal_category_summary.csv")
)

# -----------------------------
# 9. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "CLUE_NCS_candidate_integration.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "NCS_long")
writeData(wb, "NCS_long", ncs_long)

addWorksheet(wb, "Top_opposing")
writeData(wb, "Top_opposing", top_opposing)

addWorksheet(wb, "Top_similar")
writeData(wb, "Top_similar", top_similar)

addWorksheet(wb, "Raw_matches")
writeData(wb, "Raw_matches", matched_raw)

addWorksheet(wb, "Candidate_integrated")
writeData(wb, "Candidate_integrated", cand_integrated)

addWorksheet(wb, "Reversal_summary")
writeData(wb, "Reversal_summary", reversal_summary)

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

sink(file.path(log_dir, "session_info_CLUE_NCS_integration.txt"))
sessionInfo()
sink()

# -----------------------------
# 10. Console output
# -----------------------------

cat("\nCLUE NCS integration completed successfully.\n\n")

cat("Send me these files:\n")
cat(file.path(table_dir, "03_top100_opposing_negative_connectivity_perturbagens.csv"), "\n")
cat(file.path(table_dir, "05_raw_CMap_matches_to_candidate_names.csv"), "\n")
cat(file.path(table_dir, "06_candidate_table_integrated_with_CLUE_NCS.csv"), "\n")
cat(file.path(table_dir, "07_CMap_reversal_category_summary.csv"), "\n\n")

cat("Reversal summary:\n")
print(reversal_summary)
.....
.....
############################################################
# 10A_check_CLUE_folder_contents.R
# Check CLUE output folder structure before NCS integration
############################################################

packages <- c("dplyr", "readr", "stringr", "tidyr", "openxlsx")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Set project and CLUE folder
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"
clue_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION/results/CLUE"

out_dir <- file.path(project_dir, "results", "CLUE_folder_check")
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nProject directory:\n")
cat(project_dir, "\n\n")

cat("CLUE folder:\n")
cat(clue_dir, "\n\n")

if (!dir.exists(clue_dir)) {
  stop("CLUE folder not found. Please check the clue_dir path.")
}

# -----------------------------
# 2. Helper functions
# -----------------------------

format_size_mb <- function(bytes) {
  round(bytes / (1024^2), 4)
}

get_extension <- function(x) {
  ext <- tools::file_ext(x)
  ext[ext == ""] <- "no_extension"
  tolower(ext)
}

safe_read_lines <- function(file_path, n = 20) {
  out <- tryCatch(
    readLines(file_path, n = n, warn = FALSE),
    error = function(e) paste("READ_ERROR:", e$message)
  )
  return(out)
}

detect_text_file_type <- function(file_path) {
  
  ext <- tolower(tools::file_ext(file_path))
  base <- basename(file_path)
  
  if (base == "ncs.gct") return("MAIN_RESULT_NCS_GCT")
  if (ext == "gct") return("GCT_matrix")
  if (ext == "gmt") return("Gene_set_input_GMT")
  if (ext %in% c("yaml", "yml")) return("Query_config_YAML")
  if (ext %in% c("csv", "tsv", "txt")) return("Text_table_or_log")
  if (ext %in% c("json")) return("JSON_metadata")
  if (ext %in% c("pdf", "png", "jpg", "jpeg", "tiff")) return("Figure_or_document")
  if (ext %in% c("zip", "gz")) return("Compressed_file")
  
  return("Other")
}

check_gct_header <- function(file_path) {
  
  lines <- safe_read_lines(file_path, n = 8)
  
  if (length(lines) == 0) {
    return(data.frame(
      File = file_path,
      Version = NA,
      Rows = NA,
      Columns = NA,
      Row_Metadata_Columns = NA,
      Column_Metadata_Rows = NA,
      Header_Preview = NA,
      Status = "Could not read"
    ))
  }
  
  version <- trimws(lines[1])
  
  if (version == "#1.3") {
    
    dims <- suppressWarnings(as.integer(strsplit(lines[2], "\t")[[1]]))
    header <- if (length(lines) >= 3) lines[3] else NA
    
    return(data.frame(
      File = file_path,
      Version = version,
      Rows = dims[1],
      Columns = dims[2],
      Row_Metadata_Columns = dims[3],
      Column_Metadata_Rows = dims[4],
      Header_Preview = substr(header, 1, 500),
      Status = "Readable GCT 1.3"
    ))
  }
  
  if (version == "#1.2") {
    
    dims <- suppressWarnings(as.integer(strsplit(lines[2], "\t")[[1]]))
    header <- if (length(lines) >= 3) lines[3] else NA
    
    return(data.frame(
      File = file_path,
      Version = version,
      Rows = dims[1],
      Columns = dims[2],
      Row_Metadata_Columns = NA,
      Column_Metadata_Rows = NA,
      Header_Preview = substr(header, 1, 500),
      Status = "Readable GCT 1.2"
    ))
  }
  
  return(data.frame(
    File = file_path,
    Version = version,
    Rows = NA,
    Columns = NA,
    Row_Metadata_Columns = NA,
    Column_Metadata_Rows = NA,
    Header_Preview = paste(lines[1:min(length(lines), 3)], collapse = " | "),
    Status = "Not standard GCT"
  ))
}

preview_small_text_file <- function(file_path, n = 5) {
  
  lines <- safe_read_lines(file_path, n = n)
  
  data.frame(
    File = file_path,
    Line_Number = seq_along(lines),
    Preview = lines,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 3. Create full file inventory
# -----------------------------

all_files <- list.files(
  clue_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)

file_info <- file.info(all_files)

inventory <- data.frame(
  File_Path = all_files,
  Relative_Path = stringr::str_replace(all_files, fixed(clue_dir), ""),
  File_Name = basename(all_files),
  Parent_Folder = basename(dirname(all_files)),
  Extension = get_extension(all_files),
  File_Type = sapply(all_files, detect_text_file_type),
  Size_MB = format_size_mb(file_info$size),
  Modified_Time = as.character(file_info$mtime),
  stringsAsFactors = FALSE
) %>%
  arrange(Parent_Folder, File_Name)

write_csv(
  inventory,
  file.path(table_dir, "01_CLUE_folder_file_inventory.csv")
)

cat("Total files found:", nrow(inventory), "\n\n")

# -----------------------------
# 4. Folder-level summary
# -----------------------------

folder_summary <- inventory %>%
  group_by(Parent_Folder) %>%
  summarise(
    Number_of_Files = n(),
    Total_Size_MB = round(sum(Size_MB, na.rm = TRUE), 4),
    File_Types = paste(sort(unique(File_Type)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(Number_of_Files))

write_csv(
  folder_summary,
  file.path(table_dir, "02_CLUE_folder_summary.csv")
)

extension_summary <- inventory %>%
  group_by(Extension, File_Type) %>%
  summarise(
    Number_of_Files = n(),
    Total_Size_MB = round(sum(Size_MB, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Number_of_Files))

write_csv(
  extension_summary,
  file.path(table_dir, "03_CLUE_file_extension_summary.csv")
)

# -----------------------------
# 5. Check required/useful files
# -----------------------------

required_patterns <- data.frame(
  Item = c(
    "Main CMap/CLUE score matrix",
    "Query configuration",
    "UP input gene set",
    "DOWN input gene set",
    "Matrices folder",
    "ARFS folder",
    "GSEA folder"
  ),
  Expected = c(
    "ncs.gct",
    "query_config.yaml",
    "up.gmt",
    "down.gmt",
    "matrices",
    "arfs",
    "gsea"
  ),
  Importance = c(
    "Essential for integration",
    "Important for record keeping",
    "Input record only",
    "Input record only",
    "Check for extra result matrices",
    "Optional/supporting",
    "Optional/supporting"
  ),
  stringsAsFactors = FALSE
)

required_check <- required_patterns %>%
  rowwise() %>%
  mutate(
    Found = any(grepl(Expected, inventory$File_Path, ignore.case = TRUE)) ||
      dir.exists(file.path(clue_dir, Expected)),
    Matching_Files = paste(
      inventory$File_Path[grepl(Expected, inventory$File_Path, ignore.case = TRUE)],
      collapse = " || "
    )
  ) %>%
  ungroup()

write_csv(
  required_check,
  file.path(table_dir, "04_required_CLUE_files_check.csv")
)

# -----------------------------
# 6. Check GCT files
# -----------------------------

gct_files <- inventory %>%
  filter(Extension == "gct") %>%
  pull(File_Path)

if (length(gct_files) > 0) {
  gct_header_check <- bind_rows(lapply(gct_files, check_gct_header))
} else {
  gct_header_check <- data.frame(
    File = NA,
    Version = NA,
    Rows = NA,
    Columns = NA,
    Row_Metadata_Columns = NA,
    Column_Metadata_Rows = NA,
    Header_Preview = NA,
    Status = "No GCT files found"
  )
}

write_csv(
  gct_header_check,
  file.path(table_dir, "05_GCT_header_check.csv")
)

# -----------------------------
# 7. Preview GMT/YAML/small text files
# -----------------------------

preview_files <- inventory %>%
  filter(Extension %in% c("gmt", "yaml", "yml", "txt", "csv", "tsv")) %>%
  filter(Size_MB < 5) %>%
  pull(File_Path)

if (length(preview_files) > 0) {
  preview_table <- bind_rows(lapply(preview_files, preview_small_text_file, n = 5))
} else {
  preview_table <- data.frame(
    File = NA,
    Line_Number = NA,
    Preview = "No small text files available for preview"
  )
}

write_csv(
  preview_table,
  file.path(table_dir, "06_small_text_file_preview.csv")
)

# -----------------------------
# 8. Identify likely score/result files
# -----------------------------

likely_score_files <- inventory %>%
  filter(
    grepl("ncs|score|tau|connect|result|rank|sig", File_Name, ignore.case = TRUE) |
      grepl("matrix|matrices", Parent_Folder, ignore.case = TRUE)
  ) %>%
  arrange(desc(Size_MB))

write_csv(
  likely_score_files,
  file.path(table_dir, "07_likely_score_or_result_files.csv")
)

# -----------------------------
# 9. Find candidate tables from project
# -----------------------------

candidate_files <- list.files(
  project_dir,
  pattern = "01_all_61_candidates_final_clean_labels\\.csv$|01_all_61_refined_candidate_perturbagens_final\\.csv$|02_candidates_for_validation_final_clean_labels\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

candidate_check <- data.frame(
  Candidate_File = candidate_files,
  File_Name = basename(candidate_files),
  Size_MB = ifelse(length(candidate_files) > 0, format_size_mb(file.info(candidate_files)$size), NA),
  stringsAsFactors = FALSE
)

if (length(candidate_files) == 0) {
  candidate_check <- data.frame(
    Candidate_File = NA,
    File_Name = NA,
    Size_MB = NA,
    Note = "No candidate table found"
  )
}

write_csv(
  candidate_check,
  file.path(table_dir, "08_candidate_table_check.csv")
)

# -----------------------------
# 10. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "CLUE_folder_content_check.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "File_inventory")
writeData(wb, "File_inventory", inventory)

addWorksheet(wb, "Folder_summary")
writeData(wb, "Folder_summary", folder_summary)

addWorksheet(wb, "Extension_summary")
writeData(wb, "Extension_summary", extension_summary)

addWorksheet(wb, "Required_files")
writeData(wb, "Required_files", required_check)

addWorksheet(wb, "GCT_header_check")
writeData(wb, "GCT_header_check", gct_header_check)

addWorksheet(wb, "Text_file_preview")
writeData(wb, "Text_file_preview", preview_table)

addWorksheet(wb, "Likely_score_files")
writeData(wb, "Likely_score_files", likely_score_files)

addWorksheet(wb, "Candidate_table_check")
writeData(wb, "Candidate_table_check", candidate_check)

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
# 11. Save recommendation text
# -----------------------------

recommendation <- c(
  "CLUE folder check recommendation",
  "",
  "Essential file for next analysis:",
  "- ncs.gct",
  "",
  "Useful record-keeping files:",
  "- query_config.yaml",
  "- up.gmt",
  "- down.gmt",
  "",
  "Folders to inspect:",
  "- matrices: may contain additional score matrices",
  "- arfs: optional/supporting CLUE outputs",
  "- gsea: optional/supporting enrichment outputs",
  "",
  "For the next CMap/NCS integration, use ncs.gct as the main score file.",
  "GMT files are input gene sets and should not be treated as result-score files.",
  "",
  "After running this script, send these files:",
  "01_CLUE_folder_file_inventory.csv",
  "04_required_CLUE_files_check.csv",
  "05_GCT_header_check.csv",
  "07_likely_score_or_result_files.csv"
)

writeLines(
  recommendation,
  file.path(log_dir, "CLUE_folder_check_recommendation.txt")
)

sink(file.path(log_dir, "session_info_CLUE_folder_check.txt"))
sessionInfo()
sink()

# -----------------------------
# 12. Console output
# -----------------------------

cat("\nCLUE folder check completed successfully.\n\n")

cat("Folder summary:\n")
print(folder_summary)

cat("\nRequired file check:\n")
print(required_check)

cat("\nGCT header check:\n")
print(gct_header_check)

cat("\nLikely score/result files:\n")
print(likely_score_files)

cat("\nCandidate table check:\n")
print(candidate_check)

cat("\nSend me these files from:\n")
cat(table_dir, "\n\n")

cat("01_CLUE_folder_file_inventory.csv\n")
cat("04_required_CLUE_files_check.csv\n")
cat("05_GCT_header_check.csv\n")
cat("07_likely_score_or_result_files.csv\n\n")

cat("Excel workbook saved:\n")
cat(xlsx_file, "\n")
############################################################
# 10A_check_CLUE_folder_contents_FIXED.R
# Robustly check all files/folders inside CLUE output folder
############################################################

packages <- c("dplyr", "readr", "stringr", "tidyr", "openxlsx")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 1. Paths
# -----------------------------

project_dir <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION"
clue_dir    <- "C:/Users/acer/OneDrive/Desktop/Imran_MTB_DRUG-MODULATION/results/CLUE"

out_dir   <- file.path(project_dir, "results", "CLUE_folder_check_FIXED")
table_dir <- file.path(out_dir, "tables")
log_dir   <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(clue_dir)) {
  stop("CLUE folder not found. Please check clue_dir path.")
}

cat("\nChecking CLUE folder:\n")
cat(clue_dir, "\n\n")

# -----------------------------
# 2. Helper functions
# -----------------------------

size_mb <- function(x) {
  round(x / (1024^2), 4)
}

get_ext <- function(x) {
  ext <- tools::file_ext(x)
  ext[ext == ""] <- "no_extension"
  tolower(ext)
}

safe_read_lines <- function(file_path, n = 10) {
  out <- tryCatch(
    readLines(file_path, n = n, warn = FALSE),
    error = function(e) character(0)
  )
  return(out)
}

preview_file_safe <- function(file_path, n = 5) {
  
  lines <- safe_read_lines(file_path, n = n)
  
  if (length(lines) == 0) {
    return(data.frame(
      File_Path = file_path,
      File_Name = basename(file_path),
      Line_Number = 1,
      Preview = "EMPTY_FILE_OR_COULD_NOT_READ",
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    File_Path = rep(file_path, length(lines)),
    File_Name = rep(basename(file_path), length(lines)),
    Line_Number = seq_along(lines),
    Preview = lines,
    stringsAsFactors = FALSE
  )
}

detect_file_type <- function(file_path, is_dir = FALSE) {
  
  if (is_dir) return("Folder")
  
  base <- basename(file_path)
  ext <- get_ext(file_path)
  
  if (base == "ncs.gct") return("MAIN_RESULT_NCS_GCT")
  if (ext == "gct") return("GCT_matrix")
  if (ext == "gmt") return("Gene_set_input_GMT")
  if (ext %in% c("yaml", "yml")) return("Query_config_YAML")
  if (ext %in% c("csv", "tsv", "txt")) return("Text_table_or_log")
  if (ext %in% c("json")) return("JSON_metadata")
  if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "pdf")) return("Figure_or_document")
  if (ext %in% c("zip", "gz")) return("Compressed_file")
  
  return("Other")
}

check_gct_header <- function(file_path) {
  
  lines <- safe_read_lines(file_path, n = 8)
  
  if (length(lines) == 0) {
    return(data.frame(
      File_Path = file_path,
      File_Name = basename(file_path),
      Version = NA,
      Rows = NA,
      Columns = NA,
      Row_Metadata_Columns = NA,
      Column_Metadata_Rows = NA,
      Header_Preview = "EMPTY_OR_COULD_NOT_READ",
      Status = "Could not read",
      stringsAsFactors = FALSE
    ))
  }
  
  version <- trimws(lines[1])
  
  if (version == "#1.3") {
    
    dims <- suppressWarnings(as.integer(strsplit(lines[2], "\t")[[1]]))
    header <- ifelse(length(lines) >= 3, lines[3], NA)
    
    return(data.frame(
      File_Path = file_path,
      File_Name = basename(file_path),
      Version = version,
      Rows = dims[1],
      Columns = dims[2],
      Row_Metadata_Columns = dims[3],
      Column_Metadata_Rows = dims[4],
      Header_Preview = substr(header, 1, 500),
      Status = "Readable GCT 1.3",
      stringsAsFactors = FALSE
    ))
  }
  
  if (version == "#1.2") {
    
    dims <- suppressWarnings(as.integer(strsplit(lines[2], "\t")[[1]]))
    header <- ifelse(length(lines) >= 3, lines[3], NA)
    
    return(data.frame(
      File_Path = file_path,
      File_Name = basename(file_path),
      Version = version,
      Rows = dims[1],
      Columns = dims[2],
      Row_Metadata_Columns = NA,
      Column_Metadata_Rows = NA,
      Header_Preview = substr(header, 1, 500),
      Status = "Readable GCT 1.2",
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    File_Path = file_path,
    File_Name = basename(file_path),
    Version = version,
    Rows = NA,
    Columns = NA,
    Row_Metadata_Columns = NA,
    Column_Metadata_Rows = NA,
    Header_Preview = paste(lines[seq_len(min(length(lines), 3))], collapse = " | "),
    Status = "Not standard GCT",
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 3. Full folder inventory
# -----------------------------

all_items <- list.files(
  clue_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE,
  include.dirs = TRUE
)

if (length(all_items) == 0) {
  stop("No files or folders found inside CLUE folder.")
}

info <- file.info(all_items)

inventory <- data.frame(
  Full_Path = all_items,
  Relative_Path = stringr::str_remove(all_items, fixed(paste0(clue_dir, "/"))),
  Name = basename(all_items),
  Parent_Folder = basename(dirname(all_items)),
  Is_Directory = info$isdir,
  Extension = get_ext(all_items),
  Size_MB = size_mb(info$size),
  Modified_Time = as.character(info$mtime),
  stringsAsFactors = FALSE
) %>%
  mutate(
    File_Type = mapply(detect_file_type, Full_Path, Is_Directory)
  ) %>%
  arrange(desc(Is_Directory), Parent_Folder, Name)

write_csv(
  inventory,
  file.path(table_dir, "01_CLUE_full_file_and_folder_inventory.csv")
)

# -----------------------------
# 4. Folder summary
# -----------------------------

folder_summary <- inventory %>%
  filter(!Is_Directory) %>%
  group_by(Parent_Folder) %>%
  summarise(
    Number_of_Files = n(),
    Total_Size_MB = round(sum(Size_MB, na.rm = TRUE), 4),
    File_Types = paste(sort(unique(File_Type)), collapse = "; "),
    Files = paste(Name, collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(Number_of_Files))

write_csv(
  folder_summary,
  file.path(table_dir, "02_CLUE_folder_summary.csv")
)

# -----------------------------
# 5. Check expected CLUE items
# -----------------------------

expected_items <- data.frame(
  Item = c(
    "Main result score matrix",
    "Query configuration",
    "UP input gene set",
    "DOWN input gene set",
    "ARFS folder",
    "GSEA folder",
    "Matrices folder"
  ),
  Expected_Name = c(
    "ncs.gct",
    "query_config.yaml",
    "up.gmt",
    "down.gmt",
    "arfs",
    "gsea",
    "matrices"
  ),
  Role = c(
    "Essential for next integration",
    "Record of query settings",
    "Input gene-set record only",
    "Input gene-set record only",
    "Supporting output folder",
    "Supporting enrichment folder",
    "May contain extra matrix outputs"
  ),
  stringsAsFactors = FALSE
)

expected_check <- expected_items %>%
  rowwise() %>%
  mutate(
    Found = any(tolower(inventory$Name) == tolower(Expected_Name)),
    Matching_Paths = paste(
      inventory$Full_Path[tolower(inventory$Name) == tolower(Expected_Name)],
      collapse = " || "
    )
  ) %>%
  ungroup()

write_csv(
  expected_check,
  file.path(table_dir, "03_expected_CLUE_items_check.csv")
)

# -----------------------------
# 6. GCT header check
# -----------------------------

gct_files <- inventory %>%
  filter(!Is_Directory, Extension == "gct") %>%
  pull(Full_Path)

if (length(gct_files) > 0) {
  gct_check <- bind_rows(lapply(gct_files, check_gct_header))
} else {
  gct_check <- data.frame(
    File_Path = NA,
    File_Name = NA,
    Version = NA,
    Rows = NA,
    Columns = NA,
    Row_Metadata_Columns = NA,
    Column_Metadata_Rows = NA,
    Header_Preview = NA,
    Status = "No GCT files found",
    stringsAsFactors = FALSE
  )
}

write_csv(
  gct_check,
  file.path(table_dir, "04_GCT_header_check.csv")
)

# -----------------------------
# 7. Preview small readable files
# -----------------------------

preview_candidates <- inventory %>%
  filter(
    !Is_Directory,
    Extension %in% c("gmt", "yaml", "yml", "txt", "csv", "tsv"),
    Size_MB <= 5
  ) %>%
  pull(Full_Path)

if (length(preview_candidates) > 0) {
  preview_table <- bind_rows(lapply(preview_candidates, preview_file_safe, n = 5))
} else {
  preview_table <- data.frame(
    File_Path = NA,
    File_Name = NA,
    Line_Number = NA,
    Preview = "No small readable text files found",
    stringsAsFactors = FALSE
  )
}

write_csv(
  preview_table,
  file.path(table_dir, "05_small_text_file_preview.csv")
)

# -----------------------------
# 8. Check each important folder separately
# -----------------------------

important_folders <- c("arfs", "gsea", "matrices")

folder_details <- list()

for (folder_name in important_folders) {
  
  folder_path <- file.path(clue_dir, folder_name)
  
  if (dir.exists(folder_path)) {
    
    folder_files <- list.files(
      folder_path,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE,
      include.dirs = TRUE
    )
    
    if (length(folder_files) > 0) {
      
      folder_info <- file.info(folder_files)
      
      folder_details[[folder_name]] <- data.frame(
        Checked_Folder = folder_name,
        Full_Path = folder_files,
        Relative_Path = stringr::str_remove(folder_files, fixed(paste0(folder_path, "/"))),
        Name = basename(folder_files),
        Is_Directory = folder_info$isdir,
        Extension = get_ext(folder_files),
        Size_MB = size_mb(folder_info$size),
        Modified_Time = as.character(folder_info$mtime),
        stringsAsFactors = FALSE
      ) %>%
        mutate(File_Type = mapply(detect_file_type, Full_Path, Is_Directory))
      
    } else {
      
      folder_details[[folder_name]] <- data.frame(
        Checked_Folder = folder_name,
        Full_Path = folder_path,
        Relative_Path = "",
        Name = folder_name,
        Is_Directory = TRUE,
        Extension = "folder",
        Size_MB = 0,
        Modified_Time = NA,
        File_Type = "Folder exists but empty",
        stringsAsFactors = FALSE
      )
    }
    
  } else {
    
    folder_details[[folder_name]] <- data.frame(
      Checked_Folder = folder_name,
      Full_Path = folder_path,
      Relative_Path = "",
      Name = folder_name,
      Is_Directory = TRUE,
      Extension = "folder",
      Size_MB = NA,
      Modified_Time = NA,
      File_Type = "Folder not found",
      stringsAsFactors = FALSE
    )
  }
}

important_folder_inventory <- bind_rows(folder_details)

write_csv(
  important_folder_inventory,
  file.path(table_dir, "06_arfs_gsea_matrices_folder_inventory.csv")
)

# -----------------------------
# 9. Likely result / score files
# -----------------------------

likely_result_files <- inventory %>%
  filter(
    !Is_Directory,
    grepl("ncs|score|tau|connect|rank|result|sig|gsea|matrix|mat", Name, ignore.case = TRUE) |
      Parent_Folder %in% c("matrices", "arfs", "gsea")
  ) %>%
  arrange(desc(Size_MB))

write_csv(
  likely_result_files,
  file.path(table_dir, "07_likely_result_or_score_files.csv")
)

# -----------------------------
# 10. Candidate table check
# -----------------------------

candidate_files <- list.files(
  project_dir,
  pattern = "01_all_61_candidates_final_clean_labels\\.csv$|01_all_61_refined_candidate_perturbagens_final\\.csv$|02_candidates_for_validation_final_clean_labels\\.csv$|02_candidates_for_validation_final_clean_labels\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(candidate_files) > 0) {
  candidate_check <- data.frame(
    Candidate_File = candidate_files,
    File_Name = basename(candidate_files),
    Size_MB = size_mb(file.info(candidate_files)$size),
    stringsAsFactors = FALSE
  )
} else {
  candidate_check <- data.frame(
    Candidate_File = NA,
    File_Name = NA,
    Size_MB = NA,
    Note = "No candidate table found",
    stringsAsFactors = FALSE
  )
}

write_csv(
  candidate_check,
  file.path(table_dir, "08_candidate_table_check.csv")
)

# -----------------------------
# 11. Save Excel workbook
# -----------------------------

xlsx_file <- file.path(out_dir, "CLUE_folder_content_check_FIXED.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "Full_inventory")
writeData(wb, "Full_inventory", inventory)

addWorksheet(wb, "Folder_summary")
writeData(wb, "Folder_summary", folder_summary)

addWorksheet(wb, "Expected_items")
writeData(wb, "Expected_items", expected_check)

addWorksheet(wb, "GCT_header_check")
writeData(wb, "GCT_header_check", gct_check)

addWorksheet(wb, "Text_preview")
writeData(wb, "Text_preview", preview_table)

addWorksheet(wb, "ARFS_GSEA_Matrices")
writeData(wb, "ARFS_GSEA_Matrices", important_folder_inventory)

addWorksheet(wb, "Likely_result_files")
writeData(wb, "Likely_result_files", likely_result_files)

addWorksheet(wb, "Candidate_table_check")
writeData(wb, "Candidate_table_check", candidate_check)

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
# 12. Save recommendation
# -----------------------------

recommendation <- data.frame(
  Priority = c(
    "Use first",
    "Check if needed",
    "Record only",
    "Record only",
    "Record only"
  ),
  File_or_Folder = c(
    "ncs.gct",
    "matrices folder",
    "query_config.yaml",
    "up.gmt",
    "down.gmt"
  ),
  Reason = c(
    "Main CLUE normalized connectivity score matrix for candidate matching.",
    "May contain additional matrices, but not always needed if ncs.gct is readable.",
    "Stores query settings for methods/reproducibility.",
    "Input UP gene set only, not result scores.",
    "Input DOWN gene set only, not result scores."
  ),
  stringsAsFactors = FALSE
)

write_csv(
  recommendation,
  file.path(table_dir, "09_recommended_files_for_next_step.csv")
)

sink(file.path(log_dir, "session_info_CLUE_folder_check_FIXED.txt"))
sessionInfo()
sink()

# -----------------------------
# 13. Console output
# -----------------------------

cat("\nCLUE folder check completed successfully.\n\n")

cat("Folder summary:\n")
print(folder_summary)

cat("\nExpected item check:\n")
print(expected_check)

cat("\nGCT header check:\n")
print(gct_check)

cat("\nImportant folder inventory:\n")
print(important_folder_inventory)

cat("\nLikely result/score files:\n")
print(likely_result_files)

cat("\nCandidate table check:\n")
print(candidate_check)

cat("\nOutput folder:\n")
cat(table_dir, "\n\n")

cat("Please send these output files:\n")
cat("01_CLUE_full_file_and_folder_inventory.csv\n")
cat("03_expected_CLUE_items_check.csv\n")
cat("04_GCT_header_check.csv\n")
cat("06_arfs_gsea_matrices_folder_inventory.csv\n")
cat("07_likely_result_or_score_files.csv\n")
cat("09_recommended_files_for_next_step.csv\n\n")

cat("Excel workbook saved:\n")
cat(xlsx_file, "\n")