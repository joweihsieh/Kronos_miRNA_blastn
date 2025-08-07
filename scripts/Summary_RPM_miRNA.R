#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(Biostrings)

# ========== Parameter ==========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: script.R <gff3_path> <miRNA_Seq> <gff_revised> <Result_fa> <MIR_fa> <counts_path> <total_mapped_path> <colmap_path>")
}

gff3_path <- args[1]
miRNA_Seq <- args[2]
gff_revised <- args[3]
Result_fa <- args[4] #"Results.fa"
MIR_fa <- args[5]
counts_path <- args[6]
total_mapped_path <- args[7]
colmap_path <- args[8]


# ==========  ==========
gff3 <- read.csv(gff3_path, sep = "\t", header = FALSE)
colnames(gff3)[1:8] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase")
colnames(gff3)[9] <- "attributes"

siRNA_df <- gff3[grepl("_locus", gff3$type), ]
siRNA_df$miRNA_id <- sub(".*ID=([^;]+);.*", "\\1", siRNA_df$attributes)

fa_miRNA <- read.csv(miRNA_Seq, sep = " ", header = FALSE)
colnames(fa_miRNA) <- c("chr", "start", "end", "miRNA_id", "score", "strand", "sequence")
selected <- fa_miRNA[fa_miRNA$miRNA_id %in% siRNA_df$miRNA_id, ]

gff_df <- data.frame(
  seqid = selected$chr,
  source = "ShortStack",
  type = "mature_miRNA",
  start = selected$start,
  end = selected$end,
  score = ".",
  strand = selected$strand,
  phase = ".",
  attributes = paste0("ID=", selected$miRNA_id),
  stringsAsFactors = FALSE
)

combined_gff <- rbind(gff3, gff_df)
combined_gff <- combined_gff[order(combined_gff$seqid, combined_gff$start), ]

lowconf_ids <- sub("ID=", "", combined_gff$attributes[grepl("_locus", combined_gff$type)])
combined_gff$ID <- sub(".*ID=([^;]+).*", "\\1", combined_gff$attributes)

combined_gff$Parent <- ifelse(
  combined_gff$type %in% c("mature_miRNA", "miRNA-star"),
  sub("\\.(mature|star)$", "", combined_gff$ID),
  NA
)

combined_gff$attributes <- ifelse(
  !is.na(combined_gff$Parent),
  paste0(combined_gff$attributes, ";Parent=", combined_gff$Parent),
  combined_gff$attributes
)

combined_gff$attributes <- ifelse(
  combined_gff$type %in% c("siRNA21_locus", "mature_miRNA") &
    combined_gff$ID %in% lowconf_ids,
  paste0(combined_gff$attributes, ";low_confident"),
  combined_gff$attributes
)
combined_gff$attributes <- ifelse(
  combined_gff$type %in% c("siRNA22_locus", "mature_miRNA") &
    combined_gff$ID %in% lowconf_ids,
  paste0(combined_gff$attributes, ";low_confident"),
  combined_gff$attributes
)

combined_gff$type[combined_gff$type == "siRNA21_locus"] <- "MIRNA_hairpin"
combined_gff$type[combined_gff$type == "siRNA22_locus"] <- "MIRNA_hairpin"

combined_gff <- combined_gff[, !colnames(combined_gff) %in% c("ID", "Parent")]

combined_gff <- combined_gff %>%
  mutate(
    attributes = ifelse(
      type == "mature_miRNA",
      sub("ID=([^;]+)(;.*)?", "ID=\\1.mature\\2", attributes),
      attributes
    )
  )

writeLines("##gff-version 3", gff_revised)
write.table(combined_gff,
            file = gff_revised,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)


# ==========  ==========

fa <- readDNAStringSet(Result_fa)
names(fa) <- sub("::.*", "", names(fa))  # remove info after :: 

fa_cleaned <- data.frame(
  attributes = names(fa),
  sequence = as.character(fa),
  stringsAsFactors = FALSE,
  miRNA_id = sub(".*ID=([^;]+);.*", "\\1", names(fa))
)


fa_cleaned_siRNA <- fa_cleaned[fa_cleaned$miRNA_id %in% selected$miRNA_id, c("miRNA_id", "sequence")]

hairpin_df <- fa_cleaned_siRNA %>%
  select(miRNA_id, hairpin_seq = sequence)

mature_df <- fa_miRNA %>% select(miRNA_id, mature_seq = sequence)
mature_df_siRNA <- mature_df[mature_df$miRNA_id %in% selected$miRNA_id, ]

seq_tmp <- mature_df_siRNA %>%
  left_join(hairpin_df, by = "miRNA_id")

seq_wide_siRNA <- dplyr::rename(seq_tmp,
  name = miRNA_id,
  precursor_seq = hairpin_seq
)
# ==========

fa <- readDNAStringSet(MIR_fa)
names_df <- tibble(
  raw_name = names(fa),
  name = sub("\\..*", "", sub("::.*", "", names(fa))),
  type = case_when(
    grepl("\\.mature", names(fa)) ~ "mature",
    grepl("\\.star", names(fa)) ~ "star",
    TRUE ~ "precursor"
  ),
  seq = as.character(fa)
)

seq_wide <- names_df %>%
  select(name, type, seq) %>%
  pivot_wider(names_from = type, values_from = seq, values_fill = "") %>%
  dplyr::rename(mature_seq = mature, star_seq = star, precursor_seq = precursor)

full_seq_wide <- full_join(
  seq_wide_siRNA,
  seq_wide,
  by = "name",
  suffix = c("_siRNA", "_main")
) %>%
  mutate(
    mature_seq = coalesce(mature_seq_main, mature_seq_siRNA),
    precursor_seq = coalesce(precursor_seq_main, precursor_seq_siRNA),
    star_seq = star_seq
  ) %>%
  select(name, precursor_seq, mature_seq, star_seq)

counts_df <- read.table(counts_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts_annotated <- counts_df %>% left_join(full_seq_wide, by = c("Name" = "name"))

# Normalize
TotalMapped <- read.table(total_mapped_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TotalMapped$Clean_Name <- sub("_genome", "_condensed", TotalMapped$Sample_ID)
rpm_factors <- setNames(TotalMapped$genome_count / 1e6, TotalMapped$Clean_Name)

cleaned_names <- sub("^X", "", colnames(counts_annotated))
colnames(counts_annotated) <- cleaned_names
sample_cols <- intersect(colnames(counts_annotated), names(rpm_factors))

counts_rpm <- counts_annotated
for (col in sample_cols) {
  counts_rpm[[col]] <- counts_annotated[[col]] / rpm_factors[[col]]
}

colmap <- read.table(colmap_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
name_map <- setNames(colmap$Original_SB_ID, paste0(colmap$Renamed_File_Prefix, "_condensed"))

new_names <- ifelse(cleaned_names %in% names(name_map), name_map[cleaned_names], cleaned_names)
colnames(counts_annotated) <- new_names
colnames(counts_rpm) <- new_names

write.csv(counts_annotated, file = "Counts_with_sequences.csv", row.names = FALSE, quote = TRUE)
write.csv(counts_rpm, file = "RPM_with_sequences.csv", row.names = FALSE, quote = TRUE)
