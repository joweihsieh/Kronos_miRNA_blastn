#!/usr/bin/env Rscript

library(readxl)
library(Biostrings)
library(dplyr)
library(readr)
library(stringr)

# ========== 參數輸入 ==========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: script.R <xlsx_path> <sheet_name> <base_prefix> <genome_db> <genome_fa> <premiRNA_bed>")
}

xlsx_path <- args[1]
sheet_name <- args[2]
base_prefix <- args[3]
genome_db <- args[4]
genome_fa <- args[5]
premiRNA_bed <- args[6]


#xlsx_path <- "/quobyte/bcmeyersgrp/jwahsieh/Triticum_turgidum_Kronos_liftover/blast/test/Supplementary_DataSet_S1.xlsx"
#sheet_name <- "T.turgidum"
#base_prefix <- "Ttu_Svevov1_miRNA_SB_to_Kronos"
#genome_db <- "/quobyte/bcmeyersgrp/jwahsieh/genome/Triticum_turgidum_Kronos/final/Ttu_Kronos_genome_db"
#genome_fa <- "/quobyte/bcmeyersgrp/jwahsieh/genome/Triticum_turgidum_Kronos/final/Ttu_Kronos_genome_clean.fa"
#premiRNA_bed <- "/quobyte/bcmeyersgrp/jwahsieh/Triticum_turgidum_Kronos_liftover/blast/Ttu_Svevov1_premiRNA_SB_to_Kronos.bed"


out_dir <- dirname(xlsx_path)

output_fasta <- file.path(out_dir, paste0(base_prefix, ".fa"))
blast_result <- file.path(out_dir, paste0(base_prefix, "_blast.tab"))
filtered_tsv <- str_replace(blast_result, "\\.tab$", "_filtered.tab")
filtered_bed <- str_replace(blast_result, "\\.tab$", "_filtered.bed")
filtered_clean <- str_replace(blast_result, "\\.tab$", "_filtered_clean.bed")
final_bed <- str_replace(blast_result, "\\.tab$", ".bed")
precursor_bed <- premiRNA_bed
output_seq <- str_replace(final_bed, "\\.bed$", ".fa")
output_seq_up <- str_replace(final_bed, "\\.bed$", "_Up.fa")
merged_txt <- str_replace(final_bed, "\\.bed$", "_merged.txt")
summary_txt <- str_replace(final_bed, "\\.bed$", "_summary.tsv")
details_txt <- str_replace(final_bed, "\\.bed$", "_details.tsv")

# ========== 1. 建立 FASTA ==========
df <- read_excel(xlsx_path, sheet = sheet_name, skip = 1)
df2 <- df[grep("Ttu", df$Chromosome), ]
fasta <- DNAStringSet(df2$miRNA)
names(fasta) <- df2$Name
writeXStringSet(fasta, filepath = output_fasta)

# ========== 2. 執行 BLAST ==========
cmd <- sprintf("blastn -query %s -strand both -task blastn-short -db %s -out %s -word_size 15 -perc_identity 100 -no_greedy -ungapped -dust no -outfmt 6",
               output_fasta, genome_db, blast_result)
message(cmd)
system(cmd)

# ========== 3. 過濾比對 ==========
query_seqs <- readDNAStringSet(output_fasta)
df <- read_tsv(blast_result, col_names = FALSE, show_col_types = FALSE)
colnames(df)[1:12] <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send","evalue","bitscore")
qlen_map <- tibble(name = names(query_seqs), Q1_seq = as.character(query_seqs), Q1_len = width(query_seqs))
df <- left_join(df, qlen_map, by = c("qseqid" = "name"))
df$total_mismatch <- df$Q1_len - df$length
df_filtered <- df %>% filter(qstart <= 4, qend >= 18, total_mismatch <= 4)
write_tsv(df_filtered, filtered_tsv)

bed_df <- df_filtered %>%
  mutate(
    chrom = paste0("Chr", sseqid),
    start = pmin(sstart, send) - 1,
    end = pmax(sstart, send),
    strand = if_else(sstart < send, "+", "-"),
    name = qseqid
  ) %>% select(chrom, start, end, name, total_mismatch, strand)
write_tsv(bed_df, filtered_bed, col_names = FALSE)

# ========== 4. 交集 precursor ==========
system(sprintf("bedtools intersect -a %s -b %s -wa -wb > %s", filtered_bed, precursor_bed, filtered_clean))
system(sprintf("awk '$4 == $10{print $1, $2, $3, $4, $5, $6}' OFS=\"\t\" %s > %s", filtered_clean, final_bed))
file.remove(filtered_bed, filtered_clean, filtered_tsv, blast_result)

# ========== 5. 取 genome 序列並轉為大寫 ==========
system(sprintf("bedtools getfasta -fi %s -bed %s -fo %s -s", genome_fa, final_bed, output_seq))
system(sprintf("awk '{if ($0 ~ /^>/) print $0; else print toupper($0)}' %s > %s", output_seq, output_seq_up))

# ========== 6. 合併 BED 與序列 ==========
system(sprintf("awk 'FNR==NR { if ($0 ~ /^>/) { gsub(/^>/, \"\", $0); key = $0; getline; seqs[key] = $0; next } } { key = $1\":\"$2\"-\"$3\"(\"$6\")\"; seq = (key in seqs) ? seqs[key] : \"NA\"; print $1, $2+1, $3, $4, $5, $6, seq }' %s %s OFS=\"\t\" > %s", output_seq_up, final_bed, merged_txt))
