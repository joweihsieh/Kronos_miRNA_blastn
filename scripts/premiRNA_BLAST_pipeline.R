#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggforce)

# ========== 參數輸入 ==========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: script.R <base_prefix> <genome_db> <genome_fa> <premiRNA_fasta> <query_bed> <query_fai> <query_label> <target_label>")
}

base_prefix <- args[1]
genome_db <- args[2]
genome_fa <- args[3]
premiRNA_fasta <- args[4]
query_bed_path <- args[5]
query_fai_path <- args[6]
query_label <- args[7]
target_label <- args[8]

#out_dir <- dirname(premiRNA_fasta)
out_dir <- args[9]
blast_result <- file.path(out_dir, paste0(base_prefix, "_blast.tab"))
final_bed <- str_replace(blast_result, "\\.tab$", ".bed")
extracted_fa <- str_replace(final_bed, "\\.bed$", ".fa")
upper_fa <- str_replace(final_bed, "\\.bed$", "_Up.fa")

# ========== 1. 執行 BLAST ==========
cmd <- sprintf("blastn -query %s -db %s -outfmt 6 -out %s",
               premiRNA_fasta, genome_db, blast_result)
message(cmd)
system(cmd)

# ========== 2. 讀入與過濾 ==========
df <- read_tsv(blast_result, col_names = FALSE, show_col_types = FALSE)
colnames(df) <- c("query", "subject", "identity", "alignment_length", "mismatches",
                  "gap_opens", "q_start", "q_end", "s_start", "s_end",
                  "evalue", "bit_score")

bed <- read.table(query_bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bed) <- c("chrom", "start", "end", "name", "score", "strand")

# 加入染色體資訊進行比對篩選
df <- df %>%
  left_join(bed[, c("name", "chrom")], by = c("query" = "name")) %>%
  mutate(subject_chr = paste0("Chr", subject),
         Chr_same = ifelse(subject_chr == chrom, "Same", "Different"))

top_hits <- df %>%
  filter(Chr_same == "Same") %>%
  group_by(query) %>%
  arrange(desc(bit_score), desc(identity)) %>%
  slice(1) %>%
  ungroup()

# 產出 BED
write.table(top_hits %>%
              mutate(
                chrom = subject_chr,
                start = pmin(s_start, s_end) - 1,
                end = pmax(s_start, s_end),
                name = query,
                score = ".",
                strand = ifelse(s_start < s_end, "+", "-")
              ) %>%
              select(chrom, start, end, name, score, strand),
            final_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========== 2.5 提取與轉換序列 ==========
system(sprintf("bedtools getfasta -fi %s -bed %s -fo %s -s", genome_fa, final_bed, extracted_fa))
system(sprintf("awk '{if ($0 ~ /^>/) print $0; else print toupper($0)}' %s > %s", extracted_fa, upper_fa))

# ========== 3. 繪圖 ==========
query_bed <- read.table(query_bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(query_bed) <- c("chr", "start", "end", "name", "score", "strand")
query_bed$genome <- query_label
query_bed$midpoint <- (query_bed$start + query_bed$end) / 2

mapped_bed <- read.table(final_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(mapped_bed) <- c("chr", "start", "end", "name", "score", "strand")
mapped_bed$genome <- target_label
mapped_bed$midpoint <- (mapped_bed$start + mapped_bed$end) / 2

merged <- inner_join(query_bed, mapped_bed, by = "name", suffix = c("_query", "_target"))

query_fai <- read.table(query_fai_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
target_fai_path <- paste0(genome_fa, ".fai")
target_fai <- read.table(target_fai_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(query_fai) <- colnames(target_fai) <- c("chr", "length", "offset", "linebases", "linewidth")
query_fai$genome <- query_label
target_fai$genome <- target_label

used_chrs <- unique(c(merged$chr_query, merged$chr_target))
query_fai <- filter(query_fai, chr %in% used_chrs)
target_fai <- filter(target_fai, chr %in% used_chrs)
query_fai <- mutate(query_fai, y = 1.3, yend = 1.5)
target_fai <- mutate(target_fai, y = 0.5, yend = 0.7)

chr_bar_df <- bind_rows(query_fai, target_fai)

points_df <- bind_rows(
  merged %>% transmute(name, chr = chr_target, pos = midpoint_target, genome = target_label, y = 0.6),
  merged %>% transmute(name, chr = chr_query,  pos = midpoint_query,  genome = query_label,  y = 1.4)
)

lines_df <- merged %>%
  transmute(name, chr = chr_target, x = midpoint_target, xend = midpoint_query, y = 0.6, yend = 1.4)

p <- ggplot() +
  geom_rect(data = chr_bar_df,
            aes(xmin = 0, xmax = length, ymin = y, ymax = yend, fill = genome),
            color = "black", alpha = 0.3) +
  geom_point(data = points_df,
             aes(x = pos, y = y, color = genome),
             size = 2.5) +
  geom_segment(data = lines_df,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "gray40", alpha = 0.8) +
  facet_wrap(~chr, scales = "free_x") +
  scale_y_continuous(breaks = c(0.6, 1.4), labels = c(target_label, query_label)) +
  labs(x = "Genomic position", y = NULL, title = "miRNA cross-genome mapping") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

ggsave(file.path(out_dir, paste0(base_prefix, "_liftover_connections.png")), plot = p, width = 20, height = 7, dpi = 300)

identity_info <- top_hits %>% select(query, identity)
lines_df <- merged %>%
  transmute(name, chr = chr_target, x = midpoint_target, xend = midpoint_query, y = 0.6, yend = 1.4) %>%
  left_join(identity_info, by = c("name" = "query")) %>%
  mutate(identity_group = ifelse(identity == 100, "100%", "<100%"))

p2 <- ggplot() +
  geom_rect(data = chr_bar_df,
            aes(xmin = 0, xmax = length, ymin = y, ymax = yend, fill = genome),
            color = "black", alpha = 0.3) +
  geom_point(data = points_df,
             aes(x = pos, y = y),
             size = 2.5, color = "black") +
  geom_segment(data = lines_df,
               aes(x = x, xend = xend, y = y, yend = yend, color = identity_group),
               alpha = 0.9, linewidth = 0.7) +
  scale_color_manual(values = c("100%" = "black", "<100%" = "#56B4E9")) +
  facet_wrap(~chr, scales = "free_x") +
  scale_y_continuous(breaks = c(0.6, 1.4), labels = c(target_label, query_label)) +
  labs(x = "Genomic position", y = NULL, title = "miRNA cross-genome mapping", color = "Identity") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

ggsave(file.path(out_dir, paste0(base_prefix, "_liftover_connections_color.png")), plot = p2, width = 20, height = 7, dpi = 300)


missing_queries <- setdiff(unique(df$query), unique(top_hits$query))
if (length(missing_queries) > 0) {
  message("Some queries have no match on the same chromosome:")
  print(missing_queries)
} else {
  message("✅ All queries have at least one same-chromosome hit.")
}
