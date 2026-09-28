# Alpha and beta diversity for the oral virome count table.
#
# Every number printed here is written to a file. Nothing is reported that the
# pipeline did not compute.
#
# Design note on small n: there are 8 samples per group. The smallest p value a
# two-sided exact Wilcoxon rank-sum test can return at 8 vs 8 is about 0.0002,
# and the test has low power against anything but a large shift. Effect sizes
# are reported alongside p values for that reason, and the report says plainly
# what this design can and cannot support.

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file    <- args[1]
presence_file  <- args[2]
samples_file   <- args[3]
group_col      <- args[4]
out_alpha      <- args[5]
out_alpha_tests<- args[6]
out_permanova  <- args[7]
out_pcoa       <- args[8]
out_figure     <- args[9]
out_session    <- args[10]

set.seed(42)

counts_raw   <- read_tsv(counts_file, show_col_types = FALSE)
presence_raw <- read_tsv(presence_file, show_col_types = FALSE)
meta         <- read_tsv(samples_file, show_col_types = FALSE)

stopifnot(identical(counts_raw$accession, presence_raw$accession))

sample_ids <- meta$sample_id
stopifnot(all(sample_ids %in% colnames(counts_raw)))

# Mask counts by the pre-specified presence filter, then aggregate accessions
# that share a virus label so segmented genomes are not counted repeatedly.
cnt <- as.matrix(counts_raw[, sample_ids])
prs <- as.matrix(presence_raw[, sample_ids])
masked <- cnt * prs

agg <- rowsum(masked, group = counts_raw$virus)
agg <- agg[rowSums(agg) > 0, , drop = FALSE]

# vegan wants samples as rows.
mat <- t(agg)

cat("virus labels retained:", nrow(agg), "\n")
cat("samples:", nrow(mat), "\n")
cat("total masked reads:", sum(mat), "\n")

# ---------------------------------------------------------------- alpha
alpha <- tibble(
  sample_id = rownames(mat),
  observed_richness = as.integer(specnumber(mat)),
  shannon = diversity(mat, index = "shannon"),
  total_reads = as.integer(rowSums(mat))
) |>
  left_join(meta |> select(sample_id, all_of(group_col), country), by = "sample_id")

write_tsv(alpha, out_alpha)

group_vec <- factor(alpha[[group_col]])
if (nlevels(group_vec) != 2) {
  stop("expected exactly two groups, found: ", paste(levels(group_vec), collapse = ", "))
}
lv <- levels(group_vec)

# Wilcoxon rank-sum with a rank-biserial effect size.
# r_rb = 2*U/(n1*n2) - 1, reported with the direction of the second level.
wilcox_one <- function(metric) {
  x <- alpha[[metric]][group_vec == lv[1]]
  y <- alpha[[metric]][group_vec == lv[2]]
  tst <- suppressWarnings(wilcox.test(y, x, exact = TRUE, conf.int = FALSE))
  u <- unname(tst$statistic)
  r_rb <- 2 * u / (length(x) * length(y)) - 1
  tibble(
    metric = metric,
    group_1 = lv[1], n_1 = length(x),
    median_1 = median(x), iqr_1 = IQR(x),
    group_2 = lv[2], n_2 = length(y),
    median_2 = median(y), iqr_2 = IQR(y),
    test = "Wilcoxon rank-sum (exact, two-sided)",
    W = u,
    p_value = tst$p.value,
    rank_biserial_r = r_rb
  )
}

alpha_tests <- bind_rows(lapply(c("observed_richness", "shannon"), wilcox_one))
write_tsv(alpha_tests, out_alpha_tests)
print(as.data.frame(alpha_tests))

# ---------------------------------------------------------------- beta
# Bray-Curtis on relative abundance. Samples with no retained reads cannot be
# placed in the ordination and are dropped here; the count is reported.
rel <- sweep(mat, 1, pmax(rowSums(mat), 1), "/")
keep <- rowSums(mat) > 0
n_dropped <- sum(!keep)
cat("samples dropped from beta diversity (no retained viral reads):", n_dropped, "\n")

rel <- rel[keep, , drop = FALSE]
meta_beta <- alpha[match(rownames(rel), alpha$sample_id), ]

bray <- vegdist(rel, method = "bray")

perm_rows <- list()
if (nrow(rel) >= 4) {
  # Diagnosis alone, then diagnosis adjusted for country. The source paper
  # reported country explaining more variance than diagnosis, so both are run.
  f1 <- adonis2(bray ~ meta_beta[[group_col]], permutations = 9999)
  perm_rows[[1]] <- tibble(
    model = "bray ~ diagnosis",
    term = "diagnosis",
    df = f1$Df[1], sum_sq = f1$SumOfSqs[1],
    R2 = f1$R2[1], F = f1$F[1], p_value = f1$`Pr(>F)`[1],
    permutations = 9999
  )

  if (length(unique(meta_beta$country)) > 1) {
    f2 <- adonis2(
      bray ~ country + meta_beta[[group_col]],
      data = meta_beta, permutations = 9999, by = "terms"
    )
    perm_rows[[2]] <- tibble(
      model = "bray ~ country + diagnosis", term = "country",
      df = f2$Df[1], sum_sq = f2$SumOfSqs[1],
      R2 = f2$R2[1], F = f2$F[1], p_value = f2$`Pr(>F)`[1], permutations = 9999
    )
    perm_rows[[3]] <- tibble(
      model = "bray ~ country + diagnosis", term = "diagnosis",
      df = f2$Df[2], sum_sq = f2$SumOfSqs[2],
      R2 = f2$R2[2], F = f2$F[2], p_value = f2$`Pr(>F)`[2], permutations = 9999
    )
  }

  # Dispersion check. PERMANOVA can be driven by unequal spread rather than
  # by a location shift, so this is reported next to it.
  bd <- betadisper(bray, factor(meta_beta[[group_col]]))
  bd_test <- permutest(bd, permutations = 9999)
  perm_rows[[length(perm_rows) + 1]] <- tibble(
    model = "betadisper (homogeneity of dispersion)",
    term = "diagnosis",
    df = bd_test$tab$Df[1], sum_sq = bd_test$tab$`Sum Sq`[1],
    R2 = NA_real_, F = bd_test$tab$F[1],
    p_value = bd_test$tab$`Pr(>F)`[1], permutations = 9999
  )
}

permanova <- bind_rows(perm_rows)
if (nrow(permanova) == 0) {
  permanova <- tibble(
    model = NA_character_, term = NA_character_, df = NA_integer_,
    sum_sq = NA_real_, R2 = NA_real_, F = NA_real_,
    p_value = NA_real_, permutations = NA_integer_
  )
}
write_tsv(permanova, out_permanova)
print(as.data.frame(permanova))

# PCoA
pcoa_tbl <- tibble(sample_id = character(), axis1 = numeric(), axis2 = numeric())
var_exp <- c(NA_real_, NA_real_)
if (nrow(rel) >= 3) {
  pc <- cmdscale(bray, k = 2, eig = TRUE)
  pos <- pc$eig[pc$eig > 0]
  var_exp <- 100 * pc$eig[1:2] / sum(pos)
  pcoa_tbl <- tibble(
    sample_id = rownames(pc$points),
    axis1 = pc$points[, 1],
    axis2 = pc$points[, 2]
  ) |>
    left_join(meta |> select(sample_id, all_of(group_col), country), by = "sample_id") |>
    mutate(axis1_pct_var = var_exp[1], axis2_pct_var = var_exp[2])
}
write_tsv(pcoa_tbl, out_pcoa)

# ---------------------------------------------------------------- figure
pal <- c("#4C72B0", "#C44E52")
names(pal) <- lv

long_alpha <- alpha |>
  select(sample_id, all_of(group_col), observed_richness, shannon) |>
  pivot_longer(c(observed_richness, shannon), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    observed_richness = "Observed richness",
    shannon = "Shannon index"
  ))

p_alpha <- ggplot(long_alpha, aes(.data[[group_col]], value, fill = .data[[group_col]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.6) +
  geom_jitter(width = 0.14, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = NULL, title = "Alpha diversity by periodontal status") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 12, hjust = 1))

if (nrow(pcoa_tbl) > 0) {
  p_beta <- ggplot(pcoa_tbl, aes(axis1, axis2,
                                 colour = .data[[group_col]], shape = country)) +
    geom_point(size = 2.4, stroke = 0.7) +
    scale_colour_manual(values = pal) +
    labs(
      title = "PCoA, Bray-Curtis",
      x = sprintf("Axis 1 (%.1f%%)", var_exp[1]),
      y = sprintf("Axis 2 (%.1f%%)", var_exp[2]),
      colour = NULL, shape = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "right", legend.key.size = unit(0.35, "cm"))
} else {
  p_beta <- ggplot() + theme_void() +
    annotate("text", 0, 0, size = 3,
             label = "Too few samples with viral reads for an ordination")
}

ggsave(out_figure, p_alpha / p_beta, width = 7.2, height = 6.2, dpi = 200)

writeLines(capture.output(sessionInfo()), out_session)
cat("done\n")
