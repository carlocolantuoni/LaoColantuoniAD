suppressPackageStartupMessages({ library(dplyr); library(readr) })

nx_fisher_z <- function(r1, r2, n1, n2) {
  if (is.na(r1) || is.na(r2) || is.na(n1) || is.na(n2) || n1 < 4 || n2 < 4)
    return(list(z = NA, p = NA, delta_r = NA))
  r1c <- pmax(pmin(r1, 0.999), -0.999); r2c <- pmax(pmin(r2, 0.999), -0.999)
  z1 <- 0.5 * log((1 + r1c) / (1 - r1c)); z2 <- 0.5 * log((1 + r2c) / (1 - r2c))
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3)); zd <- (z1 - z2) / se
  list(z = zd, p = 2 * (1 - pnorm(abs(zd))), delta_r = r1 - r2)
}

nx_wilcox_p <- function(v1, v2) {
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) < 3 || length(v2) < 3) return(NA_real_)
  tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA_real_)
}

nx_pairwise_wilcox <- function(df, value_col, group_col, comparisons) {
  rows <- lapply(comparisons, function(pr) {
    v1 <- df[[value_col]][as.character(df[[group_col]]) == pr[1]]
    v2 <- df[[value_col]][as.character(df[[group_col]]) == pr[2]]
    v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
    p <- if (length(v1) >= 1 && length(v2) >= 1)
      tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA_real_) else NA_real_
    data.frame(group1 = pr[1], group2 = pr[2], n1 = length(v1), n2 = length(v2),
               p_wilcox = p, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_BH <- p.adjust(out$p_wilcox, method = "BH")
  out
}

INPUT_FOLDER  <- "exported_data_fixed"
OUTPUT_FOLDER <- sub("^exported_data", "analysis_output", INPUT_FOLDER)
exp_dir       <- file.path(getwd(), INPUT_FOLDER)
analysis_dir  <- file.path(getwd(), OUTPUT_FOLDER)
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
ref_dir <- file.path(getwd(), "Output4", "Doublehi")
cat(sprintf("INPUT  folder : %s\nOUTPUT folder : %s\n", exp_dir, analysis_dir))

.nx_log <- character(0)
nx_record <- function(ref_name, status) .nx_log[[length(.nx_log) + 1]] <<-
  sprintf("  %-48s %s", ref_name, status)

nx_check <- function(new_df, ref_name, ref_dir = get("ref_dir", envir = .GlobalEnv), tol = 1e-8, key = NULL) {
  f <- file.path(ref_dir, ref_name)
  if (!file.exists(f)) { cat(sprintf("[check] %-46s no existing file\n", ref_name))
    nx_record(ref_name, "no existing file"); return(invisible(NULL)) }
  old <- as.data.frame(readr::read_csv(f, show_col_types = FALSE))
  names(old) <- gsub("HiPath", "CN-Hi", gsub("LoPath", "CN-Lo", names(old)))

  if (is.null(key)) {
    cand <- intersect(c("Gene", "gene", "SampleID", "Pathway", "Feature", "ID"),
                      intersect(names(old), names(new_df)))
    key <- if (length(cand)) cand[1] else NULL
  }
  if (!is.null(key) && key %in% names(old) && key %in% names(new_df))
    old <- old[match(new_df[[key]], old[[key]]), , drop = FALSE]
  common <- intersect(names(old), names(new_df))
  num    <- common[vapply(new_df[common], is.numeric, logical(1))]
  diffs  <- character(0)
  for (cc in num) {
    a <- old[[cc]]; b <- new_df[[cc]]
    if (length(a) != length(b) || !isTRUE(all.equal(a, b, tolerance = tol))) diffs <- c(diffs, cc)
  }
  status <- if (!length(diffs)) sprintf("IDENTICAL (%d cols%s)", length(num),
                                        if (!is.null(key)) paste0(", by ", key) else "")
  else sprintf("DIFFERS: %s", paste(utils::head(diffs, 8), collapse = ", "))
  cat(sprintf("[check] %-46s %s\n", ref_name, status))
  nx_record(ref_name, status)
  invisible(list(same = !length(diffs), diffs = diffs))
}

meta    <- as.data.frame(readr::read_csv(file.path(exp_dir, "metadata_master.csv"),   show_col_types = FALSE))
prot    <- as.data.frame(readr::read_csv(file.path(exp_dir, "protein_abundance.csv"), show_col_types = FALSE))
expr_df <- as.data.frame(readr::read_csv(file.path(exp_dir, "expression_rna.csv.gz"), show_col_types = FALSE))

expr <- as.matrix(expr_df[, -1, drop = FALSE]); rownames(expr) <- expr_df[[1]]

meta$Group_Current <- factor(dplyr::recode(as.character(meta$Group_Current),
                                           "LoPath" = "CN-Lo", "HiPath" = "CN-Hi"),
                             levels = c("CN-Lo", "CN-Hi", "MCI", "AD", "Others"))

shared <- Reduce(intersect, list(as.character(meta$SampleID),
                                 as.character(prot$SampleID), colnames(expr)))
meta <- meta[match(shared, meta$SampleID), , drop = FALSE]
prot <- prot[match(shared, prot$SampleID), , drop = FALSE]
expr <- expr[, shared, drop = FALSE]
stopifnot(identical(as.character(meta$SampleID), colnames(expr)),
          identical(as.character(prot$SampleID), colnames(expr)))
cat(sprintf("Data ready: meta %d x %d | prot %d x %d | expr %d x %d (n=%d shared)\n",
            nrow(meta), ncol(meta), nrow(prot), ncol(prot), nrow(expr), ncol(expr), length(shared)))

sg_overall <- meta$SampleID
sg_control <- meta$SampleID[trimws(as.character(meta$Diagnosis_COMP)) == "Control"]
sg_cnlo    <- meta$SampleID[meta$Group_Current == "CN-Lo"]
sg_cnhi    <- meta$SampleID[meta$Group_Current == "CN-Hi"]
sg_mci     <- meta$SampleID[meta$Group_Current == "MCI"]
sg_ad      <- meta$SampleID[meta$Group_Current == "AD"]
sg_youngcon<- meta$SampleID[trimws(as.character(meta$Diagnosis_COMP)) == "YoungCon"]

prm_samps <- prot$SampleID[is.finite(as.numeric(prot$MS_NPTX2))]

ms_cols       <- grep("^MS_", colnames(prot), value = TRUE)
raw_genes     <- gsub("^MS_", "", ms_cols)
prot_genes    <- ifelse(raw_genes == "KIAA1045", "PHF24", raw_genes)
prot_data_col <- setNames(ms_cols, prot_genes)

meta$In_PRM_MS      <- as.logical(meta$In_PRM_MS)
meta$NPTX2_RNA      <- as.numeric(expr["NPTX2", meta$SampleID])
pm                  <- match(meta$SampleID, prot$SampleID)
meta$NPTX2_MS       <- as.numeric(prot$MS_NPTX2[pm])
meta$Age_Num        <- suppressWarnings(as.numeric(meta$Age))
meta$Diag_Clean     <- trimws(as.character(meta$Diagnosis_COMP))

main_df <- meta %>%
  dplyr::filter(Diag_Clean %in% c("Control", "MCI", "AD")) %>%
  dplyr::mutate(Panel = "main", Group = factor(Diag_Clean, levels = c("Control", "MCI", "AD")))

path_df <- meta %>%
  dplyr::mutate(Group = dplyr::case_when(
    Diag_Clean == "YoungCon" ~ "YoungCon",
    Group_Current == "CN-Lo" ~ "CN-Lo",
    Group_Current == "CN-Hi" ~ "CN-Hi",
    TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(Group)) %>%
  dplyr::mutate(Panel = "path", Group = factor(Group, levels = c("CN-Lo", "CN-Hi", "YoungCon")))

val_cols <- c("SampleID", "Panel", "Group", "Age_Num", "In_PRM_MS",
              "NPTX2_RNA", "NPTX2_MS")
fig1b_values <- dplyr::bind_rows(main_df[, val_cols], path_df[, val_cols])
fig1b_values$Group <- as.character(fig1b_values$Group)
readr::write_csv(fig1b_values, file.path(analysis_dir, "Fig1b_boxplot_values.csv"))

comps_main <- list(c("Control","MCI"), c("Control","AD"), c("MCI","AD"))
comps_path <- list(c("CN-Lo","CN-Hi"), c("CN-Lo","YoungCon"), c("CN-Hi","YoungCon"))

mods  <- c("NPTX2_RNA", "NPTX2_MS")
parts <- list()
for (m in mods) {
  parts[[paste0(m,"_main")]] <- cbind(Modality = m, Panel = "main",
                                      nx_pairwise_wilcox(main_df, m, "Group", comps_main))
  parts[[paste0(m,"_path")]] <- cbind(Modality = m, Panel = "path",
                                      nx_pairwise_wilcox(path_df, m, "Group", comps_path))
}
fig1b_stats <- do.call(rbind, parts); rownames(fig1b_stats) <- NULL
readr::write_csv(fig1b_stats, file.path(analysis_dir, "Fig1b_comparison_stats.csv"))
nx_check(fig1b_stats, "Boxplot_22_comparison_stats.csv")

cat("\nFigure 1b/3b analysis done -> Fig1b_boxplot_values.csv + Fig1b_comparison_stats.csv\n")

subgroups_cor <- list("Overall" = sg_overall, "Control" = sg_control,
                      "CN-Lo" = sg_cnlo, "CN-Hi" = sg_cnhi, "MCI" = sg_mci, "AD" = sg_ad)
subgroups_exp <- list("CN-Lo" = sg_cnlo, "CN-Hi" = sg_cnhi, "MCI" = sg_mci, "AD" = sg_ad)

corr_A    <- c("DPP6","GRIN2B","HOMER1","PHF24","NPTXR","RIMS1","VGF")
corr_B    <- c("ATP6V1H","GRIA4","HNRNPA2B1","NPTX1","RAB11A","RAB5A","RHEB","SYT1","VDAC1")
corr_weak <- c("LAMP1","LAMP2","SQSTM1")

nx_traj_protein_metrics <- function(gene, meta, sub_cor, sub_exp, col_map, ref_col = "MS_NPTX2") {
  vec_raw <- as.numeric(meta[[col_map[[gene]]]]); names(vec_raw) <- meta$SampleID
  ref     <- as.numeric(meta[[ref_col]]);          names(ref)     <- meta$SampleID

  fin <- vec_raw[is.finite(vec_raw)]; vmin <- min(fin); vmax <- max(fin)
  vec_sc <- if (is.finite(vmin) && is.finite(vmax) && vmax > vmin) (vec_raw - vmin) / (vmax - vmin)
  else rep(NA_real_, length(vec_raw))
  names(vec_sc) <- names(vec_raw)

  res <- list(Gene = gene); gr <- list(); gn <- list()
  for (sg in names(sub_cor)) {
    s <- sub_cor[[sg]]; v <- vec_raw[s]; r <- ref[s]; ok <- is.finite(v) & is.finite(r); n <- sum(ok)
    if (n >= 3) { rv <- cor(v[ok], r[ok], method = "pearson"); res[[paste0("Cor_", sg)]] <- rv; gr[[sg]] <- rv; gn[[sg]] <- n }
    else        { res[[paste0("Cor_", sg)]] <- NA; gr[[sg]] <- NA; gn[[sg]] <- NA }
  }
  gev <- list()
  for (sg in names(sub_exp)) {
    s <- sub_exp[[sg]]; vs <- vec_sc[s]; vs <- vs[is.finite(vs)]
    res[[paste0("Exp_", sg)]] <- if (length(vs)) mean(vs) else NA
    gev[[sg]] <- vs
  }
  fz <- function(a, b) nx_fisher_z(gr[[a]], gr[[b]], gn[[a]], gn[[b]])
  z_HiLo <- fz("CN-Hi","CN-Lo"); z_HiMci <- fz("CN-Hi","MCI"); z_HiAd <- fz("CN-Hi","AD")
  z_LoMci <- fz("CN-Lo","MCI"); z_MciAd <- fz("MCI","AD")
  res$Cor_FisherZ_Hi_vs_Lo<-z_HiLo$z;  res$Cor_FisherZ_p_Hi_vs_Lo<-z_HiLo$p;  res$Cor_DeltaR_Hi_vs_Lo<-z_HiLo$delta_r
  res$Cor_FisherZ_Hi_vs_MCI<-z_HiMci$z;res$Cor_FisherZ_p_Hi_vs_MCI<-z_HiMci$p;res$Cor_DeltaR_Hi_vs_MCI<-z_HiMci$delta_r
  res$Cor_FisherZ_Hi_vs_AD<-z_HiAd$z;  res$Cor_FisherZ_p_Hi_vs_AD<-z_HiAd$p;  res$Cor_DeltaR_Hi_vs_AD<-z_HiAd$delta_r
  res$Cor_FisherZ_Lo_vs_MCI<-z_LoMci$z;res$Cor_FisherZ_p_Lo_vs_MCI<-z_LoMci$p;res$Cor_DeltaR_Lo_vs_MCI<-z_LoMci$delta_r
  res$Cor_FisherZ_MCI_vs_AD<-z_MciAd$z;res$Cor_FisherZ_p_MCI_vs_AD<-z_MciAd$p;res$Cor_DeltaR_MCI_vs_AD<-z_MciAd$delta_r
  res$Exp_Wilcox_p_Hi_vs_Lo  <- nx_wilcox_p(gev[["CN-Hi"]], gev[["CN-Lo"]])
  res$Exp_Wilcox_p_Hi_vs_MCI <- nx_wilcox_p(gev[["CN-Hi"]], gev[["MCI"]])
  res$Exp_Wilcox_p_Hi_vs_AD  <- nx_wilcox_p(gev[["CN-Hi"]], gev[["AD"]])
  res$Exp_Wilcox_p_Lo_vs_MCI <- nx_wilcox_p(gev[["CN-Lo"]], gev[["MCI"]])
  res$Exp_Wilcox_p_MCI_vs_AD <- nx_wilcox_p(gev[["MCI"]],   gev[["AD"]])
  as.data.frame(res, check.names = FALSE)
}

nx_traj_mrna_expr <- function(gene, expr, prm, sub_exp) {
  out <- setNames(as.list(rep(NA_real_, length(sub_exp))), paste0("mRNA_", names(sub_exp)))
  if (gene %in% rownames(expr)) {
    prm_in <- intersect(prm, colnames(expr))
    v <- as.numeric(expr[gene, prm_in]); names(v) <- prm_in
    fin <- v[is.finite(v)]
    if (length(fin)) {
      vmin <- min(fin); vmax <- max(fin)
      vsc <- if (vmax > vmin) (v - vmin) / (vmax - vmin) else rep(NA_real_, length(v)); names(vsc) <- names(v)
      for (sg in names(sub_exp)) {
        s <- intersect(sub_exp[[sg]], prm_in); vv <- vsc[s]; vv <- vv[is.finite(vv)]
        out[[paste0("mRNA_", sg)]] <- if (length(vv)) mean(vv) else NA
      }
    }
  }
  data.frame(Gene = gene, out, check.names = FALSE, stringsAsFactors = FALSE)
}

df_fig1cd <- dplyr::bind_rows(lapply(prot_genes, nx_traj_protein_metrics,
                                     meta = prot, sub_cor = subgroups_cor,
                                     sub_exp = subgroups_exp, col_map = prot_data_col))

df_fig1cd <- df_fig1cd %>% dplyr::mutate(
  Cor_FisherZ_padj_Hi_vs_Lo  = p.adjust(Cor_FisherZ_p_Hi_vs_Lo,  method = "BH"),
  Cor_FisherZ_padj_Hi_vs_MCI = p.adjust(Cor_FisherZ_p_Hi_vs_MCI, method = "BH"),
  Cor_FisherZ_padj_Hi_vs_AD  = p.adjust(Cor_FisherZ_p_Hi_vs_AD,  method = "BH"),
  Cor_FisherZ_padj_Lo_vs_MCI = p.adjust(Cor_FisherZ_p_Lo_vs_MCI, method = "BH"),
  Cor_FisherZ_padj_MCI_vs_AD = p.adjust(Cor_FisherZ_p_MCI_vs_AD, method = "BH"),
  Exp_Wilcox_padj_Hi_vs_Lo   = p.adjust(Exp_Wilcox_p_Hi_vs_Lo,   method = "BH"),
  Exp_Wilcox_padj_Hi_vs_MCI  = p.adjust(Exp_Wilcox_p_Hi_vs_MCI,  method = "BH"),
  Exp_Wilcox_padj_Hi_vs_AD   = p.adjust(Exp_Wilcox_p_Hi_vs_AD,   method = "BH"),
  Exp_Wilcox_padj_Lo_vs_MCI  = p.adjust(Exp_Wilcox_p_Lo_vs_MCI,  method = "BH"),
  Exp_Wilcox_padj_MCI_vs_AD  = p.adjust(Exp_Wilcox_p_MCI_vs_AD,  method = "BH")
)

df_mrna   <- dplyr::bind_rows(lapply(prot_genes, nx_traj_mrna_expr,
                                     expr = expr, prm = prm_samps, sub_exp = subgroups_exp))
df_fig1cd <- dplyr::left_join(df_fig1cd, df_mrna, by = "Gene")

df_fig1cd <- df_fig1cd %>% dplyr::mutate(
  Corr_Traj = dplyr::case_when(
    Gene == "NPTX2"     ~ "Ref",
    Gene %in% corr_A    ~ "A",
    Gene %in% corr_B    ~ "B",
    Gene %in% corr_weak ~ "Weak Corr",
    TRUE                ~ "Weak Corr"))

readr::write_csv(df_fig1cd, file.path(analysis_dir, "Fig1cd_trajectory_metrics.csv"))
nx_check(df_fig1cd, "Table_17_Trajectory_Comparison_WithStats.csv")
nx_check(df_fig1cd, "Table_17_Trajectory_Comparison.csv")

cat("\nFigure 1c/1d analysis done -> analysis_output/Fig1cd_trajectory_metrics.csv\n")

nx_cor_to_ref <- function(expr_mat, ref_vec, samps, genes, min_n = 10) {
  s    <- intersect(intersect(samps, names(ref_vec)), colnames(expr_mat))
  M    <- t(expr_mat[genes, s, drop = FALSE])
  rv   <- as.numeric(ref_vec[s])
  n_ok <- colSums(is.finite(M) & is.finite(rv))
  r    <- suppressWarnings(as.numeric(stats::cor(M, rv, use = "pairwise.complete.obs")))
  r[n_ok < min_n] <- NA_real_
  data.frame(Gene = genes, r = r, stringsAsFactors = FALSE)
}

nx_classify_r <- function(df) {
  df <- df[is.finite(df$r), , drop = FALSE]
  df$Category <- ifelse(abs(df$r) >= 0.5, "Strong (|r|\u22650.5)",
                        ifelse(abs(df$r) >= 0.3, "Moderate (0.3\u2264|r|<0.5)", "Weak (|r|<0.3)"))
  df[order(-abs(df$r)), c("Gene", "r", "Category")]
}

grp_rna <- list("CN-Lo" = sg_cnlo, "CN-Hi" = sg_cnhi, "MCI" = sg_mci, "AD" = sg_ad,
                "YoungCon" = sg_youngcon)

nptx2_rna <- as.numeric(expr["NPTX2", ]); names(nptx2_rna) <- colnames(expr)
genes_rna <- setdiff(rownames(expr), "NPTX2")
cor_rna_overall <- nx_cor_to_ref(expr, nptx2_rna, colnames(expr), genes_rna)
cor_rna_groups  <- dplyr::bind_rows(lapply(names(grp_rna), function(g) {
  d <- nx_cor_to_ref(expr, nptx2_rna, grp_rna[[g]], genes_rna); d <- d[is.finite(d$r), ]; d$Group <- g; d
}))

nptx2_prot <- as.numeric(prot$MS_NPTX2); names(nptx2_prot) <- prot$SampleID
prm236     <- intersect(prot$SampleID[is.finite(nptx2_prot)], colnames(expr))
genes_prot <- rownames(expr)
cor_prot_overall <- nx_cor_to_ref(expr, nptx2_prot, prm236, genes_prot)
grp_prot         <- lapply(grp_rna, function(s) intersect(s, prm236))
cor_prot_groups  <- dplyr::bind_rows(lapply(names(grp_prot), function(g) {
  d <- nx_cor_to_ref(expr, nptx2_prot, grp_prot[[g]], genes_prot); d <- d[is.finite(d$r), ]; d$Group <- g; d
}))

readr::write_csv(nx_classify_r(cor_rna_overall),  file.path(analysis_dir, "Fig3c_cor_overall.csv"))
readr::write_csv(nx_classify_r(cor_prot_overall), file.path(analysis_dir, "Fig2b_cor_overall.csv"))
readr::write_csv(cor_rna_groups,  file.path(analysis_dir, "Fig3c_cor_groups.csv"))
readr::write_csv(cor_prot_groups, file.path(analysis_dir, "Fig2b_cor_groups.csv"))

nx_check(nx_classify_r(cor_rna_overall),  "Table_23.5_Transcriptome_NPTX2_RNA_Cor.csv")
nx_check(nx_classify_r(cor_prot_overall), "Table_23.6_Transcriptome_NPTX2_Protein_Cor.csv")

nx_density_ks <- function(groups_df, overall_df, modality) {
  grps <- c("CN-Lo", "CN-Hi", "MCI", "AD")
  vecs <- lapply(grps, function(g) { v <- groups_df$r[groups_df$Group == g]; v[is.finite(v)] })
  names(vecs) <- grps
  vecs[["Overall"]] <- overall_df$r[is.finite(overall_df$r)]
  combos <- utils::combn(names(vecs), 2, simplify = FALSE)
  out <- do.call(rbind, lapply(combos, function(p) {
    x <- vecs[[p[1]]]; y <- vecs[[p[2]]]
    if (length(x) < 2 || length(y) < 2) return(NULL)
    k <- suppressWarnings(ks.test(x, y))
    data.frame(Modality = modality, Group1 = p[1], Group2 = p[2],
               n1 = length(x), n2 = length(y), KS_D = unname(k$statistic),
               KS_p = k$p.value, stringsAsFactors = FALSE)
  }))
  out$KS_padj <- p.adjust(out$KS_p, method = "BH")
  out
}
fig2b3c_density_ks <- rbind(
  nx_density_ks(cor_rna_groups,  cor_rna_overall,  "NPTX2_mRNA_ref"),
  nx_density_ks(cor_prot_groups, cor_prot_overall, "NPTX2_protein_ref"))
readr::write_csv(fig2b3c_density_ks, file.path(analysis_dir, "Fig2b3c_density_KS.csv"))
cat("\nFigure 2b/3c density-correlation analysis done (+ by-group KS summary).\n")

MIN_EXPR_THRESH <- 1.0
TRI_T_STRONG    <- 0.40
TRI_T_WEAK      <- 0.25

nx_grp_ids <- function(g) intersect(meta$SampleID[meta$Group_Current == g], colnames(expr))
lo_ids  <- nx_grp_ids("CN-Lo")
hi_ids  <- nx_grp_ids("CN-Hi")
ad_ids  <- nx_grp_ids("AD")
mci_ids <- nx_grp_ids("MCI")
valid_all <- colnames(expr)
cat(sprintf("Master-class sample counts: CN-Lo=%d CN-Hi=%d AD=%d MCI=%d\n",
            length(lo_ids), length(hi_ids), length(ad_ids), length(mci_ids)))
if (length(mci_ids) < 10) stop("MCI group has <10 samples - cannot proceed.")

np_rna  <- setNames(as.numeric(expr["NPTX2", ]), colnames(expr))
np_prot <- setNames(as.numeric(prot$MS_NPTX2), prot$SampleID)

nx_cor_n_all <- function(ids, ref) {
  M  <- t(expr[, ids, drop = FALSE])
  rv <- as.numeric(ref[ids])
  n  <- colSums(is.finite(M) & is.finite(rv))
  r  <- suppressWarnings(as.numeric(stats::cor(M, rv, use = "pairwise.complete.obs")))
  r[n <= 3] <- NA_real_
  list(r = r, n = n)
}

cat("  Computing genome-wide correlations (master gene list)...\n")
union4    <- c(lo_ids, hi_ids, ad_ids, mci_ids)
mean_expr <- rowMeans(expr[, union4, drop = FALSE], na.rm = TRUE)

c_rna_lo <- nx_cor_n_all(lo_ids, np_rna);  c_rna_hi <- nx_cor_n_all(hi_ids, np_rna)
c_rna_ad <- nx_cor_n_all(ad_ids, np_rna)
c_pro_lo <- nx_cor_n_all(lo_ids, np_prot); c_pro_hi <- nx_cor_n_all(hi_ids, np_prot)
c_pro_ad <- nx_cor_n_all(ad_ids, np_prot)

R6 <- cbind(c_rna_lo$r, c_rna_hi$r, c_rna_ad$r, c_pro_lo$r, c_pro_hi$r, c_pro_ad$r)
N6 <- cbind(c_rna_lo$n, c_rna_hi$n, c_rna_ad$n, c_pro_lo$n, c_pro_hi$n, c_pro_ad$n)
valid6  <- !is.na(R6) & N6 > 3
absR    <- abs(R6); absR[!valid6] <- -Inf
best    <- max.col(absR, ties.method = "first")
has_any <- rowSums(valid6) > 0
idx     <- cbind(seq_len(nrow(R6)), best)
Max_r   <- R6[idx]; Max_r[!has_any] <- NA_real_
Max_n   <- N6[idx]; Max_n[!has_any] <- NA_real_

results_full <- tibble::tibble(Gene = rownames(expr), Mean_Expr = mean_expr,
                               Max_r = Max_r, Max_n = Max_n)

nx_pval <- function(r, n) {
  t_stat <- r * sqrt((n - 2) / pmax(1e-15, 1 - r^2))
  p <- 2 * pt(-abs(t_stat), df = pmax(1, n - 2))
  p[is.na(n) | n <= 3 | is.na(r)] <- NA
  p
}
results_qc <- results_full %>%
  dplyr::filter(Mean_Expr >= MIN_EXPR_THRESH) %>%
  dplyr::mutate(p_val_max = nx_pval(Max_r, Max_n),
                padj_max  = p.adjust(p_val_max, method = "BH"))
master_gene_list <- results_qc %>% dplyr::filter(abs(Max_r) > 0.4 & padj_max < 0.01)

pc_file <- file.path(exp_dir, "protein_coding_symbols.csv")
if (!exists("valid_pc_symbols")) {
  if (file.exists(pc_file)) {
    valid_pc_symbols <- as.character(readr::read_csv(pc_file, show_col_types = FALSE)[[1]])
    cat(sprintf("  Protein-coding symbols: loaded %d from %s\n",
                length(valid_pc_symbols), basename(pc_file)))
  } else if (requireNamespace("biomaRt", quietly = TRUE)) {
    cat("  Protein-coding symbols: fetching from Ensembl (biomaRt)...\n")
    mart <- tryCatch(biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                     error = function(e) NULL)
    if (!is.null(mart)) {
      pc_df <- biomaRt::getBM(attributes = c("hgnc_symbol", "gene_biotype"), mart = mart)
      valid_pc_symbols <- unique(pc_df$hgnc_symbol[pc_df$gene_biotype == "protein_coding"])
    } else { cat("  WARNING biomaRt unreachable - skipping PC filter.\n")
      valid_pc_symbols <- master_gene_list$Gene }
  } else { cat("  WARNING biomaRt not installed - skipping PC filter.\n")
    valid_pc_symbols <- master_gene_list$Gene }
}

master_gene_list <- master_gene_list %>%
  dplyr::filter(!grepl("^ENSG", Gene)) %>%
  dplyr::filter(Gene %in% valid_pc_symbols) %>%
  dplyr::filter(!grepl("^OR[0-9]+[A-Z]", Gene)) %>%
  dplyr::arrange(dplyr::desc(abs(Max_r)))
candidate_genes <- master_gene_list$Gene
cat(sprintf("  Master gene list: %d protein-coding genes passed QC\n", length(candidate_genes)))

nx_build_cor_matrix <- function(genes, ref_vec, subgroups) {
  mr <- matrix(NA_real_, length(genes), length(subgroups), dimnames = list(genes, names(subgroups)))
  mn <- matrix(NA_real_, length(genes), length(subgroups), dimnames = list(genes, names(subgroups)))
  for (g in genes) {
    fv <- setNames(as.numeric(expr[g, ]), colnames(expr))
    for (j in seq_along(subgroups)) {
      s <- intersect(intersect(names(fv), subgroups[[j]]), names(ref_vec))
      if (length(s) >= 3) {
        x <- fv[s]; y <- ref_vec[s]; ok <- is.finite(x) & is.finite(y); nv <- sum(ok)
        mn[g, j] <- nv
        if (nv >= 3 && sd(x[ok]) > 0 && sd(y[ok]) > 0) mr[g, j] <- cor(x[ok], y[ok])
      } else mn[g, j] <- 0
    }
  }
  list(r = mr, n = mn)
}

subgroups_ad  <- list(Overall = valid_all, "CN-Lo" = lo_ids, "CN-Hi" = hi_ids, AD  = ad_ids)
subgroups_mci <- list(Overall = valid_all, "CN-Lo" = lo_ids, "CN-Hi" = hi_ids, MCI = mci_ids)

cat("  Building correlation matrices (AD + MCI)...\n")
mat_10_rna   <- nx_build_cor_matrix(candidate_genes, np_rna,  subgroups_ad)
mat_10_prot  <- nx_build_cor_matrix(candidate_genes, np_prot, subgroups_ad)
mat_mci_rna  <- nx_build_cor_matrix(candidate_genes, np_rna,  subgroups_mci)
mat_mci_prot <- nx_build_cor_matrix(candidate_genes, np_prot, subgroups_mci)

nx_layout_tristate <- function(mat_obj, third_col = "AD",
                               t_strong = TRI_T_STRONG, t_weak = TRI_T_WEAK) {
  r_lo <- mat_obj$r[, "CN-Lo"]; r_hi <- mat_obj$r[, "CN-Hi"]; r_3 <- mat_obj$r[, third_col]
  r_lo[is.na(r_lo)] <- 0; r_hi[is.na(r_hi)] <- 0; r_3[is.na(r_3)] <- 0
  abs_lo <- abs(r_lo); abs_hi <- abs(r_hi); abs_3 <- abs(r_3)
  row_groups <- dplyr::case_when(
    (abs_lo > t_strong) & (abs_hi > t_strong) & (abs_3 > t_strong) &
      (sign(r_lo) != sign(r_hi)) & (sign(r_hi) != sign(r_3)) ~ "CN-Hi-reversed",
    (abs_lo > t_strong) & (abs_hi > t_strong) & (abs_3 < t_weak) &
      (sign(r_lo) != sign(r_hi)) ~ "CN-Hi-reversed",
    (abs_lo > t_strong) & (abs_hi > t_strong) & (sign(r_lo) == sign(r_hi)) ~ "CN-Hi-preserved",
    (abs_lo > t_strong) & (abs_hi < t_weak) & (abs_3 > t_strong) &
      (sign(r_lo) == sign(r_3)) ~ "CN-Hi-suppressed",
    (abs_lo < t_weak) & (abs_hi > t_strong) & (abs_3 < t_weak) ~ "CN-Hi-recruited",
    (abs_lo > t_strong) & (abs_hi < t_weak) & (abs_3 < t_weak) ~ "Pathology-disrupted",
    TRUE ~ "Other")
  pattern_levels <- c("CN-Hi-preserved", "CN-Hi-recruited", "CN-Hi-suppressed",
                      "CN-Hi-reversed", "Pathology-disrupted", "Other")
  row_groups <- factor(row_groups, levels = pattern_levels)
  sort_score <- dplyr::case_when(
    row_groups == "CN-Hi-preserved"    ~ r_lo,
    row_groups == "Pathology-disrupted" ~ r_lo,
    row_groups == "CN-Hi-suppressed"   ~ r_lo,
    row_groups == "CN-Hi-reversed"     ~ r_hi,
    row_groups == "CN-Hi-recruited"    ~ r_hi,
    TRUE ~ mat_obj$r[, "Overall"])
  df <- data.frame(Gene = rownames(mat_obj$r), Group = row_groups,
                   Score = sort_score, stringsAsFactors = FALSE)
  df <- df[df$Group != "Other", ]
  df <- df[order(df$Group, -df$Score), ]
  list(gene_order = df$Gene, row_groups_sorted = droplevels(df$Group))
}

layout_rna      <- nx_layout_tristate(mat_10_rna,   "AD")
layout_prot     <- nx_layout_tristate(mat_10_prot,  "AD")
layout_rna_mci  <- nx_layout_tristate(mat_mci_rna,  "MCI")
layout_prot_mci <- nx_layout_tristate(mat_mci_prot, "MCI")
cat("  Tri-state counts (RNA/AD):\n"); print(table(layout_rna$row_groups_sorted))
cat("  Tri-state counts (PROT/AD):\n"); print(table(layout_prot$row_groups_sorted))

nx_pattern_signed <- function(layout_obj, mat_obj) {
  genes <- layout_obj$gene_order; pat <- as.character(layout_obj$row_groups_sorted)
  r_lo <- mat_obj$r[genes, "CN-Lo"]; r_hi <- mat_obj$r[genes, "CN-Hi"]
  driving_r <- dplyr::case_when(
    pat == "CN-Hi-preserved"    ~ ifelse(abs(r_lo) >= abs(r_hi), r_lo, r_hi),
    pat == "CN-Hi-suppressed"   ~ r_lo,
    pat == "CN-Hi-recruited"    ~ r_hi,
    pat == "CN-Hi-reversed"     ~ r_hi,
    pat == "Pathology-disrupted" ~ r_lo,
    TRUE ~ r_lo)
  sgn <- ifelse(driving_r > 0, "+", ifelse(driving_r < 0, "-", "0"))
  setNames(paste0(pat, " (", sgn, ")"), genes)
}
ps_rna_ad  <- nx_pattern_signed(layout_rna,      mat_10_rna)
ps_rna_mci <- nx_pattern_signed(layout_rna_mci,  mat_mci_rna)
ps_prot_ad  <- nx_pattern_signed(layout_prot,     mat_10_prot)
ps_prot_mci <- nx_pattern_signed(layout_prot_mci, mat_mci_prot)
sh_rna  <- intersect(names(ps_rna_ad),  names(ps_rna_mci))
sh_prot <- intersect(names(ps_prot_ad), names(ps_prot_mci))
matched_rna_genes  <- sh_rna[ps_rna_ad[sh_rna]   == ps_rna_mci[sh_rna]]
matched_prot_genes <- sh_prot[ps_prot_ad[sh_prot] == ps_prot_mci[sh_prot]]
cat(sprintf("  RNA matched (Pattern+Sign in AD & MCI):  %d\n", length(matched_rna_genes)))
cat(sprintf("  PROT matched (Pattern+Sign in AD & MCI): %d\n", length(matched_prot_genes)))

nx_traj_from_pattern <- function(pat) dplyr::case_when(
  grepl("preserved",  pat, ignore.case = TRUE) ~ "A",
  grepl("recruited",  pat, ignore.case = TRUE) ~ "B_rec",
  grepl("suppressed", pat, ignore.case = TRUE) ~ "B_sup",
  grepl("reversed",   pat, ignore.case = TRUE) ~ "B_rev",
  grepl("Pathology",  pat, ignore.case = TRUE) ~ "C",
  TRUE ~ NA_character_)

nx_classify_full <- function(layout_obj, mat_obj, modality_tag, third_col = "AD") {
  go <- layout_obj$gene_order; rg <- as.character(layout_obj$row_groups_sorted)
  df <- data.frame(Gene = go, Pattern = rg, Traj = nx_traj_from_pattern(rg),
                   r_Overall = mat_obj$r[go, "Overall"],
                   `r_CN-Lo` = mat_obj$r[go, "CN-Lo"],
                   `r_CN-Hi` = mat_obj$r[go, "CN-Hi"],
                   check.names = FALSE, stringsAsFactors = FALSE)
  df[[paste0("r_", third_col)]] <- mat_obj$r[go, third_col]
  df <- df %>% dplyr::mutate(
    Driving_r = dplyr::case_when(
      Pattern == "CN-Hi-preserved"    ~ ifelse(abs(`r_CN-Lo`) >= abs(`r_CN-Hi`), `r_CN-Lo`, `r_CN-Hi`),
      Pattern == "CN-Hi-suppressed"   ~ `r_CN-Lo`,
      Pattern == "CN-Hi-recruited"    ~ `r_CN-Hi`,
      Pattern == "CN-Hi-reversed"     ~ `r_CN-Hi`,
      Pattern == "Pathology-disrupted" ~ `r_CN-Lo`,
      TRUE ~ `r_CN-Lo`),
    Sign = dplyr::case_when(Driving_r > 0 ~ "Positive", Driving_r < 0 ~ "Negative", TRUE ~ "Zero"),
    Pattern_Signed = paste0(Pattern, " (",
                            ifelse(Sign == "Positive", "+", ifelse(Sign == "Negative", "-", "0")), ")"))
  lvl <- c("CN-Hi-preserved", "CN-Hi-recruited", "CN-Hi-suppressed",
           "CN-Hi-reversed", "Pathology-disrupted")
  df$Pattern <- factor(df$Pattern, levels = lvl)
  df$Sign    <- factor(df$Sign,    levels = c("Positive", "Negative", "Zero"))
  df %>% dplyr::arrange(Pattern, Sign, dplyr::desc(abs(Driving_r))) %>% dplyr::mutate(Modality = modality_tag)
}

full_rna      <- nx_classify_full(layout_rna,      mat_10_rna,   "RNA",  "AD")
full_prot     <- nx_classify_full(layout_prot,     mat_10_prot,  "PROT", "AD")
full_rna_mci  <- nx_classify_full(layout_rna_mci,  mat_mci_rna,  "RNA",  "MCI")
full_prot_mci <- nx_classify_full(layout_prot_mci, mat_mci_prot, "PROT", "MCI")
full_rna$Match_MCI_AD      <- as.integer(full_rna$Gene     %in% matched_rna_genes)
full_prot$Match_MCI_AD     <- as.integer(full_prot$Gene    %in% matched_prot_genes)
full_rna_mci$Match_MCI_AD  <- as.integer(full_rna_mci$Gene %in% matched_rna_genes)
full_prot_mci$Match_MCI_AD <- as.integer(full_prot_mci$Gene %in% matched_prot_genes)

readr::write_csv(full_rna,      file.path(analysis_dir, "Fig2e3f_FullGeneList_RNA_AD.csv"))
readr::write_csv(full_prot,     file.path(analysis_dir, "Fig2e3f_FullGeneList_PROT_AD.csv"))
readr::write_csv(full_rna_mci,  file.path(analysis_dir, "Fig2e3f_FullGeneList_RNA_MCI.csv"))
readr::write_csv(full_prot_mci, file.path(analysis_dir, "Fig2e3f_FullGeneList_PROT_MCI.csv"))

nx_check(full_rna,      "Table_24_FullGeneList_RNA_AD.csv")
nx_check(full_prot,     "Table_24_FullGeneList_PROT_AD.csv")
nx_check(full_rna_mci,  "Table_24_FullGeneList_RNA_MCI.csv")
nx_check(full_prot_mci, "Table_24_FullGeneList_PROT_MCI.csv")
cat("\nFigure 2e/3f master classification done. FullGeneList is the single classification source.\n")

nx_rr_view <- function(full_df)
  data.frame(Gene = full_df$Gene,
             r_Lo = full_df[["r_CN-Lo"]], r_Hi = full_df[["r_CN-Hi"]], r_AD = full_df$r_AD,
             stringsAsFactors = FALSE)

nx_check(nx_rr_view(full_rna),  "Table_25_RvsR_RNA_TrajGenes.csv")
nx_check(nx_rr_view(full_prot), "Table_25_RvsR_PROT_TrajGenes.csv")
cat("\nR-to-R (2c,d/3d,e) confirmed: its values are columns of the FullGeneList.\n")

df_3a <- meta %>%
  dplyr::mutate(NPTX2_RNA = as.numeric(expr["NPTX2", SampleID]),
                NPTX2_MS  = as.numeric(prot$MS_NPTX2[match(SampleID, prot$SampleID)]),
                Group     = as.character(Group_Current)) %>%
  dplyr::filter(is.finite(NPTX2_RNA), is.finite(NPTX2_MS),
                Group %in% c("CN-Lo", "CN-Hi", "MCI", "AD"))
df_3a$Group <- factor(df_3a$Group, levels = c("CN-Lo", "CN-Hi", "MCI", "AD"))

nx_cor_bhp <- function(dat, grp) {
  if (nrow(dat) < 3) return(data.frame(Group = grp, r = NA_real_, p = NA_real_, n = nrow(dat)))
  ct <- cor.test(dat$NPTX2_RNA, dat$NPTX2_MS, method = "pearson")
  data.frame(Group = grp, r = unname(ct$estimate), p = ct$p.value, n = nrow(dat),
             stringsAsFactors = FALSE)
}
cor_all_3a <- dplyr::bind_rows(
  nx_cor_bhp(df_3a, "Overall"),
  dplyr::bind_rows(lapply(levels(df_3a$Group), function(g) nx_cor_bhp(df_3a[df_3a$Group == g, ], g)))
)
cor_all_3a$BHP <- p.adjust(cor_all_3a$p, method = "BH")

readr::write_csv(cor_all_3a, file.path(analysis_dir, "Fig3a_cor_stats.csv"))
readr::write_csv(df_3a[, c("SampleID", "NPTX2_RNA", "NPTX2_MS", "Group")],
                 file.path(analysis_dir, "Fig3a_points.csv"))

nx_check(cor_all_3a, "Table_23_NPTX2_mRNA_Protein_Cor.csv")
cat("\nFigure 3a scatter stats done.\n")

nx_check_ora <- function(new_df, ref_name, tol = 1e-8) {
  f <- file.path(ref_dir, ref_name)
  if (!file.exists(f)) { cat(sprintf("[check] %-46s no existing file\n", ref_name))
    nx_record(ref_name, "no existing file"); return(invisible(NULL)) }
  old <- as.data.frame(readr::read_csv(f, show_col_types = FALSE))
  names(old) <- gsub("HiPath", "CN-Hi", gsub("LoPath", "CN-Lo", names(old)))
  if ("Cluster" %in% names(old))
    old$Cluster <- gsub("HiPath", "CN-Hi", gsub("LoPath", "CN-Lo", as.character(old$Cluster)))
  k_new <- paste(new_df$Cluster, new_df$ID); k_old <- paste(old$Cluster, old$ID)
  old <- old[match(k_new, k_old), , drop = FALSE]
  common <- intersect(names(old), names(new_df))
  num <- common[vapply(new_df[common], is.numeric, logical(1))]
  diffs <- num[vapply(num, function(cc)
    !isTRUE(all.equal(old[[cc]], new_df[[cc]], tolerance = tol)), logical(1))]
  status <- if (!length(diffs)) sprintf("IDENTICAL (%d cols, by Cluster+ID)", length(num))
  else sprintf("DIFFERS: %s", paste(head(diffs, 6), collapse = ", "))
  cat(sprintf("[check] %-46s %s\n", ref_name, status)); nx_record(ref_name, status)
  invisible(list(same = !length(diffs), diffs = diffs))
}

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

  nx_ora_lists <- function(layout_obj, min_genes = 5) {
    grps <- setdiff(levels(layout_obj$row_groups_sorted), "Other")
    out <- list()
    for (g in grps) {
      genes <- setdiff(layout_obj$gene_order[layout_obj$row_groups_sorted == g], "NPTX2")
      if (length(genes) >= min_genes) out[[g]] <- genes
    }
    out
  }

  nx_ora_compute <- function(pattern_list) {
    if (length(pattern_list) == 0) return(NULL)
    cc <- clusterProfiler::compareCluster(
      geneCluster = pattern_list, fun = "enrichGO",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
      minGSSize = 15, maxGSSize = 500)
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) return(NULL)
    d <- as.data.frame(cc)
    ratio_num <- function(x) vapply(strsplit(as.character(x), "/"), function(p) {
      p <- as.numeric(p); if (length(p) == 2 && p[2] != 0) p[1] / p[2] else NA_real_ }, numeric(1))
    d$GeneRatio_num <- ratio_num(d$GeneRatio)
    d$BgRatio_num   <- ratio_num(d$BgRatio)
    d$FoldEnrichment <- d$GeneRatio_num / d$BgRatio_num
    d
  }

  nx_ora_matched <- function(d, matched) {
    if (is.null(d) || nrow(d) == 0) return(d)
    d$n_genes_in_pathway   <- vapply(strsplit(as.character(d$geneID), "/"), length, integer(1))
    d$n_matched_in_pathway <- vapply(strsplit(as.character(d$geneID), "/"),
                                     function(h) sum(h %in% matched), integer(1))
    d$frac_matched <- ifelse(d$n_genes_in_pathway > 0,
                             d$n_matched_in_pathway / d$n_genes_in_pathway, 0)
    d$is_concordant <- d$frac_matched >= 0.5
    d
  }

  ora_rna_ad   <- nx_ora_matched(nx_ora_compute(nx_ora_lists(layout_rna)),      matched_rna_genes)
  ora_prot_ad  <- nx_ora_matched(nx_ora_compute(nx_ora_lists(layout_prot)),     matched_prot_genes)
  ora_rna_mci  <- nx_ora_compute(nx_ora_lists(layout_rna_mci))
  ora_prot_mci <- nx_ora_compute(nx_ora_lists(layout_prot_mci))

  nx_save_ora <- function(d, fname, ref_name) {
    if (is.null(d)) { cat(sprintf("  ORA %s: no significant CC pathways\n", fname)); return(invisible(NULL)) }
    d <- d %>% dplyr::arrange(Cluster, dplyr::desc(FoldEnrichment))
    readr::write_csv(d, file.path(analysis_dir, fname))
    nx_check_ora(d, ref_name)
  }
  nx_save_ora(ora_rna_ad,   "Fig3g_ORA_RNA_AD_CC.csv",   "Table_11_ORA_RNA_AD_CC_Master.csv")
  nx_save_ora(ora_rna_mci,  "Fig3g_ORA_RNA_MCI_CC.csv",  "Table_11_ORA_RNA_MCI_CC_Master.csv")
  nx_save_ora(ora_prot_ad,  "Fig2f_ORA_PROT_AD_CC.csv",  "Table_11_ORA_PROT_AD_CC_Master.csv")
  nx_save_ora(ora_prot_mci, "Fig2f_ORA_PROT_MCI_CC.csv", "Table_11_ORA_PROT_MCI_CC_Master.csv")
  cat("\nFigure 2f/3g ORA done.\n")

} else {
  cat("\n[skip] Figure 2f/3g ORA: clusterProfiler / org.Hs.eg.db not installed.\n",
      "       install via BiocManager::install(c('clusterProfiler','org.Hs.eg.db')) and re-run.\n")
}

cfg_de <- list(delta_eq = 0.50, alpha_tost = 0.05, p_de_max = 0.05, fc_threshold = 0.50)

nx_tost <- function(x, y, delta) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2)
    return(list(diff = NA_real_, p_lower = NA_real_, p_upper = NA_real_,
                p_tost = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  d <- mean(x) - mean(y); vx <- var(x); vy <- var(y); nx <- length(x); ny <- length(y)
  se <- sqrt(vx / nx + vy / ny)
  if (!is.finite(se) || se == 0)
    return(list(diff = d, p_lower = NA_real_, p_upper = NA_real_,
                p_tost = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  dfw <- (vx / nx + vy / ny)^2 / ((vx / nx)^2 / (nx - 1) + (vy / ny)^2 / (ny - 1))
  p_lower <- pt((d - (-delta)) / se, df = dfw, lower.tail = FALSE)
  p_upper <- pt((d - delta) / se,    df = dfw, lower.tail = TRUE)
  t_ci <- qt(1 - 0.05, df = dfw)
  list(diff = d, p_lower = p_lower, p_upper = p_upper, p_tost = max(p_lower, p_upper),
       ci_lo = d - t_ci * se, ci_hi = d + t_ci * se)
}

nx_de_stats <- function(gene) {
  v_lo <- as.numeric(expr[gene, lo_ids]); v_hi <- as.numeric(expr[gene, hi_ids])
  v_ad <- as.numeric(expr[gene, ad_ids]); v_mci <- as.numeric(expr[gene, mci_ids])
  m_lo <- mean(v_lo, na.rm = TRUE); m_hi <- mean(v_hi, na.rm = TRUE)
  m_ad <- mean(v_ad, na.rm = TRUE); m_mci <- mean(v_mci, na.rm = TRUE)
  t_lh <- nx_tost(v_lo, v_hi, cfg_de$delta_eq)
  p_hi_ad  <- tryCatch(t.test(v_hi, v_ad)$p.value,  error = function(e) NA_real_)
  p_hi_mci <- tryCatch(t.test(v_hi, v_mci)$p.value, error = function(e) NA_real_)
  r_lo <- tryCatch(cor(v_lo, as.numeric(expr["NPTX2", lo_ids]), use = "pairwise.complete.obs"), error = function(e) NA_real_)
  r_hi <- tryCatch(cor(v_hi, as.numeric(expr["NPTX2", hi_ids]), use = "pairwise.complete.obs"), error = function(e) NA_real_)
  r_ad <- tryCatch(cor(v_ad, as.numeric(expr["NPTX2", ad_ids]), use = "pairwise.complete.obs"), error = function(e) NA_real_)
  mar <- suppressWarnings(max(abs(c(r_lo, r_hi, r_ad)), na.rm = TRUE)); if (!is.finite(mar)) mar <- NA_real_
  data.frame(Gene = gene, Mean_Lo = m_lo, Mean_Hi = m_hi, Mean_AD = m_ad, Mean_MCI = m_mci,
             Delta_Lo = m_lo - m_hi, Delta_AD = m_ad - m_hi, Delta_MCI = m_mci - m_hi,
             TOST_diff = t_lh$diff, TOST_p_lo = t_lh$p_lower, TOST_p_up = t_lh$p_upper,
             TOST_p = t_lh$p_tost, TOST_CI_lo = t_lh$ci_lo, TOST_CI_hi = t_lh$ci_hi,
             Equivalent = !is.na(t_lh$p_tost) & t_lh$p_tost < cfg_de$alpha_tost,
             p_hi_ad = p_hi_ad, p_hi_mci = p_hi_mci, Max_Abs_R = mar, stringsAsFactors = FALSE)
}

de_target <- intersect(as.character(layout_rna$gene_order), rownames(expr))
df_de <- data.frame(Gene = de_target,
                    Pattern = as.character(layout_rna$row_groups_sorted)[match(de_target, layout_rna$gene_order)],
                    stringsAsFactors = FALSE) %>%
  dplyr::left_join(dplyr::bind_rows(lapply(de_target, nx_de_stats)), by = "Gene") %>%
  dplyr::mutate(Abs_Y_AD = abs(Delta_AD), Abs_Y_MCI = abs(Delta_MCI),
                Is_Resilient_AD  = Equivalent & !is.na(p_hi_ad)  & p_hi_ad  < cfg_de$p_de_max & Abs_Y_AD > cfg_de$fc_threshold,
                Is_Resilient_MCI = Equivalent & !is.na(p_hi_mci) & p_hi_mci < cfg_de$p_de_max,
                Is_Traj_Matched  = Gene %in% matched_rna_genes,
                Is_Resilient_Both = Is_Resilient_AD & Is_Resilient_MCI & Is_Traj_Matched,
                Label_Text = paste0(Gene, " (", Pattern, ")")) %>%
  dplyr::arrange(dplyr::desc(Is_Resilient_Both), dplyr::desc(Is_Resilient_AD & Is_Resilient_MCI),
                 dplyr::desc(Is_Resilient_AD), p_hi_ad)

m_lo_all <- rowMeans(expr[, lo_ids], na.rm = TRUE); m_hi_all <- rowMeans(expr[, hi_ids], na.rm = TRUE)
m_ad_all <- rowMeans(expr[, ad_ids], na.rm = TRUE); m_mci_all <- rowMeans(expr[, mci_ids], na.rm = TRUE)
ghost_genes <- setdiff(rownames(expr), c(de_target, "NPTX2"))
df_ghost <- data.frame(Gene = ghost_genes,
                       Delta_Lo  = m_lo_all[ghost_genes]  - m_hi_all[ghost_genes],
                       Delta_AD  = m_ad_all[ghost_genes]  - m_hi_all[ghost_genes],
                       Delta_MCI = m_mci_all[ghost_genes] - m_hi_all[ghost_genes],
                       stringsAsFactors = FALSE)
df_nptx2_anchor <- data.frame(Gene = "NPTX2",
                              Delta_Lo  = m_lo_all["NPTX2"]  - m_hi_all["NPTX2"],
                              Delta_AD  = m_ad_all["NPTX2"]  - m_hi_all["NPTX2"],
                              Delta_MCI = m_mci_all["NPTX2"] - m_hi_all["NPTX2"])

readr::write_csv(df_de,          file.path(analysis_dir, "Fig4a_Double_DE.csv"))
readr::write_csv(df_ghost,       file.path(analysis_dir, "Fig4a_ghost.csv"))
readr::write_csv(df_nptx2_anchor, file.path(analysis_dir, "Fig4a_nptx2_anchor.csv"))
nx_check(df_de, "Table_14_Double_DE_RNA_AD_and_MCI.csv")
cat("\nFigure 4a Double DE done.\n")

genes_4bi <- unique(df_de$Gene[df_de$Is_Resilient_Both %in% c(TRUE, 1, "TRUE")])
genes_4bi <- genes_4bi[!is.na(genes_4bi)]
cat(sprintf("  Fig 4b-m: %d resilient-both genes -> %s\n", length(genes_4bi), paste(genes_4bi, collapse = ", ")))

nx_cor_ci <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4) return(list(r = NA, lower = NA, upper = NA, n = n, z = NA))
  r <- cor(x, y); rc <- pmax(pmin(r, 0.999), -0.999); z <- 0.5 * log((1 + rc) / (1 - rc)); se <- 1 / sqrt(n - 3)
  rci <- function(zz) (exp(2 * zz) - 1) / (exp(2 * zz) + 1)
  list(r = r, lower = rci(z - 1.96 * se), upper = rci(z + 1.96 * se), n = n, z = z)
}
nx_cor_cmp <- function(a, b) {
  if (is.na(a$z) || is.na(b$z)) return(list(p = NA, stars = "ns"))
  zd <- (a$z - b$z) / sqrt(1 / (a$n - 3) + 1 / (b$n - 3)); p <- 2 * (1 - pnorm(abs(zd)))
  list(p = p, stars = if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns")
}
nx_expr_cmp <- function(v1, v2) {
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) < 3 || length(v2) < 3) return(list(p = NA, stars = "ns"))
  p <- tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA_real_)
  list(p = p, stars = if (is.na(p)) "ns" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns")
}
nx_ci95 <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2) return(list(mean = NA, lower = NA, upper = NA))
  m <- mean(x); se <- sd(x) / sqrt(n); tc <- qt(0.975, df = n - 1)
  list(mean = m, lower = m - tc * se, upper = m + tc * se)
}
nx_minmax <- function(v) { f <- v[is.finite(v)]; lo <- min(f); hi <- max(f)
if (is.finite(lo) && is.finite(hi) && hi > lo) (v - lo) / (hi - lo) else rep(NA_real_, length(v)) }

all_s_4bm <- intersect(meta$SampleID, colnames(expr))

val_rows <- list(); cmp_rows <- list()
for (g in genes_4bi) {
  if (!g %in% rownames(expr)) { cat(sprintf("  [skip] %s not in expr\n", g)); next }
  gv <- setNames(as.numeric(expr[g, all_s_4bm]), all_s_4bm)

  s_lo  <- nx_cor_ci(gv[sg_cnlo], np_rna[sg_cnlo]); s_hi  <- nx_cor_ci(gv[sg_cnhi], np_rna[sg_cnhi])
  s_mci <- nx_cor_ci(gv[sg_mci],  np_rna[sg_mci]);  s_ad  <- nx_cor_ci(gv[sg_ad],  np_rna[sg_ad])

  gs <- nx_minmax(gv); names(gs) <- names(gv)
  e_lo  <- nx_ci95(gs[sg_cnlo]); e_hi <- nx_ci95(gs[sg_cnhi]); e_mci <- nx_ci95(gs[sg_mci]); e_ad <- nx_ci95(gs[sg_ad])
  val_rows[[g]] <- data.frame(Gene = g, State = c("CN-Lo", "CN-Hi", "MCI", "AD"),
                              Cor_r    = c(s_lo$r,     s_hi$r,     s_mci$r,     s_ad$r),
                              Cor_ymin = c(s_lo$lower, s_hi$lower, s_mci$lower, s_ad$lower),
                              Cor_ymax = c(s_lo$upper, s_hi$upper, s_mci$upper, s_ad$upper),
                              Expr_mean = c(e_lo$mean, e_hi$mean, e_mci$mean, e_ad$mean),
                              Expr_ymin = c(e_lo$lower, e_hi$lower, e_mci$lower, e_ad$lower),
                              Expr_ymax = c(e_lo$upper, e_hi$upper, e_mci$upper, e_ad$upper),
                              stringsAsFactors = FALSE)

  cc1 <- nx_cor_cmp(s_lo, s_hi); cc2 <- nx_cor_cmp(s_hi, s_mci); cc3 <- nx_cor_cmp(s_hi, s_ad)
  ec1 <- nx_expr_cmp(gs[sg_cnlo], gs[sg_cnhi]); ec2 <- nx_expr_cmp(gs[sg_cnhi], gs[sg_mci]); ec3 <- nx_expr_cmp(gs[sg_cnhi], gs[sg_ad])
  cmp_rows[[g]] <- data.frame(Gene = g, Comparison = c("CN-Lo_vs_CN-Hi", "CN-Hi_vs_MCI", "CN-Hi_vs_AD"),
                              Cor_p = c(cc1$p, cc2$p, cc3$p), Cor_stars = c(cc1$stars, cc2$stars, cc3$stars),
                              Expr_p = c(ec1$p, ec2$p, ec3$p), Expr_stars = c(ec1$stars, ec2$stars, ec3$stars), stringsAsFactors = FALSE)
}
readr::write_csv(dplyr::bind_rows(val_rows), file.path(analysis_dir, "Fig4bm_trajectory_values.csv"))
readr::write_csv(dplyr::bind_rows(cmp_rows), file.path(analysis_dir, "Fig4bm_comparison_stats.csv"))
cat(sprintf("\nFigure 4b-m trajectory stats done (%d genes, 4 groups).\n", length(genes_4bi)))

nx_check_keyed <- function(new_df, ref_name, keys, harmonize = character(0), tol = 1e-6) {
  f <- file.path(ref_dir, ref_name)
  if (!file.exists(f)) { cat(sprintf("[check] %-46s no existing file\n", ref_name))
    nx_record(ref_name, "no existing file"); return(invisible(NULL)) }
  old <- as.data.frame(readr::read_csv(f, show_col_types = FALSE))
  names(old) <- gsub("HiPath", "CN-Hi", gsub("LoPath", "CN-Lo", names(old)))
  for (cc in harmonize) if (cc %in% names(old))
    old[[cc]] <- gsub("HiPath", "CN-Hi", gsub("LoPath", "CN-Lo", as.character(old[[cc]])))
  kn <- do.call(paste, c(lapply(keys, function(k) as.character(new_df[[k]])), sep = "|"))
  ko <- do.call(paste, c(lapply(keys, function(k) as.character(old[[k]])),    sep = "|"))
  old <- old[match(kn, ko), , drop = FALSE]
  common <- intersect(names(old), names(new_df))
  num <- common[vapply(new_df[common], is.numeric, logical(1))]
  diffs <- num[vapply(num, function(cc) !isTRUE(all.equal(old[[cc]], new_df[[cc]], tolerance = tol)), logical(1))]
  status <- if (!length(diffs)) sprintf("IDENTICAL (%d cols, by %s)", length(num), paste(keys, collapse = "+"))
  else sprintf("DIFFERS: %s", paste(head(diffs, 6), collapse = ", "))
  cat(sprintf("[check] %-46s %s\n", ref_name, status)); nx_record(ref_name, status)
  invisible(list(same = !length(diffs), diffs = diffs))
}

cfg_de_supp <- list(fdr_thresh = 0.25, fc_thresh = 0.20)

protein_panel <- c("ATP6V1H", "DPP6", "GRIA4", "GRIN2B", "HNRNPA2B1", "HOMER1",
                   "KIAA1045", "LAMP1", "LAMP2", "NPTX1", "NPTX2", "NPTXR",
                   "RAB11A", "RAB5A", "RHEB", "RIMS1", "SQSTM1", "SYT1", "VDAC1", "VGF")
disp_name <- function(g) ifelse(g == "KIAA1045", "PHF24", g)

p_keep <- as.character(meta$Group_Current) %in% c("CN-Lo", "CN-Hi", "AD")
p_grp  <- factor(as.character(meta$Group_Current)[p_keep], levels = c("CN-Lo", "CN-Hi", "AD"))
p_prot <- prot[p_keep, , drop = FALSE]

nx_wilcox_one <- function(col, num_grp, den_grp) {
  vals <- suppressWarnings(as.numeric(p_prot[[col]]))
  x <- vals[p_grp == num_grp]; x <- x[is.finite(x)]
  y <- vals[p_grp == den_grp]; y <- y[is.finite(y)]
  if (length(x) < 3 || length(y) < 3)
    return(data.frame(log2FC = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                      P.Value = NA_real_, n_num = length(x), n_den = length(y)))
  ww <- tryCatch(suppressWarnings(wilcox.test(x, y, conf.int = TRUE, conf.level = 0.95)),
                 error = function(e) NULL)
  if (is.null(ww))
    return(data.frame(log2FC = mean(x) - mean(y), ci_lo = NA_real_, ci_hi = NA_real_,
                      P.Value = NA_real_, n_num = length(x), n_den = length(y)))
  data.frame(log2FC = mean(x) - mean(y), ci_lo = unname(ww$conf.int[1]),
             ci_hi = unname(ww$conf.int[2]), P.Value = ww$p.value,
             n_num = length(x), n_den = length(y), stringsAsFactors = FALSE)
}
nx_run_panel <- function(num_grp, den_grp, label) {
  out <- dplyr::bind_rows(lapply(protein_panel, function(g) {
    col <- paste0("MS_", g)
    if (!col %in% colnames(p_prot)) return(NULL)
    cbind(Gene = g, Display = disp_name(g), Contrast = label,
          nx_wilcox_one(col, num_grp, den_grp), stringsAsFactors = FALSE)
  }))
  out$adj.P.Val <- NA_real_
  ok <- is.finite(out$P.Value); out$adj.P.Val[ok] <- p.adjust(out$P.Value[ok], method = "BH")
  out
}
res_p_HL <- nx_run_panel("CN-Hi", "CN-Lo", "CN-Hi vs CN-Lo")
res_p_AH <- nx_run_panel("AD",    "CN-Hi", "AD vs CN-Hi")

pj <- dplyr::inner_join(
  res_p_HL %>% dplyr::transmute(Gene, FC_HivsLo = log2FC, FDR_HivsLo = adj.P.Val),
  res_p_AH %>% dplyr::transmute(Gene, FC_ADvsHi = log2FC, FDR_ADvsHi = adj.P.Val), by = "Gene") %>%
  dplyr::mutate(Sig_HivsLo = !is.na(FDR_HivsLo) & FDR_HivsLo < cfg_de_supp$fdr_thresh & abs(FC_HivsLo) > cfg_de_supp$fc_thresh,
                Sig_ADvsHi = !is.na(FDR_ADvsHi) & FDR_ADvsHi < cfg_de_supp$fdr_thresh & abs(FC_ADvsHi) > cfg_de_supp$fc_thresh,
                Shared = Sig_HivsLo & Sig_ADvsHi & (sign(FC_HivsLo) != sign(FC_ADvsHi)))
p_shared <- pj$Gene[pj$Shared]

nx_annotate_supp <- function(df) df %>% dplyr::mutate(
  Sig = !is.na(adj.P.Val) & adj.P.Val < cfg_de_supp$fdr_thresh & abs(log2FC) > cfg_de_supp$fc_thresh,
  Direction = dplyr::case_when(!Sig ~ "ns", log2FC > cfg_de_supp$fc_thresh ~ "up",
                               log2FC < -cfg_de_supp$fc_thresh ~ "down", TRUE ~ "ns"),
  Is_Shared = Gene %in% p_shared)
res_p_HL <- nx_annotate_supp(res_p_HL); res_p_AH <- nx_annotate_supp(res_p_AH)

protein_de <- dplyr::bind_rows(res_p_HL, res_p_AH) %>%
  dplyr::select(Gene, Display, Contrast, log2FC, ci_lo, ci_hi, P.Value, adj.P.Val,
                n_num, n_den, Direction, Is_Shared) %>%
  dplyr::arrange(Contrast, adj.P.Val)
readr::write_csv(protein_de, file.path(analysis_dir, "Supp1ab_Protein_Panel_DE.csv"))
nx_check_keyed(protein_de, "Table_33_Protein_Panel_DE.csv",
               keys = c("Contrast", "Gene"), harmonize = "Contrast")
cat("\nSupp 1a/1b protein-panel DE done.\n")

if (requireNamespace("limma", quietly = TRUE)) {
  d_keep <- as.character(meta$Group_Current) %in% c("CN-Lo", "CN-Hi", "AD")
  d_meta <- data.frame(SampleID = meta$SampleID[d_keep],
                       Group = factor(as.character(meta$Group_Current)[d_keep],
                                      levels = c("CN-Lo", "CN-Hi", "AD")))
  d_expr <- expr[, d_meta$SampleID, drop = FALSE]

  rv <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::rowVars(d_expr, na.rm = TRUE)
  else { mu <- rowMeans(d_expr, na.rm = TRUE); rowSums((d_expr - mu)^2, na.rm = TRUE) / (ncol(d_expr) - 1) }
  has_inf <- apply(d_expr, 1, function(x) any(is.infinite(x)))
  d_expr <- d_expr[is.finite(rv) & rv > 0 & !has_inf, , drop = FALSE]

  if (exists("valid_pc_symbols")) {
    keep_pc <- !grepl("^ENSG", rownames(d_expr)) &
      rownames(d_expr) %in% valid_pc_symbols & !grepl("^OR[0-9]+[A-Z]", rownames(d_expr))
    d_expr <- d_expr[keep_pc, , drop = FALSE]
  } else cat("  [warn] valid_pc_symbols missing; limma run without PC filter.\n")
  cat(sprintf("  limma DE on %d genes (CN-Lo=%d CN-Hi=%d AD=%d)\n", nrow(d_expr),
              sum(d_meta$Group == "CN-Lo"), sum(d_meta$Group == "CN-Hi"), sum(d_meta$Group == "AD")))

  design <- model.matrix(~ 0 + Group, data = d_meta)
  colnames(design) <- c("CN_Lo", "CN_Hi", "AD")
  cmat <- limma::makeContrasts(HivsLo = CN_Hi - CN_Lo, ADvsHi = AD - CN_Hi, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(d_expr, design), cmat),
                        trend = FALSE, robust = TRUE)
  nx_top <- function(coef) {
    tt <- limma::topTable(fit2, coef = coef, number = Inf, sort.by = "none")
    data.frame(Gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
               t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }
  de_HivsLo <- nx_top("HivsLo")
  de_ADvsHi <- nx_top("ADvsHi")
  flag_rna <- function(df) df %>% dplyr::mutate(
    Sig = adj.P.Val < cfg_de_supp$fdr_thresh & abs(log2FC) > cfg_de_supp$fc_thresh,
    Direction = dplyr::case_when(!Sig ~ "ns", log2FC > cfg_de_supp$fc_thresh ~ "up",
                                 log2FC < -cfg_de_supp$fc_thresh ~ "down", TRUE ~ "ns"))
  de_HivsLo <- flag_rna(de_HivsLo); de_ADvsHi <- flag_rna(de_ADvsHi)

  rj <- dplyr::inner_join(
    de_HivsLo %>% dplyr::transmute(Gene, FC_HivsLo = log2FC, FDR_HivsLo = adj.P.Val, Sig_HivsLo = Sig),
    de_ADvsHi %>% dplyr::transmute(Gene, FC_ADvsHi = log2FC, FDR_ADvsHi = adj.P.Val, Sig_ADvsHi = Sig), by = "Gene") %>%
    dplyr::mutate(Both_Sig = Sig_HivsLo & Sig_ADvsHi,
                  Shared_Concordant = Both_Sig & (sign(FC_HivsLo) != sign(FC_ADvsHi)),
                  Combined_Score = -log10(pmax(FDR_HivsLo, 1e-300)) + -log10(pmax(FDR_ADvsHi, 1e-300)))
  shared_table <- rj %>% dplyr::filter(Shared_Concordant) %>%
    dplyr::mutate(Direction_HivsLo = ifelse(FC_HivsLo > 0, "up", "down"),
                  Direction_ADvsHi = ifelse(FC_ADvsHi > 0, "up", "down"),
                  Lo_AD_relative_to_Hi = ifelse(FC_HivsLo > 0, "below Hi (gene up in Hi)", "above Hi (gene down in Hi)")) %>%
    dplyr::select(Gene, FC_HivsLo, FDR_HivsLo, Direction_HivsLo, FC_ADvsHi, FDR_ADvsHi,
                  Direction_ADvsHi, Lo_AD_relative_to_Hi, Combined_Score) %>%
    dplyr::arrange(dplyr::desc(Combined_Score))

  readr::write_csv(de_HivsLo %>% dplyr::arrange(adj.P.Val), file.path(analysis_dir, "Supp4ab_DE_HivsLo.csv"))
  readr::write_csv(de_ADvsHi %>% dplyr::arrange(adj.P.Val), file.path(analysis_dir, "Supp4ab_DE_ADvsHi.csv"))
  readr::write_csv(shared_table, file.path(analysis_dir, "Supp4ab_Shared_Concordant.csv"))
  nx_check(de_HivsLo %>% dplyr::arrange(adj.P.Val), "Table_32_DE_HivsLo.csv")
  nx_check(de_ADvsHi %>% dplyr::arrange(adj.P.Val), "Table_32_DE_ADvsHi.csv")
  nx_check(shared_table, "Table_32_Shared_Concordant.csv")
  cat("\nSupp 4a/4b RNA limma DE done.\n")
} else {
  cat("\n[skip] Supp 4a/4b limma DE: 'limma' not installed (BiocManager::install('limma')).\n")
}

rna_genes_20 <- c("ATP6V1H", "DPP6", "GRIA4", "GRIN2B", "HNRNPA2B1", "HOMER1", "PHF24",
                  "LAMP1", "LAMP2", "NPTX1", "NPTX2", "NPTXR", "RAB11A", "RAB5A",
                  "RHEB", "RIMS1", "SQSTM1", "SYT1", "VDAC1", "VGF")
rna_avail <- rna_genes_20[rna_genes_20 %in% rownames(expr)]
nx_keep_cov <- function(M) vapply(seq_len(ncol(M)), function(j) {
  v <- M[, j]; ok <- is.finite(v); sum(ok) >= 3 && sd(v[ok]) > 0 }, logical(1))

cov_rna <- cbind(PC1 = as.numeric(meta$PC1), PC2 = as.numeric(meta$PC2),
                 RIN = as.numeric(meta$RIN), PMI = as.numeric(meta$PMI),
                 Age = suppressWarnings(as.numeric(meta$Age)),
                 Braak = as.numeric(meta$SS_B), CERAD = as.numeric(meta$SS_C))
cov_rna <- cov_rna[, nx_keep_cov(cov_rna), drop = FALSE]
R_rna <- cor(t(expr[rna_avail, , drop = FALSE]), cov_rna, use = "pairwise.complete.obs")
R_rna_df <- data.frame(Gene = rownames(R_rna), R_rna, check.names = FALSE, row.names = NULL)
readr::write_csv(R_rna_df, file.path(analysis_dir, "Supp4c_mRNA_x_Covariates.csv"))
nx_check(R_rna_df, "Table_20_mRNA_x_Covariates.csv")

prm_idx  <- which(is.finite(as.numeric(prot$MS_NPTX2)))
ms_for   <- ifelse(rna_genes_20 == "PHF24", "MS_KIAA1045", paste0("MS_", rna_genes_20))
ms_ok    <- which(ms_for %in% colnames(prot))
prot_mat <- as.matrix(prot[prm_idx, ms_for[ms_ok], drop = FALSE]); colnames(prot_mat) <- rna_genes_20[ms_ok]
cov_prot <- cbind(PMI = as.numeric(meta$PMI)[prm_idx],
                  Age = suppressWarnings(as.numeric(meta$Age))[prm_idx],
                  Braak = as.numeric(meta$SS_B)[prm_idx], CERAD = as.numeric(meta$SS_C)[prm_idx])
cov_prot <- cov_prot[, nx_keep_cov(cov_prot), drop = FALSE]
R_prot <- cor(prot_mat, cov_prot, use = "pairwise.complete.obs")
R_prot_df <- data.frame(Protein = rownames(R_prot), R_prot, check.names = FALSE, row.names = NULL)
readr::write_csv(R_prot_df, file.path(analysis_dir, "Supp1c_Protein_x_Covariates.csv"))
nx_check(R_prot_df, "Table_21_Protein_x_Covariates.csv", key = "Protein")
cat("\nSupp 1c/4c covariate matrices done.\n")

pca_df <- data.frame(SampleID = meta$SampleID,
                     PC1 = as.numeric(meta$PC1), PC2 = as.numeric(meta$PC2),
                     LAMP2_expr = if ("LAMP2" %in% rownames(expr)) as.numeric(expr["LAMP2", ]) else NA_real_,
                     NPTX2_expr = as.numeric(expr["NPTX2", ]),
                     Has_PRM = is.finite(as.numeric(prot$MS_NPTX2)), stringsAsFactors = FALSE)
pca_df <- pca_df[is.finite(pca_df$PC1) & is.finite(pca_df$PC2), ]
readr::write_csv(pca_df, file.path(analysis_dir, "Supp2a_PCA_points.csv"))
cat("\nSupp 2a PCA point-data done (plot-only, nothing to verify).\n")

nx_pc_loading <- function(pc) {
  s <- setNames(as.numeric(meta[[pc]]), colnames(expr))
  r <- nx_cor_to_ref(expr, s, colnames(expr), rownames(expr))
  v <- setNames(r$r, r$Gene); sort(v[is.finite(v)], decreasing = TRUE)
}
pc1_ranks <- nx_pc_loading("PC1"); pc2_ranks <- nx_pc_loading("PC2")
readr::write_csv(data.frame(Gene = names(pc1_ranks), r = pc1_ranks), file.path(analysis_dir, "Supp2bc_PC1_loadings.csv"))
readr::write_csv(data.frame(Gene = names(pc2_ranks), r = pc2_ranks), file.path(analysis_dir, "Supp2bc_PC2_loadings.csv"))

gmt_dir <- path.expand("~/Desktop/GSEAforYuelin")
gmt_try <- file.path(gmt_dir, c("c5.go.cc.v2026.1.Hs.symbols.gmt",
                                "c5.go.cc.v2023.1.Hs.symbols.gmt",
                                "c5.go.v2023.1.Hs.symbols.gmt.txt"))
gmt_file <- gmt_try[file.exists(gmt_try)][1]
if (requireNamespace("fgsea", quietly = TRUE) && !is.na(gmt_file)) {
  gocc <- fgsea::gmtPathways(gmt_file)
  if (any(grepl("^GOBP_|^GOMF_", names(gocc)))) gocc <- gocc[grepl("^GOCC_", names(gocc))]
  set.seed(42)
  fg1 <- fgsea::fgsea(pathways = gocc, stats = pc1_ranks, minSize = 50, maxSize = 1000, nproc = 1)
  fg2 <- fgsea::fgsea(pathways = gocc, stats = pc2_ranks, minSize = 50, maxSize = 1000, nproc = 1)
  tidy_fg <- function(fg, pc) { d <- as.data.frame(fg)
  if ("leadingEdge" %in% names(d)) d$leadingEdge <- vapply(d$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  d$PC <- pc; d }
  readr::write_csv(rbind(tidy_fg(fg1, "PC1"), tidy_fg(fg2, "PC2")),
                   file.path(analysis_dir, "Supp2bc_fGSEA_GOCC.csv"))
  cat(sprintf("\nSupp 2b/2c fGSEA done (%d PC1 + %d PC2 pathways; seed=42).\n", nrow(fg1), nrow(fg2)))
} else {
  cat("\n[skip] Supp 2b/2c fGSEA: need fgsea + GO:CC GMT in ~/Desktop/GSEAforYuelin.\n",
      "       PC1/PC2 loading rank vectors were still saved.\n")
}

sex_lab <- ifelse(grepl("^f", as.character(meta$Sex), ignore.case = TRUE), "Female",
                  ifelse(grepl("^m", as.character(meta$Sex), ignore.case = TRUE), "Male", NA_character_))
sx <- data.frame(SampleID = meta$SampleID, Sex = sex_lab, Group = as.character(meta$Group_Current),
                 Age = suppressWarnings(as.numeric(meta$Age)), PMI = suppressWarnings(as.numeric(meta$PMI)),
                 RIN = suppressWarnings(as.numeric(meta$RIN)),
                 CERAD = as.numeric(meta$SS_C), Braak = as.numeric(meta$SS_B), stringsAsFactors = FALSE)
sx <- sx[sx$Sex %in% c("Female", "Male") & sx$Group %in% c("CN-Lo", "CN-Hi", "MCI", "AD"), ]
grp_lv  <- c("CN-Lo", "CN-Hi", "MCI", "AD")
sig_padj <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

nptx2_full <- setNames(as.numeric(expr["NPTX2", ]), colnames(expr))
nx_cor_dist <- function(samps, min_n = 10) {
  samps <- intersect(samps, colnames(expr))
  if (length(samps) < min_n) return(data.frame(Gene = character(0), r = numeric(0)))
  es <- expr[, samps, drop = FALSE]; ns <- nptx2_full[samps]
  r <- apply(es, 1, function(g) { ok <- is.finite(g) & is.finite(ns)
  if (sum(ok) < min_n) return(NA_real_)
  if (sd(g[ok]) == 0 || sd(ns[ok]) == 0) return(NA_real_); cor(g[ok], ns[ok]) })
  df <- data.frame(Gene = rownames(es), r = as.numeric(r), stringsAsFactors = FALSE)
  df[is.finite(df$r) & df$Gene != "NPTX2", ]
}
cor_rows <- list()
for (grp in c("Overall", grp_lv)) {
  if (grp == "Overall") { f <- sx$SampleID[sx$Sex == "Female"]; m <- sx$SampleID[sx$Sex == "Male"] }
  else { f <- sx$SampleID[sx$Group == grp & sx$Sex == "Female"]; m <- sx$SampleID[sx$Group == grp & sx$Sex == "Male"] }
  if (length(f) < 10 || length(m) < 10) next
  cor_rows[[grp]] <- rbind(cbind(nx_cor_dist(f), Sex = "Female", Group = grp),
                           cbind(nx_cor_dist(m), Sex = "Male",   Group = grp))
}
cor_all <- dplyr::bind_rows(cor_rows); cor_all$abs_r <- abs(cor_all$r)

readr::write_csv(cor_all[, c("Gene", "Group", "Sex", "r")],
                 file.path(analysis_dir, "Supp3_density_values.csv"))

n_grp <- sx %>% dplyr::group_by(Group) %>% dplyr::summarise(n_F = sum(Sex == "Female"), n_M = sum(Sex == "Male"), .groups = "drop")

stats_density <- cor_all %>% dplyr::filter(Group != "Overall") %>% dplyr::group_by(Group) %>%
  dplyr::summarise(median_F = median(r[Sex == "Female"], na.rm = TRUE),
                   median_M = median(r[Sex == "Male"],   na.rm = TRUE),
                   KS_D = unname(suppressWarnings(ks.test(r[Sex == "Female"], r[Sex == "Male"])$statistic)),
                   KS_p = suppressWarnings(ks.test(r[Sex == "Female"], r[Sex == "Male"])$p.value), .groups = "drop") %>%
  dplyr::left_join(n_grp, by = "Group") %>% dplyr::rename(Group_Current = Group)
readr::write_csv(stats_density, file.path(analysis_dir, "Supp3_Density_KS.csv"))
nx_check_keyed(stats_density, "Table_SexDiag_Density_KS.csv", keys = "Group_Current", harmonize = "Group_Current")

cov_long <- sx %>% dplyr::select(Sex, Group, Age, PMI, RIN) %>%
  tidyr::pivot_longer(c(Age, PMI, RIN), names_to = "Variable", values_to = "Value") %>% dplyr::filter(is.finite(Value))
cov_long <- dplyr::bind_rows(transform(cov_long, Group = "Overall"), cov_long)
stats_cov <- cov_long %>% dplyr::group_by(Variable, Group) %>%
  dplyr::summarise(n_F = sum(Sex == "Female"), n_M = sum(Sex == "Male"),
                   median_F = median(Value[Sex == "Female"], na.rm = TRUE), median_M = median(Value[Sex == "Male"], na.rm = TRUE),
                   Wilcox_p = if (sum(Sex == "Female") >= 3 && sum(Sex == "Male") >= 3)
                     suppressWarnings(wilcox.test(Value[Sex == "Female"], Value[Sex == "Male"])$p.value) else NA_real_,
                   .groups = "drop") %>%
  dplyr::group_by(Variable) %>% dplyr::mutate(Wilcox_padj = p.adjust(Wilcox_p, method = "BH")) %>% dplyr::ungroup() %>%
  dplyr::mutate(sig = sig_padj(Wilcox_padj))

nx_ord_stats <- function(varcol) {
  d <- data.frame(Sex = sx$Sex, Group = sx$Group, val = sx[[varcol]], stringsAsFactors = FALSE)
  d <- d[is.finite(d$val), ]; d <- rbind(transform(d, Group = "Overall"), d)
  d %>% dplyr::group_by(Group) %>%
    dplyr::summarise(n_F = sum(Sex == "Female"), n_M = sum(Sex == "Male"),
                     median_F = median(val[Sex == "Female"], na.rm = TRUE), median_M = median(val[Sex == "Male"], na.rm = TRUE),
                     Wilcox_p = if (sum(Sex == "Female") >= 3 && sum(Sex == "Male") >= 3)
                       suppressWarnings(wilcox.test(val[Sex == "Female"], val[Sex == "Male"])$p.value) else NA_real_,
                     .groups = "drop") %>%
    dplyr::mutate(Wilcox_padj = p.adjust(Wilcox_p, method = "BH"), sig = sig_padj(Wilcox_padj))
}
stats_cerad <- nx_ord_stats("CERAD"); stats_ssb <- nx_ord_stats("Braak")

stats_cor <- cor_all %>% dplyr::group_by(Group) %>%
  dplyr::summarise(n_genes_F = sum(Sex == "Female"), n_genes_M = sum(Sex == "Male"),
                   median_r_F = median(r[Sex == "Female"], na.rm = TRUE), median_r_M = median(r[Sex == "Male"], na.rm = TRUE),
                   median_absr_F = median(abs_r[Sex == "Female"], na.rm = TRUE), median_absr_M = median(abs_r[Sex == "Male"], na.rm = TRUE),
                   IQR_r_F = IQR(r[Sex == "Female"], na.rm = TRUE), IQR_r_M = IQR(r[Sex == "Male"], na.rm = TRUE),
                   SD_r_F = sd(r[Sex == "Female"], na.rm = TRUE), SD_r_M = sd(r[Sex == "Male"], na.rm = TRUE),
                   pct_F_strong_neg = mean(r[Sex == "Female"] < -0.3, na.rm = TRUE), pct_M_strong_neg = mean(r[Sex == "Male"] < -0.3, na.rm = TRUE),
                   pct_F_mild_neg = mean(r[Sex == "Female"] >= -0.3 & r[Sex == "Female"] < 0, na.rm = TRUE),
                   pct_M_mild_neg = mean(r[Sex == "Male"] >= -0.3 & r[Sex == "Male"] < 0, na.rm = TRUE),
                   pct_F_mild_pos = mean(r[Sex == "Female"] >= 0 & r[Sex == "Female"] <= 0.3, na.rm = TRUE),
                   pct_M_mild_pos = mean(r[Sex == "Male"] >= 0 & r[Sex == "Male"] <= 0.3, na.rm = TRUE),
                   pct_F_strong_pos = mean(r[Sex == "Female"] > 0.3, na.rm = TRUE), pct_M_strong_pos = mean(r[Sex == "Male"] > 0.3, na.rm = TRUE),
                   Wilcox_p_absr = suppressWarnings(wilcox.test(abs_r[Sex == "Female"], abs_r[Sex == "Male"])$p.value),
                   Wilcox_p_F_lt_M_absr = suppressWarnings(wilcox.test(abs_r[Sex == "Female"], abs_r[Sex == "Male"], alternative = "less")$p.value),
                   .groups = "drop") %>%
  dplyr::mutate(Wilcox_padj_absr = p.adjust(Wilcox_p_absr, method = "BH"),
                Wilcox_padj_F_lt_M_absr = p.adjust(Wilcox_p_F_lt_M_absr, method = "BH"),
                sig = sig_padj(Wilcox_padj_absr))

all_stats <- dplyr::bind_rows(
  stats_cov   %>% dplyr::mutate(Source = "Continuous (sample-level)"),
  stats_cerad %>% dplyr::mutate(Variable = "CERAD",      Source = "Ordinal (sample-level)"),
  stats_ssb   %>% dplyr::mutate(Variable = "SS_B_Score", Source = "Ordinal (sample-level)"),
  stats_cor   %>% dplyr::mutate(Variable = "abs_r",      Source = "Gene-level |r|") %>%
    dplyr::rename(n_F = n_genes_F, n_M = n_genes_M, median_F = median_absr_F, median_M = median_absr_M,
                  Wilcox_p = Wilcox_p_absr, Wilcox_padj = Wilcox_padj_absr)
) %>% dplyr::rename(Group_Current = Group)
readr::write_csv(all_stats, file.path(analysis_dir, "Supp3_Composite_Stats.csv"))
nx_check_keyed(all_stats, "Table_SexDiag_Composite_Stats.csv",
               keys = c("Variable", "Group_Current"), harmonize = "Group_Current")
cat("\nSupp 3 sex diagnostics done.\n")

bygroup_ks <- fig2b3c_density_ks %>% dplyr::transmute(
  Plot = "ByGroup_density",
  Modality = ifelse(Modality == "NPTX2_mRNA_ref", "Fig3c (mRNA ref)", "Fig2b (protein ref)"),
  Comparison = paste0(Group1, "_vs_", Group2), n1, n2, KS_D, KS_p, KS_padj)

sex_ks <- do.call(rbind, lapply(c("Overall", grp_lv), function(g) {
  d <- cor_all[cor_all$Group == g, ]
  x <- d$r[d$Sex == "Female"]; x <- x[is.finite(x)]
  y <- d$r[d$Sex == "Male"];   y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NULL)
  k <- suppressWarnings(ks.test(x, y))
  data.frame(Plot = "Sex_density", Modality = "gene-NPTX2 mRNA r",
             Comparison = paste0("Female_vs_Male @ ", g),
             n1 = length(x), n2 = length(y), KS_D = unname(k$statistic),
             KS_p = k$p.value, stringsAsFactors = FALSE)
}))
sex_ks$KS_padj <- p.adjust(sex_ks$KS_p, method = "BH")

density_ks_summary <- dplyr::bind_rows(bygroup_ks, sex_ks)
readr::write_csv(density_ks_summary, file.path(analysis_dir, "Density_KS_Summary.csv"))
cat(sprintf("Density KS summary written (%d by-group + %d sex comparisons).\n",
            nrow(bygroup_ks), nrow(sex_ks)))

gc <- data.frame(Group = as.character(meta$Group_Current),
                 Braak = suppressWarnings(as.numeric(meta$SS_B)),
                 CERAD = suppressWarnings(as.numeric(meta$SS_C)),
                 stringsAsFactors = FALSE)
gc <- gc[gc$Group %in% grp_lv, ]
gc$Group <- factor(gc$Group, levels = grp_lv)

braak_levels <- sort(unique(gc$Braak[is.finite(gc$Braak)]))
comp_rows <- list()
for (v in c("Braak", "CERAD")) for (g in grp_lv) {
  x <- gc[[v]][gc$Group == g]; x <- x[is.finite(x)]
  row <- data.frame(Variable = v, Group = g, n = length(x),
                    Median = median(x), IQR = IQR(x), Mean = mean(x), SD = sd(x),
                    stringsAsFactors = FALSE)
  for (lv in braak_levels)
    row[[paste0("pct_Braak_", lv)]] <- if (v == "Braak") round(mean(x == lv), 4) else NA_real_
  comp_rows[[paste(v, g)]] <- row
}
braak_cerad_composition <- bind_rows(comp_rows)
readr::write_csv(braak_cerad_composition, file.path(analysis_dir, "BraakCERAD_composition.csv"))

omni <- bind_rows(lapply(c("Braak", "CERAD"), function(v) {
  d <- gc[is.finite(gc[[v]]), ]; k <- kruskal.test(d[[v]], d$Group)
  data.frame(Variable = v, Test = "Kruskal-Wallis (across 4 groups)",
             statistic = unname(k$statistic), df = unname(k$parameter), p = k$p.value,
             stringsAsFactors = FALSE)
}))
bt <- table(gc$Group, gc$Braak); cs <- suppressWarnings(chisq.test(bt))
omni <- bind_rows(omni, data.frame(Variable = "Braak", Test = "Chi-square (Group x Braak level)",
                                   statistic = unname(cs$statistic), df = unname(cs$parameter), p = cs$p.value))
readr::write_csv(omni, file.path(analysis_dir, "BraakCERAD_omnibus.csv"))

nx_group_pairwise <- function(v) {
  out <- do.call(rbind, lapply(utils::combn(grp_lv, 2, simplify = FALSE), function(p) {
    x <- gc[[v]][gc$Group == p[1]]; x <- x[is.finite(x)]
    y <- gc[[v]][gc$Group == p[2]]; y <- y[is.finite(y)]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    ks <- suppressWarnings(ks.test(x, y)); wx <- suppressWarnings(wilcox.test(x, y))
    data.frame(Variable = v, Group1 = p[1], Group2 = p[2], n1 = length(x), n2 = length(y),
               median1 = median(x), median2 = median(y),
               KS_D = unname(ks$statistic), KS_p = ks$p.value, Wilcox_p = wx$p.value,
               stringsAsFactors = FALSE)
  }))
  out$KS_padj <- p.adjust(out$KS_p, method = "BH")
  out$Wilcox_padj <- p.adjust(out$Wilcox_p, method = "BH")
  out
}
braak_cerad_pairwise <- bind_rows(nx_group_pairwise("Braak"), nx_group_pairwise("CERAD"))
readr::write_csv(braak_cerad_pairwise, file.path(analysis_dir, "BraakCERAD_pairwise.csv"))
cat(sprintf("\nBraak/CERAD group composition done (KW + Chi-sq omnibus + %d pairwise rows).\n",
            nrow(braak_cerad_pairwise)))

cat("\n\n===================== VERIFICATION SUMMARY =====================\n")
cat(paste(.nx_log, collapse = "\n"), "\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  %d IDENTICAL  |  %d DIFFERS  |  %d new (no existing file)\n",
            sum(grepl("IDENTICAL", .nx_log)), sum(grepl("DIFFERS", .nx_log)),
            sum(grepl("no existing file", .nx_log))))
cat("================================================================\n")
