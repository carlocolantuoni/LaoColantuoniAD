suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(dplyr)
  library(ggpubr)
  library(cowplot)
  library(svglite)
})

DATA_FOLDER    <- "exported_data_fixed"
RESULTS_FOLDER <- "analysis_output_fixed"
FIG_FOLDER     <- "final_figures_fixed"
exp_dir <- file.path(getwd(), DATA_FOLDER)
in_dir  <- file.path(getwd(), RESULTS_FOLDER)
fig_dir <- file.path(getwd(), FIG_FOLDER)
png_dir <- file.path(fig_dir, "png")
svg_dir <- file.path(fig_dir, "svg")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(svg_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("RAW    : %s\nTABLES : %s\nFIGURES: %s  (png/ + svg/)\n", exp_dir, in_dir, fig_dir))

calc_y_positions <- function(y, n_brk) { y <- y[is.finite(y)]; if (!length(y)) return(rep(NA_real_, n_brk)); span <- diff(range(y)); if (span <= 0) span <- 1; base <- max(y) + 0.01 * span; base + 0.10 * span * (0:(n_brk - 1)) }

create_panel_compact <- function(data, y_col, y_lab, group_col, colors, comparisons, y_limits, show_y_title = TRUE) {
  y_vals <- data[[y_col]]
  y_visible_max <- y_limits[2]
  y_visible_min <- y_limits[1]
  span <- y_visible_max - y_visible_min
  base <- y_visible_max - 0.30 * span
  y_pos <- base + 0.08 * span * (0:(length(comparisons) - 1))
  data$Dot_Color <- ifelse(!is.na(data$Age_Num) & data$Age_Num > FIG1B3B$age_threshold, FIG1B3B$age_dot_color, "grey40")
  ggplot(data, aes(x = .data[[group_col]], y = .data[[y_col]], fill = .data[[group_col]])) +
    geom_boxplot(width = FIG1B3B$box_width, outlier.shape = NA, color = "grey20", lwd = 0.35) +
    geom_jitter(aes(color = Dot_Color), width = 0.12, size = FIG1B3B$jitter_size, alpha = FIG1B3B$jitter_alpha, shape = 16) +
    scale_color_identity() +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = function(b) { nn <- table(data[[group_col]]); sprintf("%s\n(n=%d)", b, as.integer(nn[b])) }) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = if(show_y_title) y_lab else NULL) +
    theme_classic(base_size = FIG1B3B$base_font, base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 3, 4, 3),
      axis.text.x = element_text(size = FIG1B3B$axis_x_font, margin = margin(t = 2)),
      axis.text.y = element_text(size = FIG1B3B$axis_y_font),
      axis.title.y = if(show_y_title) element_text(size = FIG1B3B$ytitle_font) else element_blank()
    ) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", p.adjust.method = "BH",
                       label = "p.adj.format", hide.ns = FALSE, vjust = -0.3,
                       size = FIG1B3B$pval_font, tip.length = 0.015, label.y = y_pos)
}

make_nptx2_boxplot_22 <- function(df_left, df_right, y_col, y_label, out_file, title_text, y_cap = NULL) {
  df_left  <- df_left  %>% filter(is.finite(.data[[y_col]]))
  df_right <- df_right %>% filter(is.finite(.data[[y_col]]))
  if(nrow(df_left) < 3 || nrow(df_right) < 3) {
    cat(sprintf("  Skipping %s: insufficient data (left=%d, right=%d)\n", title_text, nrow(df_left), nrow(df_right)))
    return(NULL)
  }
  y_all <- c(df_left[[y_col]], df_right[[y_col]])
  global_max_left  <- max(calc_y_positions(df_left[[y_col]],  length(comps_main_22)))
  global_max_right <- max(calc_y_positions(df_right[[y_col]], length(comps_path_22)))
  global_max <- max(global_max_left, global_max_right)
  upper <- if(!is.null(y_cap)) y_cap else global_max + diff(range(y_all)) * 0.05
  rng <- c(min(y_all, na.rm = TRUE), upper)
  p_left  <- create_panel_compact(df_left,  y_col, y_label, "Group", colors_main, comps_main_22, rng, TRUE)  + .transparent_theme
  p_right <- create_panel_compact(df_right, y_col, y_label, "Group", colors_path, comps_path_22, rng, FALSE) + .transparent_theme
  panels <- plot_grid(p_left, p_right, ncol = 2, rel_widths = FIG1B3B$rel_widths, align = "h")

  age_legend <- ggplot() +
    annotate("point", x = 1.00, y = 1, color = FIG1B3B$age_dot_color, size = 1.8) +
    annotate("text",  x = 1.07, y = 1, label = FIG1B3B$age_label, hjust = 0, size = FIG1B3B$age_text_size,
             color = "grey20", family = "Arial") +
    scale_x_continuous(limits = c(0.55, 3.5)) + scale_y_continuous(limits = c(0.5, 1.5)) +
    theme_void() + .transparent_theme
  final <- plot_grid(panels, age_legend, ncol = 1, rel_heights = c(1, 0.09))
  bn <- sub("\\.svg$", "", basename(out_file))
  svglite(file.path(svg_dir, paste0(bn, ".svg")), width = FIG1B3B$out_w, height = FIG1B3B$out_h, system_fonts = list(sans = "Arial"), bg = "transparent")
  print(final); dev.off()
  ggplot2::ggsave(file.path(png_dir, paste0(bn, ".png")), final, width = FIG1B3B$out_w, height = FIG1B3B$out_h, dpi = 300, bg = "transparent")
  cat(sprintf("  Saved: %s.{png,svg}\n", bn))
}

viz_read <- function(name) {
  f <- file.path(in_dir, name); if (!file.exists(f)) stop("missing table: ", name)
  as.data.frame(readr::read_csv(f, show_col_types = FALSE))
}

.transparent_theme <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))
viz_save_both <- function(plot, base, w, h, dpi = 300) {
  if (inherits(plot, "ggplot")) plot <- plot + .transparent_theme
  ggplot2::ggsave(file.path(png_dir, paste0(base, ".png")), plot, width = w, height = h,
                  dpi = dpi, units = "in", bg = "transparent")
  svglite::svglite(file.path(svg_dir, paste0(base, ".svg")), width = w, height = h,
                   system_fonts = list(sans = "Arial"), bg = "transparent")
  print(plot); dev.off()
  cat(sprintf("  Saved: %s.{png,svg}\n", base))
}
meta_aligned <- as.data.frame(readr::read_csv(file.path(exp_dir, "metadata_master.csv"), show_col_types = FALSE))

meta_aligned$Group_Current <- dplyr::recode(as.character(meta_aligned$Group_Current),
                                            "LoPath" = "CN-Lo", "HiPath" = "CN-Hi")
prot         <- as.data.frame(readr::read_csv(file.path(exp_dir, "protein_abundance.csv"), show_col_types = FALSE))
meta_aligned$MS_NPTX2 <- as.numeric(prot$MS_NPTX2[match(meta_aligned$SampleID, prot$SampleID)])
.expr_df <- as.data.frame(readr::read_csv(file.path(exp_dir, "expression_rna.csv.gz"), show_col_types = FALSE))
exprs_aligned <- as.matrix(.expr_df[, -1]); rownames(exprs_aligned) <- .expr_df[[1]]

FIG1B3B <- list(

  out_w = 5.5, out_h = 4.2, rel_widths = c(1, 0.9),

  base_font = 16, axis_x_font = 14, axis_y_font = 14, ytitle_font = 17, pval_font = 5,

  box_width = 0.45, jitter_size = 1.2, jitter_alpha = 0.55,

  ycap_rna  = 11, ycap_prot = NULL,

  colors_main = c("Control" = "cornflowerblue", "MCI" = "orange", "AD" = "brown3"),
  colors_path = c("CN-Lo" = "lightskyblue", "CN-Hi" = "royalblue", "YoungCon" = "seagreen3"),

  age_threshold = 95, age_dot_color = "red", age_label = "Age > 95 yr", age_text_size = 6,

  comps_main  = list(c("Control","MCI"), c("Control","AD"), c("MCI","AD")),
  comps_path  = list(c("CN-Lo","CN-Hi"), c("CN-Lo","YoungCon"), c("CN-Hi","YoungCon")))

colors_main   <- FIG1B3B$colors_main
colors_path   <- FIG1B3B$colors_path
comps_main_22 <- FIG1B3B$comps_main
comps_path_22 <- FIG1B3B$comps_path

df_main_22 <- meta_aligned %>%
  mutate(Diag_Clean = trimws(as.character(Diagnosis_COMP))) %>%
  filter(Diag_Clean %in% c("Control", "MCI", "AD")) %>%
  mutate(Group = factor(Diag_Clean, levels = c("Control", "MCI", "AD")))
df_path_22 <- meta_aligned %>%
  mutate(Diag_Clean = trimws(as.character(Diagnosis_COMP)),
         Group = case_when(Diag_Clean == "YoungCon" ~ "YoungCon",
                           Group_Current == "CN-Lo"  ~ "CN-Lo",
                           Group_Current == "CN-Hi"  ~ "CN-Hi",
                           TRUE ~ NA_character_)) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("CN-Lo", "CN-Hi", "YoungCon")))
age_col <- grep("^Age", colnames(meta_aligned), ignore.case = TRUE, value = TRUE)[1]
df_main_22$Age_Num <- suppressWarnings(as.numeric(meta_aligned[[age_col]][match(df_main_22$SampleID, meta_aligned$SampleID)]))
df_path_22$Age_Num <- suppressWarnings(as.numeric(meta_aligned[[age_col]][match(df_path_22$SampleID, meta_aligned$SampleID)]))

df_main_22$NPTX2_RNA <- as.numeric(exprs_aligned["NPTX2", df_main_22$SampleID])
df_path_22$NPTX2_RNA <- as.numeric(exprs_aligned["NPTX2", df_path_22$SampleID])
df_main_22$NPTX2_MS  <- meta_aligned$MS_NPTX2[match(df_main_22$SampleID, meta_aligned$SampleID)]
df_path_22$NPTX2_MS  <- meta_aligned$MS_NPTX2[match(df_path_22$SampleID, meta_aligned$SampleID)]

make_nptx2_boxplot_22(df_main_22, df_path_22, "NPTX2_RNA", "NPTX2 mRNA (Log2)",
                      file.path(fig_dir, "Fig3b_NPTX2_mRNA.svg"), "NPTX2 mRNA", y_cap = FIG1B3B$ycap_rna)
make_nptx2_boxplot_22(df_main_22, df_path_22, "NPTX2_MS", "NPTX2 Protein (PRM-MS)",
                      file.path(fig_dir, "Fig1b_NPTX2_protein.svg"), "NPTX2 Protein", y_cap = FIG1B3B$ycap_prot)

cat("\nFigure 1b / 3b ported.\n")

suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(RColorBrewer) })

SET_1cd <- list(

  tile_w_cm = 1.5, tile_h_cm = 0.70,

  dev_w = 6.0, dev_h = 8.5,

  cell_value_font = 8, row_name_font = 10, col_name_font = 11,
  row_title_font = 11, col_title_font = 13, col_names_rot = 45,

  corr_white_above = 0.55, expr_white_frac = 0.65, mrna_white_frac = 0.65,

  title_cor  = "Protein vs NPTX2 Correlation Patterns",
  title_exp  = "Protein Expression (Min-Max Scaled)",
  title_mrna = "mRNA Expression (PRM-MS samples only, n=%d)",

  cor_palette  = rev(RColorBrewer::brewer.pal(11, "RdBu")),
  expr_palette = c("#FFFFB2","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#B10026"))

df_17 <- viz_read("Fig1cd_trajectory_metrics.csv")

grp_n    <- sapply(c("CN-Lo", "CN-Hi", "MCI", "AD"), function(g)
  sum(meta_aligned$Group_Current == g & is.finite(meta_aligned$MS_NPTX2)))
col_lab4 <- sprintf("%s\n(n=%d)", c("CN-Lo", "CN-Hi", "MCI", "AD"), grp_n)

mat_cor  <- as.matrix(df_17 %>% dplyr::select(starts_with("Cor_") & !contains("FisherZ") & !contains("DeltaR")))
mat_cor  <- mat_cor[, 3:6, drop = FALSE]
mat_exp  <- as.matrix(df_17 %>% dplyr::select(starts_with("Exp_") & !contains("Wilcox")))
mat_mrna <- as.matrix(df_17 %>% dplyr::select(starts_with("mRNA_")))
rownames(mat_cor) <- df_17$Gene; rownames(mat_exp) <- df_17$Gene; rownames(mat_mrna) <- df_17$Gene
colnames(mat_cor) <- col_lab4; colnames(mat_exp) <- col_lab4; colnames(mat_mrna) <- col_lab4

traj_order <- c("A", "B", "Weak Corr", "Ref")
grp_lvls   <- c(intersect(traj_order, unique(df_17$Corr_Traj)),
                setdiff(unique(df_17$Corr_Traj), traj_order))
row_grp    <- factor(df_17$Corr_Traj, levels = grp_lvls)

lab_col <- "black"

cell_w <- unit(ncol(mat_cor) * SET_1cd$tile_w_cm, "cm")
cell_h <- unit(nrow(mat_cor) * SET_1cd$tile_h_cm, "cm")
dev_w  <- SET_1cd$dev_w; dev_h <- SET_1cd$dev_h

cor_col_fun <- colorRamp2(seq(-1, 1, length = 50), colorRampPalette(SET_1cd$cor_palette)(50))
ylOrRd_cols <- SET_1cd$expr_palette
exp_data_min <- min(mat_exp, na.rm = TRUE); exp_data_max <- max(mat_exp, na.rm = TRUE)
exp_col_fun  <- colorRamp2(seq(exp_data_min, exp_data_max, length.out = 7), ylOrRd_cols)
exp_txt_thresh <- exp_data_min + SET_1cd$expr_white_frac * (exp_data_max - exp_data_min)
mrna_data_min <- min(mat_mrna, na.rm = TRUE); mrna_data_max <- max(mat_mrna, na.rm = TRUE)
mrna_col_fun  <- colorRamp2(seq(mrna_data_min, mrna_data_max, length.out = 7), ylOrRd_cols)
mrna_txt_thresh <- mrna_data_min + SET_1cd$mrna_white_frac * (mrna_data_max - mrna_data_min)

viz_save_ht <- function(ht, base, w, h, column_title, ttl_font = 13) {
  do_draw <- function() draw(ht, column_title = column_title,
                             column_title_gp = gpar(fontsize = ttl_font, fontfamily = "Arial", fontface = "bold"),
                             padding = unit(c(5, 5, 15, 20), "mm"))
  svglite::svglite(file.path(svg_dir, paste0(base, ".svg")), width = w, height = h, system_fonts = list(sans = "Arial"), bg = "transparent")
  do_draw(); dev.off()
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file.path(png_dir, paste0(base, ".png")), width = w, height = h, units = "in", res = 300, background = "transparent")
  else grDevices::png(file.path(png_dir, paste0(base, ".png")), width = w, height = h, units = "in", res = 300, bg = "transparent")
  do_draw(); dev.off()
  cat(sprintf("  Saved: %s.{png,svg}\n", base))
}

make_traj_ht <- function(mat, name, col_fun, white_above_val) {
  Heatmap(
    mat, name = name, col = col_fun,
    width = cell_w, height = cell_h,
    row_split = row_grp, cluster_row_slices = FALSE, cluster_rows = FALSE,
    show_row_names = TRUE, row_names_side = "left",
    row_names_gp = gpar(col = lab_col, fontsize = SET_1cd$row_name_font, fontface = "bold", fontfamily = "Arial"),
    row_title_rot = 90, row_title_gp = gpar(fontsize = SET_1cd$row_title_font, fontface = "bold", fontfamily = "Arial"),
    cluster_columns = FALSE, column_title = NULL,
    column_names_gp = gpar(fontsize = SET_1cd$col_name_font, fontfamily = "Arial"), column_names_rot = SET_1cd$col_names_rot,
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat[i, j]
      if (!is.na(val)) {
        txt_col <- ifelse(val >= white_above_val, "white", "black")
        grid::grid.text(sprintf("%.2f", val), x, y, gp = grid::gpar(fontsize = SET_1cd$cell_value_font, col = txt_col, fontfamily = "Arial"))
      }
    })
}

ht_cor <- Heatmap(
  mat_cor, name = "Pearson r", col = cor_col_fun, width = cell_w, height = cell_h,
  row_split = row_grp, cluster_row_slices = FALSE, cluster_rows = FALSE,
  show_row_names = TRUE, row_names_side = "left",
  row_names_gp = gpar(col = lab_col, fontsize = SET_1cd$row_name_font, fontface = "bold", fontfamily = "Arial"),
  row_title_rot = 90, row_title_gp = gpar(fontsize = SET_1cd$row_title_font, fontface = "bold", fontfamily = "Arial"),
  cluster_columns = FALSE, column_title = NULL,
  column_names_gp = gpar(fontsize = SET_1cd$col_name_font, fontfamily = "Arial"), column_names_rot = SET_1cd$col_names_rot,
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- mat_cor[i, j]
    if (!is.na(val)) {
      txt_col <- ifelse(abs(val) >= SET_1cd$corr_white_above, "white", "black")
      grid::grid.text(sprintf("%.2f", val), x, y, gp = grid::gpar(fontsize = SET_1cd$cell_value_font, col = txt_col, fontfamily = "Arial"))
    }
  },
  heatmap_legend_param = list(title = "Pearson r", title_position = "leftcenter-rot"))
viz_save_ht(ht_cor, "Fig1c_Pattern_Correlation", dev_w, dev_h, SET_1cd$title_cor, ttl_font = SET_1cd$col_title_font)

ht_exp <- make_traj_ht(mat_exp, "Scaled\nProtein", exp_col_fun, exp_txt_thresh)
viz_save_ht(ht_exp, "Fig1d_Pattern_Expression_protein", dev_w, dev_h, SET_1cd$title_exp, ttl_font = SET_1cd$col_title_font)

prm_samps_17 <- meta_aligned$SampleID[is.finite(meta_aligned$MS_NPTX2)]
n_prm_mrna   <- length(intersect(prm_samps_17, colnames(exprs_aligned)))
ht_mrna <- make_traj_ht(mat_mrna, "Scaled\nmRNA", mrna_col_fun, mrna_txt_thresh)
viz_save_ht(ht_mrna, "Fig1d_Pattern_mRNA_PRMsamples", dev_w, dev_h,
            sprintf(SET_1cd$title_mrna, n_prm_mrna), ttl_font = SET_1cd$col_title_font)

cat("\nFigure 1c / 1d ported.\n")

cfg_traj <- list(
  font_family = "Arial", strong_value = 0.4, weak_cap = 0.25, arrow_extend = 0.20,
  pos_color = "#A32D2D", neg_color = "#185FA5", neg_alpha = 0.5, arrow_color = "black",
  line_thickness = 1.4, point_size = 3.6, arrow_size = 0.18, cap_width = 0.06,
  font_title = 18, font_axis = 17, title_angle = 270, show_subtitle = FALSE,
  y_axis_label = "Correlation with NPTX2 (r)", x_labels = c("CN-Lo", "CN-Hi", "AD/MCI"),

  plot_width = 4.6, plot_height = 11.5, bottom_panel_extra = 0.20)

make_point <- function(x, type, color) {
  s <- cfg_traj$strong_value; w <- cfg_traj$weak_cap; e <- cfg_traj$arrow_extend
  if (type == "strong+") {
    list(point = data.frame(x = x, y =  s, color = color),
         arrows = data.frame(x = x, xend = x, y =  s, yend =  s + e, color = color), caps = NULL)
  } else if (type == "strong-") {
    list(point = data.frame(x = x, y = -s, color = color),
         arrows = data.frame(x = x, xend = x, y = -s, yend = -s - e, color = color), caps = NULL)
  } else if (type == "weak") {
    list(point = data.frame(x = x, y = 0, color = color),
         arrows = data.frame(x = c(x, x), xend = c(x, x), y = c(0, 0), yend = c(w, -w), color = color),
         caps = data.frame(x = c(x - cfg_traj$cap_width, x - cfg_traj$cap_width),
                           xend = c(x + cfg_traj$cap_width, x + cfg_traj$cap_width),
                           y = c(w, -w), yend = c(w, -w), color = color))
  } else stop("Unknown point type: ", type)
}

class_defs <- list(
  list(id = "01_preserved", title = "CN-Hi-preserved", color = "#009E73",
       subtitle = "Strong in CN-Lo and CN-Hi;\nAD/MCI either weakens (case A) or stays strong (case B)",
       blue_seq = c("strong+", "strong+"), red_seq = c("strong-", "strong-"),
       branches = list(list(case = "A", linetype = "solid",  blue_end_type = "weak",    red_end_type = "weak"),
                       list(case = "B", linetype = "dashed", blue_end_type = "strong+", red_end_type = "strong-"))),
  list(id = "02_recruited", title = "CN-Hi-recruited", color = "#E69F00",
       subtitle = "Coupling appears only in CN-Hi",
       blue_seq = c("weak", "strong+", "weak"), red_seq = c("weak", "strong-", "weak"), branches = NULL),
  list(id = "03_suppressed", title = "CN-Hi-suppressed", color = "#0072B2",
       subtitle = "Coupling drops out only in CN-Hi",
       blue_seq = c("strong+", "weak", "strong+"), red_seq = c("strong-", "weak", "strong-"), branches = NULL),
  list(id = "04_reversed", title = "CN-Hi-reversed", color = "#D55E00",
       subtitle = "Sign flips at CN-Hi;\nAD/MCI either weakens (case A) or re-flips (case B)",
       blue_seq = c("strong+", "strong-"), red_seq = c("strong-", "strong+"),
       branches = list(list(case = "A", linetype = "solid",  blue_end_type = "weak",    red_end_type = "weak"),
                       list(case = "B", linetype = "dashed", blue_end_type = "strong+", red_end_type = "strong-"))),
  list(id = "05_pathology_disrupted", title = "Pathology-disrupted", color = "#CC79A7",
       subtitle = "Coupling lost from CN-Hi onward",
       blue_seq = c("strong+", "weak", "weak"), red_seq = c("strong-", "weak", "weak"), branches = NULL))

build_trajectory_schematic <- function(spec, show_x = FALSE) {
  build_line <- function(seq, color) {
    pts <- lapply(seq_along(seq), function(i) make_point(i, seq[i], color))
    list(points = bind_rows(lapply(pts, function(p) p$point)),
         arrows = bind_rows(lapply(pts, function(p) p$arrows)),
         caps   = bind_rows(lapply(pts, function(p) p$caps)))
  }
  line_segs <- function(pts, linetype = "solid") {
    if (nrow(pts) < 2) return(data.frame())
    data.frame(x = pts$x[-nrow(pts)], y = pts$y[-nrow(pts)], xend = pts$x[-1], yend = pts$y[-1],
               color = pts$color[-1], linetype = linetype)
  }
  blue <- build_line(spec$blue_seq, cfg_traj$pos_color)
  red  <- build_line(spec$red_seq,  cfg_traj$neg_color)
  segs_blue <- line_segs(blue$points, "solid"); segs_red <- line_segs(red$points, "solid")
  branch_points_blue <- data.frame(); branch_points_red <- data.frame()
  branch_arrows_blue <- data.frame(); branch_arrows_red <- data.frame()
  branch_caps_blue   <- data.frame(); branch_caps_red   <- data.frame()
  branch_segs_blue   <- data.frame(); branch_segs_red   <- data.frame(); branch_labels <- data.frame()
  if (!is.null(spec$branches)) {
    hi_blue_y <- blue$points$y[nrow(blue$points)]; hi_red_y <- red$points$y[nrow(red$points)]; ad_x <- 3
    for (b in spec$branches) {
      pt_blue <- make_point(ad_x, b$blue_end_type, cfg_traj$pos_color)
      pt_red  <- make_point(ad_x, b$red_end_type,  cfg_traj$neg_color)
      branch_points_blue <- bind_rows(branch_points_blue, pt_blue$point)
      branch_points_red  <- bind_rows(branch_points_red,  pt_red$point)
      if (!is.null(pt_blue$arrows)) branch_arrows_blue <- bind_rows(branch_arrows_blue, pt_blue$arrows)
      if (!is.null(pt_red$arrows))  branch_arrows_red  <- bind_rows(branch_arrows_red,  pt_red$arrows)
      if (!is.null(pt_blue$caps))   branch_caps_blue   <- bind_rows(branch_caps_blue,   pt_blue$caps)
      if (!is.null(pt_red$caps))    branch_caps_red    <- bind_rows(branch_caps_red,    pt_red$caps)
      branch_segs_blue <- bind_rows(branch_segs_blue, data.frame(x = 2, y = hi_blue_y, xend = ad_x, yend = pt_blue$point$y, color = cfg_traj$pos_color, linetype = b$linetype))
      branch_segs_red  <- bind_rows(branch_segs_red,  data.frame(x = 2, y = hi_red_y,  xend = ad_x, yend = pt_red$point$y,  color = cfg_traj$neg_color, linetype = b$linetype))
      branch_labels <- bind_rows(branch_labels, data.frame(x = ad_x + 0.08, y = pt_blue$point$y, label = sprintf("case %s", b$case), linetype = b$linetype))
    }
  }
  add_arrows <- function(p, df) {
    if (nrow(df) > 0) p <- p + geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend),
                                            color = cfg_traj$arrow_color, linewidth = 0.9,
                                            arrow = arrow(length = unit(cfg_traj$arrow_size, "cm"), type = "closed"))
    p
  }
  p <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5)
  if (nrow(segs_blue) > 0) p <- p + geom_segment(data = segs_blue, aes(x = x, y = y, xend = xend, yend = yend), color = cfg_traj$pos_color, linewidth = cfg_traj$line_thickness)
  if (nrow(segs_red)  > 0) p <- p + geom_segment(data = segs_red,  aes(x = x, y = y, xend = xend, yend = yend), color = cfg_traj$neg_color, alpha = cfg_traj$neg_alpha, linewidth = cfg_traj$line_thickness)
  if (nrow(branch_segs_blue) > 0) p <- p + geom_segment(data = branch_segs_blue, aes(x = x, y = y, xend = xend, yend = yend, linetype = I(linetype)), color = cfg_traj$pos_color, linewidth = cfg_traj$line_thickness)
  if (nrow(branch_segs_red)  > 0) p <- p + geom_segment(data = branch_segs_red,  aes(x = x, y = y, xend = xend, yend = yend, linetype = I(linetype)), color = cfg_traj$neg_color, alpha = cfg_traj$neg_alpha, linewidth = cfg_traj$line_thickness)
  p <- add_arrows(p, blue$arrows); p <- add_arrows(p, red$arrows)
  p <- add_arrows(p, branch_arrows_blue); p <- add_arrows(p, branch_arrows_red)
  for (df in list(blue$caps, red$caps, branch_caps_blue, branch_caps_red)) {
    if (!is.null(df) && nrow(df) > 0) p <- p + geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), color = cfg_traj$arrow_color, linewidth = 0.9, lineend = "round")
  }
  p <- p + geom_point(data = blue$points, aes(x = x, y = y), color = cfg_traj$pos_color, size = cfg_traj$point_size) +
    geom_point(data = red$points,  aes(x = x, y = y), color = cfg_traj$neg_color, alpha = cfg_traj$neg_alpha, size = cfg_traj$point_size)
  if (nrow(branch_points_blue) > 0) p <- p + geom_point(data = branch_points_blue, aes(x = x, y = y), color = cfg_traj$pos_color, size = cfg_traj$point_size)
  if (nrow(branch_points_red)  > 0) p <- p + geom_point(data = branch_points_red,  aes(x = x, y = y), color = cfg_traj$neg_color, alpha = cfg_traj$neg_alpha, size = cfg_traj$point_size)
  if (nrow(branch_labels) > 0) p <- p + geom_text(data = branch_labels, aes(x = x, y = y, label = label), hjust = 0, size = 4.0, color = "grey25", fontface = "bold", family = cfg_traj$font_family)
  x_upper <- if (!is.null(spec$branches)) 3.75 else 3.5
  p <- p +
    scale_x_continuous(breaks = 1:3, labels = cfg_traj$x_labels, limits = c(0.5, x_upper), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-0.85, 0.85), expand = c(0, 0),
                       sec.axis = dup_axis(name = spec$title, breaks = NULL)) +
    labs(subtitle = if (isTRUE(cfg_traj$show_subtitle)) spec$subtitle else NULL, x = NULL, y = NULL) +
    theme_classic(base_family = cfg_traj$font_family) +
    theme(plot.subtitle = element_text(size = cfg_traj$font_axis, color = "grey30", lineheight = 1.15, margin = margin(t = 2, b = 6)),
          axis.title.y.right = element_text(size = cfg_traj$font_title, face = "bold", angle = cfg_traj$title_angle, vjust = 1, margin = margin(l = 6), color = spec$color),
          axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y = element_text(size = cfg_traj$font_axis, color = "black"),
          axis.text.x = if (show_x) element_text(size = cfg_traj$font_axis, face = "bold", color = "black") else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey95", linetype = "dotted"),
          plot.margin = margin(2, 4, 2, 4)) +
    .transparent_theme
  p
}

n_cls    <- length(class_defs)
panels_2a <- lapply(seq_len(n_cls), function(i) build_trajectory_schematic(class_defs[[i]], show_x = (i == n_cls)))
rel_h_2a <- c(rep(1, n_cls - 1), 1 + cfg_traj$bottom_panel_extra)
stacked_2a <- cowplot::plot_grid(plotlist = panels_2a, ncol = 1, align = "v", axis = "lr", rel_heights = rel_h_2a)

fig2a <- cowplot::ggdraw() +
  cowplot::draw_plot(stacked_2a, x = 0.045, y = 0, width = 0.955, height = 1) +
  cowplot::draw_label(cfg_traj$y_axis_label, x = 0.018, y = 0.5, angle = 90,
                      fontface = "bold", size = cfg_traj$font_axis, fontfamily = cfg_traj$font_family)
viz_save_both(fig2a, "Fig2a_Pattern_Classes", cfg_traj$plot_width, cfg_traj$plot_height)
cat("\nFigure 2a ported (5 classes, vertically merged, shared x-axis).\n")

suppressPackageStartupMessages({ library(cowplot) })

SET_2b3c <- list(
  out_w = 5.5, out_h = 9.5, rel_heights = c(1, 0.75),
  base_font = 14, title_font = 15, subtitle_font = 11, axis_title_font = 13, axis_text_font = 12,
  density_aspect = 1,
  label_cutoff = 0.5, label_font = 3.4,
  ks_text_font = 4.0, ks_row_font = 4.5,
  ks_low = "#f7f4f9", ks_high = "#6a51a3",
  group_colors = c("CN-Lo" = "cornflowerblue", "CN-Hi" = "darkorchid4", "MCI" = "orange", "AD" = "brown3"),

  band_colors = c("Strong" = "#4a8c5c", "Moderate" = "#8aaa8c", "Weak" = "#b0b0b0"))

dens_group_colors <- SET_2b3c$group_colors

build_density_panel <- function(df_overall, df_groups, title_txt, subtitle_txt, label_cutoff = SET_2b3c$label_cutoff, free_range = FALSE, label_note = NULL) {
  df_overall <- df_overall %>% filter(is.finite(r)) %>%
    mutate(Category = case_when(abs(r) >= 0.5 ~ "Strong (|r|\u22650.5)",
                                abs(r) >= 0.3 ~ "Moderate (0.3\u2264|r|<0.5)",
                                TRUE          ~ "Weak (|r|<0.3)"),
           Category = factor(Category, levels = c("Strong (|r|\u22650.5)", "Moderate (0.3\u2264|r|<0.5)", "Weak (|r|<0.3)")))
  cat_colors <- c("Strong (|r|\u22650.5)" = unname(SET_2b3c$band_colors["Strong"]), "Moderate (0.3\u2264|r|<0.5)" = unname(SET_2b3c$band_colors["Moderate"]), "Weak (|r|<0.3)" = unname(SET_2b3c$band_colors["Weak"]))
  x_lo <- min(df_overall$r, na.rm = TRUE); x_hi <- max(df_overall$r, na.rm = TRUE)
  x_lim <- c(floor((x_lo - 0.05) * 20) / 20, ceiling((x_hi + 0.05) * 20) / 20)
  dens_all <- density(df_overall$r, n = 1024, from = x_lim[1], to = x_lim[2])
  dens_df <- data.frame(x = dens_all$x, y = dens_all$y) %>%
    mutate(Category = case_when(abs(x) >= 0.5 ~ "Strong (|r|\u22650.5)",
                                abs(x) >= 0.3 ~ "Moderate (0.3\u2264|r|<0.5)",
                                TRUE          ~ "Weak (|r|<0.3)"),
           Category = factor(Category, levels = c("Strong (|r|\u22650.5)", "Moderate (0.3\u2264|r|<0.5)", "Weak (|r|<0.3)")))
  df_groups <- df_groups %>% filter(is.finite(r), Group %in% names(dens_group_colors))
  df_groups$Group <- factor(df_groups$Group, levels = names(dens_group_colors))
  df_labels <- df_overall %>% filter(abs(r) >= label_cutoff) %>%
    mutate(Label = sprintf("%s (r=%.2f)", Gene, r)) %>% arrange(r)
  if (nrow(df_labels) > 0) {
    if (free_range) { pk <- max(dens_df$y, na.rm = TRUE)
    df_labels$y_pos <- pk * 1.05 + (seq_len(nrow(df_labels)) - 1) * pk * 0.07 }
    else df_labels$y_pos <- 0.8 + (seq_len(nrow(df_labels)) - 1) * 0.35
  } else df_labels <- tibble::tibble(Gene = character(), r = numeric(), Label = character(), y_pos = numeric())

  ggplot() +

    geom_area(data = dens_df, aes(x = x, y = y, fill = Category, group = 1), alpha = 0.75, color = NA) +
    geom_line(data = dens_df, aes(x = x, y = y), color = "grey20", linewidth = 0.7) +

    geom_density(data = df_groups, aes(x = r, color = Group), linewidth = 0.9, fill = NA) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#4a8c5c", linewidth = 0.5, alpha = 0.8) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "#8aaa8c", linewidth = 0.5, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "#b0b0b0", linewidth = 0.4) +
    { if (nrow(df_labels) > 0) geom_segment(data = df_labels, aes(x = r, xend = r, y = 0, yend = y_pos - 0.05), color = "grey40", linewidth = 0.3, linetype = "dotted") } +
    { if (nrow(df_labels) > 0) geom_text(data = df_labels, aes(x = r, y = y_pos, label = Label), size = SET_2b3c$label_font, color = "black", fontface = "bold.italic", angle = 45, hjust = 0, vjust = 0.5) } +
    scale_fill_manual(values = cat_colors, guide = "none") +
    scale_color_manual(values = dens_group_colors, guide = "none") +
    scale_x_continuous(limits = if (free_range) NULL else x_lim, breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = if (free_range) waiver() else 0:4,
                       limits = if (free_range) NULL else c(0, 4),
                       expand = if (free_range) expansion(mult = c(0, 0.05)) else c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = SET_2b3c$base_font, base_family = "Arial") +
    labs(title = title_txt, subtitle = subtitle_txt, x = "Pearson r", y = "Density", caption = label_note) +
    theme(plot.title = element_text(size = SET_2b3c$title_font, face = "bold"),
          plot.subtitle = element_text(size = SET_2b3c$subtitle_font, color = "grey40"),
          plot.caption = element_text(size = SET_2b3c$subtitle_font, hjust = 0, face = "italic", color = "grey25"),
          axis.title = element_text(size = SET_2b3c$axis_title_font, face = "bold"),
          axis.text = element_text(size = SET_2b3c$axis_text_font, color = "black"),
          aspect.ratio = SET_2b3c$density_aspect,
          plot.margin = margin(10, 15, 4, 10))
}

build_ks_heatmap <- function(ks_df, modality_tag) {
  grps <- c("CN-Lo", "CN-Hi", "MCI", "AD")
  d <- ks_df %>% filter(Modality == modality_tag, Group1 %in% grps, Group2 %in% grps)
  lookup <- function(a, b, col) {
    if (a == b) return(NA_real_)
    v <- d[[col]][(d$Group1 == a & d$Group2 == b) | (d$Group1 == b & d$Group2 == a)]
    if (length(v)) v[1] else NA_real_
  }
  grid <- expand.grid(row = grps, col = grps, stringsAsFactors = FALSE)
  grid$KS_D <- mapply(function(a, b) lookup(a, b, "KS_D"),   grid$row, grid$col)
  grid$padj <- mapply(function(a, b) lookup(a, b, "KS_padj"), grid$row, grid$col)
  grid$xi <- match(grid$col, grps)
  grid$yi <- (length(grps) + 1) - match(grid$row, grps)
  grid$lab <- ifelse(is.na(grid$KS_D), "", sprintf("%.2f", grid$KS_D))
  grid$lab_face <- ifelse(!is.na(grid$padj) & grid$padj < 0.05, "bold", "plain")
  col_box <- data.frame(xi = 1:4,    yi = 4.65,  grp = grps)
  row_box <- data.frame(xi = 0.30,   yi = 4:1,   grp = grps)
  row_txt <- data.frame(xi = 0.02,   yi = 4:1,   label = grps)
  ggplot() +
    geom_tile(data = grid, aes(xi, yi, fill = KS_D), color = "white", linewidth = 1.0, width = 0.94, height = 0.94) +
    geom_text(data = grid, aes(xi, yi, label = lab), fontface = grid$lab_face, size = SET_2b3c$ks_text_font, color = "grey10") +
    geom_tile(data = col_box, aes(xi, yi), fill = dens_group_colors[col_box$grp], width = 0.94, height = 0.26) +
    geom_tile(data = row_box, aes(xi, yi), fill = dens_group_colors[row_box$grp], width = 0.26, height = 0.94) +
    geom_text(data = row_txt, aes(xi, yi, label = label), hjust = 1, size = SET_2b3c$ks_row_font, fontface = "bold", family = "Arial") +
    scale_fill_gradient(low = SET_2b3c$ks_low, high = SET_2b3c$ks_high, na.value = "grey92", name = "KS D",
                        limits = c(0, max(grid$KS_D, na.rm = TRUE))) +
    coord_equal(clip = "off") +
    scale_x_continuous(limits = c(-0.9, 4.9)) +
    scale_y_continuous(limits = c(-0.4, 5.4)) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"),
          legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
          plot.margin = margin(2, 6, 6, 4))
}

render_density_with_ks <- function(cor_overall_csv, cor_groups_csv, modality_tag, title_txt, subtitle_txt, out_base, free_range = FALSE, label_note = NULL) {
  dov <- viz_read(cor_overall_csv); dgr <- viz_read(cor_groups_csv); ksd <- viz_read("Fig2b3c_density_KS.csv")
  p_dens <- build_density_panel(dov, dgr, title_txt, subtitle_txt, free_range = free_range, label_note = label_note) + .transparent_theme
  p_ks   <- build_ks_heatmap(ksd, modality_tag) + .transparent_theme
  combo  <- cowplot::plot_grid(p_dens, p_ks, ncol = 1, align = "v", axis = "l",
                               rel_heights = SET_2b3c$rel_heights)
  viz_save_both(combo, out_base, SET_2b3c$out_w, SET_2b3c$out_h)
}

n_rna  <- ncol(exprs_aligned)
n_prot <- sum(is.finite(meta_aligned$MS_NPTX2))
render_density_with_ks("Fig2b_cor_overall.csv", "Fig2b_cor_groups.csv", "NPTX2_protein_ref",
                       "Transcriptome vs. NPTX2 Protein Correlation",
                       sprintf("Overall (filled) + per-group densities | n = %d PRM-MS samples", n_prot),
                       "Fig2b_Density_Transcriptome_vs_NPTX2_Protein", free_range = TRUE)
render_density_with_ks("Fig3c_cor_overall.csv", "Fig3c_cor_groups.csv", "NPTX2_mRNA_ref",
                       "Transcriptome vs. NPTX2 mRNA Correlation",
                       sprintf("Overall (filled) + per-group densities | n = %d samples", n_rna),
                       "Fig3c_Density_Transcriptome_vs_NPTX2_mRNA",
                       label_note = "Labeled genes: |overall correlation with NPTX2| \u2265 0.5")
cat("\nFigure 2b / 3c ported.\n")

suppressPackageStartupMessages({ library(ggrepel); library(tidyr) })

cfg_rr <- list(font_family = "Arial", font_title = 16, font_subtitle = 12, font_axis_title = 15,
               font_axis_text = 13, font_legend = 12, dot_size = 2.5, dot_alpha = 0.75, label_size = 3.8,
               plot_width = 7.5, plot_height = 7.0)
traj_colors_25 <- c("A" = "#009E73", "B_rec" = "#E69F00", "B_sup" = "#0072B2",
                    "B_rev" = "#D55E00", "C" = "#CC79A7", "Ref" = "#000000")
traj_labels_25 <- c("A" = "CN-Hi-preserved", "B_rec" = "CN-Hi-recruited", "B_sup" = "CN-Hi-suppressed",
                    "B_rev" = "CN-Hi-reversed", "C" = "Pathology-disrupted", "Ref" = "NPTX2 (Ref)")

hero_genes_of <- function(fl, k = 5) {
  fl <- fl[is.finite(fl$Driving_r), ]
  out <- character(0)
  for (pat in unique(as.character(fl$Pattern))) {
    sub <- fl[as.character(fl$Pattern) == pat, , drop = FALSE]
    pos <- sub[sub$Driving_r > 0, , drop = FALSE]; pos <- pos[order(-abs(pos$Driving_r)), , drop = FALSE]
    neg <- sub[sub$Driving_r < 0, , drop = FALSE]; neg <- neg[order(-abs(neg$Driving_r)), , drop = FALSE]
    out <- c(out, head(pos$Gene, k), head(neg$Gene, k))
  }
  out
}

plot_r_vs_r <- function(df, xlab_group, ylab_group, title, subtitle, bg_df = NULL, hero_genes = NULL) {
  ok <- is.finite(df$r_x) & is.finite(df$r_y); df_ok <- df[ok, ]
  if (nrow(df_ok) < 5) return(NULL)
  df_ok$Label <- if (!is.null(hero_genes)) ifelse(df_ok$Gene %in% hero_genes, df_ok$Gene, NA_character_) else df_ok$Gene
  df_ok$Traj <- factor(df_ok$Traj, levels = c("A", "B_rec", "B_sup", "B_rev", "C", "Ref"))

  all_x <- c(df_ok$r_x, if (!is.null(bg_df)) bg_df$r_x); all_y <- c(df_ok$r_y, if (!is.null(bg_df)) bg_df$r_y)
  ok2 <- is.finite(all_x) & is.finite(all_y); ct <- suppressWarnings(cor.test(all_x[ok2], all_y[ok2]))
  rr_lab <- sprintf("all dots: r = %.2f, p %s (n = %d)", unname(ct$estimate),
                    ifelse(ct$p.value < 2.2e-16, "< 2.2e-16", sprintf("= %.2g", ct$p.value)), sum(ok2))
  ggplot(df_ok, aes(x = r_x, y = r_y, color = Traj)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_hline(yintercept = c(-0.4, 0.4), linetype = "dashed", color = "firebrick3", linewidth = 0.4, alpha = 0.6) +
    geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "firebrick3", linewidth = 0.4, alpha = 0.6) +
    geom_hline(yintercept = c(-0.25, 0.25), linetype = "dotted", color = "steelblue", linewidth = 0.4, alpha = 0.6) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted", color = "steelblue", linewidth = 0.4, alpha = 0.6) +
    { if (!is.null(bg_df) && nrow(bg_df) > 0)
      geom_point(data = bg_df, aes(x = r_x, y = r_y), inherit.aes = FALSE, color = "grey75", size = 0.5, alpha = 0.09) } +
    geom_point(size = cfg_rr$dot_size, alpha = cfg_rr$dot_alpha) +
    annotate("text", x = -0.98, y = 0.98, label = rr_lab, hjust = 0, vjust = 1,
             size = 4, fontface = "bold", color = "grey15") +
    ggrepel::geom_label_repel(aes(label = Label, color = Traj), fill = scales::alpha("white", 0.8),
                              label.size = NA, label.r = unit(0.15, "lines"), label.padding = unit(0.15, "lines"),
                              size = cfg_rr$label_size, fontface = "bold", box.padding = 0.2, point.padding = 0.15,
                              segment.color = "grey40", segment.size = 0.3, force = 0.5, force_pull = 2.0,
                              max.overlaps = Inf, min.segment.length = 0, na.rm = TRUE, show.legend = FALSE) +
    scale_color_manual(values = traj_colors_25, labels = traj_labels_25, name = "Pattern") +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_classic(base_size = 13, base_family = cfg_rr$font_family) +
    labs(title = title, subtitle = sprintf("%s (n = %d genes)", subtitle, nrow(df_ok)),
         x = sprintf("Pearson r (%s)", xlab_group), y = sprintf("Pearson r (%s)", ylab_group)) +
    theme(plot.title = element_text(size = cfg_rr$font_title, face = "bold"),
          plot.subtitle = element_text(size = cfg_rr$font_subtitle, color = "grey30"),
          axis.title = element_text(size = cfg_rr$font_axis_title, face = "bold"),
          axis.text = element_text(size = cfg_rr$font_axis_text, color = "black"),
          legend.title = element_text(size = cfg_rr$font_legend, face = "bold"),
          legend.text = element_text(size = cfg_rr$font_legend), legend.position = "right") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

render_rr <- function(cor_groups_csv, fulllist_csv, modality_label, fig_prefix) {
  cg <- viz_read(cor_groups_csv)
  wide <- as.data.frame(tidyr::pivot_wider(cg[, c("Gene", "Group", "r")], names_from = Group, values_from = r))
  fl <- viz_read(fulllist_csv)
  traj <- setNames(as.character(fl$Traj), fl$Gene)
  wide$Traj <- ifelse(wide$Gene == "NPTX2", "Ref", unname(traj[wide$Gene]))
  if (!"NPTX2" %in% wide$Gene) {
    nr <- wide[1, ]; nr[1, ] <- NA; nr$Gene <- "NPTX2"; nr$Traj <- "Ref"
    for (g in c("CN-Lo", "CN-Hi", "MCI", "AD")) if (g %in% names(nr)) nr[[g]] <- 1
    wide <- rbind(wide, nr)
  }
  traj_genes <- wide[!is.na(wide$Traj), , drop = FALSE]
  bg         <- wide[is.na(wide$Traj) & wide$Gene != "NPTX2", , drop = FALSE]
  hero_genes <- hero_genes_of(fl)
  panels <- list(
    list(x = "CN-Lo", y = "CN-Hi", t = sprintf("%s: CN-Hi vs CN-Lo", modality_label), f = "Hi_vs_Lo"),
    list(x = "AD",    y = "CN-Hi", t = sprintf("%s: CN-Hi vs AD",    modality_label), f = "Hi_vs_AD"),
    list(x = "AD",    y = "CN-Lo", t = sprintf("%s: CN-Lo vs AD",    modality_label), f = "Lo_vs_AD"))
  for (pn in panels) {
    dfp <- data.frame(Gene = traj_genes$Gene, r_x = traj_genes[[pn$x]], r_y = traj_genes[[pn$y]], Traj = traj_genes$Traj)
    bgp <- data.frame(r_x = bg[[pn$x]], r_y = bg[[pn$y]])
    p <- plot_r_vs_r(dfp, pn$x, pn$y, pn$t, "Pattern genes only", bg_df = bgp, hero_genes = hero_genes)
    if (!is.null(p)) viz_save_both(p, sprintf("%s_RvsR_%s", fig_prefix, pn$f), cfg_rr$plot_width, cfg_rr$plot_height)
  }
}

render_rr("Fig3c_cor_groups.csv", "Fig2e3f_FullGeneList_RNA_AD.csv",  "Anchor: NPTX2 mRNA",    "Fig3de")
render_rr("Fig2b_cor_groups.csv", "Fig2e3f_FullGeneList_PROT_AD.csv", "Anchor: NPTX2 protein", "Fig2cd")
cat("\nFigure 2c-d / 3d-e ported.\n")

cfg_hero <- list(
  font_family = "Arial", font_main_title = 12, font_pattern_labels = 9,
  font_row_labels = 8, font_col_labels = 8, font_cell_text = 6,
  font_legend_title = 10, font_legend_labels = 8,
  pattern_label_rot = 270, tile_height_mm = 4.0, tile_width_mm = 8.0,
  pattern_gap_mm = 4.0, canvas_margin_right_mm = 25.0, top_k = 5)
subgroup_col_fun <- colorRamp2(seq(-1, 1, length = 50),
                               colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50))
split_col_vec <- factor(c("1", "2", "2", "2"), levels = c("1", "2"))

build_hero_rows <- function(df, third_col, k = cfg_hero$top_k) {
  df <- df[is.finite(df$Driving_r), ]
  pats <- unique(as.character(df$Pattern))
  mk <- function(blk, pat) data.frame(
    gene = blk$Gene, pattern = pat, label = blk$Gene, spacer = FALSE,
    rov = blk[["r_Overall"]], rlo = blk[["r_CN-Lo"]], rhi = blk[["r_CN-Hi"]], rthird = blk[[third_col]],
    matched = (blk$Match_MCI_AD == 1), stringsAsFactors = FALSE)
  acc <- list()
  for (pat in pats) {
    sub <- df[as.character(df$Pattern) == pat, , drop = FALSE]
    pos <- sub[sub$Driving_r > 0, , drop = FALSE]; pos <- pos[order(-abs(pos$Driving_r)), , drop = FALSE]
    neg <- sub[sub$Driving_r < 0, , drop = FALSE]; neg <- neg[order(-abs(neg$Driving_r)), , drop = FALSE]
    posk <- head(pos, k); negk <- head(neg, k)
    n_hidden <- nrow(sub) - nrow(posk) - nrow(negk)
    if (nrow(posk) > 0) acc[[length(acc) + 1]] <- mk(posk, pat)
    if (n_hidden > 0)   acc[[length(acc) + 1]] <- data.frame(
      gene = paste0(".sp_", pat), pattern = pat, label = sprintf("(%d hidden)", n_hidden),
      spacer = TRUE, rov = NA, rlo = NA, rhi = NA, rthird = NA, matched = FALSE, stringsAsFactors = FALSE)
    if (nrow(negk) > 0) acc[[length(acc) + 1]] <- mk(negk, pat)
  }
  do.call(rbind, acc)
}

hero_group_n <- function(third_grp, prm_only) {
  sel  <- if (prm_only) is.finite(meta_aligned$MS_NPTX2) else rep(TRUE, nrow(meta_aligned))
  grps <- c("Overall", "CN-Lo", "CN-Hi", third_grp)
  sapply(grps, function(g) if (g == "Overall") sum(sel) else sum(meta_aligned$Group_Current == g & sel))
}

build_hero_ht <- function(df_full, third_col, third_grp, prm_only) {
  hr <- build_hero_rows(df_full, third_col)
  n_by <- hero_group_n(third_grp, prm_only)
  col_lab <- sprintf("%s\n(n=%d)", c("Overall", "CN-Lo", "CN-Hi", third_grp), n_by)
  mat_r <- as.matrix(hr[, c("rov", "rlo", "rhi", "rthird")])
  rownames(mat_r) <- hr$gene; colnames(mat_r) <- col_lab

  row_colors    <- ifelse(hr$spacer, "grey50", ifelse(hr$matched, "red", "black"))
  row_fontfaces <- ifelse(hr$spacer, "italic", "plain")
  pattern_factor <- factor(hr$pattern, levels = unique(hr$pattern))
  levels(pattern_factor) <- sub("-([a-z])", "-\n\\1", levels(pattern_factor))
  pat_colors <- sapply(levels(pattern_factor), function(lbl) {
    if (grepl("preserved",  lbl, ignore.case = TRUE)) return("#009E73")
    if (grepl("recruited",  lbl, ignore.case = TRUE)) return("#E69F00")
    if (grepl("suppressed", lbl, ignore.case = TRUE)) return("#0072B2")
    if (grepl("reversed",   lbl, ignore.case = TRUE)) return("#D55E00")
    if (grepl("Pathology",  lbl, ignore.case = TRUE)) return("#CC79A7")
    "black" })

  ht <- Heatmap(mat_r, name = "Pearson r", col = subgroup_col_fun,
                column_split = split_col_vec, cluster_column_slices = FALSE,
                row_split = pattern_factor, cluster_row_slices = FALSE, row_gap = unit(3, "mm"),
                row_title_side = "right", row_title_rot = cfg_hero$pattern_label_rot,
                row_title_gp = gpar(fontsize = cfg_hero$font_pattern_labels, fontfamily = cfg_hero$font_family,
                                    fontface = "bold", col = pat_colors),
                right_annotation = rowAnnotation(spacer = anno_empty(border = FALSE,
                                                                     width = unit(cfg_hero$pattern_gap_mm, "mm"))),
                cluster_columns = FALSE, cluster_rows = FALSE,
                show_row_names = TRUE, row_names_side = "left", row_labels = hr$label,
                row_names_gp = gpar(fontsize = cfg_hero$font_row_labels, fontfamily = cfg_hero$font_family,
                                    col = row_colors, fontface = row_fontfaces),
                column_names_gp = gpar(fontsize = cfg_hero$font_col_labels, fontfamily = cfg_hero$font_family),
                column_names_rot = 45, column_title = NULL, na_col = "grey40",
                height = nrow(mat_r) * unit(cfg_hero$tile_height_mm, "mm"),
                width  = ncol(mat_r) * unit(cfg_hero$tile_width_mm, "mm"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  v <- mat_r[i, j]
                  if (!is.na(v)) grid::grid.text(sprintf("%.2f", v), x, y,
                                                 gp = grid::gpar(fontsize = cfg_hero$font_cell_text,
                                                                 col = ifelse(abs(v) >= 0.6, "white", "black"), fontfamily = cfg_hero$font_family))
                },
                heatmap_legend_param = list(title = "Pearson r", title_position = "leftcenter-rot",
                                            title_gp = gpar(fontfamily = cfg_hero$font_family, fontface = "bold", fontsize = cfg_hero$font_legend_title),
                                            labels_gp = gpar(fontfamily = cfg_hero$font_family, fontsize = cfg_hero$font_legend_labels)))
  list(ht = ht, nrow = nrow(mat_r), ncol = ncol(mat_r), ngroups = nlevels(pattern_factor))
}

render_hero_single <- function(modality_word, prm_only, tag, third_col, out_base, show_legend = TRUE) {
  suffix <- if (prm_only) "PROT" else "RNA"
  h <- build_hero_ht(viz_read(sprintf("Fig2e3f_FullGeneList_%s_%s.csv", suffix, tag)), third_col, tag, prm_only)
  w_in <- h$ncol * (cfg_hero$tile_width_mm  / 25.4) + 3.8
  h_in <- h$nrow * (cfg_hero$tile_height_mm / 25.4) + (h$ngroups - 1) * (3 / 25.4) + 2.2
  do_draw <- function() draw(h$ht, show_heatmap_legend = show_legend,
                             column_title = sprintf("Anchor: NPTX2 %s (-%s)", modality_word, tag),
                             column_title_gp = gpar(fontsize = cfg_hero$font_main_title, fontfamily = cfg_hero$font_family, fontface = "bold"),
                             padding = unit(c(8, 12, 14, cfg_hero$canvas_margin_right_mm), "mm"))
  svglite::svglite(file.path(svg_dir, paste0(out_base, ".svg")), width = w_in, height = h_in,
                   system_fonts = list(sans = "Arial"), bg = "transparent")
  do_draw(); dev.off()
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file.path(png_dir, paste0(out_base, ".png")), width = w_in, height = h_in, units = "in", res = 300, background = "transparent")
  else grDevices::png(file.path(png_dir, paste0(out_base, ".png")), width = w_in, height = h_in, units = "in", res = 300, bg = "transparent")
  do_draw(); dev.off()
  cat(sprintf("  Saved: %s.{png,svg}\n", out_base))
}

render_hero_single("mRNA",    FALSE, "AD",  "r_AD",  "Fig3f_Hero_Heatmap_RNA_AD",  show_legend = FALSE)
render_hero_single("mRNA",    FALSE, "MCI", "r_MCI", "Fig3f_Hero_Heatmap_RNA_MCI", show_legend = TRUE)
render_hero_single("protein", TRUE,  "AD",  "r_AD",  "Fig2e_Hero_Heatmap_PROT_AD", show_legend = FALSE)
render_hero_single("protein", TRUE,  "MCI", "r_MCI", "Fig2e_Hero_Heatmap_PROT_MCI", show_legend = TRUE)
cat("\nFigure 2e / 3f ported (AD legend dropped; shown on MCI).\n")

suppressPackageStartupMessages({ library(stringr) })

cfg_ora <- list(font_family = "Arial", text_wrap_width = 25, canvas_width_in = 5.5, min_canvas_height = 3.5,
                height_scale = 0.75)

ora_read <- function(name) {
  f <- file.path(in_dir, name)
  if (file.exists(f)) as.data.frame(readr::read_csv(f, show_col_types = FALSE)) else NULL
}

run_custom_go_plot <- function(df_res, modality_prefix, classification_tag, fig_prefix, use_matched_coloring = FALSE) {
  if (is.null(df_res) || nrow(df_res) == 0) {
    cat(sprintf("  ORA %s %s: no rows, skipped\n", modality_prefix, classification_tag)); return(invisible(NULL))
  }
  df_top <- df_res %>% group_by(Cluster) %>% arrange(desc(FoldEnrichment)) %>% slice_head(n = 8) %>% ungroup()
  df_top$Description_wrap <- str_wrap(df_top$Description, width = cfg_ora$text_wrap_width)
  global_max_x <- max(df_top$GeneRatio_num, na.rm = TRUE) * 1.05
  for (pat in unique(df_top$Cluster)) {
    df_sub <- df_top %>% filter(Cluster == pat)
    df_sub$Description_wrap <- factor(df_sub$Description_wrap, levels = rev(df_sub$Description_wrap))
    if (use_matched_coloring && "is_concordant" %in% colnames(df_sub)) {
      axis_colors <- ifelse(df_sub$is_concordant[match(levels(df_sub$Description_wrap), df_sub$Description_wrap)], "red", "black")
      n_red <- sum(df_sub$is_concordant, na.rm = TRUE)
      subtitle_text <- sprintf("Pattern: %s\nRanked by Fold Enrichment | %d/%d pathways \u226550%% matched", pat, n_red, nrow(df_sub))
    } else {
      axis_colors <- rep("black", length(levels(df_sub$Description_wrap)))
      subtitle_text <- sprintf("Pattern: %s\nRanked by Fold Enrichment", pat)
    }
    mod_word <- if (modality_prefix == "PROT") "protein" else "mRNA"
    p <- ggplot(df_sub, aes(x = GeneRatio_num, y = Description_wrap)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(name = "BHP", low = "red", high = "blue",
                           labels = function(v) sprintf("%.2f", v),
                           guide = guide_colorbar(reverse = TRUE, order = 1, barwidth = unit(12, "pt"), barheight = unit(60, "pt"))) +
      scale_size_continuous(name = "Count", guide = guide_legend(order = 2)) +
      scale_x_continuous(limits = c(0, global_max_x), breaks = scales::pretty_breaks(n = 4)) +
      theme_classic(base_family = cfg_ora$font_family) +
      labs(title = sprintf("Anchor: NPTX2 %s (-%s): Cellular Component", mod_word, classification_tag),
           subtitle = subtitle_text, x = "Gene Ratio", y = NULL) +
      theme(axis.text.y = element_text(size = 11, face = "bold", color = axis_colors, lineheight = 0.8),
            axis.text.x = element_text(size = 10, face = "bold", color = "black"),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey30", lineheight = 1.1),
            plot.title.position = "plot", plot.margin = margin(15, 15, 15, 20),
            legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9),
            panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"))
    panel_height_in <- max(1.5, nrow(df_sub) * 0.5)
    canvas_height   <- cfg_ora$height_scale * max(cfg_ora$min_canvas_height, panel_height_in + 2.5)
    safe_pat <- gsub("[^A-Za-z0-9]", "_", pat)
    viz_save_both(p, sprintf("%s_ORA_%s_%s_%s", fig_prefix, modality_prefix, classification_tag, safe_pat),
                  cfg_ora$canvas_width_in, canvas_height)
  }
}

run_custom_go_plot(ora_read("Fig3g_ORA_RNA_AD_CC.csv"),   "RNA",  "AD",  "Fig3g", use_matched_coloring = TRUE)
run_custom_go_plot(ora_read("Fig3g_ORA_RNA_MCI_CC.csv"),  "RNA",  "MCI", "Fig3g", use_matched_coloring = FALSE)
run_custom_go_plot(ora_read("Fig2f_ORA_PROT_AD_CC.csv"),  "PROT", "AD",  "Fig2f", use_matched_coloring = TRUE)
run_custom_go_plot(ora_read("Fig2f_ORA_PROT_MCI_CC.csv"), "PROT", "MCI", "Fig2f", use_matched_coloring = FALSE)
cat("\nFigure 2f / 3g ported.\n")

SET_3a <- list(
  out_w = 7.5, out_h = 6,
  title = "NPTX2 mRNA vs. NPTX2 Protein (PRM-MS)",
  x_lab = "NPTX2 mRNA (Log2)", y_lab = "NPTX2 Protein (PRM-MS)",
  base_font = 16, title_font = 16, subtitle_font = 14, axis_title_font = 14,
  axis_text_font = 13, legend_font = 14, anno_font = 4,
  point_size = 2.0, point_alpha = 0.75, fit_lwd = 0.7, overall_fit_lwd = 0.9,
  colors = c("CN-Lo" = "lightskyblue", "CN-Hi" = "royalblue", "MCI" = "orange", "AD" = "brown3"))

group_colors_23 <- SET_3a$colors
fmt_bhp_23 <- function(b) if (is.na(b)) "NA" else if (b < 0.001) "<0.001" else if (b < 0.01) "<0.01" else if (b < 0.05) "<0.05" else sprintf("%.2f", b)
pts_23 <- viz_read("Fig3a_points.csv")
cor_23 <- viz_read("Fig3a_cor_stats.csv")
pts_23$Group <- factor(pts_23$Group, levels = c("CN-Lo", "CN-Hi", "MCI", "AD"))
ov_23 <- cor_23[cor_23$Group == "Overall", ]
subtitle_23 <- sprintf("Overall: r = %.2f, BHP %s (n = %d)", ov_23$r, fmt_bhp_23(ov_23$BHP), ov_23$n)
grp_23 <- cor_23[cor_23$Group != "Overall", ]
grp_anno_23 <- paste(sprintf("%s: r = %.2f, BHP %s (n=%d)", grp_23$Group, grp_23$r,
                             vapply(grp_23$BHP, fmt_bhp_23, character(1)), grp_23$n), collapse = "\n")
p_23 <- ggplot(pts_23, aes(x = NPTX2_RNA, y = NPTX2_MS, color = Group)) +
  geom_point(size = SET_3a$point_size, alpha = SET_3a$point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = SET_3a$fit_lwd, linetype = "dashed", show.legend = FALSE) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", fill = "grey80",
              linewidth = SET_3a$overall_fit_lwd, alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = group_colors_23, name = "Group") +
  annotate("text", x = Inf, y = -Inf, label = grp_anno_23, hjust = 1.05, vjust = -0.3,
           size = SET_3a$anno_font, family = "Arial", color = "grey20") +
  theme_classic(base_size = SET_3a$base_font, base_family = "Arial") +
  labs(title = SET_3a$title, subtitle = subtitle_23, x = SET_3a$x_lab, y = SET_3a$y_lab) +
  theme(plot.title = element_text(size = SET_3a$title_font, face = "bold"), plot.subtitle = element_text(size = SET_3a$subtitle_font, color = "grey30"),
        axis.title = element_text(size = SET_3a$axis_title_font, face = "bold"), axis.text = element_text(size = SET_3a$axis_text_font, color = "black"),
        legend.title = element_text(size = SET_3a$legend_font + 1, face = "bold"), legend.text = element_text(size = SET_3a$legend_font), legend.position = "right")
viz_save_both(p_23, "Fig3a_NPTX2_mRNA_vs_Protein", SET_3a$out_w, SET_3a$out_h)
cat("\nFigure 3a ported.\n")

suppressPackageStartupMessages({ library(ggrepel); library(scales) })

FIG4A <- list(width = 10.5, height = 7.0,
              delta_eq = 0.50, alpha_tost = 0.05, p_de_max = 0.05, fc_threshold = 0.50, font_family = "Arial")
cfg_de <- FIG4A
pattern_colors <- c(
  "CN-Hi-preserved"     = "#009E73",
  "CN-Hi-recruited"     = "#E69F00",
  "CN-Hi-suppressed"    = "#0072B2",
  "CN-Hi-reversed"      = "#D55E00",
  "Pathology-disrupted" = "#CC79A7")

build_de_plot <- function(df_all, df_ghost_all, df_anchor, y_col, target_flag, shared_flag,
                          third_group_name, target_legend_label, title_text, subtitle_detail) {
  df_bg     <- df_all %>% filter(!(!!sym(target_flag)))
  df_fg     <- df_all %>% filter(!!sym(target_flag))
  df_shared <- df_all %>% filter(!!sym(shared_flag))
  max_x <- max(abs(c(df_all$Delta_Lo, df_anchor$Delta_Lo)), na.rm = TRUE) * 1.05
  max_y <- max(abs(c(df_all[[y_col]], df_anchor[[y_col]])), na.rm = TRUE) * 1.05
  y_axis_lab <- bquote("\u2190" ~ "Upregulated in CN-Hi" ~ "        " ~
                         Log[2] * "FC (" * .(third_group_name) * " vs CN-Hi)" ~ "        " ~
                         "Upregulated in " * .(third_group_name) ~ "\u2192")
  x_axis_lab <- bquote("\u2190" ~ "Upregulated in CN-Hi" ~ "        " ~
                         Log[2] * "FC (CN-Lo vs CN-Hi)" ~ "        " ~ "Upregulated in CN-Lo" ~ "\u2192")
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_hline(yintercept = c(-cfg_de$fc_threshold, cfg_de$fc_threshold), linetype = "dotted", color = "grey80") +
    geom_vline(xintercept = c(-cfg_de$delta_eq, cfg_de$delta_eq), linetype = "dotted", color = "grey80") +
    geom_point(data = df_ghost_all, aes(x = Delta_Lo, y = .data[[y_col]]),
               color = "grey75", size = 0.5, alpha = 0.09) +
    geom_point(data = df_bg, aes(x = Delta_Lo, y = .data[[y_col]], color = Pattern, size = Max_Abs_R), alpha = 1.0) +
    geom_point(data = df_fg, aes(x = Delta_Lo, y = .data[[y_col]], size = Max_Abs_R, fill = target_legend_label),
               color = "black", shape = 21, stroke = 0, alpha = 1.0) +
    geom_point(data = df_shared, aes(x = Delta_Lo, y = .data[[y_col]], shape = "Resilient-Expression (AD + MCI)"),
               color = "red", size = 5.5, stroke = 0.8, fill = NA) +
    scale_fill_manual(values = setNames("black", target_legend_label), name = NULL) +
    scale_shape_manual(values = c("Resilient-Expression (AD + MCI)" = 1), name = NULL) +
    geom_point(data = df_anchor, aes(x = Delta_Lo, y = .data[[y_col]]),
               shape = 21, size = 11, color = "red", fill = "white", stroke = 1.5,
               inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = df_anchor, aes(x = Delta_Lo, y = .data[[y_col]], label = Gene),
              color = "red", fontface = "bold", size = 2.5, inherit.aes = FALSE, show.legend = FALSE) +
    ggrepel::geom_label_repel(data = df_fg,
                              aes(x = Delta_Lo, y = .data[[y_col]], label = Gene, color = Pattern),
                              fill = scales::alpha("white", 0.85), label.size = NA, label.r = unit(0.3, "lines"),
                              size = 3.0, fontface = "bold", box.padding = 0.6, point.padding = 0.5,
                              segment.color = "black", segment.size = 0.5, force = 4.0, max.overlaps = Inf, show.legend = FALSE) +
    coord_cartesian(xlim = c(-max_x, max_x), ylim = c(-max_y, max_y)) +
    theme_classic(base_family = cfg_de$font_family) +
    scale_color_manual(values = pattern_colors, name = "Pattern class", breaks = names(pattern_colors)) +
    scale_size_continuous(range = c(0.2, 6.0), limits = c(0.3, 1.0), name = "Max |NPTX2\ncorrelation|") +
    labs(title = title_text, subtitle = subtitle_detail, x = x_axis_lab, y = y_axis_lab) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey30", lineheight = 1.2),
          axis.title.x = element_text(size = 11, face = "bold"),
          axis.title.y = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 10, color = "black"),
          legend.position = "right", legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(size = 10)) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), order = 1),
           size = guide_legend(override.aes = list(color = "grey60"), order = 2),
           fill = guide_legend(override.aes = list(color = "black", fill = "black", size = 3, shape = 21, stroke = 0), order = 3),
           shape = guide_legend(override.aes = list(color = "red", size = 5.5, stroke = 0.8, fill = NA), order = 4))
}

df_de    <- viz_read("Fig4a_Double_DE.csv")
df_ghost <- viz_read("Fig4a_ghost.csv")
df_nptx2 <- viz_read("Fig4a_nptx2_anchor.csv")
for (cc in c("Is_Resilient_AD", "Is_Resilient_MCI", "Is_Resilient_Both"))
  df_de[[cc]] <- as.logical(df_de[[cc]])
n_shared <- sum(df_de$Is_Resilient_Both, na.rm = TRUE)
subtitle_ad <- sprintf(paste0(
  "Black dots (Resilient-Expression (AD)): TOST CN-Lo\u2261CN-Hi at \u0394=%.2f (\u03b1=%.2f),\n",
  "p(CN-Hi vs AD) < %.2f, |Log2FC(AD vs CN-Hi)| > %.2f.\n",
  "Red ring (AD + MCI): also resilient on MCI axis AND pattern-matched (n = %d).\n"),
  cfg_de$delta_eq, cfg_de$alpha_tost, cfg_de$p_de_max, cfg_de$fc_threshold, n_shared)
subtitle_mci <- sprintf(paste0(
  "Black dots (Resilient-Expression (MCI)): TOST CN-Lo\u2261CN-Hi at \u0394=%.2f (\u03b1=%.2f),\n",
  "p(CN-Hi vs MCI) < %.2f; no |FC| filter on MCI axis.\n",
  "Red ring (AD + MCI): also resilient on AD axis (|Log2FC|>%.2f) AND pattern-matched (n = %d).\n",
  "Red \u2605 marks NPTX2 reference."),
  cfg_de$delta_eq, cfg_de$alpha_tost, cfg_de$p_de_max, cfg_de$fc_threshold, n_shared)
p_de_ad <- build_de_plot(df_de, df_ghost, df_nptx2, "Delta_AD", "Is_Resilient_AD", "Is_Resilient_Both",
                         "AD", "Resilient-Expression (AD)", "Double DE Map: Pathology \u00d7 Progression (AD)", subtitle_ad)
p_de_mci <- build_de_plot(df_de, df_ghost, df_nptx2, "Delta_MCI", "Is_Resilient_MCI", "Is_Resilient_Both",
                          "MCI", "Resilient-Expression (MCI)", "Double DE Map: Pathology \u00d7 Symptom onset (MCI)", subtitle_mci)
viz_save_both(p_de_ad,  "Fig4a_Double_DE_AD",  FIG4A$width, FIG4A$height)
viz_save_both(p_de_mci, "Fig4a_Double_DE_MCI", FIG4A$width, FIG4A$height)
cat("\nFigure 4a ported.\n")

cfg_traj4 <- list(font_family = "Arial", font_title = 14, font_subtitle = 11, font_axis_title = 11.5,
                  font_axis_text = 11, point_size = 5.0, error_bar_width = 0.15,

                  col_width = 2.3,
                  row_h_cor = 3.0, row_h_exp = 3.0,
                  label_size = 15,

                  header_cor = "r vs NPTX2 mRNA (mean \u00b1 SE). Bars: Fisher z-transformation test (p shown)",   ylab_cor = "Pearson r",
                  header_exp = "Scaled expression 0\u20131 (box = median & IQR; dots = individual samples). Bars: Wilcoxon rank-sum test (p shown)", ylab_exp = "Scaled (0\u20131)",
                  header_font = 12, ylab_font = 12,
                  ylab_col_frac = 0.042, header_row_frac = 0.1)

calc_cor_stats_t4 <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4) return(list(r = NA, lower = NA, upper = NA, n = n, z = NA))
  r <- cor(x, y); rc <- pmax(pmin(r, 0.999), -0.999); z <- 0.5 * log((1 + rc) / (1 - rc)); se <- 1 / sqrt(n - 3)
  rci <- function(zz) (exp(2 * zz) - 1) / (exp(2 * zz) + 1)
  list(r = r, lower = rci(z - se), upper = rci(z + se), n = n, z = z)
}
fmt_p <- function(p) if (is.na(p)) "n/a" else if (p < 1e-4) sprintf("%.1e", p) else sprintf("%.3g", p)
compare_cors_t4 <- function(s1, s2) {
  if (is.na(s1$z) || is.na(s2$z)) return("n/a")
  zd <- (s1$z - s2$z) / sqrt(1 / (s1$n - 3) + 1 / (s2$n - 3)); fmt_p(2 * (1 - pnorm(abs(zd))))
}
compare_expr_t4 <- function(v1, v2) {
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) < 3 || length(v2) < 3) return("n/a")
  fmt_p(suppressWarnings(wilcox.test(v1, v2)$p.value))
}
minmax_t4 <- function(v) { fin <- v[is.finite(v)]; lo <- min(fin); hi <- max(fin)
if (is.finite(lo) && is.finite(hi) && hi > lo) (v - lo) / (hi - lo) else rep(NA, length(v)) }
ci95_t4 <- function(x) { x <- x[is.finite(x)]; n <- length(x)
if (n < 2) return(list(mean = NA, lower = NA, upper = NA))
m <- mean(x); se <- sd(x) / sqrt(n); list(mean = m, lower = m - se, upper = m + se) }

all_s_t4 <- intersect(meta_aligned$SampleID, colnames(exprs_aligned))
sg_t4 <- lapply(c("CN-Lo", "CN-Hi", "MCI", "AD"),
                function(g) intersect(meta_aligned$SampleID[meta_aligned$Group_Current == g], all_s_t4))
names(sg_t4) <- c("CN-Lo", "CN-Hi", "MCI", "AD")
ref_nptx2_t4 <- setNames(as.numeric(exprs_aligned["NPTX2", all_s_t4]), all_s_t4)
states_t4 <- factor(c("CN-Lo", "CN-Hi", "MCI", "AD"), levels = c("CN-Lo", "CN-Hi", "MCI", "AD"))

n_str_t4 <- paste(sprintf("%s=%d", names(sg_t4), lengths(sg_t4)), collapse = ", ")

traj4_theme <- function(show_x = TRUE) theme_classic(base_family = cfg_traj4$font_family) +
  theme(plot.title = element_text(size = cfg_traj4$font_title, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = if (show_x) element_text(size = cfg_traj4$font_axis_text, face = "bold", color = "black", angle = 45, hjust = 1) else element_blank(),
        axis.ticks.x = if (show_x) element_line() else element_blank(),
        axis.text.y = element_text(size = cfg_traj4$font_axis_text, color = "black"),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 6),
        panel.grid.major.y = element_line(color = "grey95", linetype = "dotted"))

generate_traj4 <- function(gene_name) {
  if (!gene_name %in% rownames(exprs_aligned)) { cat(sprintf("  traj4: %s not in matrix, skipped\n", gene_name)); return(invisible(NULL)) }
  gv <- setNames(as.numeric(exprs_aligned[gene_name, all_s_t4]), all_s_t4)

  s  <- lapply(sg_t4, function(ss) calc_cor_stats_t4(gv[ss], ref_nptx2_t4[ss]))
  df_cor <- data.frame(State = states_t4, Value = sapply(s, `[[`, "r"),
                       ymin = sapply(s, `[[`, "lower"), ymax = sapply(s, `[[`, "upper"))
  sc_lh <- compare_cors_t4(s[["CN-Lo"]], s[["CN-Hi"]])
  sc_hm <- compare_cors_t4(s[["CN-Hi"]], s[["MCI"]])
  sc_ha <- compare_cors_t4(s[["CN-Hi"]], s[["AD"]])

  gs <- minmax_t4(gv); names(gs) <- names(gv)
  ci <- lapply(sg_t4, function(ss) ci95_t4(gs[ss]))
  df_exp <- data.frame(State = states_t4, Value = sapply(ci, `[[`, "mean"),
                       ymin = sapply(ci, `[[`, "lower"), ymax = sapply(ci, `[[`, "upper"))
  se_lh <- compare_expr_t4(gs[sg_t4[["CN-Lo"]]], gs[sg_t4[["CN-Hi"]]])
  se_hm <- compare_expr_t4(gs[sg_t4[["CN-Hi"]]], gs[sg_t4[["MCI"]]])
  se_ha <- compare_expr_t4(gs[sg_t4[["CN-Hi"]]], gs[sg_t4[["AD"]]])

  yb1 <- 1.06; yb2 <- 1.40
  p_cor <- ggplot(df_cor, aes(State, Value)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.8) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = cfg_traj4$error_bar_width, linewidth = 0.8, color = "black") +
    geom_point(size = cfg_traj4$point_size, color = "black") +
    annotate("segment", x = 1.05, xend = 1.95, y = yb1, yend = yb1, linewidth = 0.8, color = "black") +
    annotate("text", x = 1.5, y = yb1 + 0.09, vjust = 0, label = sc_lh, size = 3.5, fontface = "bold", family = cfg_traj4$font_family) +
    annotate("segment", x = 2.05, xend = 2.95, y = yb1, yend = yb1, linewidth = 0.8, color = "black") +
    annotate("text", x = 2.5, y = yb1 + 0.09, vjust = 0, label = sc_hm, size = 3.5, fontface = "bold", family = cfg_traj4$font_family) +
    annotate("segment", x = 2.05, xend = 3.95, y = yb2, yend = yb2, linewidth = 0.8, color = "black") +
    annotate("text", x = 3.0, y = yb2 + 0.09, vjust = 0, label = sc_ha, size = 3.5, fontface = "bold", family = cfg_traj4$font_family) +
    scale_y_continuous(breaks = seq(-0.5, 1, 0.5)) +
    coord_cartesian(ylim = c(-0.5, 1.65)) +
    labs(title = gene_name, x = NULL, y = NULL) +
    traj4_theme(show_x = TRUE)

  eb1 <- 1.00; eb2 <- 1.26
  df_exp_pts <- do.call(rbind, lapply(names(sg_t4), function(g) {
    v <- gs[sg_t4[[g]]]; v <- as.numeric(v[is.finite(v)])
    if (length(v) == 0) NULL else data.frame(State = factor(g, levels = levels(states_t4)), Value = v)
  }))
  p_exp <- ggplot(df_exp_pts, aes(State, Value)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA, color = "firebrick3", linewidth = 0.6) +
    geom_jitter(width = 0.16, height = 0, size = 1.0, alpha = 0.40, color = "firebrick3") +
    annotate("segment", x = 1.05, xend = 1.95, y = eb1, yend = eb1, linewidth = 0.8, color = "firebrick3") +
    annotate("text", x = 1.5, y = eb1 + 0.07, vjust = 0, label = se_lh, size = 3.5, fontface = "bold", family = cfg_traj4$font_family, color = "firebrick3") +
    annotate("segment", x = 2.05, xend = 2.95, y = eb1, yend = eb1, linewidth = 0.8, color = "firebrick3") +
    annotate("text", x = 2.5, y = eb1 + 0.07, vjust = 0, label = se_hm, size = 3.5, fontface = "bold", family = cfg_traj4$font_family, color = "firebrick3") +
    annotate("segment", x = 2.05, xend = 3.95, y = eb2, yend = eb2, linewidth = 0.8, color = "firebrick3") +
    annotate("text", x = 3.0, y = eb2 + 0.07, vjust = 0, label = se_ha, size = 3.5, fontface = "bold", family = cfg_traj4$font_family, color = "firebrick3") +
    scale_y_continuous(limits = c(0, 1.50), breaks = seq(0, 1, 0.5)) +
    labs(title = gene_name, x = NULL, y = NULL) +
    traj4_theme(show_x = TRUE)

  invisible(list(cor = p_cor, exp = p_exp))
}

.f4a <- viz_read("Fig4a_Double_DE.csv")
traj4_genes <- unique(.f4a$Gene[.f4a$Is_Resilient_Both %in% c(TRUE, 1, "TRUE")])
traj4_genes <- traj4_genes[!is.na(traj4_genes)]
cat(sprintf("  Fig 4b-m: %d resilient-both genes -> %s\n", length(traj4_genes), paste(traj4_genes, collapse = ", ")))

panels_t4 <- Filter(Negate(is.null), lapply(traj4_genes, generate_traj4))
if (length(panels_t4) > 0) {
  cor_row <- lapply(panels_t4, `[[`, "cor")
  exp_row <- lapply(panels_t4, `[[`, "exp")
  n_g     <- length(panels_t4)
  cor_grid <- cowplot::plot_grid(plotlist = cor_row, nrow = 1, align = "h")
  exp_grid <- cowplot::plot_grid(plotlist = exp_row, nrow = 1, align = "h")

  row_block <- function(grid, ylab, header) {
    body <- cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(ylab, x = 0.9, angle = 90, fontface = "bold",
                                              size = cfg_traj4$ylab_font, fontfamily = cfg_traj4$font_family),
      grid, ncol = 2, rel_widths = c(cfg_traj4$ylab_col_frac, 1))
    cowplot::plot_grid(
      body,
      cowplot::ggdraw() + cowplot::draw_label(header, fontface = "bold",
                                              size = cfg_traj4$header_font, fontfamily = cfg_traj4$font_family),
      ncol = 1, rel_heights = c(1, cfg_traj4$header_row_frac))
  }
  cor_block <- row_block(cor_grid, cfg_traj4$ylab_cor, paste0(cfg_traj4$header_cor, ".  n per group: ", n_str_t4))
  exp_block <- row_block(exp_grid, cfg_traj4$ylab_exp, paste0(cfg_traj4$header_exp, ".  n per group: ", n_str_t4))
  fig_w <- cfg_traj4$col_width * n_g
  viz_save_both(cor_block, "Fig4_Resilient_Correlation", fig_w, cfg_traj4$row_h_cor)
  viz_save_both(exp_block, "Fig4_Resilient_Expression",  fig_w, cfg_traj4$row_h_exp)
}
cat("\nFigure 4b-m ported (correlation + expression as two aligned figures; SE bars; real p).\n")

write_class_counts <- function(suffix, mod_label) {
  truthy <- function(v) v %in% c(1, TRUE, "TRUE")
  rows <- list()
  for (tag in c("AD", "MCI")) {
    fl <- viz_read(sprintf("Fig2e3f_FullGeneList_%s_%s.csv", suffix, tag))
    has_match <- "Match_MCI_AD" %in% names(fl)
    other <- if (tag == "AD") "MCI" else "AD"
    for (p in unique(as.character(fl$Pattern))) {
      sub  <- fl[as.character(fl$Pattern) == p, , drop = FALSE]
      midx <- if (has_match) truthy(sub$Match_MCI_AD) else rep(FALSE, nrow(sub))
      rows[[length(rows) + 1]] <- data.frame(
        Modality = mod_label, Classified_in = tag, Recapitulated_stage = other, Class = p,
        N_genes = nrow(sub), N_recapitulated = sum(midx),
        Genes = paste(sub$Gene, collapse = "; "),
        Genes_recapitulated = paste(sub$Gene[midx], collapse = "; "),
        stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1]] <- data.frame(
      Modality = mod_label, Classified_in = tag, Recapitulated_stage = other, Class = "TOTAL",
      N_genes = nrow(fl), N_recapitulated = if (has_match) sum(truthy(fl$Match_MCI_AD)) else NA_integer_,
      Genes = paste(fl$Gene, collapse = "; "),
      Genes_recapitulated = if (has_match) paste(fl$Gene[truthy(fl$Match_MCI_AD)], collapse = "; ") else "",
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  fn <- file.path(fig_dir, sprintf("ClassCounts_%s.csv", suffix))
  readr::write_csv(out, fn)
  cat(sprintf("  Saved: %s\n", basename(fn)))
  print(out[, c("Classified_in", "Class", "N_genes", "N_recapitulated")], row.names = FALSE)
  invisible(out)
}
write_class_counts("RNA",  "Transcriptome vs NPTX2 mRNA")
write_class_counts("PROT", "Transcriptome vs NPTX2 protein")
cat("\nClassification count summaries written.\n")

suppressPackageStartupMessages({ library(ggrepel); library(ComplexHeatmap); library(circlize); library(RColorBrewer) })

supp_read <- function(name) { f <- file.path(in_dir, name); if (file.exists(f)) as.data.frame(readr::read_csv(f, show_col_types = FALSE)) else NULL }
sig_stars_supp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

SET_SUPP <- list(
  cov_tile_mm = 7,
  cov_w_pad = 3.2, cov_h_pad = 2.2,
  pca_w = 7.5, pca_h = 6.0,
  fgsea_w = 7.0, fgsea_h = 3.5,
  sexdens_overall = c(7, 5),
  sexdens_bygroup = c(8, 7),
  composite = c(16, 10))

cfg_33 <- list(fdr_thresh = 0.05, fc_thresh = 0.20, font_family = "Arial",
               point_size = 3.2, label_size = 3.5, plot_w = 7.0, plot_h = 6.5,
               col_up = "#A32D2D", col_down = "#185FA5", col_ns = "grey60")
build_arrow_xlab_33 <- function(numerator, denominator) bquote(
  "\u2190" ~ "Upregulated in " * .(denominator) ~ "      " ~
    Log[2] * " FC (" * .(numerator) * " vs " * .(denominator) * ")" ~ "      " ~
    "Upregulated in " * .(numerator) ~ "\u2192")
build_volcano_33 <- function(df, contrast_label, numerator, denominator) {

  df$Direction <- ifelse(!is.na(df$adj.P.Val) & df$adj.P.Val < cfg_33$fdr_thresh & df$log2FC >  cfg_33$fc_thresh, "up",
                         ifelse(!is.na(df$adj.P.Val) & df$adj.P.Val < cfg_33$fdr_thresh & df$log2FC < -cfg_33$fc_thresh, "down", "ns"))
  df$neglog10p <- -log10(pmax(df$adj.P.Val, 1e-300))
  x_max <- max(c(max(abs(df$log2FC), na.rm = TRUE) * 1.15, 0.5))
  y_max <- max(df$neglog10p, na.rm = TRUE) * 1.10 + 0.3
  n_up <- sum(df$Direction == "up"); n_down <- sum(df$Direction == "down")
  d_ns <- df %>% dplyr::filter(Direction == "ns")
  d_sig <- df %>% dplyr::filter(Direction != "ns")
  p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
    geom_vline(xintercept = c(-cfg_33$fc_thresh, cfg_33$fc_thresh), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_hline(yintercept = -log10(cfg_33$fdr_thresh), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3)
  if (nrow(d_ns) > 0) p <- p + geom_point(data = d_ns, color = cfg_33$col_ns, size = cfg_33$point_size, alpha = 0.55)
  if (nrow(d_sig) > 0) p <- p + geom_point(data = d_sig, aes(color = Direction), size = cfg_33$point_size, alpha = 0.95) +
    scale_color_manual(values = c("up" = cfg_33$col_up, "down" = cfg_33$col_down),
                       labels = c("up" = paste("Upregulated in", numerator), "down" = paste("Upregulated in", denominator)), name = NULL)
  df$Display <- ifelse(df$Direction != "ns", df$Display, "")
  p <- p + ggrepel::geom_text_repel(aes(label = Display), size = cfg_33$label_size, fontface = "bold.italic",
                                    family = cfg_33$font_family, color = "black", segment.size = 0.25, segment.alpha = 0.6,
                                    box.padding = 0.40, point.padding = 0.30, min.segment.length = 0, max.overlaps = Inf, force = 1.8)
  p + coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +
    theme_classic(base_family = cfg_33$font_family) +
    labs(title = contrast_label,
         subtitle = sprintf("PRM-MS panel (Wilcoxon, n=%d) | BHP<%.2f, |log2FC|>%.2f | DE up=%d  DE down=%d",
                            nrow(df), cfg_33$fdr_thresh, cfg_33$fc_thresh, n_up, n_down),
         x = build_arrow_xlab_33(numerator, denominator), y = expression(-log[10]("BH-adjusted p"))) +
    theme(plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 9, color = "grey30"),
          axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
          axis.title.y = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10, color = "black"),
          legend.position = "top", legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(color = "grey95", linetype = "dotted"))
}
fde <- supp_read("Supp1ab_Protein_Panel_DE.csv")
if (!is.null(fde)) {
  d_hl  <- fde[fde$Contrast == "CN-Hi vs CN-Lo", ]
  d_ah  <- fde[fde$Contrast == "AD vs CN-Hi", ]
  viz_save_both(build_volcano_33(d_hl, "CN-Hi vs CN-Lo", "CN-Hi", "CN-Lo") + .transparent_theme,
                "Supp1a_Protein_Volcano_HivsLo", cfg_33$plot_w, cfg_33$plot_h)
  viz_save_both(build_volcano_33(d_ah, "AD vs CN-Hi",    "AD",    "CN-Hi") + .transparent_theme,
                "Supp1b_Protein_Volcano_ADvsHi", cfg_33$plot_w, cfg_33$plot_h)
}

cov_col_fun <- colorRamp2(seq(-1, 1, length.out = 11), rev(RColorBrewer::brewer.pal(11, "RdBu")))
make_cov_hm <- function(csv, id_col, out_base, title_txt) {
  d <- supp_read(csv); if (is.null(d)) return(invisible(NULL))
  M <- as.matrix(d[, setdiff(names(d), id_col), drop = FALSE]); rownames(M) <- d[[id_col]]
  ht <- Heatmap(M, name = "Pearson r", col = cov_col_fun,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
                column_names_gp = gpar(fontsize = 9, fontfamily = "Arial", col = "red"),
                column_names_rot = 55,
                width = ncol(M) * unit(7, "mm"), height = nrow(M) * unit(7, "mm"),
                cell_fun = function(j, i, x, y, w, h, fill) grid::grid.text(sprintf("%.2f", M[i, j]), x, y,
                                                                            gp = grid::gpar(fontsize = 6, fontfamily = "Arial", col = ifelse(abs(M[i, j]) >= 0.6, "white", "black"))),
                heatmap_legend_param = list(title = "Pearson r",
                                            title_gp = gpar(fontfamily = "Arial", fontface = "bold"), labels_gp = gpar(fontfamily = "Arial")))
  viz_save_ht(ht, out_base, ncol(M) * (SET_SUPP$cov_tile_mm / 25.4) + SET_SUPP$cov_w_pad, nrow(M) * (SET_SUPP$cov_tile_mm / 25.4) + SET_SUPP$cov_h_pad, title_txt)
}
make_cov_hm("Supp1c_Protein_x_Covariates.csv", "Protein", "Supp1c_Protein_x_Covariates", "Protein")
make_cov_hm("Supp4c_mRNA_x_Covariates.csv",    "Gene",    "Supp4c_mRNA_x_Covariates",    "mRNA")

pca <- supp_read("Supp2a_PCA_points.csv")
if (!is.null(pca)) {
  pca_theme <- theme_classic(base_family = "Arial") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, color = "black"), legend.position = "right")
  plot_pca_grad <- function(col, lab, title_txt, out) {
    d <- pca[order(pca[[col]], na.last = FALSE), ]
    p <- ggplot(d, aes(PC1, PC2)) + geom_point(aes(color = .data[[col]]), size = 2.4, alpha = 0.80) +
      scale_color_gradientn(colours = c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"),
                            name = lab, na.value = "grey80") +
      labs(title = title_txt, subtitle = sprintf("All %d samples", nrow(pca)), x = "PC1", y = "PC2") + pca_theme
    viz_save_both(p, out, SET_SUPP$pca_w, SET_SUPP$pca_h)
  }
  if (any(is.finite(pca$LAMP2_expr))) plot_pca_grad("LAMP2_expr", "LAMP2\n(Log2 expr)", "PCA: LAMP2 Expression Gradient", "Supp2a_PCA_LAMP2")
  plot_pca_grad("NPTX2_expr", "NPTX2\n(Log2 expr)", "PCA: NPTX2 Expression Gradient", "Supp2a_PCA_NPTX2")
  prm <- pca$Has_PRM %in% c(TRUE, "TRUE")
  p_prm <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(color = "grey", size = 2.0, alpha = 0.6) +
    geom_point(data = pca[prm, ], aes(shape = "PRM-MS"), color = "red", size = 4.0, stroke = 0.8, fill = NA) +
    scale_shape_manual(values = c("PRM-MS" = 1), name = NULL) +
    labs(title = "PCA: PRM-MS Samples Highlighted",
         subtitle = sprintf("%d PRM-MS samples (red rings) among %d total", sum(prm), nrow(pca)), x = "PC1", y = "PC2") + pca_theme
  viz_save_both(p_prm, "Supp2a_PCA_PRM_overlay", SET_SUPP$pca_w, SET_SUPP$pca_h)
}

fg <- supp_read("Supp2bc_fGSEA_GOCC.csv")
if (!is.null(fg) && all(c("pathway", "NES", "padj", "size", "PC") %in% names(fg))) {
  pc_colors_26 <- c("PC1" = "#D55E00", "PC2" = "#0B7285")
  prep_top10_26 <- function(pc_label, direction) {
    d <- fg[fg$PC == pc_label & is.finite(fg$NES), ]
    d <- if (direction == "pos") d[d$NES > 0, ] else d[d$NES < 0, ]
    d <- head(d[order(d$padj), ], 10)
    if (nrow(d) == 0) return(d)
    d$Pretty <- vapply(d$pathway, function(x) {
      x <- sub("^GOCC_", "", x); x <- gsub("_", " ", x); x <- tools::toTitleCase(tolower(x))
      paste(strwrap(x, width = 30), collapse = "\n") }, character(1))
    d
  }
  plot_fgsea_compact <- function(df_sub, pc_name, direction, fill_color, out_base) {
    if (nrow(df_sub) == 0) { cat(sprintf("  [skip] no %s pathways for %s\n", direction, pc_name)); return(invisible(NULL)) }
    if (direction == "pos") { df_sub$Pretty <- factor(df_sub$Pretty, levels = df_sub$Pretty[order(df_sub$NES)]);  text_x <- 0.05;  text_hjust <- 0 }
    else                    { df_sub$Pretty <- factor(df_sub$Pretty, levels = df_sub$Pretty[order(-df_sub$NES)]); text_x <- -0.05; text_hjust <- 1 }
    dir_label <- ifelse(direction == "pos", "Positively", "Negatively")
    p <- ggplot(df_sub, aes(x = NES, y = Pretty)) +
      geom_col(fill = fill_color, color = "black", linewidth = 0.25, width = 0.5) +
      geom_text(aes(label = sprintf("BHP=%.1e (n=%d)", padj, size)), x = text_x, hjust = text_hjust,
                size = 2.8, family = "Arial", color = "white", fontface = "bold") +
      theme_classic(base_size = 10, base_family = "Arial") +
      labs(title = sprintf("%s: %s Enriched GO:CC", pc_name, dir_label),
           subtitle = "Top 10 by BHP | Large gene sets (min 50)", x = "Normalized Enrichment Score (NES)", y = NULL) +
      theme(plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9, color = "grey40"),
            axis.text.y = element_text(size = 8, face = "bold", color = "black", lineheight = 0.85),
            axis.text.x = element_text(size = 9, color = "black"), axis.title.x = element_text(size = 10, face = "bold"),
            panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"), plot.margin = margin(8, 12, 8, 8))
    viz_save_both(p, out_base, SET_SUPP$fgsea_w, SET_SUPP$fgsea_h)
  }
  plot_fgsea_compact(prep_top10_26("PC1", "pos"), "PC1", "pos", pc_colors_26["PC1"], "Supp2b_PC1_pos_CC_Identity")
  plot_fgsea_compact(prep_top10_26("PC1", "neg"), "PC1", "neg", pc_colors_26["PC1"], "Supp2b_PC1_neg_CC_Identity")
  plot_fgsea_compact(prep_top10_26("PC2", "pos"), "PC2", "pos", pc_colors_26["PC2"], "Supp2c_PC2_pos_CC_Identity")
  plot_fgsea_compact(prep_top10_26("PC2", "neg"), "PC2", "neg", pc_colors_26["PC2"], "Supp2c_PC2_neg_CC_Identity")
} else {
  cat("  [skip] Supp 2b/2c: fGSEA table not found (fgsea/GMT unavailable in 01_analysis).\n")
}

sex_colors <- c(Female = "#D81B60", Male = "#1E88E5")
build_density_plot <- function(df, title_txt, ks_obj, n_F, n_M) {
  med_F <- median(df$r[df$Sex == "Female"], na.rm = TRUE); med_M <- median(df$r[df$Sex == "Male"], na.rm = TRUE)
  subtitle_txt <- sprintf("n_F = %d | n_M = %d | KS D = %.3f, p = %.2g | median r: F = %.3f, M = %.3f",
                          n_F, n_M, ks_obj$statistic, ks_obj$p.value, med_F, med_M)
  ggplot(df, aes(x = r, color = Sex, fill = Sex)) +
    geom_density(linewidth = 0.9, alpha = 0.20) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = med_F, linetype = "dotted", color = sex_colors["Female"], linewidth = 0.6) +
    geom_vline(xintercept = med_M, linetype = "dotted", color = sex_colors["Male"],   linewidth = 0.6) +
    scale_color_manual(values = sex_colors) + scale_fill_manual(values = sex_colors) +
    scale_x_continuous(breaks = seq(-1, 1, 0.25), limits = c(-1, 1)) +
    theme_classic(base_size = 12) +
    labs(title = title_txt, subtitle = subtitle_txt, x = "Pearson r (gene mRNA vs NPTX2 mRNA)", y = "Density") +
    theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 10, color = "grey30"),
          axis.title = element_text(face = "bold", size = 12), axis.text = element_text(color = "black", size = 11),
          legend.position = "top", legend.title = element_blank())
}
build_facet_density <- function(df, title_txt = NULL, panel_tags = c("a", "b", "c", "d")) {
  grp_levels <- levels(df$Group)
  ann_all <- df %>% dplyr::distinct(Group, n_F, n_M, ks_D, ks_p, med_F, med_M) %>%
    dplyr::mutate(lbl = sprintf("n_F=%d  n_M=%d\nKS D=%.3f, p=%.2g\nmed_F=%.3f  med_M=%.3f", n_F, n_M, ks_D, ks_p, med_F, med_M))
  one_panel <- function(g) {
    dfg  <- df[df$Group == g, ]
    anng <- ann_all[ann_all$Group == g, ]
    ggplot(dfg, aes(x = r, color = Sex, fill = Sex)) +
      geom_density(linewidth = 0.9, alpha = 0.20) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_text(data = anng, aes(x = 0.98, y = Inf, label = lbl), inherit.aes = FALSE, hjust = 1, vjust = 1.15,
                size = 3, color = "grey25", lineheight = 0.95) +
      scale_color_manual(values = sex_colors) + scale_fill_manual(values = sex_colors) +
      scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1)) +
      theme_classic(base_size = 12) +
      labs(title = g, x = "Pearson r (gene mRNA vs NPTX2 mRNA)", y = "Density") +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            axis.title = element_text(face = "bold", size = 11), legend.title = element_blank()) +
      .transparent_theme
  }

  leg    <- cowplot::get_legend(one_panel(grp_levels[1]) + theme(legend.position = "top"))
  panels <- lapply(grp_levels, function(g) one_panel(g) + theme(legend.position = "none"))
  grid   <- cowplot::plot_grid(plotlist = panels, ncol = 2, align = "hv")
  cowplot::plot_grid(leg, grid, ncol = 1, rel_heights = c(0.06, 1))
}
sd_vals <- supp_read("Supp3_density_values.csv")
sd_ks   <- supp_read("Supp3_Density_KS.csv")
if (!is.null(sd_vals)) {
  ov <- sd_vals[sd_vals$Group == "Overall", ]
  if (nrow(ov) > 0) {
    ks_ov <- suppressWarnings(ks.test(ov$r[ov$Sex == "Female"], ov$r[ov$Sex == "Male"]))
    p_ov <- build_density_plot(ov, NULL,
                               ks_ov, sum(ov$Sex == "Female"), sum(ov$Sex == "Male")) + .transparent_theme
    viz_save_both(p_ov, "Supp3_SexDensity_Overall", SET_SUPP$sexdens_overall[1], SET_SUPP$sexdens_overall[2])
  }
  if (!is.null(sd_ks)) {
    pg  <- sd_vals[sd_vals$Group != "Overall", ]
    ann <- data.frame(Group = sd_ks$Group_Current, n_F = sd_ks$n_F, n_M = sd_ks$n_M,
                      ks_D = sd_ks$KS_D, ks_p = sd_ks$KS_p, med_F = sd_ks$median_F, med_M = sd_ks$median_M)
    pg  <- merge(pg, ann, by = "Group", all.x = TRUE)
    pg$Group <- factor(pg$Group, levels = c("CN-Lo", "CN-Hi", "MCI", "AD"))
    p_facet <- build_facet_density(pg, NULL) + .transparent_theme
    viz_save_both(p_facet, "Supp3_SexDensity_byGroup", SET_SUPP$sexdens_bygroup[1], SET_SUPP$sexdens_bygroup[2])
  }
}

sex_lab_all <- ifelse(grepl("^f", as.character(meta_aligned$Sex), ignore.case = TRUE), "Female",
                      ifelse(grepl("^m", as.character(meta_aligned$Sex), ignore.case = TRUE), "Male", NA_character_))
sex_meta_ej <- data.frame(SampleID = meta_aligned$SampleID, Sex_Label = sex_lab_all,
                          Group_Current = as.character(meta_aligned$Group_Current),
                          Age = suppressWarnings(as.numeric(meta_aligned$Age)), PMI = suppressWarnings(as.numeric(meta_aligned$PMI)),
                          RIN = suppressWarnings(as.numeric(meta_aligned$RIN)), CERAD = suppressWarnings(as.numeric(meta_aligned$SS_C)),
                          Braak = suppressWarnings(as.numeric(meta_aligned$SS_B)), stringsAsFactors = FALSE)
sex_meta_ej <- sex_meta_ej[sex_meta_ej$Sex_Label %in% c("Female", "Male") &
                             sex_meta_ej$Group_Current %in% c("CN-Lo", "CN-Hi", "MCI", "AD"), ]
gpos_ej <- c("Overall", "CN-Lo", "CN-Hi", "MCI", "AD")
sex_colors_ej <- c("Female" = "#CC79A7", "Male" = "#0072B2")
sig_from_padj <- function(p) dplyr::case_when(is.na(p) ~ "n/a", p < 1e-4 ~ sprintf("%.1e", p), TRUE ~ sprintf("%.3g", p))

cov_long_ej <- sex_meta_ej %>% dplyr::select(Sex_Label, Group_Current, Age, PMI, RIN) %>%
  tidyr::pivot_longer(c(Age, PMI, RIN), names_to = "Variable", values_to = "Value") %>% dplyr::filter(is.finite(Value))
cov_long_ej <- dplyr::bind_rows(dplyr::mutate(cov_long_ej, Group_Current = "Overall"), cov_long_ej) %>%
  dplyr::mutate(Group_Current = factor(Group_Current, levels = gpos_ej),
                Sex_Label = factor(Sex_Label, levels = c("Female", "Male")),
                Variable = factor(Variable, levels = c("Age", "PMI", "RIN")))
build_ordinal_ej <- function(value_col) {
  d <- sex_meta_ej %>% dplyr::select(Sex_Label, Group_Current, val = dplyr::all_of(value_col)) %>% dplyr::filter(is.finite(val))
  dplyr::bind_rows(dplyr::mutate(d, Group_Current = "Overall"), d) %>%
    dplyr::mutate(Group_Current = factor(Group_Current, levels = gpos_ej),
                  Sex_Label = factor(Sex_Label, levels = c("Female", "Male")),
                  val_factor = factor(val, levels = sort(unique(d$val))))
}
df_cerad_ej <- build_ordinal_ej("CERAD"); df_braak_ej <- build_ordinal_ej("Braak")

cor_ej <- supp_read("Supp3_density_values.csv")
if (!is.null(cor_ej)) cor_ej <- cor_ej %>%
  dplyr::transmute(Group_Current = factor(Group, levels = gpos_ej),
                   Sex_Label = factor(Sex, levels = c("Female", "Male")), r = r, abs_r = abs(r))

stats_cov_ej <- cov_long_ej %>% dplyr::group_by(Variable, Group_Current) %>%
  dplyr::summarise(n_F = sum(Sex_Label == "Female"), n_M = sum(Sex_Label == "Male"),
                   Wilcox_p = if (n_F >= 3 && n_M >= 3) suppressWarnings(stats::wilcox.test(Value[Sex_Label == "Female"], Value[Sex_Label == "Male"])$p.value) else NA_real_,
                   .groups = "drop") %>% dplyr::group_by(Variable) %>%
  dplyr::mutate(sig = sig_from_padj(p.adjust(Wilcox_p, method = "BH"))) %>% dplyr::ungroup()
ordinal_stats_ej <- function(df_all) df_all %>% dplyr::group_by(Group_Current) %>%
  dplyr::summarise(n_F = sum(Sex_Label == "Female"), n_M = sum(Sex_Label == "Male"),
                   Wilcox_p = if (n_F >= 3 && n_M >= 3) suppressWarnings(stats::wilcox.test(val[Sex_Label == "Female"], val[Sex_Label == "Male"])$p.value) else NA_real_,
                   .groups = "drop") %>% dplyr::mutate(sig = sig_from_padj(p.adjust(Wilcox_p, method = "BH")))
stats_cerad_ej <- ordinal_stats_ej(df_cerad_ej); stats_braak_ej <- ordinal_stats_ej(df_braak_ej)

make_cov_panel_ej <- function(var_name, y_lab) {
  df <- cov_long_ej %>% dplyr::filter(Variable == var_name)
  st <- stats_cov_ej %>% dplyr::filter(Variable == var_name)
  ymin <- min(df$Value, na.rm = TRUE); ymax <- max(df$Value, na.rm = TRUE); span <- ymax - ymin; y_pos <- ymax + 0.08 * span
  ggplot(df, aes(x = Group_Current, y = Value, fill = Sex_Label)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.5, position = position_dodge(width = 0.75)) +
    geom_point(aes(color = Sex_Label), position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.75), size = 0.5, alpha = 0.55, show.legend = FALSE) +
    geom_text(data = st, aes(x = Group_Current, y = y_pos, label = sig), inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_fill_manual(values = sex_colors_ej) + scale_color_manual(values = sex_colors_ej) +
    coord_cartesian(ylim = c(ymin, y_pos + 0.05 * span)) + theme_classic(base_size = 11) +
    labs(title = var_name, x = NULL, y = y_lab) +
    theme(plot.title = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 25, hjust = 1, size = 9.5, color = "black"), axis.text.y = element_text(size = 9, color = "black"),
          legend.position = "none")
}
make_corr_panel_ej <- function() {
  df <- cor_ej; st <- stats_cor_ej
  ymax <- max(df$abs_r, na.rm = TRUE); span <- ymax; y_pos <- ymax + 0.08 * span
  ggplot(df, aes(x = Group_Current, y = abs_r, fill = Sex_Label)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.65, position = position_dodge(width = 0.75)) +
    geom_text(data = st, aes(x = Group_Current, y = y_pos, label = sig), inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_fill_manual(values = sex_colors_ej) + coord_cartesian(ylim = c(0, y_pos + 0.05 * span)) + theme_classic(base_size = 11) +
    labs(title = "|r| correlation (gene mRNA vs NPTX2)", x = NULL, y = "|Pearson r|") +
    theme(plot.title = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 25, hjust = 1, size = 9.5, color = "black"), axis.text.y = element_text(size = 9, color = "black"),
          legend.position = "none")
}
make_stacked_panel_ej <- function(df_all, stats_df, title_text, palette_name) {
  cell_totals  <- df_all %>% dplyr::group_by(Group_Current, Sex_Label) %>% dplyr::summarise(total_n = dplyr::n(), .groups = "drop")
  stage_counts <- df_all %>% dplyr::group_by(Group_Current, Sex_Label, val_factor) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop")
  df_long <- stage_counts %>% dplyr::left_join(cell_totals, by = c("Group_Current", "Sex_Label")) %>%
    dplyr::mutate(prop = n / total_n) %>% dplyr::filter(total_n > 0)
  present <- df_long %>% dplyr::distinct(Group_Current, Sex_Label) %>% dplyr::arrange(Group_Current, Sex_Label) %>%
    dplyr::mutate(x_pos_label = paste(Group_Current, ifelse(Sex_Label == "Female", "F", "M")))
  df_long$x_pos <- factor(paste(df_long$Group_Current, ifelse(df_long$Sex_Label == "Female", "F", "M")), levels = present$x_pos_label)
  pg <- present %>% dplyr::group_by(Group_Current) %>% dplyr::summarise(n_sex = dplyr::n(), .groups = "drop") %>% dplyr::filter(n_sex == 2)
  centers <- present %>% dplyr::mutate(x_idx = dplyr::row_number()) %>% dplyr::filter(Group_Current %in% pg$Group_Current) %>%
    dplyr::group_by(Group_Current) %>% dplyr::summarise(x_center = mean(x_idx), .groups = "drop")
  st <- stats_df %>% dplyr::inner_join(centers, by = "Group_Current")
  n_stages <- length(levels(df_long$val_factor))
  pal <- if (n_stages <= 9) RColorBrewer::brewer.pal(max(3, n_stages), palette_name)[seq_len(n_stages)]
  else colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))(n_stages)
  names(pal) <- levels(df_long$val_factor)
  ggplot(df_long, aes(x = x_pos, y = prop, fill = val_factor)) +
    geom_col(position = "stack", width = 0.8, color = "white", linewidth = 0.25) +
    geom_text(data = st, aes(x = x_center, y = 1.10, label = sig), inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.18), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    scale_fill_manual(values = pal, name = "Stage", drop = FALSE) + theme_classic(base_size = 11) +
    labs(title = title_text, x = NULL, y = "Proportion within Sex") +
    theme(plot.title = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 10),
          axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1), axis.text.y = element_text(size = 9, color = "black"),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
}

if (!is.null(cor_ej)) {
  stats_cor_ej <- cor_ej %>% dplyr::group_by(Group_Current) %>%
    dplyr::summarise(Wilcox_p = suppressWarnings(stats::wilcox.test(abs_r[Sex_Label == "Female"], abs_r[Sex_Label == "Male"])$p.value), .groups = "drop") %>%
    dplyr::mutate(sig = sig_from_padj(p.adjust(Wilcox_p, method = "BH")))
  p_age   <- make_cov_panel_ej("Age", "Age (years)") + .transparent_theme
  p_pmi   <- make_cov_panel_ej("PMI", "PMI (h)")     + .transparent_theme
  p_cerad <- make_stacked_panel_ej(df_cerad_ej, stats_cerad_ej, "C score (CERAD)", "Blues") + .transparent_theme
  p_braak <- make_stacked_panel_ej(df_braak_ej, stats_braak_ej, "B score (Braak) ", "OrRd") + .transparent_theme
  p_rin   <- make_cov_panel_ej("RIN", "RIN") + .transparent_theme
  p_corr  <- make_corr_panel_ej() + .transparent_theme
  composite_grid <- cowplot::plot_grid(p_age, p_pmi, p_cerad, p_braak, p_rin, p_corr,
                                       ncol = 3, align = "vh", axis = "tblr", labels = c("e", "f", "g", "h", "i", "j"), label_size = 14, label_fontface = "bold")
  comp_title <- cowplot::ggdraw() + cowplot::draw_label("Sex differences across covariates and NPTX2 coupling", fontface = "bold", size = 14, hjust = 0.5)
  comp_sub   <- cowplot::ggdraw() + cowplot::draw_label("Female vs Male per group; two-sided Wilcoxon rank-sum test, BHP-adjusted within panel (numbers shown are the adjusted p). n per group on each panel. Panels g,h: ordinal stacked proportions; panel j: gene-level |r|.",
                                                        fontface = "italic", size = 9.0, color = "grey25", hjust = 0.5)
  composite_final <- cowplot::plot_grid(comp_title, comp_sub, composite_grid, ncol = 1, rel_heights = c(0.05, 0.04, 1))
  viz_save_both(composite_final, "Supp3_SexComposite_6Panel", SET_SUPP$composite[1], SET_SUPP$composite[2])
}

cfg_32 <- list(font_family = "Arial", fdr_thresh = 0.01, fc_thresh = 0.5, point_size = 1.8,
               col_ns = "grey70", col_up = "firebrick3", col_down = "steelblue3", label_size = 3,
               plot_w = 7, plot_h = 6)
build_arrow_xlab <- function(num, den) sprintf("log\u2082 FC  (%s \u2192 %s)", den, num)
build_volcano <- function(df, contrast_label, numerator, denominator) {
  df <- df %>% dplyr::mutate(
    Direction = dplyr::case_when(
      !is.na(adj.P.Val) & adj.P.Val < cfg_32$fdr_thresh & log2FC >  cfg_32$fc_thresh ~ "up",
      !is.na(adj.P.Val) & adj.P.Val < cfg_32$fdr_thresh & log2FC < -cfg_32$fc_thresh ~ "down",
      TRUE ~ "ns"),
    neglog10p = -log10(pmax(adj.P.Val, 1e-300)),
    Is_Label = Gene %in% top_labels)
  x_max <- max(c(max(abs(df$log2FC), na.rm = TRUE) * 1.05, 1.0))
  n_up <- sum(df$Direction == "up"); n_down <- sum(df$Direction == "down")
  p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
    geom_vline(xintercept = c(-cfg_32$fc_thresh, cfg_32$fc_thresh), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_hline(yintercept = -log10(cfg_32$fdr_thresh), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_point(data = df %>% dplyr::filter(Direction == "ns"), color = cfg_32$col_ns, size = cfg_32$point_size, alpha = 0.35)
  d_sig <- df %>% dplyr::filter(Direction != "ns")
  if (nrow(d_sig) > 0) p <- p + geom_point(data = d_sig, aes(color = Direction), size = cfg_32$point_size, alpha = 0.75) +
    scale_color_manual(values = c("up" = cfg_32$col_up, "down" = cfg_32$col_down),
                       labels = c("up" = paste("Upregulated in", numerator), "down" = paste("Upregulated in", denominator)), name = NULL)
  d_lab <- df %>% dplyr::filter(Is_Label)
  if (nrow(d_lab) > 0) p <- p + ggrepel::geom_text_repel(data = d_lab, aes(label = Gene),
                                                         size = cfg_32$label_size, fontface = "bold.italic", family = cfg_32$font_family, color = "black",
                                                         segment.size = 0.25, segment.alpha = 0.6, box.padding = 0.35, point.padding = 0.25,
                                                         min.segment.length = 0, max.overlaps = Inf, force = 1.5)
  p + coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, max(df$neglog10p, na.rm = TRUE) * 1.05 + 0.1)) +
    theme_classic(base_family = cfg_32$font_family) +
    labs(title = contrast_label,
         subtitle = sprintf("Protein-coding only (n=%d) | BHP<%.2f, |log2FC|>%.2f | DE up=%d  DE down=%d",
                            nrow(df), cfg_32$fdr_thresh, cfg_32$fc_thresh, n_up, n_down),
         x = build_arrow_xlab(numerator, denominator), y = expression(-log[10] * "(BH-adjusted p)")) +
    theme(plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 9, color = "grey30"),
          axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
          axis.title.y = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10, color = "black"),
          legend.position = "top", legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(color = "grey95", linetype = "dotted"))
}
de_hl <- supp_read("Supp4ab_DE_HivsLo.csv"); de_ah <- supp_read("Supp4ab_DE_ADvsHi.csv")
if (!is.null(de_hl) && !is.null(de_ah)) {
  top_labels   <- "NPTX2"
  viz_save_both(build_volcano(de_hl, "CN-Hi vs CN-Lo", "CN-Hi", "CN-Lo") + .transparent_theme, "Supp4a_Volcano_HivsLo", cfg_32$plot_w, cfg_32$plot_h)
  viz_save_both(build_volcano(de_ah, "AD vs CN-Hi",    "AD",    "CN-Hi") + .transparent_theme, "Supp4b_Volcano_ADvsHi", cfg_32$plot_w, cfg_32$plot_h)
}
cat("\nSupplementary figures ported.\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  cat("  [skip] install.packages('openxlsx') to build the summary workbook.\n")
} else {
  meta_full <- as.data.frame(readr::read_csv(file.path(exp_dir, "metadata_master.csv"), show_col_types = FALSE))

  grp6_levels <- c("CN-Lo", "CN-Hi", "YoungCon", "MCI", "AD", "Others")
  g_cur <- as.character(meta_full$Group_Current)
  dx    <- as.character(meta_full$Diagnosis_COMP)
  meta_full$Group6 <- factor(
    ifelse(g_cur %in% c("CN-Lo", "CN-Hi", "MCI", "AD"), g_cur,
           ifelse(dx %in% "YoungCon", "YoungCon", "Others")), levels = grp6_levels)

  tab <- addmargins(table(Source = meta_full$Source, Group = meta_full$Group6))
  tab_counts <- as.data.frame.matrix(tab)
  names(tab_counts)[names(tab_counts) == "Sum"] <- "Total"
  tab_counts <- cbind(Source = rownames(tab_counts), tab_counts, row.names = NULL)
  tab_counts$Source[tab_counts$Source == "Sum"] <- "TOTAL"

  id_cols  <- c("SampleID", "Source_ID", "TGID", "Source", "Group_Current", "Group6", "PC1", "PC2")
  num_cols <- setdiff(names(meta_full)[vapply(meta_full, is.numeric, logical(1))], id_cols)
  grp_list <- c(list(Overall = rep(TRUE, nrow(meta_full))),
                setNames(lapply(grp6_levels, function(g) !is.na(meta_full$Group6) & meta_full$Group6 == g), grp6_levels))
  summ_one <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  if (n == 0) return(c(N = 0, Mean = NA, SD = NA, Median = NA, Min = NA, Max = NA))
  c(N = n, Mean = mean(x), SD = stats::sd(x), Median = stats::median(x), Min = min(x), Max = max(x)) }
  cont_rows <- list()
  for (v in num_cols) for (gname in names(grp_list)) {
    s <- summ_one(meta_full[[v]][grp_list[[gname]]])
    cont_rows[[length(cont_rows) + 1]] <- data.frame(
      Variable = v, Group = gname, N = s[["N"]],
      Mean = round(s[["Mean"]], 3), SD = round(s[["SD"]], 3), Median = round(s[["Median"]], 3),
      Min = round(s[["Min"]], 3), Max = round(s[["Max"]], 3),
      MeanSD = ifelse(is.na(s[["Mean"]]), NA, sprintf("%.2f \u00b1 %.2f", s[["Mean"]], s[["SD"]])),
      `Range(min-max)` = ifelse(is.na(s[["Min"]]), NA, sprintf("%.2f \u2013 %.2f", s[["Min"]], s[["Max"]])),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  summary_cont <- do.call(rbind, cont_rows)
  names(summary_cont)[names(summary_cont) == "MeanSD"] <- "Mean\u00b1SD"

  cat_cols <- intersect(c("Sex", "APOE", "Diagnosis_COMP", "PathDx", "LivingDx",
                          "NrtPlaquesCERAD_COMPhl", "TauBraakNFT_COMPhl", "In_PRM_MS"), names(meta_full))
  cat_rows <- list()
  for (v in cat_cols) {
    lev <- as.character(meta_full[[v]]); lev[is.na(lev)] <- "(NA)"
    ct  <- as.data.frame.matrix(table(Level = lev, Group = meta_full$Group6))
    ct$Overall <- rowSums(ct)
    cat_rows[[v]] <- cbind(Variable = v, Level = rownames(ct), ct, row.names = NULL)
  }
  summary_cat <- do.call(rbind, cat_rows)

  def_map <- c(
    SampleID = "Unique RNA-seq sample identifier (analysis primary key).",
    Source_ID = "Original subject/sample identifier from the source cohort.",
    TGID = "Source cohort genomic/transcriptomic identifier.",
    Source = "Cohort / brain bank of origin.",
    Group_Current = "Analysis group: CN-Lo, CN-Hi, MCI, AD, or Others.",
    Diagnosis_COMP = "Composite diagnosis (Control, MCI, AD, YoungCon, OtherPath, DementNoAD, NA).",
    PathDx = "Neuropathological diagnosis.",
    LivingDx = "Clinical (antemortem) diagnosis.",
    Sex = "Biological sex.",
    APOE = "APOE genotype.",
    Age = "Age at death (years).",
    PMI = "Post-mortem interval (hours).",
    RIN = "RNA integrity number.",
    NrtPlaquesCERAD_COMPhl = "CERAD neuritic-plaque burden dichotomized (CERAD.Hi vs CERAD.Lo).",
    SS_C = "NIA-AA CERAD 'C' score (0-3); C>=2 = high plaque burden.",
    TauBraakNFT_COMPhl = "Braak NFT stage dichotomized (Braak.Hi vs Braak.Lo).",
    SS_B = "NIA-AA Braak 'B' score (0-3); B>=2 = Braak stage III-VI.",
    In_PRM_MS = "TRUE if the sample is in the PRM-MS targeted-proteomics subset.")
  base_cols <- setdiff(names(meta_full), c("Group6", "PC1", "PC2"))
  defs <- data.frame(
    Column = base_cols,
    Definition = ifelse(base_cols %in% names(def_map), def_map[base_cols], "(definition - please fill in)"),
    Type = vapply(meta_full[base_cols], function(c) class(c)[1], character(1)),
    N_unique  = vapply(meta_full[base_cols], function(c) length(unique(c)), integer(1)),
    N_missing = vapply(meta_full[base_cols], function(c) sum(is.na(c)), integer(1)),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

  out_xlsx <- file.path(fig_dir, "Cohort_Metadata_Summary.xlsx")
  openxlsx::write.xlsx(
    list(Group_x_Source      = tab_counts,
         Continuous_Summary  = summary_cont,
         Categorical_Summary = summary_cat,
         Column_Definitions  = defs,
         Metadata            = meta_full[, base_cols, drop = FALSE]),
    file = out_xlsx, overwrite = TRUE, asTable = TRUE)
  cat(sprintf("\nWrote cohort + metadata summary workbook -> %s\n", out_xlsx))
}
