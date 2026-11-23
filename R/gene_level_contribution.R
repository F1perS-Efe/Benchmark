#' Leave-one-out gene contribution (delta-AUC) and immune-axis annotation
#'
#' Evaluates each matched signature gene across cohorts by computing the full
#' signature AUC (`full_auc`), the leave-one-out AUC after removing that gene
#' (`loo_auc`), `delta_auc = full_auc - loo_auc` (`> 0` indicates a positive
#' driver gene whose removal degrades discrimination; `< 0` indicates noise),
#' and `single_gene_auc`. Attaches a curated `immune_axis` annotation.
#'
#' @param datasets Named list of SigBenchR datasets (`validate_dataset()`).
#' @param signature A signature definition (character vector, named sign
#'   vector, or list with `genes` / `genes_up` / `genes_down`).
#' @param method Scoring method passed to [score_signature()]. Default `"zscore"`.
#' @param sd_handling Response policy (`"exclude"`, `"as_0"`, `"as_1"`; `"na"`
#'   is accepted as an alias for `"exclude"`).
#' @param direction Default ROC direction (`"<"` or `">"`). Overridden by
#'   `signature$direction` when present.
#' @param apply_log2 Passed to [prepare_expression()]. Default `TRUE`.
#' @param min_genes Minimum matched genes required in the full signature per
#'   cohort. Default `2L`.
#'
#' @return A list with `per_cohort` (`data.frame`) and `summary` (`data.frame`
#'   aggregated across cohorts with `mean_delta_auc`, `median_delta_auc`,
#'   `sd_delta_auc`, `n_cohorts`, `immune_axis`).
#' @export
analyze_gene_contributions <- function(datasets,
                                       signature,
                                       method = "zscore",
                                       sd_handling = "exclude",
                                       direction = "<",
                                       apply_log2 = TRUE,
                                       min_genes = 2L) {
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` must be a named list.", call. = FALSE)
  }
  if (identical(tolower(as.character(sd_handling)[1]), "na")) {
    sd_handling <- "exclude"
  }
  sig <- .as_signature(signature)
  sig_dir <- .normalize_roc_direction(sig$direction, default = direction)
  parsed <- .parse_signature_genes(sig$genes, sig$genes_up, sig$genes_down, weights = sig$weights)

  rows <- list()
  row_i <- 0L

  for (ds_name in names(datasets)) {
    prep <- .prepare_cohort(
      datasets[[ds_name]],
      name = ds_name,
      apply_log2 = apply_log2,
      sd_handling = sd_handling,
      filter_response = TRUE
    )
    expr_full <- prep$expr
    y_kept <- prep$y_kept

    eff_min_genes <- .resolve_min_genes(sig, default = min_genes)
    full_scored <- score_signature(
      expr_full,
      genes = sig$genes,
      genes_up = sig$genes_up,
      genes_down = sig$genes_down,
      weights = sig$weights,
      method = method,
      min_genes = eff_min_genes
    )
    if (isTRUE(full_scored$skipped) || is.null(full_scored$scores)) {
      next
    }

    s_ids <- intersect(names(full_scored$scores), names(y_kept))
    resp_full <- y_kept[s_ids]
    sc_full <- full_scored$scores[s_ids]
    keep_full <- is.finite(sc_full) & !is.na(resp_full)
    if (length(unique(resp_full[keep_full])) < 2L) {
      next
    }
    full_auc <- auc_mw(
      orient_score(sc_full[keep_full], direction = sig_dir),
      resp_full[keep_full]
    )
    if (!is.finite(full_auc)) {
      next
    }

    cov_up <- if (length(parsed$genes_up) > 0L) match_genes(expr_full, parsed$genes_up)$found else character(0)
    cov_dn <- if (length(parsed$genes_down) > 0L) match_genes(expr_full, parsed$genes_down)$found else character(0)

    arm_genes <- c(cov_up, cov_dn)
    arm_labels <- c(rep("up", length(cov_up)), rep("down", length(cov_dn)))
    if (length(arm_genes) < 2L) {
      next
    }

    w_all <- if (isTRUE(parsed$weighted)) {
      .resolve_found_weights(arm_genes, parsed$weights, full_scored$coverage$alias_rescued)
    } else {
      NULL
    }

    for (idx in seq_along(arm_genes)) {
      g <- arm_genes[idx]
      arm <- arm_labels[idx]
      up_rem <- if (arm == "up") setdiff(cov_up, g) else cov_up
      dn_rem <- if (arm == "down") setdiff(cov_dn, g) else cov_dn

      loo_scored <- if (isTRUE(parsed$weighted)) {
        rem_all <- c(up_rem, dn_rem)
        if (length(rem_all) > 0L) {
          score_signature(expr_full, weights = w_all[rem_all], method = method, min_genes = 1L)
        } else {
          list(skipped = TRUE)
        }
      } else if (length(up_rem) > 0L && length(dn_rem) > 0L) {
        score_signature(expr_full, genes_up = up_rem, genes_down = dn_rem, method = method, min_genes = 1L)
      } else if (length(up_rem) > 0L) {
        score_signature(expr_full, genes = up_rem, method = method, min_genes = 1L)
      } else if (length(dn_rem) > 0L) {
        score_signature(expr_full, genes_up = character(0), genes_down = dn_rem, method = method, min_genes = 1L)
      } else {
        list(skipped = TRUE)
      }

      loo_auc <- NA_real_
      if (!isTRUE(loo_scored$skipped) && !is.null(loo_scored$scores)) {
        l_ids <- intersect(names(loo_scored$scores), names(y_kept))
        resp_l <- y_kept[l_ids]
        sc_l <- loo_scored$scores[l_ids]
        keep_l <- is.finite(sc_l) & !is.na(resp_l)
        if (length(unique(resp_l[keep_l])) >= 2L) {
          loo_auc <- auc_mw(
            orient_score(sc_l[keep_l], direction = sig_dir),
            resp_l[keep_l]
          )
        }
      }

      g_vals <- as.numeric(expr_full[g, names(y_kept)])
      if (identical(arm, "down")) {
        g_vals <- -g_vals
      }
      keep_g <- is.finite(g_vals) & !is.na(y_kept)
      single_auc <- if (length(unique(y_kept[keep_g])) >= 2L) {
        auc_mw(orient_score(g_vals[keep_g], direction = sig_dir), y_kept[keep_g])
      } else {
        NA_real_
      }

      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        cohort = ds_name,
        gene = g,
        arm = arm,
        immune_axis = .annotate_immune_axis(g),
        full_auc = full_auc,
        loo_auc = loo_auc,
        delta_auc = if (is.finite(loo_auc)) full_auc - loo_auc else NA_real_,
        single_gene_auc = single_auc,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    empty_pc <- data.frame(
      cohort = character(), gene = character(), arm = character(),
      immune_axis = character(), full_auc = numeric(), loo_auc = numeric(),
      delta_auc = numeric(), single_gene_auc = numeric(),
      stringsAsFactors = FALSE
    )
    empty_sm <- data.frame(
      gene = character(), arm = character(), immune_axis = character(),
      n_cohorts = integer(), mean_delta_auc = numeric(),
      median_delta_auc = numeric(), sd_delta_auc = numeric(),
      mean_single_gene_auc = numeric(), stringsAsFactors = FALSE
    )
    return(list(per_cohort = empty_pc, summary = empty_sm))
  }

  pc <- do.call(rbind, rows)
  rownames(pc) <- NULL

  grp_key <- paste(pc$gene, pc$arm, sep = "\r")
  groups <- split(pc, grp_key)
  sm_rows <- lapply(groups, function(df_g) {
    d <- df_g$delta_auc[is.finite(df_g$delta_auc)]
    sg <- df_g$single_gene_auc[is.finite(df_g$single_gene_auc)]
    data.frame(
      gene = df_g$gene[1],
      arm = df_g$arm[1],
      immune_axis = df_g$immune_axis[1],
      n_cohorts = as.integer(length(d)),
      mean_delta_auc = if (length(d) > 0L) mean(d) else NA_real_,
      median_delta_auc = if (length(d) > 0L) stats::median(d) else NA_real_,
      sd_delta_auc = if (length(d) > 1L) stats::sd(d) else NA_real_,
      mean_single_gene_auc = if (length(sg) > 0L) mean(sg) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  sm <- do.call(rbind, sm_rows)
  sm <- sm[order(!is.finite(sm$mean_delta_auc), -sm$mean_delta_auc, sm$gene), , drop = FALSE]
  rownames(sm) <- NULL

  list(per_cohort = pc, summary = sm)
}

#' Plot leave-one-out gene contributions (delta-AUC lollipop chart)
#'
#' @param contrib Output of [analyze_gene_contributions()] or its `summary`
#'   `data.frame`.
#' @param top_n Maximum number of genes to display (ordered by absolute
#'   `mean_delta_auc`). Default `25L`.
#'
#' @return A `ggplot` object.
#' @export
plot_gene_contributions <- function(contrib, top_n = 25L) {
  df <- if (is.list(contrib) && !is.data.frame(contrib) && "summary" %in% names(contrib)) {
    contrib$summary
  } else {
    contrib
  }
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("No gene contribution rows available to plot.", call. = FALSE)
  }
  df <- df[is.finite(df$mean_delta_auc), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No finite mean_delta_auc values to plot.", call. = FALSE)
  }
  if (nrow(df) > top_n) {
    df <- df[order(-abs(df$mean_delta_auc)), , drop = FALSE]
    df <- df[seq_len(top_n), , drop = FALSE]
  }
  df$direction_effect <- ifelse(df$mean_delta_auc > 0, "Driver (delta > 0)", "Noise / Redundant (delta <= 0)")
  df$gene_label <- sprintf("%s [%s]", df$gene, df$immune_axis)
  df$gene_label <- factor(df$gene_label, levels = df$gene_label[order(df$mean_delta_auc)])

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$mean_delta_auc, y = .data$gene_label, color = .data$direction_effect)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$mean_delta_auc, y = .data$gene_label, yend = .data$gene_label),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = c("Driver (delta > 0)" = "#0F3D5E", "Noise / Redundant (delta <= 0)" = "#B91C1C"),
      name = "Contribution"
    ) +
    ggplot2::labs(
      title = "Leave-One-Out Gene Contribution (Mean delta-AUC)",
      x = "Mean delta-AUC (full_auc - loo_auc)",
      y = NULL,
      caption = "Positive delta-AUC: removing gene lowers AUC (signal driver). Negative: removing gene improves AUC."
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Curated immuno-oncology gene-to-axis annotation lookup
#'
#' @noRd
.annotate_immune_axis <- function(gene) {
  g <- toupper(trimws(as.character(gene)[1]))
  axes <- list(
    IFN_gamma_JAK_STAT = c("IFNG", "STAT1", "IRF1", "JAK1", "JAK2", "IDO1", "GBP1", "GBP5", "SOCS1", "ISG15"),
    Cytolytic_Effector = c("GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "GNLY", "FASLG", "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "IL2RG", "TBX21", "EOMES"),
    Chemokine_Recruitment = c("CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CCR5", "CXCR3", "CXCR6"),
    Immune_Checkpoint = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR", "CD276", "VTCN1", "BTLA", "ICOS", "ICOSLG", "TNFRSF4", "TNFSF4", "TNFRSF9", "TNFSF9", "TNFRSF18", "CD27", "CD28", "CD80", "CD86", "CD40", "CD40LG"),
    Antigen_Presentation = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-E", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "PSMB10", "CIITA", "NLRC5"),
    TGFb_Stromal = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "FAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "VEGFA", "MMP2", "MMP9"),
    Myeloid_Suppression = c("ARG1", "ARG2", "IL10", "CD68", "CD163", "MRC1", "CSF1R", "CCR2", "S100A8", "S100A9", "NOS2", "PTGS2")
  )
  for (ax in names(axes)) {
    if (g %in% axes[[ax]]) {
      return(ax)
    }
  }
  "Other"
}
