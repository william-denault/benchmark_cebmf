library(flashier)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(data.table)
library(ica)
library(softImpute)


#------------------------------------------------------------
# Configuration
#------------------------------------------------------------

sim_path <- "./cEBMF_sim_and_fit"
plot_path <- file.path(sim_path, "plots")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

n_simulations <- 100L
simulation_ids <- 0:(n_simulations - 1L)

max_oracle_rank <- 30L
softimpute_lambda <- 1
ica_maxit <- 1000L
max_components_to_plot <- 10L


#------------------------------------------------------------
# General helpers
#------------------------------------------------------------

read_numeric_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }

  out <- as.matrix(
    data.table::fread(
      path,
      header = FALSE,
      sep = "\t"
    )
  )

  if (!is.numeric(out) || anyNA(out)) {
    stop("Matrix must be numeric and contain no missing values: ", path)
  }

  out
}


validate_pair <- function(Z, X_true) {
  stopifnot(identical(dim(Z), dim(X_true)))
  stopifnot(is.numeric(Z), is.numeric(X_true))
  stopifnot(!anyNA(Z), !anyNA(X_true))
}


resolve_centered_rank <- function(max_rank, n, p) {
  feasible_rank <- max(0L, min(n - 1L, p))

  if (is.null(max_rank)) {
    return(feasible_rank)
  }

  if (
    length(max_rank) != 1L ||
    is.na(max_rank) ||
    !is.finite(max_rank) ||
    max_rank < 0 ||
    max_rank != as.integer(max_rank)
  ) {
    stop("max_rank must be one nonnegative integer")
  }

  min(as.integer(max_rank), feasible_rank)
}


column_mean_matrix <- function(center, n) {
  matrix(
    center,
    nrow = n,
    ncol = length(center),
    byrow = TRUE
  )
}


add_column_center <- function(X_centered, center) {
  sweep(X_centered, 2, center, FUN = "+")
}


matrix_rmse <- function(X_hat, X_true) {
  stopifnot(identical(dim(X_hat), dim(X_true)))
  sqrt(mean((X_hat - X_true)^2))
}


reconstruction_rmse <- function(fit, X_true) {
  if (!is.null(fit$X_hat)) {
    X_hat <- fit$X_hat
  } else {
    X_hat <- fit$L_pm %*% t(fit$F_pm)
  }

  matrix_rmse(X_hat, X_true)
}


weighted_svd_reconstruction <- function(fit) {
  # softImpute drops matrix dimensions from u and v for a rank-one fit.
  # Restore the documented matrix representation before reconstruction.
  u <- fit$u
  v <- fit$v
  if (is.null(dim(u))) {
    u <- matrix(u, ncol = 1L)
  }
  if (is.null(dim(v))) {
    v <- matrix(v, ncol = 1L)
  }

  keep <- which(fit$d > 0)

  if (length(keep) == 0L) {
    return(
      matrix(
        0,
        nrow = nrow(u),
        ncol = nrow(v)
      )
    )
  }

  scores <- sweep(
    u[, keep, drop = FALSE],
    2,
    fit$d[keep],
    FUN = "*"
  )

  tcrossprod(scores, v[, keep, drop = FALSE])
}


#------------------------------------------------------------
# FLASH
#------------------------------------------------------------

fit_flash_model <- function(Z, factor_prior) {
  flash(
    Z,
    greedy_Kmax = 10,
    ebnm_fn = c(
      ebnm_generalized_binary,
      factor_prior
    ),
    backfit = TRUE,
    verbose = 0
  )
}


#------------------------------------------------------------
# Oracle PCA
#------------------------------------------------------------

fit_oracle_pca <- function(Z, X_true, max_rank = 30L) {
  validate_pair(Z, X_true)

  n <- nrow(Z)
  p <- ncol(Z)
  max_rank <- resolve_centered_rank(max_rank, n, p)

  center <- colMeans(Z)
  mean_fit <- column_mean_matrix(center, n)
  ranks <- 0:max_rank

  if (max_rank == 0L) {
    current_rmse <- matrix_rmse(mean_fit, X_true)

    return(
      list(
        method = "PCA (oracle)",
        L_pm = matrix(numeric(0), nrow = n, ncol = 0L),
        F_pm = matrix(numeric(0), nrow = p, ncol = 0L),
        center = center,
        X_hat = mean_fit,
        selected_rank = 0L,
        effective_rank = 0L,
        rmse = current_rmse,
        rmse_by_rank = data.frame(
          rank = 0L,
          effective_rank = 0L,
          rmse = current_rmse,
          converged = NA
        ),
        fit = NULL
      )
    )
  }

  # prcomp uses an SVD of the centered, unscaled matrix. Its x field contains
  # component scores and rotation contains the variable loadings.
  pca_fit <- stats::prcomp(
    Z,
    center = TRUE,
    scale. = FALSE,
    rank. = max_rank,
    retx = TRUE
  )

  reconstruct_rank <- function(k) {
    if (k == 0L) {
      return(mean_fit)
    }

    selected <- seq_len(k)
    X_centered_hat <- tcrossprod(
      pca_fit$x[, selected, drop = FALSE],
      pca_fit$rotation[, selected, drop = FALSE]
    )
    add_column_center(X_centered_hat, center)
  }

  rmse_by_rank <- vapply(
    ranks,
    function(k) matrix_rmse(reconstruct_rank(k), X_true),
    numeric(1)
  )

  best_index <- which.min(rmse_by_rank)
  selected_rank <- ranks[best_index]
  X_hat <- reconstruct_rank(selected_rank)

  if (selected_rank == 0L) {
    scores <- matrix(numeric(0), nrow = n, ncol = 0L)
    loadings <- matrix(numeric(0), nrow = p, ncol = 0L)
  } else {
    selected <- seq_len(selected_rank)
    scores <- pca_fit$x[, selected, drop = FALSE]
    loadings <- pca_fit$rotation[, selected, drop = FALSE]
  }

  list(
    method = "PCA (oracle)",
    L_pm = scores,
    F_pm = loadings,
    center = center,
    X_hat = X_hat,
    selected_rank = selected_rank,
    effective_rank = selected_rank,
    rmse = rmse_by_rank[best_index],
    rmse_by_rank = data.frame(
      rank = ranks,
      effective_rank = ranks,
      rmse = rmse_by_rank,
      converged = NA
    ),
    fit = pca_fit
  )
}


#------------------------------------------------------------
# Oracle ICA
#------------------------------------------------------------

fit_oracle_ica <- function(
    Z,
    X_true,
    max_rank = 30L,
    seed = 1L,
    maxit = 1000L
) {
  validate_pair(Z, X_true)

  n <- nrow(Z)
  p <- ncol(Z)
  max_rank <- resolve_centered_rank(max_rank, n, p)

  center <- colMeans(Z)
  Z_centered <- sweep(Z, 2, center, FUN = "-")
  mean_fit <- column_mean_matrix(center, n)
  ranks <- 0:max_rank

  rmse_by_rank <- rep(NA_real_, length(ranks))
  converged_by_rank <- rep(NA, length(ranks))
  rmse_by_rank[1L] <- matrix_rmse(mean_fit, X_true)

  best_rmse <- rmse_by_rank[1L]
  selected_rank <- 0L
  best_fit <- NULL
  best_mixing <- NULL
  best_X_hat <- mean_fit

  if (max_rank > 0L) {
    for (k in seq_len(max_rank)) {
      set.seed(seed + k)

      initial_rotation <- qr.Q(
        qr(matrix(rnorm(k^2), nrow = k, ncol = k))
      )

      current_fit <- tryCatch(
        ica::ica(
          Z_centered,
          nc = k,
          method = "fast",
          center = FALSE,
          maxit = maxit,
          Rmat = initial_rotation
        ),
        error = function(error) {
          warning(
            "ICA failed at rank ", k,
            ": ", conditionMessage(error)
          )
          NULL
        }
      )

      if (is.null(current_fit)) {
        next
      }

      # In ica 1.0-3, the rank-one FastICA special case returns M as 1 x p,
      # whereas ranks above one return the documented p x k mixing matrix.
      # Normalize that special case before applying X = tcrossprod(S, M).
      mixing <- current_fit$M
      if (identical(dim(mixing), c(k, p))) {
        mixing <- t(mixing)
      }
      if (!identical(dim(mixing), c(p, k))) {
        stop(
          "Unexpected ICA mixing-matrix dimensions at rank ", k,
          ": ", paste(dim(current_fit$M), collapse = " x ")
        )
      }

      # The documented ICA model is X = tcrossprod(S, M) + E.
      X_centered_hat <- tcrossprod(current_fit$S, mixing)
      X_hat <- add_column_center(X_centered_hat, center)
      current_rmse <- matrix_rmse(X_hat, X_true)

      curve_index <- match(k, ranks)
      rmse_by_rank[curve_index] <- current_rmse
      converged_by_rank[curve_index] <- current_fit$converged

      if (current_rmse < best_rmse) {
        best_rmse <- current_rmse
        selected_rank <- k
        best_fit <- current_fit
        best_mixing <- mixing
        best_X_hat <- X_hat
      }
    }
  }

  if (selected_rank == 0L) {
    scores <- matrix(numeric(0), nrow = n, ncol = 0L)
    mixing <- matrix(numeric(0), nrow = p, ncol = 0L)
    selected_converged <- NA
  } else {
    scores <- best_fit$S
    mixing <- best_mixing
    selected_converged <- best_fit$converged
  }

  if (isFALSE(selected_converged)) {
    warning(
      "Selected ICA fit did not converge for seed ", seed,
      " at rank ", selected_rank
    )
  }

  list(
    method = "ICA (oracle)",
    L_pm = scores,
    F_pm = mixing,
    center = center,
    X_hat = best_X_hat,
    selected_rank = selected_rank,
    effective_rank = selected_rank,
    rmse = best_rmse,
    rmse_by_rank = data.frame(
      rank = ranks,
      effective_rank = ranks,
      rmse = rmse_by_rank,
      converged = converged_by_rank
    ),
    fit = best_fit,
    converged = selected_converged
  )
}


#------------------------------------------------------------
# Oracle-rank SoftImpute
#------------------------------------------------------------

fit_oracle_softimpute <- function(
    Z,
    X_true,
    max_rank = 30L,
    lambda = 1.1
) {
  validate_pair(Z, X_true)

  if (
    length(lambda) != 1L ||
    is.na(lambda) ||
    !is.finite(lambda) ||
    lambda < 0
  ) {
    stop("lambda must be one finite nonnegative number")
  }

  n <- nrow(Z)
  p <- ncol(Z)
  max_rank <- resolve_centered_rank(max_rank, n, p)

  center <- colMeans(Z)
  Z_centered <- sweep(Z, 2, center, FUN = "-")
  mean_fit <- column_mean_matrix(center, n)
  ranks <- 0:max_rank

  rmse_by_rank <- rep(NA_real_, length(ranks))
  effective_rank_by_rank <- integer(length(ranks))
  rmse_by_rank[1L] <- matrix_rmse(mean_fit, X_true)

  best_rmse <- rmse_by_rank[1L]
  selected_rank <- 0L
  best_effective_rank <- 0L
  best_fit <- NULL
  best_X_hat <- mean_fit

  if (max_rank > 0L) {
    for (k in seq_len(max_rank)) {
      current_fit <- softImpute::softImpute(
        Z_centered,
        rank.max = k,
        lambda = lambda,
        type = "svd",
        thresh = 1e-7,
        maxit = 1000L,
        trace.it = FALSE
      )

      keep <- which(current_fit$d > 0)
      effective_rank <- length(keep)
      X_centered_hat <- weighted_svd_reconstruction(current_fit)
      X_hat <- add_column_center(X_centered_hat, center)
      current_rmse <- matrix_rmse(X_hat, X_true)

      curve_index <- match(k, ranks)
      rmse_by_rank[curve_index] <- current_rmse
      effective_rank_by_rank[curve_index] <- effective_rank

      if (current_rmse < best_rmse) {
        best_rmse <- current_rmse
        selected_rank <- k
        best_effective_rank <- effective_rank
        best_fit <- current_fit
        best_X_hat <- X_hat
      }
    }
  }

  if (selected_rank == 0L || best_effective_rank == 0L) {
    scores <- matrix(numeric(0), nrow = n, ncol = 0L)
    loadings <- matrix(numeric(0), nrow = p, ncol = 0L)
  } else {
    best_u <- best_fit$u
    best_v <- best_fit$v
    if (is.null(dim(best_u))) {
      best_u <- matrix(best_u, ncol = 1L)
    }
    if (is.null(dim(best_v))) {
      best_v <- matrix(best_v, ncol = 1L)
    }

    keep <- which(best_fit$d > 0)
    scores <- sweep(
      best_u[, keep, drop = FALSE],
      2,
      best_fit$d[keep],
      FUN = "*"
    )
    loadings <- best_v[, keep, drop = FALSE]
  }

  list(
    method = "SoftImpute (oracle)",
    L_pm = scores,
    F_pm = loadings,
    center = center,
    X_hat = best_X_hat,
    selected_rank = selected_rank,
    effective_rank = best_effective_rank,
    rmse = best_rmse,
    rmse_by_rank = data.frame(
      rank = ranks,
      effective_rank = effective_rank_by_rank,
      rmse = rmse_by_rank,
      converged = NA
    ),
    fit = best_fit,
    lambda = lambda
  )
}


#------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------

empty_component_plot <- function(title, subtitle) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Column-mean reconstruction"
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void() +
    labs(
      title = title,
      subtitle = subtitle
    )
}


make_component_heatmap <- function(
    scores,
    simulation,
    method_name,
    rmse,
    rank_text,
    max_components = 10L
) {
  title <- paste0("Simulation ", simulation, ": ", method_name)
  subtitle <- paste0(
    rank_text,
    "; RMSE = ", format(round(rmse, 4), nsmall = 4)
  )

  if (is.null(scores) || ncol(scores) == 0L) {
    return(empty_component_plot(title, subtitle))
  }

  shown <- seq_len(min(ncol(scores), max_components))
  score_df <- reshape2::melt(scores[, shown, drop = FALSE])

  if (ncol(scores) > length(shown)) {
    subtitle <- paste0(
      subtitle,
      "; showing first ", length(shown), " components"
    )
  }

  ggplot(
    score_df,
    aes(
      x = Var2,
      y = Var1,
      fill = value
    )
  ) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_y_reverse() +
    theme_minimal(base_size = 9) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Component",
      y = "Sample",
      fill = "Score"
    ) +
    theme(
      plot.title = element_text(size = 11),
      plot.subtitle = element_text(size = 8.5),
      axis.text.x = element_text(size = 7),
      legend.position = "right"
    )
}


make_ica_line_plot <- function(
    fit,
    simulation,
    max_components = 10L
) {
  title <- paste0("Simulation ", simulation, ": selected ICA sources")
  subtitle <- paste0(
    "Oracle rank = ", fit$selected_rank,
    "; RMSE = ", format(round(fit$rmse, 4), nsmall = 4)
  )

  if (ncol(fit$L_pm) == 0L) {
    return(empty_component_plot(title, subtitle))
  }

  shown <- seq_len(min(ncol(fit$L_pm), max_components))
  component_df <- as.data.frame(fit$L_pm[, shown, drop = FALSE])
  names(component_df) <- paste0("IC", shown)
  component_df$sample <- seq_len(nrow(component_df))
  component_df <- reshape2::melt(
    component_df,
    id.vars = "sample",
    variable.name = "component",
    value.name = "score"
  )

  if (ncol(fit$L_pm) > length(shown)) {
    subtitle <- paste0(
      subtitle,
      "; showing first ", length(shown), " components"
    )
  }

  ggplot(
    component_df,
    aes(
      x = sample,
      y = score,
      group = component
    )
  ) +
    geom_line(linewidth = 0.35, colour = "#2C7FB8") +
    facet_wrap(~component, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 10) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample",
      y = "Source score"
    ) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}


#------------------------------------------------------------
# Priors and storage
#------------------------------------------------------------

factor_priors <- list(
  Laplace = ebnm_point_laplace
)

results <- list()
all_fits <- setNames(
  vector("list", n_simulations),
  paste0("simulation_", simulation_ids)
)
simulation_panels <- setNames(
  vector("list", n_simulations),
  paste0("simulation_", simulation_ids)
)
oracle_rank_curves <- list()

ptime <- proc.time()


#------------------------------------------------------------
# Fit all simulations
#------------------------------------------------------------

for (i in simulation_ids) {
  message("Processing simulation ", i)

  Z <- read_numeric_matrix(
    file.path(sim_path, paste0("Z_matrix_", i, ".txt"))
  )
  X_true <- read_numeric_matrix(
    file.path(sim_path, paste0("X_matrix_", i, ".txt"))
  )
  L_cebmf <- read_numeric_matrix(
    file.path(sim_path, paste0("L_cebmf_", i, ".txt"))
  )
  fit_cebmf <- read_numeric_matrix(
    file.path(sim_path, paste0("fit_cebmf_", i, ".txt"))
  )

  validate_pair(Z, X_true)
  stopifnot(identical(dim(fit_cebmf), dim(X_true)))
  stopifnot(nrow(L_cebmf) == nrow(Z))

  set.seed(i)

  simulation_fits <- list()
  simulation_plots <- list()

  # Raw noisy-matrix baseline: necessary to measure whether denoising helped.
  raw_rmse <- matrix_rmse(Z, X_true)
  results[[length(results) + 1L]] <- data.frame(
    simulation = i,
    method = "No denoising (Z)",
    n_components = NA_integer_,
    selected_rank = NA_integer_,
    effective_rank = NA_integer_,
    rmse = raw_rmse,
    elbo = NA_real_,
    lambda = NA_real_
  )

  for (method_name in names(factor_priors)) {
    message("  Fitting FLASH: ", method_name)

    fit <- fit_flash_model(
      Z = Z,
      factor_prior = factor_priors[[method_name]]
    )
    current_rmse <- reconstruction_rmse(fit, X_true)
    current_rank <- ncol(fit$L_pm)

    simulation_fits[[method_name]] <- fit
    results[[length(results) + 1L]] <- data.frame(
      simulation = i,
      method = method_name,
      n_components = current_rank,
      selected_rank = current_rank,
      effective_rank = current_rank,
      rmse = current_rmse,
      elbo = tail(fit$elbo, 1),
      lambda = NA_real_
    )

    simulation_plots[[method_name]] <- make_component_heatmap(
      scores = fit$L_pm,
      simulation = i,
      method_name = paste0("GB + ", method_name),
      rmse = current_rmse,
      rank_text = paste0("Fitted factors = ", current_rank),
      max_components = max_components_to_plot
    )
  }

  # Imported cEBMF result.
  cebmf_rmse <- matrix_rmse(fit_cebmf, X_true)
  cebmf_rank <- ncol(L_cebmf)
  simulation_fits[["cEBMF= cGB + pt Normal"]] <- list(
    L_pm = L_cebmf,
    X_hat = fit_cebmf
  )
  results[[length(results) + 1L]] <- data.frame(
    simulation = i,
    method = "cEBMF= cGB + pt Normal",
    n_components = cebmf_rank,
    selected_rank = cebmf_rank,
    effective_rank = cebmf_rank,
    rmse = cebmf_rmse,
    elbo = NA_real_,
    lambda = NA_real_
  )
  simulation_plots[["cEBMF= cGB + pt Normal"]] <- make_component_heatmap(
    scores = L_cebmf[,1:6],
    simulation = i,
    method_name = "cEBMF= cGB + pt Normal",
    rmse = cebmf_rmse,
    rank_text = paste0("Fitted factors = ", 4),
    max_components = 6
  )

  message("  Fitting oracle PCA")
  pca_fit <- fit_oracle_pca(
    Z = Z,
    X_true = X_true,
    max_rank = max_oracle_rank
  )
  simulation_fits[["PCA (oracle)"]] <- pca_fit
  results[[length(results) + 1L]] <- data.frame(
    simulation = i,
    method = "PCA (oracle)",
    n_components = pca_fit$effective_rank,
    selected_rank = pca_fit$selected_rank,
    effective_rank = pca_fit$effective_rank,
    rmse = pca_fit$rmse,
    elbo = NA_real_,
    lambda = NA_real_
  )
  simulation_plots[["PCA (oracle)"]] <- make_component_heatmap(
    scores = pca_fit$L_pm,
    simulation = i,
    method_name = "oracle PCA",
    rmse = pca_fit$rmse,
    rank_text = paste0("Oracle rank = ", pca_fit$selected_rank),
    max_components = max_components_to_plot
  )
  oracle_rank_curves[[length(oracle_rank_curves) + 1L]] <-
    pca_fit$rmse_by_rank |>
    mutate(
      simulation = i,
      method = "PCA (oracle)",
      lambda = NA_real_,
      .before = 1
    )

  message("  Fitting oracle ICA")
  ica_fit <- fit_oracle_ica(
    Z = Z,
    X_true = X_true,
    max_rank = max_oracle_rank,
    seed = 1000L + i,
    maxit = ica_maxit
  )
  simulation_fits[["ICA (oracle)"]] <- ica_fit
  results[[length(results) + 1L]] <- data.frame(
    simulation = i,
    method = "ICA (oracle)",
    n_components = ica_fit$effective_rank,
    selected_rank = ica_fit$selected_rank,
    effective_rank = ica_fit$effective_rank,
    rmse = ica_fit$rmse,
    elbo = NA_real_,
    lambda = NA_real_
  )
  simulation_plots[["ICA (oracle)"]] <- make_component_heatmap(
    scores = ica_fit$L_pm,
    simulation = i,
    method_name = "oracle ICA",
    rmse = ica_fit$rmse,
    rank_text = paste0("Oracle rank = ", ica_fit$selected_rank),
    max_components = max_components_to_plot
  )
  oracle_rank_curves[[length(oracle_rank_curves) + 1L]] <-
    ica_fit$rmse_by_rank |>
    mutate(
      simulation = i,
      method = "ICA (oracle)",
      lambda = NA_real_,
      .before = 1
    )

  ica_line_plot <- make_ica_line_plot(
    fit = ica_fit,
    simulation = i,
    max_components = max_components_to_plot
  )
  ggsave(
    filename = file.path(
      plot_path,
      paste0("ica_components_simulation_", i, ".png")
    ),
    plot = ica_line_plot,
    width = 11,
    height = max(
      4.5,
      2.3 * ceiling(
        min(ica_fit$effective_rank, max_components_to_plot) / 2
      )
    ),
    dpi = 300
  )

  message(
    "  Fitting oracle-rank SoftImpute (lambda = ",
    softimpute_lambda,
    ")"
  )
  soft_fit <- fit_oracle_softimpute(
    Z = Z,
    X_true = X_true,
    max_rank = max_oracle_rank,
    lambda = softimpute_lambda
  )
  simulation_fits[["SoftImpute (oracle)"]] <- soft_fit
  results[[length(results) + 1L]] <- data.frame(
    simulation = i,
    method = "SoftImpute (oracle)",
    n_components = soft_fit$effective_rank,
    selected_rank = soft_fit$selected_rank,
    effective_rank = soft_fit$effective_rank,
    rmse = soft_fit$rmse,
    elbo = NA_real_,
    lambda = soft_fit$lambda
  )
  simulation_plots[["SoftImpute (oracle)"]] <- make_component_heatmap(
    scores = soft_fit$L_pm,
    simulation = i,
    method_name = "oracle-rank SoftImpute",
    rmse = soft_fit$rmse,
    rank_text = paste0(
      "Oracle rank.max = ", soft_fit$selected_rank,
      "; effective rank = ", soft_fit$effective_rank,
      "; lambda = ", soft_fit$lambda
    ),
    max_components = max_components_to_plot
  )
  oracle_rank_curves[[length(oracle_rank_curves) + 1L]] <-
    soft_fit$rmse_by_rank |>
    mutate(
      simulation = i,
      method = "SoftImpute (oracle)",
      lambda = soft_fit$lambda,
      .before = 1
    )

  all_fits[[paste0("simulation_", i)]] <- simulation_fits

  # Five component panels: 3 columns x 2 rows is readable at this size.
  panel_ncol <- 3L
  panel_nrow <- ceiling(length(simulation_plots) / panel_ncol)
  combined_panel <- gridExtra::arrangeGrob(
    grobs = simulation_plots,
    ncol = panel_ncol,
    top = grid::textGrob(
      paste0("Denoising factor estimates: simulation ", i),
      gp = grid::gpar(
        fontsize = 15,
        fontface = "bold"
      )
    )
  )

  simulation_panels[[paste0("simulation_", i)]] <- combined_panel

  grid::grid.newpage()
  grid::grid.draw(combined_panel)

  ggsave(
    filename = file.path(
      plot_path,
      paste0("component_heatmaps_simulation_", i, ".png")
    ),
    plot = combined_panel,
    width = 18,
    height = 5 * panel_nrow,
    dpi = 300
  )
}


#------------------------------------------------------------
# Numerical output
#------------------------------------------------------------

elapsed_time <- proc.time() - ptime
print(elapsed_time)

rmse_results <- dplyr::bind_rows(results)


# Create example factor
 rmse_results$method[which( rmse_results$method=="Laplace")]="EBMF= GB+pt Laplace"

raw_results <- rmse_results |>
  filter(method == "No denoising (Z)") |>
  select(simulation, raw_rmse = rmse)

rmse_results <- rmse_results |>
  left_join(raw_results, by = "simulation") |>
  mutate(
    rmse_gain = raw_rmse - rmse,
    relative_rmse_gain = if_else(
      raw_rmse > 0,
      (raw_rmse - rmse) / raw_rmse,
      NA_real_
    )
  )

write.csv(
  rmse_results,
  file = file.path(plot_path, "rmse_results.csv"),
  row.names = FALSE
)
print(rmse_results)




rmse_summary <- rmse_results |>
  group_by(method) |>
  summarise(
    mean_rmse = mean(rmse),
    sd_rmse = sd(rmse),
    mean_rmse_gain = mean(rmse_gain),
    mean_n_components = if (all(is.na(n_components))) {
      NA_real_
    } else {
      mean(n_components, na.rm = TRUE)
    },
    .groups = "drop"
  ) |>
  arrange(mean_rmse)
print(rmse_summary)

oracle_choices <- rmse_results |>
  filter(grepl("oracle", method, ignore.case = TRUE)) |>
  select(
    simulation,
    method,
    selected_rank,
    effective_rank,
    lambda,
    rmse
  )
print(oracle_choices)

oracle_rank_results <- dplyr::bind_rows(oracle_rank_curves)
write.csv(
  oracle_rank_results,
  file = file.path(plot_path, "oracle_rmse_by_rank.csv"),
  row.names = FALSE
)


#------------------------------------------------------------
# Readable RMSE comparison plot
#------------------------------------------------------------

method_order <- c(
  "No denoising (Z)",
  names(factor_priors),
  "cEBMF= cGB + pt Normal",
  "PCA (oracle)",
  "ICA (oracle)",
  "SoftImpute (oracle)"
)

method_order <- method_order[method_order %in% unique(rmse_results$method)]

rmse_plot_data <- rmse_results |>
  mutate(method = factor(method, levels = method_order))

if (anyNA(rmse_plot_data$method)) {
  stop(
    "Some result methods are missing from method_order: ",
    paste(
      setdiff(unique(rmse_results$method), method_order),
      collapse = ", "
    )
  )
}


rmse_comparison_plot <- ggplot(
  rmse_results[ -which (rmse_results$method =="No denoising (Z)"),],
  aes(
    x = simulation,
    y = rmse,
    colour = method,
    linetype = method,
    shape = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = simulation_ids) +

  theme_minimal(base_size = 11) +
  labs(
    title = "Denoising RMSE by simulation",
    subtitle = paste0(
      "RMSE no denoising = ",floor(mean (rmse_plot_data$rmse[   which (rmse_plot_data$method =="No denoising (Z)") ])*1e4)*1e-4
    ),
    x = "Simulation",
    y = "Reconstruction RMSE",
    colour = "Method",
    linetype = "Method",
    shape = "Method"
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

print(rmse_comparison_plot)
ggsave(
  filename = file.path(plot_path, "rmse_comparison.png"),
  plot = rmse_comparison_plot,
  width = 11,
  height = 6.5,
  dpi = 300
)


#------------------------------------------------------------
# Oracle RMSE-versus-rank diagnostic plot
#------------------------------------------------------------

rank_method_order <- c(
  "PCA (oracle)",
  "ICA (oracle)",
  "SoftImpute (oracle)"
)

oracle_rank_plot_data <- oracle_rank_results |>
  mutate(method = factor(method, levels = rank_method_order))

oracle_rank_plot <- ggplot(
  oracle_rank_plot_data,
  aes(
    x = rank,
    y = rmse,
    colour = method,
    group = method
  )
) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  geom_point(size = 1.2, na.rm = TRUE) +
  facet_wrap(~simulation, ncol = 5, scales = "free_y") +
  scale_colour_manual(
    values = method_colours[rank_method_order],
    drop = FALSE
  ) +
  theme_minimal(base_size = 9) +
  labs(
    title = "Oracle reconstruction RMSE by candidate rank",
    subtitle = paste0(
      "SoftImpute uses fixed lambda = ", softimpute_lambda
    ),
    x = "Candidate rank",
    y = "RMSE",
    colour = "Method"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(oracle_rank_plot)
ggsave(
  filename = file.path(plot_path, "oracle_rmse_by_rank.png"),
  plot = oracle_rank_plot,
  width = 14,
  height = 7,
  dpi = 300
)

