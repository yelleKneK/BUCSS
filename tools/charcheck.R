# Characterization harness (developer tool, not shipped).
# Runs a broad grid of VALID inputs through every ss_buc_* planner and records
# (the size row, ncp_adjusted). Used to prove that refactors leave
# valid-input outputs bit-for-bit unchanged, beyond the single-example oracles.
#
# Usage:
#   Rscript tools/charcheck.R save  /tmp/bucss_baseline.rds
#   Rscript tools/charcheck.R check /tmp/bucss_baseline.rds

suppressMessages(pkgload::load_all(".", quiet = TRUE))

grids <- list(
  ss_buc_independent_t = function() {
    # t = 30 exercises the uncapped root (adjusted parameter far above the 1.x
    # grid bound of 100).
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6, 30), n = c(15, 20, 40),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8, .9), KEEP.OUT.ATTRS = FALSE)
    rows <- lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
    # add a total-N slice (odd N exercises the N path's own two-branch rounding)
    g2 <- expand.grid(t_observed = c(3, 4), N = c(41, 80),
                      alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                      desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    c(rows, lapply(seq_len(nrow(g2)), function(i) as.list(g2[i, ])))
  },
  ss_buc_paired_t = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6), N = c(20, 40, 80),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8, .9), KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_one_way_anova = function() {
    # F = 60 exercises the uncapped root; N = 121 does not divide evenly, so
    # the round-down branch of the conservative two-sided rounding can bind.
    g <- expand.grid(F_observed = c(4, 6, 10, 60), N = c(60, 120, 121),
                     levels_A = c(3, 4),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8), KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_factorial_anova = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(120, 180),
                     levels_A = 2, levels_B = 3,
                     effect = c("factor_A", "factor_B", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_factorial_anova_general = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = 120, cells = 6, df_numerator = 2,
                     df_denominator = c(114, 104, 90),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_rm_anova = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(40, 60),
                     levels_A = 2, levels_B = 3,
                     effect = c("factor_A", "factor_B", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_rm_anova_general = function() {
    g <- expand.grid(F_observed = c(4, 6.5, 10), N = c(60, 80), df_numerator = c(1, 2),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_mixed_anova = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(60, 90),
                     levels_between = 2, levels_within = 3,
                     effect = c("between", "within", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_mixed_anova_general = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = 90, df_numerator = 2, num_groups = 3,
                     effect = "between_only", alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    rows <- lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
    # add a between_within slice that needs df_num_within
    g2 <- expand.grid(F_observed = c(4, 6, 10), N = 90, df_numerator = 2, num_groups = 3,
                      effect = "between_within", df_num_within = 3,
                      alpha_prior = c(.05, 1), assurance = c(.5, .8, .9), desired_power = .8,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    # and a within_only slice, so all three effect arms are characterized
    g3 <- expand.grid(F_observed = c(4, 6, 10), N = 90, df_numerator = 2, num_groups = 3,
                      effect = "within_only", alpha_prior = c(.05, 1),
                      assurance = c(.5, .8, .9), desired_power = .8,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    c(rows,
      lapply(seq_len(nrow(g2)), function(i) as.list(g2[i, ])),
      lapply(seq_len(nrow(g3)), function(i) as.list(g3[i, ])))
  },
  ss_buc_reg_coef = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6), N = c(80, 150), p = c(2, 3, 5),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_R2 = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(100, 150), p = c(3, 4, 6),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  },
  ss_buc_reg_joint = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(120, 150), p = c(4, 5),
                     p_joint = c(2, 3), alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
  }
)

run_grid <- function() {
  out <- list()
  for (fn_name in names(grids)) {
    fn <- get(fn_name)
    cases <- grids[[fn_name]]()
    for (args in cases) {
      key <- paste0(fn_name, "(", paste(names(args), unlist(args), sep = "=",
                                        collapse = ","), ")")
      res <- tryCatch(
        suppressWarnings(do.call(fn, args)),
        error = function(e) structure("error", msg = conditionMessage(e)))
      if (inherits(res, "bucss_power")) {
        out[[key]] <- c(n = res$value[res$term %in% .BUCSS_SIZE_TERMS][1],
                        ncp = res$value[res$term == "ncp_adjusted"])
      } else {
        out[[key]] <- c(n = NA_real_, ncp = NA_real_)
      }
    }
  }
  out
}

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; path <- args[2]
res <- run_grid()
cat("grid size:", length(res), "calls\n")

if (mode == "save") {
  saveRDS(res, path)
  cat("baseline saved to", path, "\n")
  cat("valid (non-error) cases:", sum(!vapply(res, function(x) is.na(x[["n"]]), logical(1))), "\n")
} else if (mode == "check") {
  base <- readRDS(path)
  stopifnot(identical(names(base), names(res)))
  diffs <- 0
  for (k in names(base)) {
    a <- base[[k]]; b <- res[[k]]
    if (!isTRUE(all.equal(a, b, tolerance = 0))) {
      diffs <- diffs + 1
      cat("DIFF:", k, "\n  was:", a["n"], a["ncp"], "\n  now:", b["n"], b["ncp"], "\n")
    }
  }
  if (diffs == 0) cat("OK: all", length(res), "grid outputs identical to baseline.\n") else
    cat("FAIL:", diffs, "differing cases.\n")
}
