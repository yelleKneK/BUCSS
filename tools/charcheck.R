# Characterization harness (developer tool, not shipped: tools/ is excluded by
# .Rbuildignore). Runs a broad grid of inputs through every ss_buc_* planner
# and records, for each call, every number the result reports: the size row,
# total_N, actual_power, ncp_adjusted, and the assurance_ceiling attribute. A
# call the planner refuses is recorded as well, by the kind of refusal (the
# zero noncentrality parameter, a nonsignificant prior, or another error) and
# the ceiling its message names, since a changed refusal is a changed output.
# Used to prove that a refactor leaves the outputs bit for bit unchanged, and,
# when an output is meant to change, to list exactly which cases moved. It is
# a characterization, never evidence of correctness: a wrong number recorded
# in the baseline stays wrong in every later comparison.
#
# Usage, from the package root:
#   Rscript tools/charcheck.R save  BASELINE.rds [--lib=PATH]
#   Rscript tools/charcheck.R check BASELINE.rds [--lib=PATH] [--out=DIFFS.csv]
#
# Without --lib the source tree is loaded with pkgload; with --lib the BUCSS
# installed in PATH is used (so a baseline can be saved from an earlier commit
# installed into its own library). 'check' prints every differing case with
# the fields that moved, writes them to DIFFS.csv when --out is given, and
# exits with status 1 on any difference or on a grid that no longer matches
# the baseline's; 'save' exits 0 once the baseline is written.

args <- commandArgs(trailingOnly = TRUE)
flag_value <- function(flag) {
  v <- sub(paste0("^--", flag, "="), "", grep(paste0("^--", flag, "="), args, value = TRUE))
  if (length(v)) v[length(v)] else NA_character_
}
positional <- args[!grepl("^--", args)]
mode <- positional[1]
path <- positional[2]
if (length(positional) < 2L || !mode %in% c("save", "check")) {
  cat("usage: Rscript tools/charcheck.R save|check BASELINE.rds [--lib=PATH] [--out=DIFFS.csv]\n")
  quit(status = 2L)
}
lib <- flag_value("lib")
out_csv <- flag_value("out")

if (!is.na(lib)) {
  suppressPackageStartupMessages(library(BUCSS, lib.loc = normalizePath(lib, mustWork = TRUE)))
  if (!identical(normalizePath(dirname(getNamespaceInfo("BUCSS", "path"))), normalizePath(lib))) {
    cat("BUCSS did not load from ", lib, "\n", sep = "")
    quit(status = 2L)
  }
} else {
  if (!file.exists("DESCRIPTION") ||
      !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]), "BUCSS")) {
    cat("Run from the BUCSS package root, or give --lib=PATH.\n")
    quit(status = 2L)
  }
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
}
bucss_ns <- asNamespace("BUCSS")
size_terms <- get(".BUCSS_SIZE_TERMS", envir = bucss_ns)
provenance <- sprintf("BUCSS %s from %s; %s", as.character(utils::packageVersion("BUCSS")),
                      getNamespaceInfo("BUCSS", "path"), R.version.string)
cat(provenance, "\n")

as_cases <- function(g) lapply(seq_len(nrow(g)), function(i) as.list(g[i, , drop = FALSE]))

grids <- list(
  ss_buc_independent_t = function() {
    # t = 30 exercises the uncapped root (adjusted parameter far above the 1.x
    # grid bound of 100).
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6, 30), n = c(15, 20, 40),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8, .9), KEEP.OUT.ATTRS = FALSE)
    # a total-N slice (odd N exercises the N path's own two-branch rounding)
    g2 <- expand.grid(t_observed = c(3, 4), N = c(41, 80),
                      alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                      desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    # unequal prior group sizes given as a vector: harmonic means that are
    # whole (10 and 15 give 12; 20 and 30 give 24; 30 and 70 give 42), equal
    # sizes written as a pair, and harmonic means that are not whole
    g3 <- expand.grid(t_observed = c(2.5, 3, 4), pair = 1:6,
                      alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                      desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    pairs <- list(c(10, 15), c(20, 30), c(30, 70), c(49, 49), c(14, 21), c(25, 40))
    c3 <- lapply(as_cases(g3), function(a) {
      n <- pairs[[a$pair]]
      a$pair <- NULL
      c(a[1L], list(n = n), a[-1L])
    })
    c(as_cases(g), as_cases(g2), c3)
  },
  ss_buc_paired_t = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6), N = c(20, 40, 80),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8, .9), KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_one_way_anova = function() {
    # F = 60 exercises the uncapped root; N = 121 does not divide evenly, so
    # the round-down branch of the conservative two-sided rounding can bind.
    g <- expand.grid(F_observed = c(4, 6, 10, 60), N = c(60, 120, 121),
                     levels_A = c(3, 4),
                     alpha_prior = c(.05, .10, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8), KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_factorial_anova = function() {
    # N = 125 over 6 cells is an unbalanced prior study (125/6 is not whole),
    # so both branches of the two-sided rounding are exercised.
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(120, 125, 180),
                     levels_A = 2, levels_B = 3,
                     effect = c("factor_A", "factor_B", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_factorial_anova_general = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = 120, cells = 6, df_numerator = 2,
                     df_denominator = c(114, 104, 90),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    # an unbalanced prior study: 125 over 6 cells, error df 125 - 6
    g2 <- expand.grid(F_observed = c(4, 6, 10), N = 125, cells = 6, df_numerator = 2,
                      df_denominator = 119,
                      alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                      desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    c(as_cases(g), as_cases(g2))
  },
  ss_buc_rm_anova = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(40, 60),
                     levels_A = 2, levels_B = 3,
                     effect = c("factor_A", "factor_B", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_rm_anova_general = function() {
    g <- expand.grid(F_observed = c(4, 6.5, 10), N = c(60, 80), df_numerator = c(1, 2),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_mixed_anova = function() {
    # N = 61 and 91 do not divide over the 2 groups, so the split-plot's own
    # two-sided rounding of an unbalanced prior study is exercised.
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(60, 61, 90, 91),
                     levels_between = 2, levels_within = 3,
                     effect = c("between", "within", "interaction"),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_mixed_anova_general = function() {
    # All three effect arms, at N = 90 (balanced over 3 groups) and at N = 61
    # and 91 (unbalanced). between_within needs df_num_within, and its
    # df_numerator is the between df (at most num_groups - 1) times
    # df_num_within: 2 = 2 x 1 and 6 = 2 x 3 are designs that exist.
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(61, 90, 91), df_numerator = 2, num_groups = 3,
                     effect = "between_only", alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    g2 <- rbind(
      expand.grid(F_observed = c(4, 6, 10), N = c(61, 90, 91), df_numerator = 2, num_groups = 3,
                  effect = "between_within", df_num_within = 1,
                  alpha_prior = c(.05, 1), assurance = c(.5, .8, .9), desired_power = .8,
                  stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE),
      expand.grid(F_observed = c(3, 4, 6), N = c(61, 90, 91), df_numerator = 6, num_groups = 3,
                  effect = "between_within", df_num_within = 3,
                  alpha_prior = c(.05, 1), assurance = c(.5, .8, .9), desired_power = .8,
                  stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
    # The combination the grid carried before 2.0.0's guard (df_numerator 2
    # with df_num_within 3, a design that cannot exist), kept at N = 90 so the
    # comparison with an older baseline shows the guard's effect.
    g2b <- expand.grid(F_observed = c(4, 6, 10), N = 90, df_numerator = 2, num_groups = 3,
                       effect = "between_within", df_num_within = 3,
                       alpha_prior = c(.05, 1), assurance = c(.5, .8, .9), desired_power = .8,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    g3 <- expand.grid(F_observed = c(4, 6, 10), N = c(61, 90, 91), df_numerator = 2, num_groups = 3,
                      effect = "within_only", alpha_prior = c(.05, 1),
                      assurance = c(.5, .8, .9), desired_power = .8,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    c(as_cases(g), as_cases(g2), as_cases(g2b), as_cases(g3))
  },
  ss_buc_reg_coef = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6), N = c(80, 150), p = c(2, 3, 5),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_R2 = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(100, 150), p = c(3, 4, 6),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_correlation = function() {
    g <- expand.grid(r_observed = c(.2, .35, .5, .7), N = c(60, 100, 151),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_one_sample_t = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 4, 6), N = c(20, 40, 81),
                     alpha_prior = c(.05, 1), assurance = c(.5, .8, .9),
                     desired_power = c(.8, .9), KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_ancova = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(120, 121), cells = c(3, 4),
                     n_covariates = c(0, 1, 3), alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_manova = function() {
    g <- expand.grid(F_observed = c(3, 4.5, 8), N = c(80, 81, 120),
                     p_variables = c(1, 3, 5), alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_chisq_diff = function() {
    g <- expand.grid(chisq_observed = c(5, 9.5, 20), N = c(150, 250, 400),
                     df_full = 42, df_restricted = 42 + c(1, 2, 5),
                     alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  },
  ss_buc_welch_t = function() {
    g <- expand.grid(t_observed = c(2.5, 3, 5), n_1 = c(30, 40), n_2 = c(40, 55),
                     sd_ratio = c(1, 1.5, 2.5), alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8,
                     KEEP.OUT.ATTRS = FALSE)
    # large effects from small prior studies, where the plan is near or below
    # the planner's floor of 10 per group
    g2 <- expand.grid(t_observed = c(6, 8), n_1 = c(12, 20), n_2 = 20,
                      sd_ratio = c(1, 2), alpha_prior = c(.05, 1),
                      assurance = c(.5, .8), desired_power = .8,
                      KEEP.OUT.ATTRS = FALSE)
    c(as_cases(g), as_cases(g2))
  },
  ss_buc_reg_joint = function() {
    g <- expand.grid(F_observed = c(4, 6, 10), N = c(120, 150), p = c(4, 5),
                     p_joint = c(2, 3), alpha_prior = c(.05, 1),
                     assurance = c(.5, .8, .9), desired_power = .8, KEEP.OUT.ATTRS = FALSE)
    as_cases(g)
  }
)

# A readable, stable key for one call: a vector argument is written c(a,b).
case_key <- function(fn_name, a) {
  vals <- vapply(a, function(v) {
    if (length(v) == 1L) as.character(v) else paste0("c(", paste(v, collapse = ","), ")")
  }, character(1))
  paste0(fn_name, "(", paste(names(a), vals, sep = "=", collapse = ","), ")")
}

# The kind of a refusal, read from the message in words both the original and
# the current wording carry, and the ceiling the message names, if any.
refusal_kind <- function(msg) {
  if (grepl("noncentrality parameter is\\s+zero", msg)) "zero_ncp"
  else if (grepl("nonsignificant based on your specified", msg) ||
           grepl("significant at the prior study's own degrees", msg)) "nonsignificant"
  else "other_error"
}
message_ceiling <- function(msg) {
  m <- regmatches(msg, regexec("at most about (\\.[0-9]+)|is about (\\.[0-9]+)", msg))[[1L]]
  if (!length(m)) return(NA_real_)
  as.numeric(m[nzchar(m)][2L])
}

fields <- c("n", "total_N", "actual_power", "ncp", "ceiling", "refusal", "refusal_ceiling")

run_grid <- function() {
  rows <- list()
  for (fn_name in names(grids)) {
    fn <- get(fn_name, envir = bucss_ns)
    for (a in grids[[fn_name]]()) {
      res <- tryCatch(suppressWarnings(do.call(fn, a)), error = function(e) e)
      pick <- function(term) {
        v <- res$value[res$term == term]
        if (length(v)) v[1L] else NA_real_
      }
      row <- if (inherits(res, "bucss_plan")) {
        ceil <- attr(res, "assurance_ceiling")
        data.frame(planner = fn_name, key = case_key(fn_name, a),
                   n = res$value[res$term %in% size_terms][1L],
                   total_N = pick("total_N"), actual_power = pick("actual_power"),
                   ncp = pick("ncp_adjusted"),
                   ceiling = if (is.null(ceil)) NA_real_ else ceil,
                   refusal = "", refusal_ceiling = NA_real_, stringsAsFactors = FALSE)
      } else {
        msg <- if (inherits(res, "error")) conditionMessage(res) else "not a bucss_plan"
        data.frame(planner = fn_name, key = case_key(fn_name, a),
                   n = NA_real_, total_N = NA_real_, actual_power = NA_real_,
                   ncp = NA_real_, ceiling = NA_real_,
                   refusal = refusal_kind(msg), refusal_ceiling = message_ceiling(msg),
                   stringsAsFactors = FALSE)
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}

t0 <- Sys.time()
res <- run_grid()
cat(sprintf("grid size: %d calls over %d planners (%d plans, %d refusals); %.0f s\n",
            nrow(res), length(unique(res$planner)), sum(res$refusal == ""),
            sum(res$refusal != ""), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
if (anyDuplicated(res$key)) {
  cat("duplicate case keys in the grid:", head(res$key[duplicated(res$key)]), "\n")
  quit(status = 2L)
}

if (mode == "save") {
  attr(res, "provenance") <- provenance
  saveRDS(res, path)
  cat("baseline saved to", path, "\n")
  quit(status = 0L)
}

base <- readRDS(path)
cat("baseline:", if (is.null(attr(base, "provenance"))) "(no provenance)" else attr(base, "provenance"), "\n")
if (!is.data.frame(base) || !all(c("key", fields) %in% names(base))) {
  cat("FAIL: the baseline was not saved by this version of charcheck.R (it lacks the recorded fields).\n")
  quit(status = 1L)
}
only_base <- setdiff(base$key, res$key)
only_now <- setdiff(res$key, base$key)
if (length(only_base) || length(only_now)) {
  cat(sprintf("FAIL: the grid differs from the baseline's (%d cases only in the baseline, %d only now).\n",
              length(only_base), length(only_now)))
  for (k in head(only_base, 10)) cat("  only in baseline:", k, "\n")
  for (k in head(only_now, 10)) cat("  only now:", k, "\n")
  quit(status = 1L)
}
base <- base[match(res$key, base$key), ]
same <- function(a, b) identical(a, b) || (is.na(a) && is.na(b) && is.numeric(a) && is.numeric(b))
diff_rows <- list()
for (i in seq_len(nrow(res))) {
  moved <- fields[!vapply(fields, function(f) same(base[[f]][i], res[[f]][i]), logical(1))]
  if (!length(moved)) next
  for (f in moved) {
    diff_rows[[length(diff_rows) + 1L]] <- data.frame(
      planner = res$planner[i], key = res$key[i], field = f,
      was = format(base[[f]][i], digits = 17), now = format(res[[f]][i], digits = 17),
      stringsAsFactors = FALSE)
  }
}
diffs <- if (length(diff_rows)) do.call(rbind, diff_rows) else
  data.frame(planner = character(), key = character(), field = character(),
             was = character(), now = character())
if (!is.na(out_csv)) {
  utils::write.csv(diffs, out_csv, row.names = FALSE)
  cat("differences written to", out_csv, "\n")
}
n_cases <- length(unique(diffs$key))
if (n_cases == 0L) {
  cat("OK: all", nrow(res), "grid outputs identical to baseline, in every recorded field.\n")
  quit(status = 0L)
}
for (k in unique(diffs$key)) {
  d <- diffs[diffs$key == k, ]
  cat("DIFF:", k, "\n")
  for (j in seq_len(nrow(d))) cat(sprintf("  %-15s was %s  now %s\n", d$field[j], d$was[j], d$now[j]))
}
cat(sprintf("FAIL: %d differing cases (%d field differences); by planner: %s\n", n_cases, nrow(diffs),
            paste(sprintf("%s %d", names(table(diffs$planner[!duplicated(diffs$key)])),
                          table(diffs$planner[!duplicated(diffs$key)])), collapse = ", ")))
quit(status = 1L)
