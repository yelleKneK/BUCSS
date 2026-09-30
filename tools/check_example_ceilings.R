# Developer check (not shipped: tools/ is excluded by .Rbuildignore). Every
# planner's @examples block ends with a try() example that requests more
# assurance than the prior result supports, under a comment that states the
# ceiling the error will report ("here near .90"). R CMD check cannot verify
# that comment: try() swallows the error, and a missing error is silent, so an
# example can claim a ceiling the message never reports, or fail to error at
# all and illustrate nothing. Both happened during development.
#
# For each of the 18 planner files (R/ss_buc_*.R) this script
#   1. collapses the roxygen @examples block and parses it as R, so a comment
#      that wraps across lines ("(here near" / "# .90, ...") is read whole;
#   2. for every try() example whose preceding comment states a ceiling,
#      evaluates the call and confirms that it stops with the classed
#      zero-noncentrality refusal (bucss_zero_ncp) and that the ceiling the
#      message reports ("about .90") equals the stated one;
#   3. re-runs the same call at the reported ceiling, which the message says
#      yields a usable plan, and confirms that it does and that the plan's
#      assurance_ceiling attribute, floored to two decimals, is the same value.
# It then asserts that exactly one ceiling comment was checked in each planner
# file, 18 in all (every planner carries one, ss_buc_one_sample_t included),
# so a reworded comment, a comment detached from its try(), or a new planner
# without one cannot pass silently by lowering the count.
#
# Usage, from the package root:
#   Rscript tools/check_example_ceilings.R              # loads the source tree (pkgload)
#   Rscript tools/check_example_ceilings.R --lib=PATH   # uses the BUCSS installed in PATH
#   Rscript tools/check_example_ceilings.R --pkg=DIR    # reads and loads another copy
#
# Exit status 0 when every comment matches and the count is right, 1
# otherwise. This is a developer Rscript, never sourced by users, so it quits
# with a status for tools/release_gate.R to read. Run it after editing any
# planner's example, any prior-study value in one, or any ceiling code.

expected_count <- 18L

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  v <- sub(paste0("^--", flag, "="), "", grep(paste0("^--", flag, "="), args, value = TRUE))
  if (length(v)) v[length(v)] else NA_character_
}
pkg_dir <- arg_value("pkg")
if (is.na(pkg_dir)) pkg_dir <- "."
pkg_dir <- normalizePath(pkg_dir, mustWork = TRUE)
lib <- arg_value("lib")

desc_file <- file.path(pkg_dir, "DESCRIPTION")
if (!file.exists(desc_file) ||
    !identical(unname(read.dcf(desc_file, fields = "Package")[1L, 1L]), "BUCSS")) {
  cat("Run from the BUCSS package root, or give --pkg=DIR.\n")
  quit(status = 1L)
}

if (!is.na(lib)) {
  suppressPackageStartupMessages(library(BUCSS, lib.loc = normalizePath(lib, mustWork = TRUE)))
  source_of <- dirname(getNamespaceInfo("BUCSS", "path"))
  if (!identical(normalizePath(source_of), normalizePath(lib))) {
    cat("BUCSS loaded from ", source_of, ", not from ", lib, "\n", sep = "")
    quit(status = 1L)
  }
} else {
  # export_all = FALSE: the examples see only what a user's session sees.
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = FALSE))
}
cat("BUCSS ", as.character(utils::packageVersion("BUCSS")), " from ",
    getNamespaceInfo("BUCSS", "path"), "\n", sep = "")

# The @examples block of a roxygen source file, with the roxygen markers
# removed: every line from "#' @examples" to the next tag or the end of the
# roxygen block.
examples_block <- function(lines) {
  start <- grep("^#'\\s*@examples\\s*$", lines)
  if (length(start) != 1L) return(NULL)
  i <- start + 1L
  out <- character()
  while (i <= length(lines) && grepl("^#'", lines[i]) && !grepl("^#'\\s*@[A-Za-z]", lines[i])) {
    out <- c(out, sub("^#' ?", "", lines[i]))
    i <- i + 1L
  }
  out
}

ceiling_rx <- "here\\s+near\\s+(\\.[0-9]+)"
fmt2 <- function(p) sub("^0", "", sprintf("%.2f", p))

# Evaluate one example call in a fresh environment, returning its value or its
# error; warnings are collected (in a small environment, not by superassignment)
# and muffled.
run_call <- function(cl) {
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character()
  value <- withCallingHandlers(
    tryCatch(eval(cl, envir = new.env(parent = globalenv())),
             error = function(e) e),
    warning = function(w) {
      seen$warnings <- c(seen$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = seen$warnings)
}

files <- sort(list.files(file.path(pkg_dir, "R"), pattern = "^ss_buc_.*[.]R$", full.names = TRUE))
rows <- list()
problems <- character()
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  block <- examples_block(lines)
  if (is.null(block)) {
    problems <- c(problems, sprintf("%s: no single @examples block", basename(f)))
    next
  }
  exprs <- tryCatch(parse(text = block, keep.source = TRUE), error = function(e) e)
  if (inherits(exprs, "error")) {
    problems <- c(problems, sprintf("%s: examples do not parse: %s", basename(f), conditionMessage(exprs)))
    next
  }
  is_comment <- grepl("^\\s*#", block)
  # Every ceiling claim in the block's comments, counted on the comment text
  # joined across lines, so a claim that wraps is still one claim.
  comment_text <- paste(sub("^\\s*#+\\s*", "", block[is_comment]), collapse = " ")
  claims_in_block <- lengths(regmatches(comment_text, gregexpr("here\\s+near", comment_text)))
  checked_here <- 0L
  srcrefs <- attr(exprs, "srcref")
  for (k in seq_along(exprs)) {
    e <- exprs[[k]]
    if (!(is.call(e) && identical(e[[1L]], as.name("try")))) next
    first <- srcrefs[[k]][1L]
    j <- first - 1L
    above <- character()
    while (j >= 1L && is_comment[j]) { above <- c(block[j], above); j <- j - 1L }
    claim <- paste(sub("^\\s*#+\\s*", "", above), collapse = " ")
    m <- regmatches(claim, regexec(ceiling_rx, claim))[[1L]]
    if (!length(m)) next
    stated <- m[2L]
    checked_here <- checked_here + 1L
    inner <- e[[2L]]
    live_run <- run_call(inner)
    v <- live_run$value
    cls <- if (inherits(v, "error")) class(v)[1L] else "(no error)"
    live <- if (inherits(v, "error")) {
      a <- regmatches(conditionMessage(v), regexec("about (\\.[0-9]+)", conditionMessage(v)))[[1L]]
      if (length(a)) a[2L] else "(no ceiling in message)"
    } else "(no error)"
    ok <- inherits(v, "bucss_zero_ncp") && identical(stated, live)
    rerun <- "not run"
    if (ok) {
      cl2 <- inner
      cl2$assurance <- as.numeric(live)
      r2 <- run_call(cl2)$value
      attr_c <- if (inherits(r2, "bucss_plan")) attr(r2, "assurance_ceiling") else NULL
      attr_s <- if (is.null(attr_c)) "(none)" else fmt2(floor(attr_c * 100 - 1e-9) / 100)
      rerun <- if (inherits(r2, "bucss_plan")) {
        sprintf("re-run at %s plans n = %s, attribute %s", live,
                format(r2$value[1L]), attr_s)
      } else if (inherits(r2, "error")) {
        sprintf("re-run at %s FAILS: %s", live, conditionMessage(r2))
      } else "re-run returned no plan"
      ok <- inherits(r2, "bucss_plan") && identical(attr_s, live)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      file = basename(f), stated = stated, live = live, class = cls,
      rerun = rerun, ok = ok, stringsAsFactors = FALSE)
  }
  if (checked_here != 1L || claims_in_block != checked_here) {
    problems <- c(problems, sprintf(
      "%s: %d ceiling claim(s) in the comments, %d attached to a try() example (expected exactly 1 of each)",
      basename(f), claims_in_block, checked_here))
  }
}

res <- if (length(rows)) do.call(rbind, rows) else
  data.frame(file = character(), stated = character(), live = character(),
             class = character(), rerun = character(), ok = logical())
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-8s %-34s doc=%s live=%s [%s]; %s\n",
              if (res$ok[i]) "OK" else "MISMATCH", res$file[i], res$stated[i],
              res$live[i], res$class[i], res$rerun[i]))
}
for (p in problems) cat("PROBLEM ", p, "\n", sep = "")
n_bad <- sum(!res$ok)
count_ok <- nrow(res) == expected_count && length(files) == expected_count
cat(sprintf("ceiling comments checked: %d of %d expected (%d planner files) | mismatches: %d | problems: %d\n",
            nrow(res), expected_count, length(files), n_bad, length(problems)))
if (n_bad > 0L || !count_ok || length(problems) > 0L) quit(status = 1L)
