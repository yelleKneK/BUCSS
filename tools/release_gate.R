# The release gate for BUCSS: the check CRAN runs, run here first, in full.
#
# No tarball of BUCSS leaves this machine without a full-mode run of this
# script printing ALL GATES PASSED and writing its record (the commit and the
# tarball's sha256) to the package's quality control folder. Modeled on the
# gates of the author's other packages; it encodes the portfolio's release
# rules (the quality control standard, Sections 6.2 and 10.2). Nothing here is
# optional, nothing is skipped for speed in full mode, and the script is never
# edited to make a failing run pass.
#
# Usage, from the package root:
#   Rscript tools/release_gate.R --quick         # the static half, after any documentation edit
#   caffeinate -i Rscript tools/release_gate.R   # full mode, before any tarball leaves
# Options:
#   --record-dir=DIR   where full mode writes its record (default: the folder
#                      "Quality Control/BUCSS" beside the package's folder)
#   --dest=DIR         where a certified tarball is copied (default: the
#                      folder that holds the package's folder)
#   --replace-stale    let full mode run although a tarball of this version
#                      already sits at the hand-off directory or in the
#                      package root; they are listed in the output and the
#                      record, never deleted, and only the one at the
#                      hand-off path is overwritten, by the certified
#                      tarball, after every gate has passed
#
# Quick mode certifies nothing. It reports the state of the working tree and
# the commit, runs every detector's positive control, and runs the static
# half: the version, and the release date in CITATION.cff and .zenodo.json,
# in lockstep with DESCRIPTION; stale citation phrases; stored rounding;
# constant seed defaults and the example seed rule; the patterns CRAN forbids
# in package code, in R/ and inst/tools/; Rd markup inside math; example
# blocks that CRAN's pretest runs twice or never; non-ASCII characters and em
# and en dashes; the letter to CRAN (cran-comments.md: no dashes, no stale
# citation phrase, no placeholder, and a Status line in its check results,
# which quick mode reports as not yet verified); roxygen regenerated into a
# copy with the pinned roxygen2 version (the gate refuses a stale man/ or
# NAMESPACE); the PDF manual rendered through LaTeX; the spelling check; the
# pkgdown reference index; and tools/check_example_ceilings.R, which must
# exit 0.
#
# Full mode refuses to start on macOS unless this process runs under
# caffeinate -i (the timing gates read wall-clock times, and a machine that
# sleeps mid-run records minutes of elapsed time against seconds of CPU), and
# refuses a working tree that is not clean (the gate certifies a commit, not a
# draft). It fails at once, and runs on so the record shows everything else,
# if the hand-off directory cannot be written or already holds a tarball of
# this version (as may the package root), unless --replace-stale is given.
# It then runs the static half, and adds: the DESCRIPTION Date within three
# days and a version above the one on CRAN; the tarball built from git
# archive HEAD in a temporary directory that no cloud service synchronizes;
# R CMD check --as-cran with both manuals and the remote checks on, which
# must exit 0 and leave a complete log (a Status line that agrees with its
# steps, a plausible number of steps, and 00check.log) before any gate on
# its output can pass; the letter's quoted check results against that log;
# the example time budget read from the check's timings file; the tarball
# installed into a scratch library, with the full test suite (NOT_CRAN=true),
# the example ceilings, and the shipped simulation harness in full mode run
# against it; and sweeps of the tarball and of every ref that can be pushed.
# Only when every gate has passed is the tarball copied out, and the copy's
# sha256 must equal the certified one. The record is written once full mode
# is past its two refusals, pass or fail, after the hand-off, and carries the
# hand-off path, R CMD check's output, and its 00check.log.
#
# R CMD check results. No ERROR and no WARNING passes. A NOTE passes only if
# it is on this explicit allowlist:
#   - "HTML Tidy not recent enough" (the Tidy that ships with macOS is from
#     2006): reported as SKIPPED LOCALLY, never as passed; win-builder runs it.
#   - "unable to verify current time" (the check could not reach its time
#     server): reported as SKIPPED LOCALLY; environmental.
#   - In the CRAN incoming feasibility NOTE, the Maintainer line; "New
#     submission" only while the package is not on CRAN; possibly misspelled
#     words only if each is declared in inst/WORDLIST; and a URL whose 5xx
#     status is transient (probed again, and returns 200). Every other line of
#     that NOTE fails the gate, including a 404 for the package's own site,
#     which must fail until the site is live.
#
# Positive controls: every detector is shown a known-bad input (and a
# known-good one) on every run, through the same function the gate itself
# calls, on planted input: a temporary .Rd file, a fake check log built from
# the strings R's own check code prints, a planted _pkgdown.yml, a small
# package with planted misspellings, a temporary git repository, and a copy
# of the committed sources (git archive HEAD, never the working tree, so a
# defect in the tree is never blamed on the detector). A detector that does
# not flag its bad input fails the gate, so no detector passes vacuously. If
# the committed sources fail the ceiling tool unaltered, that control has no
# known-good case and says so: it fails, except in quick mode on an
# uncommitted tree (whose commit may predate the tool), where it is reported
# as skipped, never as passed. Quick mode also lists any tarball of this
# version at the hand-off directory or in the package root, on which full
# mode fails unless --replace-stale is given.
#
# This file is not shipped: tools/ is excluded by .Rbuildignore. Every string
# it polices (the run-time idioms CRAN forbids, the example block tags, the
# stale citation phrases, the letter's placeholders, the instruction file's
# name, and the names of AI assistants) is assembled from character codes, so
# a search of the repository for any of them finds real offenders only, never
# this file. Do not edit this file while a run is in progress: R reads it
# incrementally.

args <- commandArgs(trailingOnly = TRUE)
known_flags <- c("--quick", "--record-dir", "--dest", "--replace-stale")
unknown <- args[!sub("=.*", "", args) %in% known_flags]
if (length(unknown)) {
  cat("unknown argument(s):", unknown, "\n",
      "usage: Rscript tools/release_gate.R [--quick] [--record-dir=DIR] [--dest=DIR] [--replace-stale]\n")
  quit(status = 2L)
}
quick <- "--quick" %in% args
replace_stale <- "--replace-stale" %in% args
flag_value <- function(flag, default) {
  v <- sub(paste0("^", flag, "="), "", grep(paste0("^", flag, "="), args, value = TRUE))
  if (length(v)) v[length(v)] else default
}

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  cat("Run from the package root.\n")
  quit(status = 2L)
}
desc <- read.dcf(file.path(root, "DESCRIPTION"))
pkg_name <- unname(desc[1L, "Package"])
if (!identical(pkg_name, "BUCSS")) {
  cat("This is the BUCSS gate; DESCRIPTION names", pkg_name, "\n")
  quit(status = 2L)
}
pkg_version <- unname(desc[1L, "Version"])
tarball_name <- sprintf("%s_%s.tar.gz", pkg_name, pkg_version)
record_dir <- flag_value("--record-dir", file.path(dirname(root), "Quality Control", pkg_name))
dest_dir <- normalizePath(flag_value("--dest", dirname(root)), mustWork = FALSE)

## ---- Budgets and expected counts ------------------------------------------
# CRAN flags an example page over 5 s on its farm, and the slowest farm
# flavor runs about three times slower than this machine, so a page may take
# 1.7 s here. The incoming pretest returns a package whose whole check takes
# over 10 minutes on its slowest flavor; a third of that is the local budget
# for the check with the network-bound feasibility step excluded. After the
# release is accepted, measure its check here and on the farm and tighten.
page_ceiling_s <- 1.7
check_budget_s <- 200
# The shipped harness and the example ceiling tool each state how many checks
# they run; a changed count fails the gate until the count here is updated on
# purpose, so a check cannot silently disappear.
expected_harness_checks <- 24L
expected_ceiling_comments <- 18L
# Help-page examples seed with 113, and only help-page examples do.
example_seed <- 113L

## ---- Results ------------------------------------------------------------------
# Results live in a small environment, so no function here assigns outside its
# own frame.
state <- new.env(parent = emptyenv())
state$results <- list()
state$facts <- character()
state$test_lines <- "(not run)"
state$harness_lines <- "(not run)"
record_result <- function(name, status, detail) {
  cat(sprintf("[%s] %s\n", status, name))
  for (d in detail) cat("       ", d, "\n", sep = "")
  state$results[[length(state$results) + 1L]] <- list(name = name, status = status, detail = detail)
}
gate <- function(name, ok, detail = character()) {
  record_result(name, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  invisible(isTRUE(ok))
}
# A step that did not run is recorded as skipped, with the reason; it is
# never counted as a pass.
gate_skip <- function(name, reason, label = "SKIPPED LOCALLY") {
  record_result(name, label, reason)
  invisible(NA)
}
fact <- function(...) {
  line <- paste0(...)
  cat("        ", line, "\n", sep = "")
  state$facts <- c(state$facts, line)
}
statuses <- function() vapply(state$results, `[[`, "", "status")
# Everything neither passed nor failed: steps skipped, and results not yet
# verified (quick mode's reading of the letter's check results).
skip_summary <- function() {
  other <- Filter(function(r) !r$status %in% c("PASS", "FAIL"), state$results)
  if (!length(other)) return("none skipped")
  kind <- vapply(other, function(r) if (startsWith(r$status, "SKIPPED")) "skipped" else tolower(r$status), "")
  nm <- vapply(other, `[[`, "", "name")
  paste(vapply(unique(kind), function(k) sprintf("%d %s [%s]", sum(kind == k), k, paste(nm[kind == k], collapse = "; ")), ""),
        collapse = "; ")
}

work <- file.path(tempdir(), "bucss_release_gate")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
t0 <- Sys.time()
cat(sprintf("BUCSS release gate, %s mode, %s\n", if (quick) "quick" else "full", format(t0, "%Y-%m-%d %H:%M:%S")))

## ---- Tools: this R, a git that runs, pandoc, and TeX ----------------------------
r_bin <- file.path(R.home("bin"), "R")
rscript_bin <- file.path(R.home("bin"), "Rscript")
runs_ok <- function(bin, arg = "--version") {
  nzchar(bin) && file.exists(bin) && identical(suppressWarnings(tryCatch(
    system2(bin, arg, stdout = FALSE, stderr = FALSE), error = function(e) 1L)), 0L)
}
resolve_tool <- function(name, dirs) {
  cands <- unique(c(unname(Sys.which(name)), file.path(dirs, name)))
  cands <- cands[nzchar(cands)]
  ok <- vapply(cands, runs_ok, logical(1))
  if (any(ok)) cands[ok][1L] else NA_character_
}
# A git binary that exits nonzero (a stale one early on the path) would make
# every git check pass on empty output, so each candidate is run first.
git_bin <- resolve_tool("git", c("/usr/bin", "/opt/homebrew/bin", "/usr/local/bin"))
pandoc_bin <- resolve_tool("pandoc", c(
  "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
  "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64",
  "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
  "/opt/homebrew/bin", "/usr/local/bin"))
tex_bin <- resolve_tool("pdflatex", c("/Library/TeX/texbin", "/opt/homebrew/bin", "/usr/local/bin"))
if (!is.na(pandoc_bin)) Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_bin))
Sys.setenv(PATH = paste(unique(c(dirname(na.omit(c(pandoc_bin, tex_bin))),
                                 strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]])),
                        collapse = .Platform$path.sep))
git_in <- function(dir, ...) suppressWarnings(system2(git_bin, c("-C", shQuote(dir), ...), stdout = TRUE, stderr = FALSE))
git <- function(...) git_in(root, ...)
fact(R.version.string, "; ", R.home())
fact("git: ", git_bin, "; pandoc: ", pandoc_bin, "; pdflatex: ", tex_bin)
gate("tools resolved (this R, a git that runs, pandoc, pdflatex)",
     !is.na(git_bin) && !is.na(pandoc_bin) && !is.na(tex_bin),
     c(if (is.na(git_bin)) "no git ran (Sys.which, /usr/bin, /opt/homebrew/bin, /usr/local/bin)",
       if (is.na(pandoc_bin)) "no pandoc ran (Sys.which and RStudio's bundled copy)",
       if (is.na(tex_bin)) "no pdflatex ran (Sys.which, /Library/TeX/texbin)"))

## ---- Strings the gate polices, assembled from character codes ----------------
chr <- function(...) intToUtf8(c(...))
rx_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
s_rng      <- chr(46, 82, 97, 110, 100, 111, 109, 46, 115, 101, 101, 100)   # the saved RNG state
s_genv     <- chr(46, 71, 108, 111, 98, 97, 108, 69, 110, 118)             # the global environment
s_genv_fn  <- chr(103, 108, 111, 98, 97, 108, 101, 110, 118)               # and its accessor
s_super    <- chr(60, 60, 45)                                              # superassignment
s_super_r  <- chr(45, 62, 62)                                              # and its right form
s_quit     <- chr(113, 117, 105, 116)
s_q        <- chr(113)
s_seed_fn  <- chr(115, 101, 116, 46, 115, 101, 101, 100)                   # the base seed setter
s_opts     <- chr(111, 112, 116, 105, 111, 110, 115)
s_warn     <- chr(119, 97, 114, 110)
tag_dt     <- paste0("\\", chr(100, 111, 110, 116, 116, 101, 115, 116))   # the example block run twice
tag_dr     <- paste0("\\", chr(100, 111, 110, 116, 114, 117, 110))        # and the one never run
conv_file  <- chr(67, 76, 65, 85, 68, 69, 46, 109, 100)                   # the instruction file
conv_dir   <- chr(46, 99, 108, 97, 117, 100, 101)                         # and its settings folder
ai_rx <- paste(c(chr(99, 108, 97, 117, 100, 101), chr(97, 110, 116, 104, 114, 111, 112, 105, 99),
                 chr(99, 111, 45, 97, 117, 116, 104, 111, 114, 101, 100, 45, 98, 121),
                 chr(99, 104, 97, 116, 103, 112, 116), chr(111, 112, 101, 110, 97, 105)),
               collapse = "|")
stale_words <- c(chr(105, 110, 32, 112, 114, 101, 115, 115), chr(115, 117, 98, 109, 105, 116, 116, 101, 100),
                 chr(117, 110, 100, 101, 114, 32, 114, 101, 118, 105, 101, 119),
                 chr(102, 111, 114, 116, 104, 99, 111, 109, 105, 110, 103),
                 chr(105, 110, 32, 112, 114, 101, 112, 97, 114, 97, 116, 105, 111, 110))
# A phrase may wrap across lines, and across a roxygen or comment marker.
stale_rx <- sprintf("(?i)\\b(%s)\\b", paste(gsub(" ", "[\\\\s#']+", stale_words), collapse = "|"))
dash_rx <- paste0("[", chr(0x2014, 0x2013), "]")
# Placeholders a letter to CRAN must not carry: the to-do tag, the triple X,
# the to-be-determined abbreviation, the phrase for a check not yet run, and
# the portfolio's fill marker (two brackets, then FILL).
s_todo <- chr(84, 79, 68, 79)
s_xxx  <- strrep(chr(88), 3L)
s_tbd  <- chr(84, 66, 68)
s_tbr  <- chr(116, 111, 32, 98, 101, 32, 114, 117, 110)
s_fill <- paste0(strrep("[", 2L), chr(70, 73, 76, 76))
placeholder_rx <- sprintf("\\b(%s|%s|%s)\\b|(?i:\\b%s\\b)|%s", s_todo, s_xxx, s_tbd,
                          gsub(" ", "\\\\s+", s_tbr), rx_escape(s_fill))
# Phrases of R CMD check's own output. A letter that states one of them as a
# result of the local check must find it in the local check's log.
check_phrases <- c("Possibly misspelled words", "invalid URLs", "invalid DOIs", "New submission",
                   "Skipping checking", "Days since last update", "Package was archived",
                   "Version contains large components")
check_status_rx <- "^Status: (OK|[0-9]+ (ERROR|WARNING|NOTE)s?(, [0-9]+ (ERROR|WARNING|NOTE)s?)*)$"
# Private working documents, by file name: the instruction file and its
# folder, Word documents, handoffs, reviews, audits, evaluations, decision
# documents, notes, and console history.
private_doc_rx <- paste0("(^|/)", rx_escape(conv_file), "$|(^|/)", rx_escape(conv_dir), "(/|$)|",
                         "[.]docx$|HANDOFF|GAP-ANALYSIS|RELEASE_QC|(^|[/_])QC_|_QC([_.]|$)|",
                         "audit|evaluation|decision|(^|/)notes/|(^|/)[.]Rhistory$")
# And what else never belongs in a tarball.
build_junk_rx <- "(^|/)[.]DS_Store$|cran-comments|[.]Rproj$|(^|/)[.]Rproj[.]user(/|$)|[.]Rcheck(/|$)|[.]tar[.]gz$|Rplots[.]pdf$"
std_surfaces <- c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R", "build", "data", "demo",
                  "exec", "inst", "man", "po", "src", "tests", "vignettes", "LICENSE", "LICENSE.md")

## ---- Detectors ------------------------------------------------------------------
read_lines <- function(f) suppressWarnings(readLines(f, warn = FALSE, encoding = "UTF-8"))
files_in <- function(dir, pattern, recursive = FALSE) {
  d <- file.path(root, dir)
  if (!dir.exists(d)) return(character())
  list.files(d, pattern = pattern, full.names = TRUE, recursive = recursive)
}
rel <- function(f) sub(paste0("^", rx_escape(root), "/"), "", f)
grep_lines <- function(pattern, files, ...) {
  hits <- character()
  for (f in files) {
    m <- grep(pattern, read_lines(f), ...)
    if (length(m)) hits <- c(hits, sprintf("%s:%d", rel(f), m))
  }
  hits
}

# The roxygen @examples block of each source file, as lines with the roxygen
# marker removed, and the line numbers they came from.
roxygen_examples <- function(lines) {
  out <- integer(); inside <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^#'\\s*@examples\\s*$", lines[i])) { inside <- TRUE; next }
    if (inside && (!grepl("^#'", lines[i]) || grepl("^#'\\s*@[A-Za-z]", lines[i]))) inside <- FALSE
    if (inside) out <- c(out, i)
  }
  out
}

# A pattern in text, matched across line breaks; returns the line number on
# which each match starts.
rx_line_hits <- function(lines, rx) {
  txt <- paste(lines, collapse = "\n")
  m <- gregexpr(rx, txt, perl = TRUE)[[1L]]
  if (m[1L] < 0) return(integer())
  vapply(m, function(p) 1L + lengths(regmatches(substr(txt, 1L, p), gregexpr("\n", substr(txt, 1L, p)))), 1L)
}
# Stale citation phrases.
stale_hits <- function(lines) rx_line_hits(lines, stale_rx)

# Example blocks CRAN's pretest runs twice or never, in source and Rd files.
example_block_hits <- function(files) c(grep_lines(tag_dt, files, fixed = TRUE), grep_lines(tag_dr, files, fixed = TRUE))
# The Rd tags R's own parser assigns in a parsed Rd file.
rd_tags <- function(x) unique(c(attr(x, "Rd_tag"), if (is.list(x)) unlist(lapply(x, rd_tags))))
# Non-ASCII bytes in files.
nonascii_hits <- function(files) {
  out <- character()
  for (f in files) {
    m <- grep("[^\\x01-\\x7F]", read_lines(f), perl = TRUE, useBytes = TRUE)
    if (length(m)) out <- c(out, sprintf("%s:%d", rel(f), m))
  }
  out
}

# The patterns CRAN forbids in package code, from the parse data, so comments
# and ordinary words never match: the RNG state and the global environment
# (as names or in strings), the global environment's accessor, superassignment
# in either direction, a warn entry in options() or withr's option helpers,
# and quitting the session.
cran_forbidden <- function(text, label) {
  pd <- tryCatch(utils::getParseData(parse(text = text, keep.source = TRUE)), error = function(e) e)
  if (inherits(pd, "error")) return(sprintf("%s: does not parse (%s)", label, conditionMessage(pd)))
  if (is.null(pd) || !nrow(pd)) return(character())
  tok <- pd[pd$terminal, ]
  hit <- function(rows, what) if (any(rows)) sprintf("%s:%d: %s", label, tok$line1[rows], what) else character()
  out <- c(
    hit(tok$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL") & tok$text %in% c(s_rng, s_genv), "global state by name"),
    hit(tok$token == "STR_CONST" & (grepl(s_rng, tok$text, fixed = TRUE) | grepl(s_genv, tok$text, fixed = TRUE)),
        "global state by string"),
    hit(tok$token == "SYMBOL_FUNCTION_CALL" & tok$text == s_genv_fn, "global environment accessor"),
    hit(tok$token %in% c("LEFT_ASSIGN", "RIGHT_ASSIGN") & tok$text %in% c(s_super, s_super_r), "superassignment"),
    hit(tok$token == "SYMBOL_FUNCTION_CALL" & tok$text %in% c(s_q, s_quit), "quits the session"))
  parent <- stats::setNames(pd$parent, pd$id)
  opt_fns <- c(s_opts, paste0("local_", s_opts), paste0("with_", s_opts))
  opt_calls <- unname(parent[as.character(tok$parent[tok$token == "SYMBOL_FUNCTION_CALL" & tok$text %in% opt_fns])])
  warn_subs <- which(tok$token == "SYMBOL_SUB" & tok$text == s_warn)
  for (w in warn_subs) {
    a <- tok$parent[w]
    while (!is.na(a) && a != 0) {
      if (a %in% opt_calls) { out <- c(out, sprintf("%s:%d: a %s entry in %s()", label, tok$line1[w], s_warn, s_opts)); break }
      a <- unname(parent[as.character(a)])
    }
  }
  out
}

# Stored rounding: round() or signif() in code that builds a value. Allowed in
# the print and format helpers (by name), in an argument that only builds a
# string (a message, a label), and in an integer comparison.
string_ctx <- c("paste", "paste0", "sprintf", "format", "formatC", "prettyNum", "cat", "message",
                "warning", "stop", "errorCondition", "warningCondition", "simpleError",
                "simpleWarning", "simpleCondition", "gettextf", "sQuote", "dQuote", "toString",
                "print", "as.character", "==", "!=")
display_fn_rx <- "^(print|format|knit_print|summary|toString|as[.]character)[.]|^[.](format|fmt)|^planning_sentence$"
call_name <- function(e) {
  f <- e[[1L]]
  if (is.name(f)) return(as.character(f))
  if (is.call(f) && is.name(f[[1L]]) && as.character(f[[1L]]) %in% c("::", ":::")) return(as.character(f[[3L]]))
  ""
}
# Whether element i of a list of call arguments or formals is the empty
# symbol (an argument with no default, or x[, 1]). The element is tested in
# place, never bound to a variable, since a variable bound to the empty
# symbol cannot be read.
empty_at <- function(lst, i) { el <- lst[i]; is.name(el[[1L]]) && !nzchar(as.character(el[[1L]])) }
walk_args <- function(e, f) {
  a <- as.list(e)
  for (i in seq_along(a)) if (!empty_at(a, i) && !is.null(a[[i]])) f(a[[i]])
}
stored_rounding <- function(exprs, label) {
  acc <- new.env(parent = emptyenv()); acc$hits <- character()
  walk <- function(e, fname, in_str) {
    if (is.pairlist(e)) { walk_args(e, function(a) walk(a, fname, in_str)); return(invisible()) }
    if (!is.call(e)) return(invisible())
    nm <- call_name(e)
    if (nm %in% c("round", "signif") && !in_str && !grepl(display_fn_rx, fname))
      acc$hits <- c(acc$hits, sprintf("%s: %s(): %s", label, fname, paste(deparse(e, width.cutoff = 200L)[1L])))
    walk_args(e[-1L], function(a) walk(a, fname, in_str || nm %in% string_ctx))
    invisible()
  }
  for (e in exprs) {
    if (is.call(e) && call_name(e) %in% c("<-", "=") && is.name(e[[2L]])) walk(e[[3L]], as.character(e[[2L]]), FALSE)
    else walk(e, "(top level)", FALSE)
  }
  acc$hits
}

# Seeds in signatures and calls: a function whose seed argument has a
# constant default, and any call to the base seed setter (package code seeds
# through withr::local_seed(), which restores the caller's state).
seed_problems <- function(exprs, label) {
  acc <- new.env(parent = emptyenv()); acc$hits <- character()
  walk <- function(e) {
    if (is.pairlist(e)) { walk_args(e, walk); return(invisible()) }
    if (!is.call(e)) return(invisible())
    nm <- call_name(e)
    if (identical(nm, "function")) {
      fml <- as.list(e[[2L]])
      for (i in seq_along(fml)) {
        if (empty_at(fml, i)) next
        if (grepl("seed", names(fml)[i], ignore.case = TRUE) && is.numeric(fml[[i]]))
          acc$hits <- c(acc$hits, sprintf("%s: argument '%s' has the constant default %s", label, names(fml)[i], format(fml[[i]])))
      }
    }
    if (identical(nm, s_seed_fn)) acc$hits <- c(acc$hits, sprintf("%s: calls %s()", label, s_seed_fn))
    walk_args(e[-1L], walk)
    invisible()
  }
  for (e in exprs) walk(e)
  acc$hits
}

# Literal seeds in lines: the base seed setter called with a number, withr's
# local_seed(n) and with_seed(n), and seed = n. Returns the values, one per match, named by line number.
seed_rx <- sprintf("(%s|local_seed|with_seed)\\(\\s*([0-9]+)|\\bseed\\s*=\\s*([0-9]+)", rx_escape(s_seed_fn))
seed_literals <- function(lines) {
  out <- integer()
  for (i in seq_along(lines)) {
    m <- regmatches(lines[i], gregexpr(seed_rx, lines[i], perl = TRUE))[[1L]]
    for (x in m) out <- c(out, stats::setNames(as.integer(sub(".*?([0-9]+)\\s*$", "\\1", x, perl = TRUE)), i))
  }
  out
}
example_seed_problems <- function(example_lines, other_lines, label_ex, label_other) {
  ex <- seed_literals(example_lines); ot <- seed_literals(other_lines)
  c(if (any(ex != example_seed)) sprintf("%s line %s: example seed %d, not %d", label_ex,
                                          names(ex)[ex != example_seed], ex[ex != example_seed], example_seed),
    if (any(ot == example_seed)) sprintf("%s line %s: seed %d outside a help-page example", label_other,
                                          names(ot)[ot == example_seed], example_seed))
}
# The same rule on one roxygen source file: its @examples lines against the
# rest, each read from a copy of the file with the other lines blanked, so
# that line numbers in the messages are the file's own.
r_file_seed_problems <- function(l, label) {
  ex <- roxygen_examples(l)
  ex_only <- rep("", length(l)); ex_only[ex] <- l[ex]
  rest <- l; rest[ex] <- ""
  example_seed_problems(ex_only, rest, label, label)
}

# Rd markup inside \eqn{} or \deqn{} (the PDF manual stops on it), parsed.
math_markup <- "\\\\(code|emph|link|strong|pkg|bold|verb|var|samp|file|env|dQuote|sQuote|href|url)\\{"
math_nodes <- function(x) {
  tag <- attr(x, "Rd_tag")
  here <- if (!is.null(tag) && tag %in% c("\\eqn", "\\deqn")) paste(unlist(x), collapse = "") else character()
  c(here, if (is.list(x)) unlist(lapply(x, math_nodes)) else character())
}
rd_math_markup <- function(rd) any(grepl(math_markup, math_nodes(rd)))
rd_examples_text <- function(rd) {
  ex <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\examples"), rd)
  if (!length(ex)) return(character())
  strsplit(paste(unlist(ex), collapse = ""), "\n", fixed = TRUE)[[1L]]
}

# Version lockstep among DESCRIPTION, the top NEWS heading, inst/CITATION,
# CITATION.cff, and .zenodo.json (the last two when present), and the
# release date of the last two against the DESCRIPTION Date.
lockstep_problems <- function(ver, date, news_lines, citation_text, citation_src, cff_lines, zenodo_lines) {
  out <- character()
  head <- grep("^#\\s", news_lines, value = TRUE)[1L]
  news_ver <- if (is.na(head)) NA_character_ else sub("^#\\s+BUCSS\\s+(\\S+)\\s*$", "\\1", head)
  if (is.na(head) || !identical(news_ver, ver)) out <- c(out, sprintf("NEWS.md top heading: %s", head))
  if (!is.null(citation_text)) {
    if (!any(grepl(ver, citation_text, fixed = TRUE)))
      out <- c(out, sprintf("inst/CITATION: no entry carries version %s", ver))
    # A version literal stands alone: not inside a DOI or other path (after
    # a slash) and not part of a longer dotted number.
    lit <- unique(unlist(regmatches(citation_src, gregexpr("(?<![/.0-9])[0-9]+[.][0-9]+[.][0-9]+(?![.0-9])", citation_src, perl = TRUE))))
    if (length(setdiff(lit, ver))) out <- c(out, sprintf("inst/CITATION hard-codes version %s", setdiff(lit, ver)))
  }
  cff_field <- function(field) sub(sprintf('^%s:\\s*"?([^"[:space:]]+)"?\\s*$', field), "\\1",
                                   grep(sprintf("^%s:", field), cff_lines, value = TRUE))
  if (!is.null(cff_lines)) {
    v <- cff_field("version")
    if (!identical(v, ver)) out <- c(out, sprintf("CITATION.cff version: %s", if (length(v)) v else "(none)"))
    d <- cff_field("date-released")
    if (!identical(d, date)) out <- c(out, sprintf("CITATION.cff date-released: %s, not the DESCRIPTION Date %s",
                                                   if (length(d)) d else "(none)", date))
  }
  zen_field <- function(field) {
    txt <- paste(zenodo_lines, collapse = " ")
    m <- regmatches(txt, regexpr(sprintf('"%s"\\s*:\\s*"[^"]*"', field), txt))
    if (length(m)) sub('.*"([^"]*)"$', "\\1", m) else character()
  }
  if (!is.null(zenodo_lines)) {
    v <- zen_field("version")
    if (!identical(v, ver)) out <- c(out, sprintf(".zenodo.json version: %s", if (length(v)) v else "(none)"))
    d <- zen_field("publication_date")
    if (!identical(d, date)) out <- c(out, sprintf(".zenodo.json publication_date: %s, not the DESCRIPTION Date %s",
                                                   if (length(d)) d else "(none)", date))
  }
  out
}

# Every file under man/ and NAMESPACE and DESCRIPTION, compared by content
# between the tree and a regenerated copy.
doc_differences <- function(a, b) {
  files_of <- function(d) c(if (dir.exists(file.path(d, "man"))) file.path("man", list.files(file.path(d, "man"), pattern = "[.]Rd$")),
                            "NAMESPACE", "DESCRIPTION")
  all <- union(files_of(a), files_of(b))
  md5 <- function(d, f) { p <- file.path(d, f); if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_ }
  all[vapply(all, function(f) !identical(md5(a, f), md5(b, f)), logical(1))]
}

# This process kept awake by caffeinate -i. On current macOS, caffeinate runs
# the command in its own process and asserts PreventUserIdleSystemSleep "on
# behalf of" that pid, so the assertion's pid must be this process or one of
# its ancestors; a caffeinate -i ancestor (older behavior) also counts.
process_chain <- function(pid = Sys.getpid()) {
  pids <- integer(); cmds <- character()
  for (i in seq_len(64L)) {
    line <- suppressWarnings(tryCatch(system2("ps", c("-o", "ppid=,args=", "-p", pid), stdout = TRUE, stderr = FALSE),
                                      error = function(e) character()))
    if (!length(line)) break
    line <- trimws(line[1L])
    pids <- c(pids, as.integer(pid)); cmds <- c(cmds, sub("^[0-9]+\\s+", "", line))
    ppid <- suppressWarnings(as.integer(sub("\\s.*$", "", line)))
    if (is.na(ppid) || ppid <= 1L) break
    pid <- ppid
  }
  list(pids = pids, cmds = cmds)
}
kept_awake <- function(assertions, chain) {
  heads <- grep("[(]caffeinate[)]:", assertions, useBytes = TRUE)
  idle_pids <- integer()
  for (h in heads) {
    if (!grepl("PreventUserIdleSystemSleep", assertions[h], useBytes = TRUE)) next
    d <- assertions[h + 1L]
    if (!is.na(d) && grepl("on behalf of", d, useBytes = TRUE))
      idle_pids <- c(idle_pids, as.integer(sub(".*[(]pid ([0-9]+)[)].*", "\\1", d, useBytes = TRUE)))
  }
  ancestor <- any(grepl("^(\\S*/)?caffeinate\\s+(-[A-Za-z]+\\s+)*-[A-Za-z]*i[A-Za-z]*(\\s|$)", chain$cmds))
  any(idle_pids %in% chain$pids) || ancestor
}

# R CMD check log: split into steps, each with its status and its lines.
parse_check_log <- function(log) {
  starts <- grep("^\\* ", log)
  steps <- lapply(seq_along(starts), function(k) {
    s <- starts[k]
    e <- if (k < length(starts)) starts[k + 1L] - 1L else length(log)
    body <- if (e > s) log[(s + 1L):e] else character()
    status <- regmatches(log[s], regexpr("(OK|NOTE|WARNING|ERROR)$", log[s]))
    if (!length(status)) {
      j <- grep("^\\s*(\\[[0-9]+s/[0-9]+s\\]\\s*)?(OK|NOTE|WARNING|ERROR)$", body)[1L]
      status <- if (is.na(j)) "" else sub("^.*(OK|NOTE|WARNING|ERROR)$", "\\1", trimws(body[j]))
      if (!is.na(j)) body <- body[-j]
    }
    if (status == "" && grepl("^\\* DONE$", log[s])) status <- "DONE"
    list(head = log[s], status = status, body = body)
  })
  status_line <- grep("^Status: ", log, value = TRUE)
  count <- function(what) {
    m <- regmatches(status_line, regexpr(sprintf("[0-9]+ %ss?", what), status_line))
    if (length(m)) as.integer(sub(" .*", "", m)) else 0L
  }
  list(steps = steps, status_line = if (length(status_line)) status_line[1L] else NA_character_,
       counts = c(ERROR = count("ERROR"), WARNING = count("WARNING"), NOTE = count("NOTE")))
}
tidy_rx <- c("^Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy[.]$",
             "^Please obtain a recent version of HTML Tidy by downloading a binary$",
             "^release or compiling the source code from <https://www[.]html-tidy[.]org/>[.]$")
time_rx <- "^unable to verify current time$"
# Classify the NOTEs: a list of failures (each with its lines), a list of
# steps skipped locally, and the URLs whose 5xx was transient.
classify_notes <- function(steps, on_cran, wordlist, probe) {
  fails <- list(); skips <- list(); transient <- character()
  norm <- function(x) gsub("[\u2018\u2019]", "'", x)
  for (st in steps) {
    if (st$status != "NOTE") next
    body <- norm(trimws(st$body)); body <- body[nzchar(body)]
    head <- st$head
    if (grepl("checking HTML version of manual", head)) {
      if (length(body) && all(Reduce(`|`, lapply(tidy_rx, grepl, x = body)))) {
        skips[[length(skips) + 1L]] <- list(name = "HTML manual: validation by HTML Tidy (win-builder runs it)", lines = body)
      } else fails[[length(fails) + 1L]] <- list(head = head, lines = st$body)
      next
    }
    if (grepl("checking for future file timestamps", head)) {
      if (length(body) && all(grepl(time_rx, body))) {
        skips[[length(skips) + 1L]] <- list(name = "future file timestamps (the time server did not answer)", lines = body)
      } else fails[[length(fails) + 1L]] <- list(head = head, lines = st$body)
      next
    }
    if (grepl("checking CRAN incoming feasibility", head)) {
      bad <- character(); section <- ""; url <- NA_character_
      raw <- norm(st$body)
      url_status <- function(i) {
        k <- i + 1L
        while (k <= length(raw) && grepl("^\\s{3,}", raw[k])) {
          if (grepl("^\\s*Status:", raw[k])) return(trimws(sub("^\\s*Status:", "", raw[k])))
          k <- k + 1L
        }
        "?"
      }
      for (i in seq_along(raw)) {
        l <- raw[i]; tl <- trimws(l)
        if (!nzchar(tl)) next
        if (!grepl("^\\s", l)) {
          section <- ""
          if (grepl("^Maintainer:", tl)) next
          if (grepl("^New submission$", tl) && !on_cran) next
          if (grepl("^Possibly misspelled words in DESCRIPTION:$", tl)) { section <- "words"; next }
          if (grepl("^Found the following \\(possibly\\) invalid URLs:$", tl)) { section <- "urls"; next }
          bad <- c(bad, l); next
        }
        if (section == "words") {
          w <- gsub("^'|'.*$|\\s.*$", "", tl)
          if (!w %in% wordlist) bad <- c(bad, sprintf("undeclared word: %s", tl))
          next
        }
        if (section == "urls") {
          if (grepl("^URL: ", tl)) {
            u <- sub("^URL: (\\S+).*$", "\\1", tl); s <- url_status(i)
            if (grepl("^5", s) && identical(probe(u), "200")) { transient <- c(transient, sprintf("%s [%s at check time, 200 on re-probe]", u, s)); next }
            bad <- c(bad, sprintf("%s [%s]", tl, s))
          }
          next
        }
        bad <- c(bad, l)
      }
      if (length(bad)) fails[[length(fails) + 1L]] <- list(head = head, lines = c(bad, "full block:", st$body))
      next
    }
    fails[[length(fails) + 1L]] <- list(head = head, lines = st$body)
  }
  list(fails = fails, skips = skips, transient = transient)
}
probe_url <- function(u) {
  tryCatch(suppressWarnings(system2("curl", c("-s", "-o", "/dev/null", "-w", "%{http_code}", "-L", "--max-time", "20", shQuote(u)),
                                    stdout = TRUE, stderr = FALSE))[1L], error = function(e) "000")
}

# Every verdict on R CMD check's output, as a named list of gates (name, ok,
# detail), plus the NOTE steps skipped locally and the transient URLs. The
# first gate asks whether the check ran to the end: exit status 0, a Status
# line whose counts agree with the steps, more than min_steps steps, a
# "* DONE" line, and a 00check.log with the same Status line. Unless it
# passes, no other gate here can pass, so a check that never started, or was
# cut off, is never recorded as clean.
dt_pass <- paste0("checking examples with --run-", substring(tag_dt, 2L))
check_verdicts <- function(lines, status, log00, on_cran, wordlist, probe, min_steps = 40L) {
  pl <- parse_check_log(lines)
  st <- vapply(pl$steps, `[[`, "", "status")
  found <- c(ERROR = sum(st == "ERROR"), WARNING = sum(st == "WARNING"), NOTE = sum(st == "NOTE"))
  status <- suppressWarnings(as.integer(status))
  s00 <- grep("^Status: ", log00, value = TRUE)[1L]
  run_problems <- c(
    if (!identical(status, 0L)) sprintf("exit status %s", if (length(status) && !is.na(status)) status else "(none: the check did not run)"),
    if (is.na(pl$status_line)) "no Status line",
    if (length(pl$steps) <= min_steps) sprintf("%d steps, more than %d expected", length(pl$steps), min_steps),
    if (!any(grepl("^\\* DONE$", lines))) "no '* DONE' line: the check did not finish",
    if (!is.na(pl$status_line) && !identical(unname(found), unname(pl$counts)))
      sprintf("its steps (%s) disagree with its Status line", paste(names(found), found, collapse = ", ")),
    if (!length(log00)) "no 00check.log in the check directory"
    else if (!identical(s00, pl$status_line)) sprintf("00check.log's Status line (%s) differs from the output's", s00))
  ran <- !length(run_problems)
  no_log <- "R CMD check left no complete log (see the first check gate), so this cannot pass"
  after_run <- function(name, ok, detail = character())
    list(name = name, ok = ran && isTRUE(ok), detail = c(if (!ran) no_log, detail))
  ew <- Filter(function(s) s$status %in% c("ERROR", "WARNING"), pl$steps)
  notes <- classify_notes(pl$steps, on_cran = on_cran, wordlist = wordlist, probe = probe)
  math_skip <- grep("Skipping checking math rendering", lines, value = TRUE)
  other_skip <- setdiff(grep("^\\s*Skipping checking", lines, value = TRUE),
                        c(math_skip, grep("HTML validation", lines, value = TRUE)))
  manuals <- grep("version of manual", lines, value = TRUE)
  second <- grep(dt_pass, lines, value = TRUE, fixed = TRUE)
  list(parsed = pl, ran = ran, skips = notes$skips, transient = notes$transient, gates = list(
    ran = list(name = sprintf("R CMD check ran to the end (exit status %s, %d steps, a Status line that agrees with them, and 00check.log)",
                              if (length(status) && !is.na(status)) status else "none", length(pl$steps)),
               ok = ran, detail = c(run_problems, if (!ran) tail(lines, 8L))),
    errors = after_run("R CMD check --as-cran: no ERROR and no WARNING", !length(ew),
                       unlist(lapply(ew, function(s) c(s$head, paste0("  ", s$body))))),
    notes = after_run("R CMD check: no NOTE outside the allowlist (the package's own site must resolve)", !length(notes$fails),
                      unlist(lapply(notes$fails, function(f) c(f$head, paste0("  ", f$lines))))),
    math = after_run("HTML manual math rendering checked (V8 on the check's library path)", !length(math_skip),
                     c(math_skip, if (length(math_skip)) "install V8 into the library the check uses")),
    other_skip = after_run("no other check step skipped", !length(other_skip), other_skip),
    manuals = after_run("PDF and HTML manuals built by the check",
                        any(grepl("checking PDF version of manual \\.\\.\\. .*OK$", lines)) &&
                          any(grepl("checking HTML version of manual", lines)), manuals),
    second_pass = after_run(sprintf("no second example pass (no %s{} block)", tag_dt), !length(second), second)))
}

# The Rd2pdf log: exit status, the PDF written, and no LaTeX error.
rd2pdf_problems <- function(status, out_exists, lines) {
  c(if (!identical(suppressWarnings(as.integer(status)), 0L)) sprintf("Rd2pdf exit status %s", status),
    if (!isTRUE(out_exists)) "no PDF written",
    head(grep("^! |LaTeX errors|Missing \\$ inserted", lines, value = TRUE), 10L))
}
# The pkgdown reference index: TRUE, or pkgdown's complaint.
pkgdown_problem <- function(dir) tryCatch({ utils::capture.output(suppressMessages(pkgdown::check_pkgdown(dir))); TRUE },
                                          error = function(e) conditionMessage(e))
# The spelling check: a data frame of unknown words, or an error.
spelling_result <- function(dir) tryCatch(spelling::spell_check_package(dir), error = function(e) e)

# Tarballs of this version already present in the given directories.
stale_tarballs <- function(dirs, name) {
  p <- unique(file.path(normalizePath(dirs[dir.exists(dirs)]), name))
  p[file.exists(p)]
}
describe_file <- function(p) sprintf("%s (%.1f KB, modified %s, sha256 %s)", p, file.info(p)$size / 1024,
                                     format(file.info(p)$mtime, "%Y-%m-%d %H:%M"), unname(tools::sha256sum(p)))
# The hand-off: copy the certified tarball, then read the copy back. The copy
# counts only if file.copy() says so and the file now at the path has the
# certified sha256; otherwise the reason is returned, and any file left at
# the path is named as not the certified tarball.
hand_off <- function(from, to_dir, name, sha) {
  to <- file.path(to_dir, name)
  acc <- new.env(parent = emptyenv()); acc$w <- character()
  copied <- dir.exists(to_dir) && withCallingHandlers(isTRUE(file.copy(from, to, overwrite = TRUE)),
    warning = function(w) { acc$w <- c(acc$w, conditionMessage(w)); invokeRestart("muffleWarning") })
  got <- if (file.exists(to)) unname(tools::sha256sum(to)) else NA_character_
  problems <- c(
    if (!dir.exists(to_dir)) sprintf("the hand-off directory %s does not exist", to_dir),
    if (dir.exists(to_dir) && !copied) c("file.copy() returned FALSE", acc$w),
    if (is.na(got)) sprintf("no file at %s", to)
    else if (!identical(got, sha)) sprintf("the file at %s has sha256 %s, not the certified %s: it is NOT the certified tarball; do not upload it",
                                           to, got, sha))
  list(ok = !length(problems), path = to, sha = got, problems = problems)
}

# The letter to CRAN. Its static problems: dashes, stale citation phrases,
# and placeholders, each by line.
letter_static_problems <- function(lines) {
  at <- function(n, what) if (length(n)) sprintf("cran-comments.md:%d: %s: %s", n, what, trimws(lines[n])) else character()
  c(at(grep(dash_rx, lines), "an em or en dash"),
    at(unique(stale_hits(lines)), "a stale citation status"),
    at(unique(rx_line_hits(lines, placeholder_rx)), "a placeholder"))
}
# The line numbers of the letter's R CMD check results section: from its
# heading to the next heading at the same level or above, without any
# subsection about another machine's check (win-builder, R-hub,
# mac-builder), whose log this gate never sees.
letter_check_section <- function(lines) {
  h <- grep("^#{1,6}\\s", lines)
  lvl <- function(i) nchar(sub("\\s.*$", "", lines[i]))
  s <- grep("^#{1,6}\\s.*R CMD check results", lines, ignore.case = TRUE)[1L]
  if (is.na(s)) return(integer())
  end_of <- function(i) { nx <- h[h > i & lvl(h) <= lvl(i)]; if (length(nx)) nx[1L] - 1L else length(lines) }
  idx <- if (end_of(s) > s) seq.int(s + 1L, end_of(s)) else integer()
  for (r in h[h > s & h <= end_of(s) & grepl("win-?builder|r-?hub|mac-?builder", lines[h], ignore.case = TRUE)])
    idx <- setdiff(idx, seq.int(r, end_of(r)))
  idx
}
# In those lines: the lines quoted as check output (Markdown code, indented
# or fenced), the quotations in the prose ("..." or curly quotes, each split
# at an elision), and the prose itself, all with their line numbers.
letter_quoted <- function(lines, idx) {
  code <- integer(); prose <- integer()
  in_fence <- FALSE; in_code <- FALSE; in_list <- FALSE; prev_blank <- TRUE
  for (i in idx) {
    l <- lines[i]
    if (grepl("^\\s*```", l)) { in_fence <- !in_fence; prev_blank <- FALSE; next }
    if (in_fence) { code <- c(code, i); next }
    if (!nzchar(trimws(l))) { prev_blank <- TRUE; next }
    ind <- nchar(sub("\\S.*$", "", gsub("\t", "    ", l)))
    if (ind >= 4L && (in_code || (prev_blank && !in_list))) { in_code <- TRUE; code <- c(code, i); prev_blank <- FALSE; next }
    in_code <- FALSE
    if (ind < 4L && grepl("^\\s*([*+-]|[0-9]+[.)])\\s", l)) in_list <- TRUE
    else if (ind == 0L && prev_blank) in_list <- FALSE
    prose <- c(prose, i); prev_blank <- FALSE
  }
  txt <- paste(lines[prose], collapse = "\n")
  line_at <- function(p) vapply(p, function(pp) prose[1L + lengths(regmatches(substr(txt, 1L, pp), gregexpr("\n", substr(txt, 1L, pp))))], 1L)
  m <- gregexpr("\"[^\"]+\"|\u201c[^\u201d]+\u201d", txt, perl = TRUE)[[1L]]
  quotes <- if (m[1L] < 0) list() else lapply(seq_along(m), function(k) {
    q <- substr(txt, m[k] + 1L, m[k] + attr(m, "match.length")[k] - 2L)
    list(line = line_at(m[k]), text = q)
  })
  list(code = code, prose = prose, prose_text = txt, line_at = line_at, quotes = quotes)
}
norm_text <- function(x) gsub("\\s+", " ", trimws(gsub("[\u201c\u201d]", "\"", gsub("[\u2018\u2019]", "'", x))))
split_elision <- function(x) { f <- trimws(strsplit(x, "\\s*(\\.\\.\\.|\u2026)\\s*")[[1L]]); f[nzchar(f)] }
# The letter's Status lines, in its check results section.
letter_status_lines <- function(lines, idx) { s <- trimws(lines[idx]); s[grepl(check_status_rx, s)] }
# The letter against a check's log: its check results section must quote a
# Status line, every one equal to the check's; every quoted line and every
# quotation (each fragment of an elided one) must appear in the log; and a
# phrase of the check's own output stated in the prose must appear there.
letter_log_problems <- function(lines, log_lines, status_line) {
  if (length(status_line) != 1L || is.na(status_line)) return("no Status line from a check to compare the letter with")
  idx <- letter_check_section(lines)
  if (!length(idx)) return("cran-comments.md has no 'R CMD check results' section")
  at <- function(i, what) sprintf("cran-comments.md:%d: %s", i, what)
  log_text <- norm_text(paste(log_lines, collapse = " "))
  log_line_set <- norm_text(log_lines)
  in_log <- function(frag) { f <- norm_text(frag)
    grepl(f, log_text, fixed = TRUE) || grepl(sub("[.,;:]$", "", f), log_text, fixed = TRUE) }
  q <- letter_quoted(lines, idx)
  st <- idx[grepl(check_status_rx, trimws(lines[idx]))]
  out <- c(
    if (!length(st)) "cran-comments.md: its check results section quotes no Status line",
    vapply(st[trimws(lines[st]) != status_line], function(i)
      at(i, sprintf("Status line '%s' differs from the check's '%s'", trimws(lines[i]), status_line)), ""))
  for (i in setdiff(q$code, st)) {
    l <- lines[i]
    if (!nzchar(trimws(l))) next
    ok <- if (grepl("\\.\\.\\.|\u2026", l)) all(vapply(split_elision(l), in_log, NA)) else norm_text(l) %in% log_line_set
    if (!ok) out <- c(out, at(i, sprintf("quoted as check output, not in the log: %s", trimws(l))))
  }
  for (x in q$quotes) {
    miss <- Filter(Negate(in_log), split_elision(norm_text(x$text)))
    if (length(miss)) out <- c(out, at(x$line, sprintf("quoted as check output, not in the log: \"%s\"", paste(miss, collapse = " ... "))))
  }
  for (ph in check_phrases) {
    m <- gregexpr(gsub(" ", "\\\\s+", ph), q$prose_text, perl = TRUE)[[1L]]
    if (m[1L] > 0 && !grepl(ph, log_text, fixed = TRUE))
      out <- c(out, at(unique(q$line_at(m)), sprintf("states check output the log does not contain: '%s'", ph)))
  }
  unique(out)
}

# Text files anywhere under a directory that name an AI assistant.
ai_content_hits <- function(dir) {
  out <- character()
  for (f in list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)) {
    con <- file(f, "rb"); head_bytes <- readBin(con, "raw", 8000L); close(con)
    if (any(head_bytes == as.raw(0L))) next
    m <- grep(ai_rx, suppressWarnings(readLines(f, warn = FALSE)), ignore.case = TRUE, useBytes = TRUE)
    if (length(m)) out <- c(out, sprintf("%s:%d", sub(paste0("^", rx_escape(dir), "/"), "", f), m))
  }
  out
}

# The refs of a repository that can be pushed, and what their history holds:
# every file name ever tracked, every commit message line, and every author
# and committer. A generated site branch (gh-pages) and a local-only private
# branch (named *-private, with no upstream) are not swept.
ref_sweep <- function(repo) {
  refs <- git_in(repo, "for-each-ref", shQuote("--format=%(refname)"), "refs/heads", "refs/remotes", "refs/tags")
  site_refs <- grep("gh-pages$", refs, value = TRUE)
  private_refs <- grep("^refs/heads/.*-private$", refs, value = TRUE)
  tracked_private <- private_refs[vapply(private_refs, function(b)
    length(git_in(repo, "rev-parse", "--abbrev-ref", shQuote(paste0(sub("^refs/heads/", "", b), "@{upstream}")))) > 0L, logical(1))]
  sweep_refs <- setdiff(refs, c(site_refs, setdiff(private_refs, tracked_private)))
  ever <- if (length(sweep_refs)) unique(git_in(repo, "log", shQuote(sweep_refs), "--name-only", shQuote("--pretty=format:"))) else character()
  ever <- ever[nzchar(ever)]
  msgs <- if (length(sweep_refs)) git_in(repo, "log", shQuote(sweep_refs), shQuote("--format=%B")) else character()
  ids <- if (length(sweep_refs)) unique(git_in(repo, "log", shQuote(sweep_refs), shQuote("--format=%an <%ae>%n%cn <%ce>"))) else character()
  list(site_refs = site_refs, local_private = setdiff(private_refs, tracked_private), sweep_refs = sweep_refs,
       ever = ever, private_ever = grep(private_doc_rx, ever, value = TRUE, ignore.case = TRUE, perl = TRUE),
       msgs = msgs, msg_hits = grep(ai_rx, msgs, value = TRUE, ignore.case = TRUE), ids = ids)
}

# The strings R's own check code prints, read from its source, so that a
# fake check log carries what a real one would, not what this file assumes.
r_check_strings <- function() {
  acc <- new.env(parent = emptyenv()); acc$s <- character()
  walk <- function(e) {
    if (is.character(e)) { acc$s <- c(acc$s, e); return(invisible()) }
    if (is.call(e) || is.pairlist(e)) walk_args(e, walk)
    invisible()
  }
  f <- tryCatch(get(".check_packages", envir = asNamespace("tools")), error = function(e) NULL)
  if (is.function(f)) walk(body(f))
  unique(sub("\n$", "", acc$s))
}

# Example timings: pages over the ceiling, on the larger of CPU and elapsed.
slow_pages <- function(et, ceiling_s) {
  et$max <- pmax(et$user + et$system, et$elapsed)
  et[et$max > ceiling_s, c("name", "max")]
}

# Tarball file list: top-level entries that are not package surfaces, and
# any file anywhere that is a private document or build debris.
tarball_problems <- function(paths, pkg) {
  inner <- sub(paste0("^", pkg, "/"), "", paths)
  inner <- inner[nzchar(inner)]
  top <- unique(sub("/.*$", "", inner))
  c(sprintf("not a package surface: %s", setdiff(top, std_surfaces)),
    sprintf("private or build file: %s", grep(paste0(private_doc_rx, "|", build_junk_rx), inner,
                                                value = TRUE, ignore.case = TRUE, perl = TRUE)))
}
# Commit identities: each must name a person in Authors@R (or be GitHub's own
# committer for a change made on the web), and none may name an AI assistant.
identity_problems <- function(ids, people) {
  norm <- function(x) trimws(gsub("\\s+", " ", tolower(gsub("[^A-Za-z]", " ", x))))
  person_ok <- function(id) {
    nm <- norm(sub("\\s*<.*$", "", id))
    any(vapply(people, function(p) {
      fam <- norm(p$family); giv <- strsplit(norm(p$given), " ")[[1L]][1L]
      grepl(sprintf("\\b%s\\b", fam), nm) && grepl(sprintf("\\b%s\\b", giv), nm)
    }, logical(1)))
  }
  web <- ids == "GitHub <noreply@github.com>"
  ids[(!vapply(ids, person_ok, logical(1)) & !web) | grepl(ai_rx, ids, ignore.case = TRUE)]
}
# The shipped harness's closing lines, parsed.
harness_verdict <- function(lines, status, lib) {
  from <- sub("^BUCSS \\S+ from ", "", grep("^BUCSS \\S+ from ", lines, value = TRUE))
  mode <- trimws(grep("^mode: ", lines, value = TRUE))
  summ <- regmatches(lines, regexpr("^[0-9]+ checks, [0-9]+ failed[.]$", lines))
  n <- if (length(summ)) as.integer(sub(" checks.*", "", summ)) else NA_integer_
  nf <- if (length(summ)) as.integer(sub(".*, ([0-9]+) failed[.]$", "\\1", summ)) else NA_integer_
  passes <- sum(grepl("^\\s*\\[PASS\\]", lines)); fails <- sum(grepl("^\\s*\\[FAIL\\]", lines))
  problems <- c(
    if (!identical(as.integer(status), 0L)) sprintf("exit status %s", status),
    if (length(from) != 1L || !startsWith(normalizePath(from, mustWork = FALSE), normalizePath(lib, mustWork = FALSE)))
      sprintf("tested BUCSS from %s, not from the scratch library", if (length(from)) from[1L] else "(no provenance line)"),
    if (!identical(mode, "mode: full")) "not run in full mode",
    if (is.na(n)) "no 'N checks, M failed.' line (an empty result is a failure of the tool)",
    if (!is.na(nf) && nf != 0L) sprintf("%d failed", nf),
    if (!is.na(n) && n != expected_harness_checks) sprintf("%d checks run, %d expected", n, expected_harness_checks),
    if (!is.na(n) && passes + fails != n) sprintf("%d result lines for %d checks", passes + fails, n))
  list(ok = !length(problems), problems = problems, summary = if (length(summ)) summ else "(none)")
}
# The test suite's one-line summary, parsed.
tests_verdict <- function(lines) {
  l <- grep("^TESTS ", lines, value = TRUE)
  if (length(l) != 1L) return(list(ok = FALSE, summary = "no test summary line (an empty result is a failure of the tool)"))
  num <- function(k) as.integer(sub(sprintf(".*\\b%s=([0-9]+).*", k), "\\1", l))
  ok <- num("expectations") > 0L && num("failed") == 0L && num("errors") == 0L && num("warnings") == 0L
  list(ok = ok, summary = sub("^TESTS ", "", l), skipped = num("skipped"))
}

## ---- 0. Positive controls: every detector flags its known-bad case ----------------
# The commit and the state of the working tree are read first: the example
# ceiling control runs on the committed sources.
head_sha <- if (is.na(git_bin)) character() else git("rev-parse", "HEAD")
sha <- if (length(head_sha) == 1L) substr(head_sha, 1L, 7L) else "(no commit)"
dirty <- if (is.na(git_bin)) "(no git)" else git("status", "--porcelain")

# Every control runs the detector's own function on input planted here, in
# a scratch folder of the gate's temporary directory.
ctrl_dir <- file.path(work, "controls")
dir.create(ctrl_dir)
planted <- function(name, lines) {
  p <- file.path(ctrl_dir, name)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, p, useBytes = TRUE)
  p
}
# A control that errors has failed.
safe <- function(expr) tryCatch(isTRUE(expr), error = function(e) FALSE)

# The example ceiling tool on copies of the committed sources (git archive
# HEAD, never the working tree). Unaltered, the copy must exit 0 with every
# comment checked: the known-good case. With one comment's value altered the
# tool must report that file, and with one comment reworded it must count
# one fewer. If the unaltered copy fails, the control has no known-good case
# and says so (NA), rather than blaming the tool for a defect in the sources.
ceiling_ctrl <- local({
  arch <- file.path(ctrl_dir, "head_sources.tar")
  arch_st <- if (is.na(git_bin)) 1L else suppressWarnings(system2(git_bin, c("-C", shQuote(root), "archive", "--format=tar",
    paste0("--output=", shQuote(arch)), "HEAD", "DESCRIPTION", "NAMESPACE", "R"), stdout = FALSE, stderr = FALSE))
  run <- function(label, mutate = NULL) {
    d <- file.path(ctrl_dir, paste0("ceiling_", label)); dir.create(d)
    utils::untar(arch, exdir = d)
    f <- file.path(d, "R", "ss_buc_independent_t.R")
    if (!file.exists(f)) return(list(status = NA, lines = "the archive holds no R/ss_buc_independent_t.R"))
    if (!is.null(mutate)) {
      txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
      new <- mutate(txt)
      if (identical(new, txt)) return(list(status = NA, lines = "the control's alteration did not apply"))
      writeLines(new, f)
    }
    log <- file.path(d, "control.log")
    st <- system2(rscript_bin, c(shQuote(file.path(root, "tools", "check_example_ceilings.R")), paste0("--pkg=", shQuote(d))),
                  stdout = log, stderr = log)
    list(status = st, lines = readLines(log, warn = FALSE))
  }
  if (!identical(as.integer(arch_st), 0L) || !file.exists(arch))
    return(list(ok = NA, why = "git archive HEAD failed, so there is no copy of the committed sources to plant in"))
  base <- run("baseline")
  if (!identical(as.integer(base$status), 0L) ||
      !any(grepl(sprintf("checked: %d of %d expected", expected_ceiling_comments, expected_ceiling_comments), base$lines, fixed = TRUE)))
    return(list(ok = NA, why = c("the committed sources fail the ceiling tool unaltered, so the control has no known-good case:",
                                 head(grep("^(MISMATCH|PROBLEM)|^Error|checked:", base$lines, value = TRUE), 5L))))
  a <- run("value", function(t) sub("(here[[:space:]#']+near[[:space:]#']+)[.]([0-9]+)", "\\1.07", t, perl = TRUE))
  b <- run("wording", function(t) sub("here([[:space:]#']+)near", "here\\1about", t, perl = TRUE))
  ok <- identical(as.integer(a$status), 1L) && any(grepl("^MISMATCH\\s+ss_buc_independent_t[.]R", a$lines)) &&
    identical(as.integer(b$status), 1L) &&
    any(grepl(sprintf("checked: %d of %d", expected_ceiling_comments - 1L, expected_ceiling_comments), b$lines, fixed = TRUE))
  list(ok = ok, why = if (!ok) c("value altered:", tail(a$lines, 2L), "comment reworded:", tail(b$lines, 2L)))
})

controls <- local({
  bad_44 <- c(sprintf("f <- function() get('%s', envir = %s)", s_rng, s_genv),
              sprintf("f <- function() exists('%s')", s_rng),
              sprintf("f <- function() assign('z', 1, envir = %s())", s_genv_fn),
              sprintf("f <- function() { k <- 0; g <- function() k %s k + 1 }", s_super),
              sprintf("f <- function() %s(%s = -1)", s_opts, s_warn),
              sprintf("f <- function() withr::local_%s(list(%s = -1))", s_opts, s_warn),
              sprintf("f <- function() %s(status = 1)", s_quit),
              sprintf("f <- function() %s('no')", s_q))
  good_44 <- paste("f <- function(x) { withr::local_seed(1); old <- options(digits = 3)",
                   "on.exit(options(old)); y <- list(q = 1); y$q } # the global environment", sep = "; ")
  rnd_bad <- parse(text = c("f <- function(x) { y <- round(x * 2, 2); y }", "h <- function(x) base::signif(x, 3)"))
  rnd_good <- parse(text = c("g <- function(x) stop('value ', signif(x, 3))", ".format_x <- function(x) round(x, 2)",
                             "print.foo <- function(x, ...) cat(round(x, 2))", "k <- function(n) n == round(n)"))
  seed_bad <- parse(text = c("f <- function(x, seed = 2017) x", sprintf("g <- function() %s(1)", s_seed_fn)))
  seed_good <- parse(text = c("f <- function(x, seed = NULL) withr::local_seed(seed)"))
  rd_bad <- tools::parse_Rd(textConnection(c("\\name{x}", "\\alias{x}", "\\title{x}",
                                             "\\description{The region is \\eqn{(-\\code{delta}, \\code{delta})}.}")))
  rd_good <- tools::parse_Rd(textConnection(c("\\name{x}", "\\alias{x}", "\\title{x}",
                                              "\\description{The region is \\eqn{(-\\delta, \\delta)}{(-delta, delta)}.}")))
  rd_file <- function(name, examples) planted(paste0(name, ".Rd"), c(sprintf("\\name{%s}", name), sprintf("\\alias{%s}", name),
                                                                    "\\title{x}", "\\description{x}", "\\examples{", examples, "}"))
  log_lines <- c(
    "* using log directory '/x/BUCSS.Rcheck'",
    "* checking CRAN incoming feasibility ... [4s/16s] NOTE",
    "Maintainer: 'A B <a@b.c>'", "",
    "Found the following (possibly) invalid URLs:",
    "  URL: https://example.invalid/site", "    From: DESCRIPTION", "    Status: 404", "    Message: Not Found",
    "  URL: https://example.invalid/busy", "    From: man/x.Rd", "    Status: 503", "    Message: Service Unavailable",
    "* checking examples ... [1s/1s] OK",
    "* checking tests ...", "  Running 'testthat.R' [5s/5s]", " [5s/5s] WARNING",
    "Running the tests in 'tests/testthat.R' produced warnings",
    "* checking HTML version of manual ... NOTE",
    "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.",
    "Please obtain a recent version of HTML Tidy by downloading a binary",
    "release or compiling the source code from <https://www.html-tidy.org/>.",
    "* checking Rd cross-references ... NOTE", "Undeclared package 'foo' in Rd xrefs",
    "* DONE", "Status: 1 WARNING, 3 NOTEs")
  pl <- parse_check_log(log_lines)
  n_found <- table(factor(vapply(pl$steps, `[[`, "", "status"), levels = c("ERROR", "WARNING", "NOTE")))
  cl_200 <- classify_notes(pl$steps, on_cran = TRUE, wordlist = character(), probe = function(u) "200")
  cl_503 <- classify_notes(pl$steps, on_cran = TRUE, wordlist = character(), probe = function(u) "503")
  fail_heads <- function(cl) sort(sub(" \\.\\.\\..*$", "", sub("^\\* checking ", "", vapply(cl$fails, `[[`, "", "head"))))

  # Fake check logs for the verdicts on R CMD check's output, carrying the
  # lines R's own check code prints for the second example pass and for the
  # skipped math rendering. min_steps is lowered to fit a short log.
  rs <- r_check_strings()
  second_steps <- grep("^examples with --run-", rs, value = TRUE)
  math_line <- grep("^Skipping checking math rendering", rs, value = TRUE)
  tidy_lines <- log_lines[20:22]
  fake_check <- function(steps = character(), html = "* checking HTML version of manual ... [0s/1s] OK",
                         pdf = "* checking PDF version of manual ... [4s/5s] OK", status = "Status: OK")
    c("* using log directory '/x/BUCSS.Rcheck'", "* checking for file 'BUCSS/DESCRIPTION' ... OK",
      "* checking examples ... [1s/1s] OK", steps, pdf, html, "* DONE", "", status)
  cv <- function(lines, status = 0L, log00 = lines, min_steps = 3L)
    check_verdicts(lines, status, log00, on_cran = TRUE, wordlist = character(), probe = function(u) "200", min_steps = min_steps)
  oks <- function(v) vapply(v$gates, `[[`, NA, "ok")
  html_note <- "* checking HTML version of manual ... [0s/1s] NOTE"
  cv_good <- cv(fake_check())
  cv_second <- cv(fake_check(steps = sprintf("* checking %s ... [2s/2s] OK", second_steps)))
  cv_math <- cv(fake_check(html = c(html_note, math_line), status = "Status: 1 NOTE"))
  cv_other <- cv(fake_check(html = c(html_note, tidy_lines, "Skipping checking the planted step: no reason"), status = "Status: 1 NOTE"))
  cv_tidy <- cv(fake_check(html = c(html_note, tidy_lines), status = "Status: 1 NOTE"))
  cv_nopdf <- cv(fake_check(pdf = character()))
  # Output that never reached its end: none, cut off, an error exit, no 00check.log.
  cv_dead <- list(cv(character(), NA, character(), min_steps = 40L), cv(character(), 127L, character(), min_steps = 40L),
                  cv(head(fake_check(), 3L), 137L, character()), cv(fake_check(), 1L), cv(fake_check(), 0L, character()),
                  cv(fake_check(), 0L, c("* DONE", "Status: 1 NOTE")), cv(fake_check()[-6L], 0L), cv(fake_check(), 0L, min_steps = 40L))

  git_ctrl <- !is.na(git_bin) && tryCatch(local({
    d <- file.path(work, "git_control"); dir.create(d)
    system2(git_bin, c("-C", shQuote(d), "init", "-q"), stdout = FALSE, stderr = FALSE)
    writeLines("x", file.path(d, "untracked.txt"))
    out <- suppressWarnings(system2(git_bin, c("-C", shQuote(d), "status", "--porcelain"), stdout = TRUE, stderr = FALSE))
    identical(out, "?? untracked.txt")
  }), error = function(e) FALSE)
  ids_bad <- c("Ken Kelley <kkelley@nd.edu>", "GitHub <noreply@github.com>", "Samantha F Anderson <s@a.edu>",
               sprintf("Somebody <noreply@%s.com>", chr(97, 110, 116, 104, 114, 111, 112, 105, 99)),
               "Anonymous Author <anon@example.com>")
  people <- list(list(given = "Ken", family = "Kelley"), list(given = "Samantha F.", family = "Anderson"))
  ai_email <- sprintf("noreply@%s.com", chr(97, 110, 116, 104, 114, 111, 112, 105, 99))
  ai_trailer <- sprintf("%s: Somebody <%s>", chr(67, 111, 45, 65, 117, 116, 104, 111, 114, 101, 100, 45, 66, 121), ai_email)
  asserts <- c("   pid 900(caffeinate): [0x1] 00:00:00 PreventUserIdleSystemSleep named: \"caffeinate command-line tool\"",
               "\tDetails: caffeinate asserting on behalf of 'Rscript' (pid 4242)")
  asserts_d <- sub("PreventUserIdleSystemSleep", "PreventUserIdleDisplaySleep", asserts)
  et <- data.frame(name = c("fast", "slow"), user = c(0.1, 1.5), system = c(0, 0.3), elapsed = c(0.1, 1.6))

  # A small package, for the spelling check and the pkgdown index.
  mini_pkg <- function(label, rd_description = "Returns a value.", topics = c("a", "b")) {
    d <- file.path(ctrl_dir, label); dir.create(file.path(d, "man"), recursive = TRUE)
    writeLines(c("Package: pkctrl", "Version: 0.0.1", "Title: A Control Package", "Description: A control for the release gate.",
                 "License: GPL-3", "Authors@R: person('A', 'B', email = 'a@b.c', role = c('aut', 'cre'))",
                 "Encoding: UTF-8", "Language: en-US", "URL: https://example.invalid/pkctrl/"), file.path(d, "DESCRIPTION"))
    writeLines(c("export(a)", "export(b)"), file.path(d, "NAMESPACE"))
    for (x in c("a", "b")) writeLines(c(sprintf("\\name{%s}", x), sprintf("\\alias{%s}", x), sprintf("\\title{The %s Function}", x),
                                        sprintf("\\description{%s}", if (x == "a") rd_description else "Returns a value."),
                                        sprintf("\\usage{%s()}", x)), file.path(d, "man", paste0(x, ".Rd")))
    writeLines(c("url: https://example.invalid/pkctrl/", "reference:", "- title: Functions", "  contents:", paste0("  - ", topics)),
               file.path(d, "_pkgdown.yml"))
    d
  }

  # A letter and a log for the letter gates.
  letter <- function(status, code_extra = character(), prose_extra = character())
    c("# Letter", "", "## R CMD check results", "", "### Local", "", "On the local machine:", "", paste0("    ", status), code_extra, "",
      "* \"Found the following (possibly) invalid URLs: URL:", "  https://example.invalid/site ... Status: 404\". The site.",
      "* \"Skipping checking HTML validation: 'tidy' doesn't look like recent enough", "  HTML Tidy.\" An old Tidy.", "",
      prose_extra, "", "### win-builder", "", "Possibly misspelled words in DESCRIPTION: Welch.", "",
      "## Reverse dependencies", "", "None.")
  letter_log <- c(log_lines[2:9], log_lines[19:22], "* DONE", "", "Status: 2 NOTEs")
  ll <- function(...) letter_log_problems(letter(...), letter_log, "Status: 2 NOTEs")

  c(
    cran_forbidden = all(vapply(seq_along(bad_44), function(i) length(cran_forbidden(bad_44[i], "bad")) > 0L, logical(1))) &&
      !length(cran_forbidden(good_44, "good")),
    stored_rounding = length(stored_rounding(rnd_bad, "bad")) == 2L && !length(stored_rounding(rnd_good, "good")),
    seed_defaults = length(seed_problems(seed_bad, "bad")) == 2L && !length(seed_problems(seed_good, "good")),
    example_seed_rule = safe(
      length(example_seed_problems(c(sprintf("%s(2017)", s_seed_fn), "x <- 1"), c("s <- f(seed = 113)"), "ex", "other")) == 2L &&
        !length(example_seed_problems(c(sprintf("%s(113)", s_seed_fn), "ss_buc_sensitivity(p, 2, seed = 113)"),
                                      c("withr::local_seed(2017)", sprintf("%s(99)", s_seed_fn)), "ex", "other")) &&
        length(r_file_seed_problems(c("#' @examples", sprintf("#' %s(2017)", s_seed_fn), "#' @export",
                                      "f <- function() withr::local_seed(113)"), "planted.R")) == 2L &&
        !length(r_file_seed_problems(c("#' @examples", sprintf("#' %s(113)", s_seed_fn), "#' @export",
                                       "f <- function(seed = NULL) withr::local_seed(2017)"), "planted.R")) &&
        length(example_seed_problems(rd_examples_text(tools::parse_Rd(rd_file("seed_bad", sprintf("%s(2017)", s_seed_fn)))),
                                     character(), "seed_bad.Rd", "")) == 1L &&
        !length(example_seed_problems(rd_examples_text(tools::parse_Rd(rd_file("seed_ok", sprintf("%s(113)", s_seed_fn)))),
                                      character(), "seed_ok.Rd", ""))),
    stale_citation = length(stale_hits(c("Anderson (2024,", paste0("#' ", strsplit(stale_words[1L], " ")[[1L]][1L]),
                                         paste0("#' ", strsplit(stale_words[1L], " ")[[1L]][2L], ").")))) == 1L &&
      !length(stale_hits("Anderson and Kelley (2024), Psychological Methods, 29, 844-867.")),
    rd_math_markup = rd_math_markup(rd_bad) && !rd_math_markup(rd_good),
    # R's own Rd parser must read each planted tag as that tag, so a
    # misspelled tag fails here rather than finding nothing in man/.
    example_blocks = safe({
      f_dt <- rd_file("blk_dt", c(paste0(tag_dt, "{"), "x <- 1", "}"))
      f_dr <- rd_file("blk_dr", c(paste0(tag_dr, "{"), "x <- 1", "}"))
      f_ok <- rd_file("blk_ok", "x <- 1")
      tag_dt %in% rd_tags(suppressWarnings(tools::parse_Rd(f_dt))) && tag_dr %in% rd_tags(suppressWarnings(tools::parse_Rd(f_dr))) &&
        length(example_block_hits(c(f_dt, f_dr, f_ok))) == 2L && !length(example_block_hits(f_ok))
    }),
    non_ascii = safe(length(nonascii_hits(planted("nonascii_bad.R", c("x <- 1", paste0("y <- '", chr(0xE9), "'"))))) == 1L &&
                       !length(nonascii_hits(planted("nonascii_ok.R", "x <- 'plain'")))),
    dashes = safe(length(grep_lines(dash_rx, planted("dashes_bad.md", c(paste0("a ", chr(0x2014), " b"), paste0("1", chr(0x2013), "2"),
                                                                        "1--2 and a-b")))) == 2L &&
                    !length(grep_lines(dash_rx, planted("dashes_ok.md", "1--2 and a-b")))),
    lockstep = length(lockstep_problems("2.0.0", "2026-09-29", c("# BUCSS 1.9.9", ""), "R package version 1.9.9", 'x <- "1.9.9"',
                                        c("version: 1.9.9", 'date-released: "2019-01-01"'),
                                        '{"version": "1.9.9", "publication_date": "2019-01-01"}')) == 7L &&
      length(lockstep_problems("2.0.0", "2026-09-29", "# BUCSS 2.0.0", "R package version 2.0.0", "meta$Version",
                               c("version: 2.0.0", "date-released: 2019-01-01"),
                               c("{", '  "version": "2.0.0",', '  "publication_date": "2019-01-01"', "}"))) == 2L &&
      !length(lockstep_problems("2.0.0", "2026-09-29", "# BUCSS 2.0.0", "R package version 2.0.0",
                                c("meta$Version", 'doi = "10.1080/00273171.2017.1289361"'),
                                c("version: 2.0.0", 'date-released: "2026-09-29"'),
                                c("{", '  "version": "2.0.0",', '  "publication_date": "2026-09-29"', "}"))),
    doc_staleness = local({
      a <- file.path(work, "doc_a"); b <- file.path(work, "doc_b")
      for (d in c(a, b)) { dir.create(file.path(d, "man"), recursive = TRUE); writeLines("x", file.path(d, "NAMESPACE"))
                           writeLines("x", file.path(d, "DESCRIPTION")); writeLines("\\name{x}", file.path(d, "man", "x.Rd")) }
      same <- doc_differences(a, b)
      writeLines("\\name{y}", file.path(b, "man", "x.Rd")); writeLines("\\name{z}", file.path(b, "man", "z.Rd"))
      !length(same) && setequal(doc_differences(a, b), c("man/x.Rd", "man/z.Rd"))
    }),
    rd2pdf = length(rd2pdf_problems(0L, TRUE, c("This is pdfTeX", "! Missing $ inserted.", "l.42 x"))) == 1L &&
      length(rd2pdf_problems(1L, FALSE, c("Error in running tools::texi2pdf()", "LaTeX errors:", "! Undefined control sequence."))) == 4L &&
      !length(rd2pdf_problems(0L, TRUE, c("This is pdfTeX", "Output written on Rd2.pdf (30 pages, 200000 bytes)."))),
    spelling = safe({
      good <- spelling_result(mini_pkg("spell_good"))
      bad <- spelling_result(mini_pkg("spell_bad", rd_description = "The planer recieves teh prior study."))
      is.data.frame(good) && !nrow(good) && is.data.frame(bad) && all(c("recieves", "teh") %in% bad$word)
    }),
    pkgdown = safe(isTRUE(pkgdown_problem(mini_pkg("pkgdown_good"))) &&
                     is.character(pkgdown_problem(mini_pkg("pkgdown_bad", topics = "a")))),
    dirty_tree = isTRUE(git_ctrl),
    awake = kept_awake(asserts, list(pids = c(4242L, 1L), cmds = c("Rscript x", "zsh"))) &&
      !kept_awake(asserts, list(pids = c(5151L, 1L), cmds = c("Rscript x", "zsh"))) &&
      !kept_awake(asserts_d, list(pids = c(4242L, 1L), cmds = c("Rscript x", "zsh"))) &&
      kept_awake(character(), list(pids = 1:2, cmds = c("Rscript x", "/usr/bin/caffeinate -dims Rscript x"))) &&
      !kept_awake(character(), list(pids = 1:2, cmds = c("Rscript x", "caffeinate -d Rscript x"))),
    check_log_parser = identical(as.integer(n_found), unname(as.integer(pl$counts))) &&
      identical(unname(pl$counts), c(0L, 1L, 3L)),
    note_allowlist = identical(fail_heads(cl_200), c("CRAN incoming feasibility", "Rd cross-references")) &&
      any(grepl("example.invalid/site", unlist(lapply(cl_200$fails, `[[`, "lines")), fixed = TRUE)) &&
      length(cl_200$skips) == 1L && length(cl_200$transient) == 1L &&
      any(grepl("example.invalid/busy", unlist(lapply(cl_503$fails, `[[`, "lines")), fixed = TRUE)) &&
      !length(classify_notes(parse_check_log(c("* checking CRAN incoming feasibility ... NOTE", "Maintainer: 'A'",
                                               "New submission", "* DONE", "Status: 1 NOTE"))$steps,
                             on_cran = FALSE, wordlist = character(), probe = function(u) "200")$fails),
    # A check that never reached its end fails every verdict on its output.
    check_ran = safe(all(oks(cv_good)) && all(vapply(cv_dead, function(v) !any(oks(v)), NA))),
    second_example_pass = safe(length(second_steps) >= 1L && !oks(cv_second)[["second_pass"]] &&
                                 all(oks(cv_second)[names(oks(cv_second)) != "second_pass"]) &&
                                 length(cv_second$gates$second_pass$detail) == length(second_steps)),
    math_skip = safe(length(math_line) == 1L && !oks(cv_math)[["math"]] && oks(cv_math)[["ran"]] && oks(cv_math)[["other_skip"]]),
    other_skip = safe(!oks(cv_other)[["other_skip"]] && oks(cv_other)[["math"]] &&
                        identical(cv_other$gates$other_skip$detail, "Skipping checking the planted step: no reason") &&
                        oks(cv_tidy)[["other_skip"]] && oks(cv_tidy)[["math"]] && length(cv_tidy$skips) == 1L),
    manuals = safe(!oks(cv_nopdf)[["manuals"]] && oks(cv_good)[["manuals"]]),
    example_timings = identical(slow_pages(et, page_ceiling_s)$name, "slow"),
    tarball_list = setequal(tarball_problems(c("BUCSS/", "BUCSS/DESCRIPTION", "BUCSS/R/a.R", "BUCSS/tools/x.R",
                                               paste0("BUCSS/inst/", conv_file), "BUCSS/inst/doc/b_QC_2026.docx"), "BUCSS"),
                            c("not a package surface: tools", paste0("private or build file: inst/", conv_file),
                              "private or build file: inst/doc/b_QC_2026.docx")) &&
      !length(tarball_problems(c("BUCSS/DESCRIPTION", "BUCSS/inst/tools/verify_by_simulation.R", "BUCSS/build/vignette.rds"), "BUCSS")),
    ai_content = safe(local({
      d <- file.path(ctrl_dir, "ai_content"); dir.create(d)
      writeLines(c("plain", ai_trailer), file.path(d, "a.txt"))
      con <- file(file.path(d, "b.bin"), "wb"); writeBin(c(as.raw(0L), charToRaw(ai_trailer)), con); close(con)
      writeLines("Fix the Welch floor", file.path(d, "c.txt"))
      identical(ai_content_hits(d), "a.txt:2")
    })),
    history_names = all(grepl(private_doc_rx, c(conv_file, paste0(conv_dir, "/settings.json"), "a/BUCSS_QC_deep.md",
                                                "notes/x.md", "Decision_memo.docx"), ignore.case = TRUE, perl = TRUE)) &&
      !any(grepl(private_doc_rx, c("R/ss_buc_R2.R", "cran-comments.md", "inst/tools/verify_by_simulation.R"),
                 ignore.case = TRUE, perl = TRUE)),
    ai_attribution = grepl(ai_rx, paste0(chr(67, 111, 45, 65, 117, 116, 104, 111, 114, 101, 100, 45, 66, 121), ": X <x@y.z>"),
                           ignore.case = TRUE) &&
      grepl(ai_rx, paste0("Author <noreply@", chr(97, 110, 116, 104, 114, 111, 112, 105, 99), ".com>"), ignore.case = TRUE) &&
      !grepl(ai_rx, "Fix the Welch floor", ignore.case = TRUE),
    identities = identical(identity_problems(ids_bad, people), ids_bad[4:5]),
    # A repository with a private file and an AI trailer and identity on a
    # branch that can be pushed, and a private file on a local-only private
    # branch, which is not swept.
    ref_sweep = !is.na(git_bin) && safe(local({
      d <- file.path(ctrl_dir, "sweep_repo"); dir.create(d)
      g <- function(...) git_in(d, "-c", "commit.gpgsign=false", ...)
      g("init", "-q")
      writeLines("x", file.path(d, "a.R")); g("add", "a.R")
      g("-c", shQuote("user.name=Ken Kelley"), "-c", "user.email=kkelley@nd.edu", "commit", "-q", "-m", "First")
      writeLines("x", file.path(d, conv_file)); g("add", conv_file)
      g("-c", "user.name=Somebody", "-c", paste0("user.email=", ai_email), "commit", "-q", "-m", "Second", "-m", shQuote(ai_trailer))
      g("checkout", "-q", "-b", "draft-private")
      dir.create(file.path(d, "notes")); writeLines("x", file.path(d, "notes", "plan.md")); g("add", "notes/plan.md")
      g("-c", shQuote("user.name=Ken Kelley"), "-c", "user.email=kkelley@nd.edu", "commit", "-q", "-m", "Third")
      s <- ref_sweep(d)
      identical(s$private_ever, conv_file) && length(s$msg_hits) == 1L && length(s$local_private) == 1L &&
        identical(identity_problems(s$ids, people), sprintf("Somebody <%s>", ai_email))
    })),
    harness_parser = !harness_verdict(c("BUCSS 2.0.0 from /L/BUCSS", "mode: full", "  [PASS] a", "  [FAIL] b", "2 checks, 1 failed."), 1L, "/L")$ok &&
      !harness_verdict(character(), 0L, "/L")$ok &&
      harness_verdict(c("BUCSS 2.0.0 from /L/BUCSS", "mode: full", rep("  [PASS] a", expected_harness_checks),
                        sprintf("%d checks, 0 failed.", expected_harness_checks)), 0L, "/L")$ok,
    tests_parser = !tests_verdict("TESTS files=1 tests=2 expectations=5 failed=1 errors=0 warnings=0 skipped=0")$ok &&
      !tests_verdict(character())$ok &&
      tests_verdict("TESTS files=1 tests=2 expectations=5 failed=0 errors=0 warnings=0 skipped=0")$ok,
    stale_tarballs = safe(local({
      d1 <- file.path(ctrl_dir, "stale_1"); d2 <- file.path(ctrl_dir, "stale_2"); dir.create(d1); dir.create(d2)
      writeLines("old", file.path(d1, "BUCSS_9.9.9.tar.gz")); writeLines("other", file.path(d2, "BUCSS_9.9.8.tar.gz"))
      s <- stale_tarballs(c(d1, d2, d1, file.path(ctrl_dir, "no_such_dir")), "BUCSS_9.9.9.tar.gz")
      length(s) == 1L && basename(dirname(s)) == "stale_1" && !length(stale_tarballs(d2, "BUCSS_9.9.9.tar.gz"))
    })),
    # The hand-off to a directory that does not exist, onto a read-only stale
    # file, and to a fresh directory: only the last may count.
    hand_off = safe(local({
      from <- planted("handoff/src/BUCSS_9.9.9.tar.gz", "certified bytes")
      want <- unname(tools::sha256sum(from))
      fresh <- file.path(ctrl_dir, "handoff", "fresh"); dir.create(fresh)
      ro <- file.path(ctrl_dir, "handoff", "readonly"); dir.create(ro)
      writeLines("stale bytes", file.path(ro, "BUCSS_9.9.9.tar.gz")); Sys.chmod(file.path(ro, "BUCSS_9.9.9.tar.gz"), "0444")
      gone <- hand_off(from, file.path(ctrl_dir, "handoff", "no", "such"), "BUCSS_9.9.9.tar.gz", want)
      blocked <- hand_off(from, ro, "BUCSS_9.9.9.tar.gz", want)
      good <- hand_off(from, fresh, "BUCSS_9.9.9.tar.gz", want)
      !gone$ok && !blocked$ok && any(grepl("NOT the certified tarball", blocked$problems, fixed = TRUE)) &&
        good$ok && identical(good$sha, want)
    })),
    letter_static = safe(length(letter_static_problems(c(
      "# Letter", paste0("A dash ", chr(0x2014), " here."), sprintf("Anderson (2024, %s).", stale_words[1L]),
      sprintf("Results: %s", s_tbd), sprintf("Win-builder is yet %s.", s_tbr), sprintf("%s: x]]", s_fill),
      sprintf("[%s] check", s_todo), sprintf("Size: %s KB.", s_xxx)))) == 7L &&
        !length(letter_static_problems(c("# Letter", "Plain text, pages 1--2, a well-known a-b.", "    Status: OK")))),
    letter_log = safe(!length(ll("Status: 2 NOTEs")) &&
                        any(grepl("differs from the check's", ll("Status: 1 NOTE"), fixed = TRUE)) &&
                        any(grepl("not in the log: * checking tests", ll("Status: 2 NOTEs", code_extra = "    * checking tests ... OK"), fixed = TRUE)) &&
                        any(grepl("invalid DOIs", ll("Status: 2 NOTEs", prose_extra = "* \"Found the following (possibly) invalid DOIs: DOI: 10.1/x\"."))) &&
                        any(grepl("'Possibly misspelled words'", ll("Status: 2 NOTEs", prose_extra = "Possibly misspelled words in DESCRIPTION are A and B."),
                                  fixed = TRUE)) &&
                        any(grepl("no Status line", letter_log_problems(c("## R CMD check results", "", "All clean."), letter_log,
                                                                        "Status: 2 NOTEs"), fixed = TRUE))),
    example_ceilings = ceiling_ctrl$ok
  )
})
ctrl_why <- list(example_ceilings = ceiling_ctrl$why)
ctrl_failed <- names(controls)[!is.na(controls) & !controls]
ctrl_inconclusive <- names(controls)[is.na(controls)]
# Quick mode on an uncommitted tree may leave a control on the committed
# sources inconclusive (they may predate the tools that check them); it is
# then reported as skipped, never as passed. Anywhere else it fails.
ctrl_may_skip <- quick && length(dirty) > 0L
gate(sprintf("detector positive controls (%d detectors%s, each run on planted input: it flags its known-bad case and passes its known-good one)",
             length(controls), if (ctrl_may_skip && length(ctrl_inconclusive)) sprintf(", %d skipped below", length(ctrl_inconclusive)) else ""),
     !length(ctrl_failed) && (!length(ctrl_inconclusive) || ctrl_may_skip),
     c(unlist(lapply(ctrl_failed, function(n) c(sprintf("control failed: %s", n), if (length(ctrl_why[[n]])) paste0("  ", ctrl_why[[n]])))),
       if (!ctrl_may_skip) unlist(lapply(ctrl_inconclusive, function(n) c(sprintf("control inconclusive: %s", n), paste0("  ", ctrl_why[[n]]))))))
if (ctrl_may_skip) for (n in ctrl_inconclusive)
  gate_skip(sprintf("positive control %s, on the committed sources", n),
            c(ctrl_why[[n]], "quick mode on an uncommitted tree: commit, and the control runs"), label = "SKIPPED")

## ---- 1. The working tree and the commit; full mode's refusals ------------------------
fact("commit ", if (length(head_sha)) head_sha else "(none)", "; working tree ",
     if (length(dirty)) sprintf("has %d uncommitted change(s)", length(dirty)) else "clean")
stale_at_start <- stale_tarballs(c(dest_dir, root), tarball_name)
if (!quick) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    assertions <- tryCatch(suppressWarnings(system2("pmset", c("-g", "assertions"), stdout = TRUE, stderr = FALSE)),
                           error = function(e) character())
    awake <- kept_awake(assertions, process_chain())
    gate("this run is kept awake by caffeinate -i (the timing gates read wall-clock times)", awake,
         if (!awake) "start it with: caffeinate -i Rscript tools/release_gate.R")
    if (!awake) { cat("\nGATE REFUSED: full mode runs only under caffeinate -i; nothing recorded\n"); quit(status = 1L) }
  }
  clean <- !is.na(git_bin) && length(head_sha) == 1L && !length(dirty) && is.null(attr(dirty, "status"))
  gate("working tree clean (the gate certifies a commit, not a draft)", clean, head(dirty, 40L))
  if (!clean) { cat("\nGATE REFUSED: commit or discard every change first; nothing recorded\n"); quit(status = 1L) }
  # The hand-off, checked before the long run rather than after it.
  dest_ok <- dir.exists(dest_dir) && file.access(dest_dir, 2L) == 0L
  gate(sprintf("hand-off directory %s exists and can be written", dest_dir), dest_ok,
       if (!dest_ok) "give an existing, writable directory with --dest=DIR")
  gate(if (length(stale_at_start) && replace_stale)
         sprintf("%d stale %s reported (--replace-stale given)", length(stale_at_start), tarball_name)
       else sprintf("no %s already at the hand-off directory or in the package root", tarball_name),
       !length(stale_at_start) || replace_stale,
       c(if (length(stale_at_start)) describe_file(stale_at_start),
         if (length(stale_at_start) && replace_stale)
           c("none is deleted; on a pass the certified tarball overwrites the one at the hand-off path, if any,",
             "and any other stays where it is, for you to delete")
         else if (length(stale_at_start))
           "an uncertified file with the certified file's name invites the wrong upload: delete or move it, or rerun with --replace-stale"))
} else if (length(stale_at_start)) {
  fact("already present (full mode fails on it unless --replace-stale is given): ", paste(stale_at_start, collapse = ", "))
}

## ---- 2. Static invariants ---------------------------------------------------------------
r_files <- files_in("R", "[.][Rr]$")
tool_files <- files_in("inst/tools", "[.][Rr]$")
rd_files <- files_in("man", "[.]Rd$")
vig_files <- files_in("vignettes", "[.]Rmd$")
test_files <- files_in("tests", "[.][Rr]$", recursive = TRUE)
opt_file <- function(f) { p <- file.path(root, f); if (file.exists(p)) p else character() }
parsed <- lapply(stats::setNames(c(r_files, tool_files), rel(c(r_files, tool_files))),
                 function(f) tryCatch(parse(f, keep.source = FALSE), error = function(e) e))
unparsable <- names(Filter(function(x) inherits(x, "error"), parsed))
gate("every file in R/ and inst/tools/ parses", !length(unparsable), unparsable)
parsed <- Filter(function(x) !inherits(x, "error"), parsed)

# Version lockstep.
news <- read_lines(file.path(root, "NEWS.md"))
cit_path <- file.path(root, "inst", "CITATION")
cit_text <- if (file.exists(cit_path)) tryCatch(
  unlist(lapply(utils::readCitationFile(cit_path, meta = as.list(desc[1L, ])), function(b) format(b, style = "text"))),
  error = function(e) sprintf("inst/CITATION does not read: %s", conditionMessage(e))) else NULL
cff <- if (length(opt_file("CITATION.cff"))) read_lines(opt_file("CITATION.cff")) else NULL
zen <- if (length(opt_file(".zenodo.json"))) read_lines(opt_file(".zenodo.json")) else NULL
desc_date_str <- if ("Date" %in% colnames(desc)) unname(desc[1L, "Date"]) else NA_character_
lock <- lockstep_problems(pkg_version, desc_date_str, news, cit_text, if (file.exists(cit_path)) read_lines(cit_path) else NULL, cff, zen)
gate(sprintf("version %s in lockstep (DESCRIPTION, NEWS.md, inst/CITATION%s%s)%s", pkg_version,
             if (!is.null(cff)) ", CITATION.cff" else "", if (!is.null(zen)) ", .zenodo.json" else "",
             if (!is.null(cff) || !is.null(zen)) sprintf(", and the release date %s in %s", desc_date_str,
                                                         paste(c(if (!is.null(cff)) "CITATION.cff", if (!is.null(zen)) ".zenodo.json"),
                                                               collapse = " and ")) else ""),
     !length(lock), lock)

# Stale citation phrases in every prose surface.
prose <- c(r_files, rd_files, vig_files, opt_file("DESCRIPTION"), opt_file("README.md"), opt_file("NEWS.md"),
           opt_file("inst/CITATION"), opt_file("CITATION.cff"), opt_file(".zenodo.json"))
stale <- unlist(lapply(prose, function(f) { h <- stale_hits(read_lines(f)); if (length(h)) sprintf("%s:%d", rel(f), h) }))
gate(sprintf("no stale citation status (%d phrases, matched across line breaks) in %d prose files",
             length(stale_words), length(prose)), !length(stale), stale)

# Stored rounding in package code.
rnd <- unlist(lapply(names(parsed)[startsWith(names(parsed), "R/")], function(n) stored_rounding(parsed[[n]], n)))
gate("no stored rounding: round()/signif() only in print and format helpers, messages, and integer comparisons",
     !length(rnd), rnd)

# Seeds: no constant default, no base seed setter in package code or the
# shipped harness, and the example seed only in help-page examples.
sp <- unlist(lapply(names(parsed), function(n) seed_problems(parsed[[n]], n)))
gate(sprintf("no constant seed default and no %s() call in R/ or inst/tools/", s_seed_fn), !length(sp), sp)
seed_rule <- character()
for (f in r_files) seed_rule <- c(seed_rule, r_file_seed_problems(read_lines(f), rel(f)))
rd_db <- tools::Rd_db(dir = root)
for (n in names(rd_db)) seed_rule <- c(seed_rule, example_seed_problems(rd_examples_text(rd_db[[n]]), character(), paste0("man/", n), ""))
for (f in c(test_files, vig_files, tool_files, opt_file("README.md")))
  seed_rule <- c(seed_rule, example_seed_problems(character(), read_lines(f), "", rel(f)))
gate(sprintf("help-page examples seed with %d, and nothing else does", example_seed), !length(seed_rule), seed_rule)

# The patterns CRAN forbids, in package code and in the shipped script.
forbidden <- unlist(lapply(c(r_files, tool_files), function(f) cran_forbidden(read_lines(f), rel(f))))
gate(sprintf("R/ and inst/tools/ (%d files) free of the patterns CRAN forbids (global state, superassignment, %s(%s), quitting)",
             length(r_files) + length(tool_files), s_opts, s_warn), !length(forbidden), forbidden)

# Rd markup inside math, parsed.
hits <- names(Filter(rd_math_markup, rd_db))
gate(sprintf("no Rd markup inside \\eqn{} or \\deqn{} (%d pages, parsed)", length(rd_db)), !length(hits), hits)

# Example blocks CRAN's pretest runs twice or never.
blocks <- example_block_hits(c(r_files, rd_files))
gate(sprintf("no %s{} or %s{} in R/ or man/", tag_dt, tag_dr), !length(blocks), blocks)

# Characters.
nonascii <- nonascii_hits(c(r_files, tool_files, opt_file("DESCRIPTION"), opt_file("NAMESPACE"), opt_file("inst/CITATION")))
gate("no non-ASCII characters in R/, inst/tools/, DESCRIPTION, NAMESPACE, or inst/CITATION", !length(nonascii), nonascii)
dashes <- grep_lines(dash_rx, prose)
gate("no em or en dashes in any prose file (house style)", !length(dashes), dashes)

# The letter to CRAN: no dash, no stale citation status, no placeholder, and
# a check results section that quotes a Status line. Whether that line and
# the rest of what it quotes match the check's log is full mode's to say.
letter_path <- file.path(root, "cran-comments.md")
letter_lines <- if (file.exists(letter_path)) read_lines(letter_path) else NULL
letter_idx <- if (is.null(letter_lines)) integer() else letter_check_section(letter_lines)
letter_status <- if (is.null(letter_lines)) character() else letter_status_lines(letter_lines, letter_idx)
letter_static <- if (is.null(letter_lines)) "cran-comments.md not found" else c(
  letter_static_problems(letter_lines),
  if (!length(letter_idx)) "cran-comments.md has no 'R CMD check results' section"
  else if (!length(letter_status)) "cran-comments.md: its check results section quotes no Status line")
gate("letter to CRAN (cran-comments.md): no em or en dash, stale citation status, or placeholder; a Status line in its check results",
     !length(letter_static), letter_static)

## ---- 3. Documentation current: roxygen into a copy, with the pinned version ------------------
rox_pin <- if ("Config/roxygen2/version" %in% colnames(desc)) unname(desc[1L, "Config/roxygen2/version"]) else NA_character_
rox_have <- tryCatch(as.character(utils::packageVersion("roxygen2")), error = function(e) NA_character_)
if (is.na(rox_pin) || !identical(rox_pin, rox_have)) {
  gate(sprintf("roxygen2 %s installed, as DESCRIPTION pins", rox_pin), FALSE,
       sprintf("installed: %s; a different roxygen2 rewrites the Rd files, so the comparison would mean nothing", rox_have))
} else {
  rox_copy <- file.path(work, "roxygen_copy")
  dir.create(rox_copy)
  for (x in c("DESCRIPTION", "NAMESPACE", "R", "man", "man-roxygen", "inst"))
    if (file.exists(file.path(root, x))) file.copy(file.path(root, x), rox_copy, recursive = TRUE)
  rox_log <- file.path(work, "roxygen.log")
  rox_st <- system2(rscript_bin, c("-e", shQuote(sprintf("roxygen2::roxygenise('%s')", rox_copy))),
                    stdout = rox_log, stderr = rox_log)
  stale_docs <- doc_differences(root, rox_copy)
  gate(sprintf("documentation current (roxygen2 %s regenerates man/, NAMESPACE, and DESCRIPTION unchanged)", rox_pin),
       identical(as.integer(rox_st), 0L) && !length(stale_docs),
       c(if (!identical(as.integer(rox_st), 0L)) c(sprintf("roxygenise exit status %s", rox_st), tail(readLines(rox_log), 8L)),
         if (length(stale_docs)) c("would change (regenerate with roxygen2 and commit):", stale_docs)))
}

## ---- 4. PDF manual, spelling, pkgdown ------------------------------------------------------
pdf_out <- file.path(work, "manual.pdf")
pdf_log <- file.path(work, "Rd2pdf.log")
old_wd <- setwd(work)   # Rd2pdf writes its scratch directory into the working directory
pdf_st <- system2(r_bin, c("CMD", "Rd2pdf", "--no-preview", "--force", paste0("--output=", shQuote(pdf_out)), shQuote(root)),
                  stdout = pdf_log, stderr = pdf_log)
setwd(old_wd)
pdf_problems <- rd2pdf_problems(pdf_st, file.exists(pdf_out), readLines(pdf_log, warn = FALSE))
gate("PDF manual renders (R CMD Rd2pdf, a full LaTeX pass)", !length(pdf_problems), pdf_problems)
spell <- spelling_result(root)
if (inherits(spell, "error")) {
  gate("spelling (inst/WORDLIST)", FALSE, conditionMessage(spell))
} else {
  gate(sprintf("spelling (inst/WORDLIST): no unknown words in DESCRIPTION, man/, and vignettes/"), nrow(spell) == 0L,
       if (nrow(spell)) sprintf("%s: %s", spell$word, vapply(spell$found, paste, "", collapse = ", ")))
}
pk <- pkgdown_problem(root)
gate("pkgdown reference index complete (every help topic listed in _pkgdown.yml)", isTRUE(pk), if (!isTRUE(pk)) pk)

## ---- 5. The example ceilings --------------------------------------------------------------
ceil_log <- file.path(work, "ceilings.log")
ceil_st <- system2(rscript_bin, c(shQuote(file.path(root, "tools", "check_example_ceilings.R")), paste0("--pkg=", shQuote(root))),
                   stdout = ceil_log, stderr = ceil_log)
ceil_lines <- readLines(ceil_log, warn = FALSE)
ceil_summary <- grep("^ceiling comments checked:", ceil_lines, value = TRUE)
gate(sprintf("example ceilings (tools/check_example_ceilings.R exits 0): %s",
             if (length(ceil_summary)) sub("^ceiling comments checked: ", "", ceil_summary) else "no summary line"),
     identical(as.integer(ceil_st), 0L) && length(ceil_summary) == 1L &&
       grepl(sprintf("checked: %d of %d expected", expected_ceiling_comments, expected_ceiling_comments), ceil_summary, fixed = TRUE),
     grep("^(MISMATCH|PROBLEM)|^Error", ceil_lines, value = TRUE))

if (quick) {
  # What the letter quotes from the check can be compared only with a check.
  if (!is.null(letter_lines))
    gate_skip(sprintf("letter's check results match the check's log (it quotes %s)",
                      if (length(letter_status)) paste(sprintf("'%s'", letter_status), collapse = ", ") else "no Status line"),
              "quick mode runs no check; full mode compares the letter's Status line and quoted output with the check's log",
              label = "NOT YET VERIFIED")
  n_fail <- sum(statuses() == "FAIL")
  cat(sprintf("\nquick mode: %d gates, %d failed, %s; commit %s%s; %.0f s\n", length(state$results), n_fail,
              skip_summary(), sha, if (length(dirty)) sprintf(" plus %d uncommitted change(s)", length(dirty)) else "",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (n_fail) { cat("QUICK GATES FAILED\n"); quit(status = 1L) }
  cat("QUICK GATES PASSED (not a release: run full mode before any tarball leaves)\n")
  quit(status = 0L)
}

## ---- 6. A version and a date that can be released ---------------------------------------------
cran_ver <- tryCatch({
  ap <- utils::available.packages(repos = "https://cloud.r-project.org")
  if (pkg_name %in% rownames(ap)) unname(ap[pkg_name, "Version"]) else "none"
}, error = function(e) NA_character_)
on_cran <- !is.na(cran_ver) && cran_ver != "none"
ver_ok <- !is.na(cran_ver) && (cran_ver == "none" || package_version(pkg_version) > package_version(cran_ver))
gate(sprintf("version %s is above the version on CRAN (%s)", pkg_version,
             if (is.na(cran_ver)) "lookup failed" else cran_ver), ver_ok,
     if (is.na(cran_ver)) "CRAN could not be reached (the URL checks need the network too)"
     else if (!ver_ok) "raise Version in DESCRIPTION and open its NEWS heading")
desc_date <- tryCatch(as.Date(unname(desc[1L, "Date"])), error = function(e) as.Date(NA))
date_age <- as.numeric(Sys.Date() - desc_date)
gate(sprintf("DESCRIPTION Date %s is within three days of certification", format(desc_date)),
     !is.na(date_age) && date_age >= 0 && date_age <= 3,
     if (is.na(date_age)) "DESCRIPTION has no parseable Date"
     else if (date_age < 0 || date_age > 3) sprintf("%s day(s) old; date the release the day it is certified", date_age))

## ---- 7. The tarball, built from the committed tree outside any synchronized folder ------------
home <- normalizePath("~", mustWork = FALSE)
sync_roots <- file.path(home, c("Dropbox", "Library/CloudStorage", "Library/Mobile Documents", "Google Drive", "OneDrive"))
sync_roots <- unique(c(sync_roots, normalizePath(sync_roots[file.exists(sync_roots)])))
in_sync <- any(vapply(sync_roots, function(s) startsWith(normalizePath(work), paste0(s, "/")), logical(1)))
src <- file.path(work, "src")
dir.create(file.path(src, pkg_name), recursive = TRUE)
arch_st <- system(sprintf("%s -C %s archive --format=tar HEAD | tar -x -C %s", shQuote(git_bin), shQuote(root),
                          shQuote(file.path(src, pkg_name))))
bld <- local({ old <- setwd(src); on.exit(setwd(old)); log <- file.path(work, "build.log")
  st <- system2(r_bin, c("CMD", "build", pkg_name), stdout = log, stderr = log); list(status = st, lines = readLines(log, warn = FALSE)) })
tarball <- file.path(src, tarball_name)
built <- identical(as.integer(arch_st), 0L) && identical(as.integer(bld$status), 0L) && file.exists(tarball)
gate(sprintf("tarball built from commit %s (git archive HEAD, in %s)", sha, work), built && !in_sync,
     c(if (in_sync) "the work directory is inside a synchronized folder",
       if (!identical(as.integer(arch_st), 0L)) "git archive failed",
       if (!identical(as.integer(bld$status), 0L)) tail(bld$lines, 10L)))
if (built) fact(sprintf("tarball %s, %.1f KB, sha256 %s", tarball_name, file.info(tarball)$size / 1024,
                        unname(tools::sha256sum(tarball))))

## ---- 8. R CMD check --as-cran, both manuals, remote checks on ----------------------------------
# The check's output (stdout and stderr) and its own 00check.log, read from
# the check directory before the temporary directory is removed at exit.
chk <- list(status = NA_integer_, lines = character(), log00 = character())
chk_wall_s <- NA_real_
check_log_path <- file.path(src, paste0(pkg_name, ".Rcheck"), "00check.log")
if (built) {
  chk_env <- c(`_R_CHECK_CRAN_INCOMING_REMOTE_` = "true", `_R_CHECK_FORCE_SUGGESTS_` = "true",
               `_R_CHECK_TIMINGS_` = "0", `_R_CHECK_XREFS_PKGS_ARE_DECLARED_` = "true")
  old_env <- Sys.getenv(c(names(chk_env), "NOT_CRAN"), unset = NA)
  do.call(Sys.setenv, as.list(chk_env)); Sys.unsetenv("NOT_CRAN")
  chk_t0 <- Sys.time()
  chk <- local({ old <- setwd(src); on.exit(setwd(old)); log <- file.path(work, "check.log")
    st <- system2(r_bin, c("CMD", "check", "--as-cran", "--timings", tarball_name), stdout = log, stderr = log)
    list(status = st, lines = readLines(log, warn = FALSE),
         log00 = if (file.exists(check_log_path)) readLines(check_log_path, warn = FALSE) else character()) })
  chk_wall_s <- as.numeric(difftime(Sys.time(), chk_t0, units = "secs"))
  for (n in names(old_env)) if (is.na(old_env[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old_env[[n]]), n))
}
wordlist <- if (file.exists(file.path(root, "inst", "WORDLIST"))) read_lines(file.path(root, "inst", "WORDLIST")) else character()
cv <- check_verdicts(chk$lines, chk$status, chk$log00, on_cran = on_cran, wordlist = wordlist, probe = probe_url)
parsed_log <- cv$parsed
check_ran <- cv$ran
fact("check ", if (is.na(parsed_log$status_line)) "(no Status line)" else parsed_log$status_line,
     sprintf("; exit status %s; %.0f s wall clock", chk$status, chk_wall_s))
emit <- function(k) gate(cv$gates[[k]]$name, cv$gates[[k]]$ok, cv$gates[[k]]$detail)
emit("ran"); emit("errors")
for (t in cv$transient) fact("transient at check time, re-probed: ", t)
emit("notes")
for (s in cv$skips) gate_skip(s$name, s$lines)
for (k in c("math", "other_skip", "manuals", "second_pass")) emit(k)
# Every later gate that reads the check's output or files passes only if the
# check ran to the end.
no_check <- "R CMD check left no complete log (see the first check gate), so this cannot pass"
check_gate <- function(name, ok, detail = character()) gate(name, check_ran && isTRUE(ok), c(if (!check_ran) no_check, detail))

# The letter against the check's log: its Status line equal to the check's,
# and every line it quotes as check output present in 00check.log.
letter_log <- if (length(chk$log00)) chk$log00 else chk$lines
letter_vs_log <- if (is.null(letter_lines)) "cran-comments.md not found" else
  letter_log_problems(letter_lines, letter_log, parsed_log$status_line)
check_gate(sprintf("letter's check results match the check's log (Status line %s; every line it quotes as check output is in 00check.log)",
                   if (length(letter_status)) paste(sprintf("'%s'", letter_status), collapse = ", ") else "(none)"),
           !length(letter_vs_log), letter_vs_log)
step_times <- grep("\\[[0-9]+s/[0-9]+s\\]", chk$lines, value = TRUE)
for (l in grep("checking (examples|tests|re-building of vignette outputs|PDF version|CRAN incoming)|Running", step_times, value = TRUE))
  fact(trimws(l))

## ---- 9. Time budgets, from the check's own timings file ----------------------------------------
ex_timings <- file.path(src, paste0(pkg_name, ".Rcheck"), paste0(pkg_name, "-Ex.timings"))
et <- if (file.exists(ex_timings)) tryCatch(utils::read.table(ex_timings, header = TRUE), error = function(e) NULL) else NULL
if (is.null(et) || !all(c("name", "user", "system", "elapsed") %in% names(et)) || !nrow(et)) {
  check_gate("example timings file read after the check", FALSE, ex_timings)
} else {
  et$max <- pmax(et$user + et$system, et$elapsed)
  et <- et[order(-et$max), ]
  over <- slow_pages(et, page_ceiling_s)
  check_gate(sprintf("no help page over %.1f s (slowest %.2f s, %s; %d pages, %.1f s in all)", page_ceiling_s,
               et$max[1L], et$name[1L], nrow(et), sum(et$elapsed)), !nrow(over),
       sprintf("%s %.2f s", over$name, over$max))
  cpu_ex <- sum(et$user + et$system); el_ex <- sum(et$elapsed)
  check_gate(sprintf("examples ran awake and unshared (%.1f s elapsed against %.1f s CPU)", el_ex, cpu_ex),
       el_ex <= 1.5 * cpu_ex + 5, if (el_ex > 1.5 * cpu_ex + 5) "elapsed far above CPU: the machine slept or was busy; rerun")
}
feas_line <- grep("checking CRAN incoming feasibility", chk$lines, value = TRUE)[1L]
feas_s <- if (is.na(feas_line)) 0 else as.numeric(sub(".*\\[[0-9]+s/([0-9]+)s\\].*", "\\1", feas_line))
if (is.na(feas_s)) feas_s <- 0
check_gate(sprintf("whole check %.0f s (the network-bound feasibility step excluded), within the %d s budget",
                   chk_wall_s - feas_s, check_budget_s), !is.na(chk_wall_s) && chk_wall_s - feas_s <= check_budget_s)

## ---- 10. The tarball in a scratch library: tests, ceilings, and the shipped harness --------------
gate_lib <- file.path(work, "lib")
dir.create(gate_lib)
inst <- if (built) local({ log <- file.path(work, "install.log")
  st <- system2(r_bin, c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(gate_lib)), shQuote(tarball)),
                stdout = log, stderr = log); list(status = st, lines = readLines(log, warn = FALSE)) }) else list(status = 1L, lines = "no tarball")
installed <- identical(as.integer(inst$status), 0L) && dir.exists(file.path(gate_lib, pkg_name))
gate("tarball installs into a scratch library", installed, if (!installed) tail(inst$lines, 8L))
lib_env <- sprintf("R_LIBS=%s", shQuote(paste(c(gate_lib, .libPaths()), collapse = .Platform$path.sep)))

if (installed) {
  test_script <- file.path(work, "run_tests.R")
  writeLines(c(
    sprintf("tests_dir <- %s", deparse(file.path(src, pkg_name, "tests", "testthat"))),
    "res <- testthat::test_dir(tests_dir, package = 'BUCSS', load_package = 'installed', reporter = 'silent',",
    "                          stop_on_failure = FALSE, stop_on_warning = FALSE)",
    "df <- as.data.frame(res)",
    "cat('BUCSS from', find.package('BUCSS'), '\\n')",
    "bad <- df[df$failed > 0 | df$error | df$warning > 0, c('file', 'test'), drop = FALSE]",
    "for (i in seq_len(nrow(bad))) cat('PROBLEM', bad$file[i], '::', bad$test[i], '\\n')",
    "sk <- df[df$skipped, c('file', 'test'), drop = FALSE]",
    "for (i in seq_len(nrow(sk))) cat('SKIPPED', sk$file[i], '::', sk$test[i], '\\n')",
    "cat(sprintf('TESTS files=%d tests=%d expectations=%d failed=%d errors=%d warnings=%d skipped=%d\\n',",
    "            length(unique(df$file)), nrow(df), sum(df$nb), sum(df$failed), sum(df$error), sum(df$warning), sum(df$skipped)))"),
    test_script)
  tst_log <- file.path(work, "tests.log")
  tst_st <- local({ old <- setwd(work); on.exit(setwd(old))
    system2(rscript_bin, shQuote(test_script), stdout = tst_log, stderr = tst_log, env = c(lib_env, "NOT_CRAN=true")) })
  tst_lines <- readLines(tst_log, warn = FALSE)
  state$test_lines <- tst_lines
  tv <- tests_verdict(tst_lines)
  from_lib <- any(startsWith(normalizePath(sub("^BUCSS from (\\S+).*$", "\\1", grep("^BUCSS from ", tst_lines, value = TRUE)),
                                           mustWork = FALSE), normalizePath(gate_lib)))
  gate(sprintf("full test suite against the tarball (NOT_CRAN=true): %s", tv$summary),
       tv$ok && from_lib && identical(as.integer(tst_st), 0L),
       c(if (!from_lib) "the suite did not load BUCSS from the scratch library",
         grep("^PROBLEM|^Error", tst_lines, value = TRUE), if (!tv$ok) tail(tst_lines, 10L)))
  if (isTRUE(tv$skipped > 0L)) gate_skip(sprintf("%d test(s) skipped with NOT_CRAN=true", tv$skipped),
                                         grep("^SKIPPED", tst_lines, value = TRUE), label = "SKIPPED")

  ceil2_log <- file.path(work, "ceilings_tarball.log")
  ceil2_st <- system2(rscript_bin, c(shQuote(file.path(root, "tools", "check_example_ceilings.R")),
                                     paste0("--pkg=", shQuote(file.path(src, pkg_name))), paste0("--lib=", shQuote(gate_lib))),
                      stdout = ceil2_log, stderr = ceil2_log)
  ceil2 <- readLines(ceil2_log, warn = FALSE)
  gate(sprintf("example ceilings against the tarball: %s",
               sub("^ceiling comments checked: ", "", c(grep("^ceiling comments checked:", ceil2, value = TRUE), "no summary line")[1L])),
       identical(as.integer(ceil2_st), 0L), grep("^(MISMATCH|PROBLEM)|^Error", ceil2, value = TRUE))

  hv_dir <- file.path(work, "harness"); dir.create(hv_dir)   # no DESCRIPTION here, so the harness uses the installed BUCSS
  hv_log <- file.path(work, "harness.log")
  hv_st <- local({ old <- setwd(hv_dir); on.exit(setwd(old))
    system2(rscript_bin, c("-e", shQuote("source(system.file('tools', 'verify_by_simulation.R', package = 'BUCSS', mustWork = TRUE))")),
            stdout = hv_log, stderr = hv_log, env = c(lib_env, "BUCSS_VERIFY=full")) })
  hv_lines <- readLines(hv_log, warn = FALSE)
  state$harness_lines <- hv_lines
  hv <- harness_verdict(hv_lines, hv_st, gate_lib)
  gate(sprintf("shipped harness, full mode, against the tarball: %s", hv$summary), hv$ok,
       c(hv$problems, grep("\\[FAIL\\]|^Error", hv_lines, value = TRUE)))
}

## ---- 11. Nothing in the tarball or the history that should not be -----------------------------------
tar_paths <- if (built) utils::untar(tarball, list = TRUE) else character()
tp <- tarball_problems(tar_paths, pkg_name)
gate(sprintf("tarball holds only package surfaces and no private or build files (%d entries)", length(tar_paths)),
     built && length(tar_paths) > 0L && !length(tp), tp)
content_hits <- character()
if (built) {
  x_dir <- file.path(work, "tarball_contents"); dir.create(x_dir)
  utils::untar(tarball, exdir = x_dir)
  content_hits <- ai_content_hits(x_dir)
}
gate("no AI attribution anywhere in the tarball's files", built && !length(content_hits), content_hits)

sw <- ref_sweep(root)
fact("swept refs: ", paste(sw$sweep_refs, collapse = ", "))
if (length(sw$site_refs)) fact("generated site branch, not swept: ", paste(sw$site_refs, collapse = ", "))
if (length(sw$local_private)) fact("local-only private branch, not swept: ", paste(sw$local_private, collapse = ", "))
gate(sprintf("no private document ever tracked on a ref that can be pushed (%d refs, %d file names)", length(sw$sweep_refs), length(sw$ever)),
     length(sw$sweep_refs) > 0L && length(sw$ever) > 0L && !length(sw$private_ever), sw$private_ever)
gate(sprintf("no AI attribution in any commit message (%d message lines)", length(sw$msgs)),
     length(sw$msgs) > 0L && !length(sw$msg_hits), unique(sw$msg_hits))
ids <- sw$ids
auth <- eval(parse(text = desc[1L, "Authors@R"]))
people <- lapply(auth, function(p) list(given = paste(p$given, collapse = " "), family = p$family))
bad_ids <- identity_problems(ids, people)
fact("commit identities: ", paste(ids, collapse = "; "))
gate("every author and committer is a person in Authors@R (or GitHub's web committer), and none is an AI assistant",
     length(ids) > 0L && !length(bad_ids), bad_ids)

## ---- Hand-off, and record ---------------------------------------------------------------------
# Only when every gate has passed is the tarball copied out, and the copy
# counts only if file.copy() says so and the file at the hand-off path then
# has the certified sha256. The one file the gate may overwrite is the one at
# that path; no other tarball is touched. The record is written after, so it
# carries the hand-off.
sha256 <- if (built) unname(tools::sha256sum(tarball)) else NA_character_
handoff <- NULL
if (!sum(statuses() == "FAIL")) {
  handoff <- hand_off(tarball, dest_dir, tarball_name, sha256)
  gate(sprintf("certified tarball handed off to %s, its sha256 read back and verified", handoff$path), handoff$ok, handoff$problems)
}
n_fail <- sum(statuses() == "FAIL")
ok_all <- n_fail == 0L
stale_now <- setdiff(stale_tarballs(c(dest_dir, root), tarball_name), if (ok_all) handoff$path)
elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
record <- c(
  sprintf("# %s release gate, %s, commit %s", pkg_name, format(Sys.time(), "%Y-%m-%d %H:%M"), sha), "",
  sprintf("Result: %s (%d gates: %d passed, %d failed, %s; %.1f minutes)",
          if (ok_all) "ALL GATES PASSED" else "FAILED", length(state$results),
          sum(statuses() == "PASS"), n_fail, skip_summary(), elapsed_min), "",
  sprintf("Commit: %s", if (length(head_sha)) head_sha else "(none)"),
  sprintf("Tarball: %s (%s)", tarball_name, if (built) sprintf("%.1f KB", file.info(tarball)$size / 1024) else "not built"),
  sprintf("sha256: %s", sha256), sprintf("Check: %s (exit status %s)", parsed_log$status_line, chk$status),
  sprintf("Hand-off: %s", if (is.null(handoff)) "none (a gate failed, so nothing was copied)"
          else if (handoff$ok) sprintf("%s (sha256 read back: %s)", handoff$path, handoff$sha)
          else sprintf("FAILED to %s: %s", handoff$path, paste(handoff$problems, collapse = "; "))),
  sprintf("Tarballs of this version present at the start: %s%s",
          if (length(stale_at_start)) paste(stale_at_start, collapse = ", ") else "none",
          if (length(stale_at_start) && replace_stale) " (--replace-stale given)" else ""),
  sprintf("Other tarballs of this version present now: %s", if (length(stale_now)) paste(stale_now, collapse = ", ") else "none"), "",
  "## Gates", "",
  unlist(lapply(state$results, function(r) c(sprintf("- [%s] %s", r$status, r$name),
                                             if (length(r$detail)) paste0("    ", r$detail)))),
  "", "## Facts", "", paste0("- ", state$facts), "",
  "## 00check.log, verbatim", "",
  sprintf("Read from %s before the gate's temporary directory was removed.", check_log_path), "",
  if (length(chk$log00)) paste0("    ", chk$log00) else "    (no 00check.log was written)", "",
  "## R CMD check output (stdout and stderr, as captured)", "",
  if (length(chk$lines)) paste0("    ", chk$lines) else "    (no output)", "",
  "## The test suite against the tarball", "", paste0("    ", state$test_lines), "",
  "## The shipped harness, full mode, against the tarball", "", paste0("    ", state$harness_lines))
dir.create(record_dir, showWarnings = FALSE, recursive = TRUE)
rec_path <- file.path(record_dir, sprintf("%s_release_gate_%s_%s.md", pkg_name, stamp, sha))
writeLines(record, rec_path)
cat(sprintf("\nrecord: %s\n", rec_path))
if (ok_all) {
  cat(sprintf("tarball: %s\nsha256:  %s (read back from the hand-off path)\n", handoff$path, handoff$sha))
  cat(sprintf("\nALL GATES PASSED (%d gates, %s; %.1f min)\n", length(state$results), skip_summary(), elapsed_min))
  cat("next: (1) R CMD INSTALL the certified tarball, so the console and the submission are the same\n",
      "package; (2) run win-builder on R-devel and R-release with this tarball (HTML validation and the\n",
      "DESCRIPTION spell check run only there); (3) upload exactly the file above, after checking its\n",
      "sha256 against this record, and delete every other tarball of this version",
      if (length(stale_now)) sprintf(":\n%s\n", paste0("    ", stale_now, collapse = "\n")) else ".\n", sep = "")
  quit(status = 0L)
}
cat(sprintf("\nGATE FAILED: %d of %d gates failed (%s); %s\n", n_fail, length(state$results), skip_summary(),
            if (!is.null(handoff) && !handoff$ok) "the hand-off failed, so no certified tarball was delivered (see above)"
            else "no tarball handed off"))
quit(status = 1L)
