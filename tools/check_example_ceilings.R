# Developer check (not shipped): every planner's zero-NCP @examples block
# claims a ceiling ("here near .86"). R CMD check cannot verify it, because the
# try() swallows the error and a MISSING error is silent, so an example can
# claim a ceiling that the message never reports, or fail to error at all and
# illustrate nothing. Both happened during development. This script evaluates
# each try() example and compares the stated value against the ceiling the live
# message reports.
#
# Usage:
#   Rscript tools/check_example_ceilings.R
#
# Run it after editing any planner's example or any prior-study value in one.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
files <- list.files("R", pattern = "^ss_buc_.*[.]R$", full.names = TRUE)
bad <- 0; checked <- 0
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  i <- grep("here near ", txt, fixed = TRUE)
  if (!length(i)) next
  stated <- regmatches(txt[i[1]], regexpr("[.][0-9]+", txt[i[1]]))
  j <- grep("try(", txt, fixed = TRUE); j <- j[j > i[1]][1]
  if (is.na(j)) next
  code <- character(0); k <- j
  repeat {
    code <- c(code, sub("^#' ", "", txt[k]))
    if (grepl("))", txt[k], fixed = TRUE)) break
    k <- k + 1; if (k > j + 6) break
  }
  expr <- paste(code, collapse = " ")
  inner <- sub("^try\\(", "", sub("\\)\\s*$", "", expr))
  msg <- tryCatch({ eval(parse(text = inner)); "NO ERROR" },
                  error = function(e) conditionMessage(e))
  live <- if (grepl("about", msg, fixed = TRUE)) {
    regmatches(msg, regexpr("[.][0-9]+", msg))
  } else if (identical(msg, "NO ERROR")) "NO ERROR" else "(mode 2, no ceiling)"
  checked <- checked + 1
  if (!identical(stated, live)) {
    bad <- bad + 1
    cat(sprintf("MISMATCH %-34s doc=%s live=%s\n", basename(f), stated, live))
  }
}
cat("ceiling comments checked:", checked, "| mismatches:", bad, "\n")
