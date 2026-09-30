# Shared by the tests of the planned sample size search. A regression that
# brought back an unbounded loop would otherwise hang the suite; the limit
# turns it into a failure. setTimeLimit() persists for the rest of the
# top-level computation, so it is always cleared on exit.
with_time_limit <- function(code, seconds = 30) {
  setTimeLimit(elapsed = seconds, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  force(code)
}
