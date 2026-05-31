#' Defunct functions in BUCSS
#'
#' The dotted-name functions from BUCSS 1.x are defunct as of BUCSS 2.0.0. They
#' have been renamed (and, in the case of the between-subjects ANOVA, split into
#' separate one-way and two-way planners) under the \code{ss_buc_*} family, which
#' returns a tidy \code{\link{bucss_power}} object rather than a two-element
#' list. Calling a defunct function now raises an error that names its
#' replacement; update your scripts to the function listed below.
#'
#' \tabular{ll}{
#'   \strong{Defunct} \tab \strong{Replacement} \cr
#'   \code{ss.power.it} \tab \code{\link{ss_buc_indep_t}} (alias \code{\link{ss_buc_smd}}) \cr
#'   \code{ss.power.dt} \tab \code{\link{ss_buc_paired_t}} (alias \code{\link{ss_buc_smd_paired}}) \cr
#'   \code{ss.power.ba} \tab \code{\link{ss_buc_one_way_anova}}, \code{\link{ss_buc_factorial_anova}} \cr
#'   \code{ss.power.ba.general} \tab \code{\link{ss_buc_factorial_anova_general}} \cr
#'   \code{ss.power.wa} \tab \code{\link{ss_buc_rm_anova}} \cr
#'   \code{ss.power.wa.general} \tab \code{\link{ss_buc_rm_anova_general}} \cr
#'   \code{ss.power.spa} \tab \code{\link{ss_buc_mixed_anova}} \cr
#'   \code{ss.power.spa.general} \tab \code{\link{ss_buc_mixed_anova_general}} \cr
#'   \code{ss.power.reg1} \tab \code{\link{ss_buc_reg_coef}} \cr
#'   \code{ss.power.reg.all} \tab \code{\link{ss_buc_R2}} \cr
#'   \code{ss.power.reg.joint} \tab \code{\link{ss_buc_reg_joint}} \cr
#' }
#'
#' BUCSS 2.0.0 is not backward compatible with the 1.x API: function names,
#' argument names, and the returned object all changed. To run scripts written
#' for BUCSS 1.x without modification, install version 1.2.1 from the CRAN
#' archive at \url{https://cran.r-project.org/src/contrib/Archive/BUCSS/}.
#'
#' @param ... Arguments from the BUCSS 1.x calling convention. They are ignored:
#'   these functions are defunct and only raise an error naming their
#'   replacement.
#'
#' @return These functions do not return a value; each stops with an error.
#'
#' @name bucss-defunct
#' @keywords internal
NULL

.bucss_defunct <- function(old, replacement) {
  .Defunct(
    msg = paste0(
      "'", old, "' is defunct as of BUCSS 2.0.0. ", replacement,
      " BUCSS 2.0.0 renamed and restructured the function family, and the 1.x ",
      "dotted-name API is not backward compatible. To run scripts written for ",
      "BUCSS 1.x unchanged, install version 1.2.1 from the CRAN archive: ",
      "https://cran.r-project.org/src/contrib/Archive/BUCSS/."
    )
  )
}

#' @rdname bucss-defunct
#' @export
ss.power.it <- function(...) {
  .bucss_defunct("ss.power.it",
                 "Use 'ss_buc_indep_t' (or its effect-size alias 'ss_buc_smd') instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.dt <- function(...) {
  .bucss_defunct("ss.power.dt",
                 "Use 'ss_buc_paired_t' (or its effect-size alias 'ss_buc_smd_paired') instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.ba <- function(...) {
  .bucss_defunct("ss.power.ba",
                 paste("Use 'ss_buc_one_way_anova' for a single between-subjects factor",
                       "or 'ss_buc_factorial_anova' for a two-factor design instead."))
}

#' @rdname bucss-defunct
#' @export
ss.power.ba.general <- function(...) {
  .bucss_defunct("ss.power.ba.general",
                 "Use 'ss_buc_factorial_anova_general' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.wa <- function(...) {
  .bucss_defunct("ss.power.wa", "Use 'ss_buc_rm_anova' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.wa.general <- function(...) {
  .bucss_defunct("ss.power.wa.general", "Use 'ss_buc_rm_anova_general' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.spa <- function(...) {
  .bucss_defunct("ss.power.spa", "Use 'ss_buc_mixed_anova' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.spa.general <- function(...) {
  .bucss_defunct("ss.power.spa.general", "Use 'ss_buc_mixed_anova_general' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.reg1 <- function(...) {
  .bucss_defunct("ss.power.reg1", "Use 'ss_buc_reg_coef' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.reg.all <- function(...) {
  .bucss_defunct("ss.power.reg.all", "Use 'ss_buc_R2' instead.")
}

#' @rdname bucss-defunct
#' @export
ss.power.reg.joint <- function(...) {
  .bucss_defunct("ss.power.reg.joint", "Use 'ss_buc_reg_joint' instead.")
}
