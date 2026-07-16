# Canonical generator for the BUCSS hex logo (man/figures/logo.svg + PNGs).
# The design is LOCKED in the family master (KenKelleyHexFamily): sodalite
# ground (#283a5e) + copper border (#b87333). A replication forest of thin
# confidence intervals straddles the dashed null; the one significant CI
# (copper, right of 0, italic p) is tied to its small n_Published (a copper dot
# grid), and BUCSS plans the much larger n_Planned (an aqua #46c2bf dot grid --
# the family's one accent, used only here). No arrows. Art code lives in
# build_BUCSS_inner() in the family builder; this script renders it into the
# package.
#   Rscript dev/make_BUCSS_hex.R      (from the BUCSS package root)
args0 <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", args0[grepl("^--file=", args0)])
HERE <- if (length(fa)) normalizePath(dirname(fa)) else getwd()   # BUCSS/dev
PKG <- normalizePath(file.path(HERE, ".."))                       # BUCSS
FAM <- normalizePath(file.path(PKG, "..", "KenKelleyHexFamily"))  # the family
source(file.path(FAM, "_render.R"))
source(file.path(FAM, "make_family_hexes.R"))   # provides build_BUCSS_inner(), svg_of(), KK_HEX

fig <- file.path(PKG, "man", "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
svg <- svg_of(build_BUCSS_inner(KK_HEX[["BUCSS"]]))   # viewBox only, matches DMAR's logo.svg
writeLines(svg, file.path(fig, "logo.svg"))
# Chrome sizes the render window to the SVG's width/height, so rasterize from a
# dimensioned copy (the shipped logo.svg stays viewBox-only like DMAR's).
dimsvg <- sub('<svg ', '<svg width="200" height="210" ', svg, fixed = TRUE)
tmp <- file.path(fig, ".logo_render.svg"); writeLines(dimsvg, tmp)
render_svg(tmp, file.path(fig, "logo.png"), 200, 210, 6)          # 1200 x 1260, like DMAR
render_svg(tmp, file.path(fig, "logo-hires.png"), 200, 210, 20)   # 4000 x 4200, like DMAR
unlink(tmp)
cat("wrote", file.path(fig, "logo.svg"), "+ logo.png + logo-hires.png\n")
