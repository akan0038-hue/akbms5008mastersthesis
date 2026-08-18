# =============================================================================
# lcms_traces.R
#
# Reads a Shimadzu LabSolutions Excel export and draws the extracted ion
# chromatograms, one panel per sample.
#
#   Set the two paths below, open in RStudio and click Source.
#
#   install.packages(c("tidyverse", "readxl", "patchwork"))
#
# Structure the parser expects, which is what LabSolutions writes:
#   [Header]                     application, version, data file
#   [MS Chromatogram]            one block per extracted ion trace
#     m/z            1-1MS(E+)m/z 828.0000
#     Interval(sec)  2
#     # of Points    960
#     Start Time(min), End Time(min)
#     R.Time(min) | Absolute Intensity | Relative Intensity
#     ... data ...
#   [MS Chromatogram]            further traces follow in the same sheet
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})

# ---- SET THESE ------------------------------------------------------------
BASE <- if (.Platform$OS.type == "windows") {
  "H:/path/to/your/Experimental"
} else {
  "/path/to/your/Experimental"
}
XLSX <- file.path(BASE, "lcmsstuff.xlsx")  # LabSolutions export
OUT  <- file.path(BASE, "lcms_output")        # folder for the figures
# ---------------------------------------------------------------------------

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(XLSX)) stop("Cannot find ", normalizePath(XLSX, mustWork = FALSE))

# ------------------------------------------------------------------ parsing --
# Each sheet holds one or more chromatogram blocks. Rather than assume fixed row
# offsets, the parser finds each block marker and then searches forward for the
# data header, so it still works if LabSolutions adds or removes a metadata line.
read_one_sheet <- function(path, sheet) {
  x <- suppressMessages(read_excel(path, sheet = sheet, col_names = FALSE,
                                   .name_repair = "minimal"))
  names(x) <- paste0("c", seq_along(x))
  key <- as.character(x$c1)
  starts <- which(key == "[MS Chromatogram]")
  if (!length(starts)) return(NULL)

  map_dfr(starts, function(i) {
    stop_at <- c(starts[starts > i], nrow(x) + 1)[1] - 1
    win     <- i:stop_at
    grab    <- function(label) {
      j <- win[which(key[win] == label)][1]
      if (is.na(j)) NA_character_ else as.character(x$c2[j])
    }
    mz_label <- grab("m/z")
    n_points <- suppressWarnings(as.integer(grab("# of Points")))
    hdr <- win[which(str_starts(key[win], "R.Time"))][1]
    if (is.na(hdr)) return(NULL)

    last <- if (is.na(n_points)) stop_at else min(hdr + n_points, stop_at)
    dat  <- x[(hdr + 1):last, 1:3]
    names(dat) <- c("rt", "intensity", "relative")
    dat %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      filter(!is.na(rt), !is.na(intensity)) %>%
      mutate(sample   = sheet,
             mz       = suppressWarnings(as.numeric(
                          str_extract(mz_label, "(?<=m/z )[0-9.]+"))),
             polarity = if_else(str_detect(mz_label, fixed("E+")), "+", "-"),
             trace    = paste0("m/z ", sprintf("%.1f", mz), " (", polarity, ")"))
  })
}

read_lcms <- function(path) {
  sheets <- excel_sheets(path)
  out <- map_dfr(sheets, ~ read_one_sheet(path, .x))
  out %>% mutate(sample = factor(sample, levels = sheets))
}

message("Reading ", XLSX)
lcms <- read_lcms(XLSX)
message("  ", nrow(lcms), " points across ",
        n_distinct(lcms$sample), " samples and ",
        n_distinct(paste(lcms$sample, lcms$mz)), " traces")
message("  m/z present in the export: ",
        paste(sort(unique(round(lcms$mz, 1))), collapse = ", "))

# ------------------------------------------------------------------- theme --
OKABE <- c("#0072B2", "#E69F00", "#009E73", "#D55E00",
           "#56B4E9", "#CC79A7", "#F0E442", "#000000")
theme_lcms <- function(base = 9.5) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(linewidth = 0.8, colour = "black"),
          axis.ticks = element_line(linewidth = 0.65, colour = "black"),
          axis.title = element_text(size = base + 0.5, face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top", legend.title = element_blank(),
          legend.key.size = unit(3.4, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(size = base, face = "bold", hjust = 0),
          plot.title = element_blank())
}
theme_set(theme_lcms())

# ------------------------------------------------------------------ plotting --
# normalise = "each"   every trace scaled to its own maximum, so weak traces stay
#                      visible. Use when comparing retention times.
# normalise = "sample" all traces in a sample scaled to the largest in that
#                      sample, preserving relative abundance within a run.
# normalise = "none"   raw counts. Only sensible when comparing runs acquired
#                      under identical conditions.
plot_traces <- function(d, normalise = c("each", "sample", "none"),
                        samples = NULL, mz_keep = NULL, rt_range = NULL,
                        ncol = 1) {
  normalise <- match.arg(normalise)
  if (!is.null(samples)) d <- filter(d, sample %in% samples)
  if (!is.null(mz_keep)) d <- filter(d, round(mz) %in% round(mz_keep))
  if (!is.null(rt_range)) d <- filter(d, rt >= rt_range[1], rt <= rt_range[2])

  d <- d %>%
    group_by(sample, mz) %>%
    mutate(y_each = 100 * intensity / max(intensity, na.rm = TRUE)) %>%
    group_by(sample) %>%
    mutate(y_sample = 100 * intensity / max(intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(y = switch(normalise, each = y_each, sample = y_sample,
                      none = intensity))

  ylab <- switch(normalise,
                 each   = "Relative Intensity (%)",
                 sample = "Relative Intensity (%)",
                 none   = "Absolute Intensity")

  ggplot(d, aes(rt, y, colour = trace)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ sample, ncol = ncol, scales = if (normalise == "none")
      "free_y" else "fixed") +
    scale_colour_manual(values = OKABE) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Retention Time (min)", y = ylab)
}

save_lcms <- function(p, name, w = 18, h = 24) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         units = "cm", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         units = "cm", dpi = 600)
  if (requireNamespace("devEMF", quietly = TRUE)) {
    devEMF::emf(file.path(OUT, paste0(name, ".emf")), width = w / 2.54,
                height = h / 2.54, coordDPI = 300, emfPlus = TRUE)
    print(p); grDevices::dev.off()
  }
  message("  wrote ", name)
}

# --------------------------------------------------------------- peak table --
# Local maxima above a threshold, for annotating retention times or for a
# supplementary table. Not a substitute for integration in LabSolutions.
peak_table <- function(d, min_rel = 20, min_gap = 0.2) {
  d %>%
    group_by(sample, trace, mz) %>%
    arrange(rt, .by_group = TRUE) %>%
    mutate(rel = 100 * intensity / max(intensity, na.rm = TRUE),
           is_max = intensity > lag(intensity, default = 0) &
                    intensity >= lead(intensity, default = 0)) %>%
    filter(is_max, rel >= min_rel) %>%
    mutate(keep = c(TRUE, diff(rt) > min_gap)) %>%
    filter(keep) %>%
    summarise(peaks = paste(sprintf("%.2f", rt), collapse = ", "),
              n_peaks = n(),
              top_rt = rt[which.max(intensity)],
              top_intensity = max(intensity), .groups = "drop")
}

# ------------------------------------------------------------------- output --
# Following the supervisor's direction, the thesis shows traces for the best
# system only. Everything else is summarised in Table 3.7 and supported by an
# appendix figure covering all seven syntheses.

SYNTH <- c("spps manual"          = "Manual",
           "LB normal coupling 1" = "Pris1",
           "LB normal coupling 2" = "Pris2",
           "LB normal coupling 3" = "Pris3",
           "LB new coupling 1"    = "Pris4",
           "LB new coupling 2"    = "Pris5",
           "LB new coupling 3"    = "Pris6")
COND  <- c(Manual = "manual",            Pris1 = "DIC/Oxyma, double 50 \u00b0C",
           Pris2  = "DIC/Oxyma, single 50 \u00b0C", Pris3 = "DIC/Oxyma, single 65 \u00b0C",
           Pris4  = "DIPEA/HCTU, single 50 \u00b0C", Pris5 = "DIPEA/HCTU, double 50 \u00b0C",
           Pris6  = "DIPEA/HCTU, single 65 \u00b0C")
TARGET_COL <- "#D55E00"; DELETION_COL <- "#0072B2"

# ---- ion definitions -------------------------------------------------------
# Corrected per Julien: one decimal to match the theoretical average masses.
#   target linear hydrazide   average 827.9   [M+H]+ 828.9   [M-H]- 826.9
#   des-Phe deletion          average 680.8   [M+H]+ 681.8   [M-H]- 679.8
# The earlier 828 / 682 / 680 were each a Dalton out. Matching is done on a
# tolerance rather than round(), because round(828.9) is 829 and the old
# round(mz) == 828 test silently drops every trace in the re-extracted export.
MZ_TARGET_POS   <- 828.9
MZ_TARGET_NEG   <- 826.9
MZ_DELETION_POS <- 681.8
MZ_DELETION_NEG <- 679.8
MZ_TOL          <- 0.7

LAB_TARGET   <- sprintf("m/z %.1f, target",  MZ_TARGET_POS)
LAB_DELETION <- sprintf("m/z %.1f, des-Phe", MZ_DELETION_POS)

syn <- lcms %>%
  filter(sample %in% names(SYNTH)) %>%
  mutate(synthesis = factor(unname(SYNTH[as.character(sample)]),
                            levels = unname(SYNTH)))

one_trace <- function(d, mz_keep, pol_keep, tol = MZ_TOL) {
  out <- d %>% filter(abs(mz - mz_keep) <= tol, polarity == pol_keep) %>% arrange(rt)
  if (!nrow(out)) {
    warning("No trace within ", tol, " of m/z ", mz_keep, " (", pol_keep,
            "). m/z present: ",
            paste(sort(unique(d$mz)), collapse = ", "), call. = FALSE)
  }
  # if the export holds both an old and a re-extracted trace near the same mass,
  # keep whichever is closest to the requested value
  if (n_distinct(out$mz) > 1) {
    keep <- out$mz[which.min(abs(out$mz - mz_keep))]
    out <- filter(out, mz == keep)
  }
  out
}

# ---- Figure 3.7: the best system only, target beside deletion ---------------
best <- filter(syn, synthesis == "Pris2")
f37 <- bind_rows(
  one_trace(best, MZ_TARGET_POS,   "+") %>% mutate(ion = LAB_TARGET),
  one_trace(best, MZ_DELETION_POS, "+") %>% mutate(ion = LAB_DELETION)) %>%
  mutate(ion = factor(ion, levels = c(LAB_TARGET, LAB_DELETION)))

p37 <- ggplot(f37, aes(rt, intensity, colour = ion)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ ion, nrow = 1) +
  scale_colour_manual(values = c(TARGET_COL, DELETION_COL), guide = "none") +
  scale_y_continuous(labels = function(v) ifelse(v == 0, "0",
                       paste0(v / 1e6, "M")),
                     expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Retention Time (min)", y = "Absolute Intensity")
save_lcms(p37, "Figure_3.7", 17, 6.5)

# ---- Figure A.3: all seven syntheses, target and deletion -------------------
# Both ions are now available for every synthesis, so the appendix figure shows
# the pair. Scales are free per panel because the routes differ by two orders of
# magnitude; a shared scale would flatten the DIPEA/HCTU rows to the baseline.
fa3 <- bind_rows(
  one_trace(syn, MZ_TARGET_POS,   "+") %>% mutate(ion = LAB_TARGET),
  one_trace(syn, MZ_DELETION_POS, "+") %>% mutate(ion = LAB_DELETION)) %>%
  mutate(ion = factor(ion, levels = c(LAB_TARGET, LAB_DELETION)),
         panel = factor(paste0(synthesis, "\n",
                               unname(COND[as.character(synthesis)])),
                        levels = paste0(levels(syn$synthesis), "\n",
                                        unname(COND))))
pa3 <- ggplot(fa3, aes(rt, intensity, colour = ion)) +
  geom_line(linewidth = 0.45) +
  facet_grid(panel ~ ion, scales = "free_y", switch = "y") +
  scale_colour_manual(values = c(TARGET_COL, DELETION_COL), guide = "none") +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Retention Time (min)", y = "Absolute Intensity") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5),
        strip.placement = "outside")
save_lcms(pa3, "Figure_A3_all_syntheses", 17, 24)

# ---- the numbers behind Table 3.7 ------------------------------------------
syn %>%
  group_by(synthesis, trace) %>%
  summarise(apex = max(intensity), rt_at_apex = rt[which.max(intensity)],
            .groups = "drop") %>%
  arrange(synthesis, trace) %>%
  write_csv(file.path(OUT, "Table_3.7_source.csv"))

message("Done. Figures and tables in ", normalizePath(OUT))

# =============================================================================
# After running, check these four things against the thesis
# -----------------------------------------------------------------------------
#  1. The startup message lists the m/z present. If it still shows 828, 682 and
#     680 rather than 828.9, 681.8 and 679.8, the workbook is the old export and
#     the traces need re-extracting in LabSolutions before this script is useful.
#
#  2. Table_3.7_source.csv now holds the apex intensities behind Table 3.8 and
#     Table F.1. Both tables must be rebuilt from it. The old Table F.1 values
#     came from the previous extraction and no longer apply.
#
#  3. Section 3.7 states that in the best system the target and the deletion
#     apexed within 0.06 min of one another. Recompute that from the new apexes:
#
#       syn %>% filter(synthesis == "Pris2") %>%
#         group_by(trace) %>% slice_max(intensity, n = 1) %>%
#         summarise(rt = rt) %>% pull(rt) %>% diff() %>% abs()
#
#  4. Section 3.7 says the manual synthesis gave the target only at trace level
#     and that the preparative fraction did likewise. Confirm both against the
#     new numbers rather than carrying the old wording across.
#
# Nothing anywhere in the thesis should still read 828, 826, 682 or 680.
# =============================================================================
