# =============================================================================
# te_analysis.R
#
# Every figure and every statistic for the AlphaFold3 results chapter.
#
# TWO WAYS TO RUN IT
#
#   Set the four file paths below, then in RStudio open this file and click
#   Source. Nothing else is needed: no setwd(), no read.csv().
#
#   "Rscript te_analysis.R" is a terminal command, not an R command. Typing it
#   at the R prompt gives "Error: unexpected symbol".
#
# Uses the four te_autoscan long files:
#   ramocasenucleophile_distances_long.csv      (Sections 3.1, 3.2)
#   ramonucleophile_distances_long.csv          (Sections 3.1, 3.3)
#   prisnucleophile_distances_long.csv          (Sections 3.1, 3.4, 3.6)
#   prisanalognucleophile_distances_long.csv    (Sections 3.1, 3.5)
# They may live in different folders.
#
# Install once:
#   install.packages(c("tidyverse", "rstatix", "patchwork", "pROC", "scales"))
#   install.packages("devEMF")   # optional, writes vector figures for Word
#
# Tested against tidyverse 2.0, rstatix 0.7, pROC 1.18.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(patchwork)
  library(pROC)
  library(scales)
})

# ---- SET YOUR FILE PATHS HERE ---------------------------------------------
# Forward slashes throughout; correct on macOS and on Windows.
# You do not need to read these files yourself; the script does it, with proper
# column headers. If you already ran read.csv(..., header = FALSE), clear the
# environment first (broom icon in the Environment pane) so nothing clashes.

# Everything below hangs off BASE, so moving the project is a one-line change.
BASE <- "/path/to/your/Experimental"

F_RAMOCASE  <- file.path(BASE, "Bioinformatics/ramocase/ramocasenucleophile_distances_long.csv")
F_RAMOSHORT <- file.path(BASE, "Bioinformatics/ramoshort/ramonucleophile_distances_long.csv")
F_PRIS      <- file.path(BASE, "Bioinformatics/prisnucleophile_distances_long.csv")
F_ANALOG    <- file.path(BASE, "Bioinformatics/prisanalogues/prisanalognucleophile_distances_long.csv")

OUT <- file.path(BASE, "Bioinformatics/R_output")
# ---------------------------------------------------------------------------

paths <- c(ramocase = F_RAMOCASE, ramoshort = F_RAMOSHORT,
           pris = F_PRIS, analogues = F_ANALOG)
missing_files <- paths[!file.exists(paths)]
if (length(missing_files))
  stop("Cannot find:\n  ", paste(missing_files, collapse = "\n  "),
       "\nCheck the paths at the top of the script.")
# Create the output folders. Google Drive paths sometimes refuse a deep
# recursive create, so each level is made in turn and then checked.
for (d in c(OUT, file.path(OUT, "figures"), file.path(OUT, "stats")))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
bad <- c(OUT, file.path(OUT, "figures"), file.path(OUT, "stats"))
bad <- bad[!dir.exists(bad)]
if (length(bad))
  stop("Could not create:\n  ", paste(bad, collapse = "\n  "),
       "\nIf OUT is on Google Drive, try a local folder instead, for example\n",
       "  OUT <- path.expand(\"~/te_results\")")
message("Reading four files\nWriting to ", normalizePath(OUT))

# ---------------------------------------------------------------- constants --
# A near-attack conformation requires all three at once. Each criterion alone
# passes large numbers of unproductive poses; only the conjunction discriminates.
NE2_CUT   <- 3.5    # nucleophile to the catalytic nitrogen, hydrogen bonding
ATK_CUT   <- 4.0    # nucleophile to the electrophilic carbon, bond forming
ANG_LO    <- 95     # Burgi-Dunitz window, lower
ANG_HI    <- 115    # Burgi-Dunitz window, upper

# The ring-closing nucleophile sits at residue 3 in the ramoplanin-family
# peptides, because the lipid cap is counted as residue 1, and at residue 2 in
# the streptogramins, which have no cap.
FILES <- tribble(
  ~key,        ~ring_residue, ~label,
  "ramocase",   3, "Ramoplanin Family",
  "ramoshort",  3, "Ramoplanin (Truncated)",
  "pris",       2, "Pristinamycin (Native Sequence)",
  "analogues",  2, "Pristinamycin Variants (Position 2)"
) %>% mutate(file = unname(paths[key]))

# Wong (2011) palette: eight hues distinguishable under all common forms of
# colour vision deficiency, and in greyscale.
OKABE <- c(orange   = "#E69F00", skyblue = "#56B4E9", green  = "#009E73",
           yellow   = "#F0E442", blue    = "#0072B2", vermillion = "#D55E00",
           purple   = "#CC79A7", black   = "#000000")

PAL_DATASET <- unname(OKABE[c("blue", "skyblue", "vermillion", "orange")])
PAL_PAIRING <- c(Native = unname(OKABE["blue"]), `Cross-Domain` = unname(OKABE["skyblue"]))
PAL_CLASS   <- c(Pristinamycin = unname(OKABE["blue"]),
                 Virginiamycin = unname(OKABE["vermillion"]))
PAL_CHEM    <- c(Hydroxyl = unname(OKABE["blue"]), Amine = unname(OKABE["orange"]))
PAL_SCOPE   <- c(Pooled = unname(OKABE["blue"]),
                 `Pristinamycin 1` = unname(OKABE["skyblue"]),
                 `Pristinamycin 2` = unname(OKABE["green"]))
BAND_NEAR <- "#B7E4C7"; BAND_INT <- "#FFE8A3"

LW  <- 0.8    # data line weight
AXW <- 0.8    # axis line weight

theme_thesis <- function(base = 9.5) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.line       = element_line(linewidth = AXW, colour = "black",
                                     lineend = "round"),
      axis.ticks      = element_line(linewidth = AXW * 0.8, colour = "black"),
      axis.ticks.length = unit(1.1, "mm"),
      # Explicit margins keep the titles close to their axes; the defaults leave
      # a visible gap once the tick length and a two-line title are added.
      axis.title.x    = element_text(size = base + 0.5, colour = "black",
                                     face = "bold", margin = margin(t = 3)),
      axis.title.y    = element_text(size = base + 0.5, colour = "black",
                                     face = "bold", margin = margin(r = 3)),
      axis.text.x     = element_text(size = base - 0.5, colour = "black",
                                     margin = margin(t = 1.5)),
      axis.text.y     = element_text(size = base - 0.5, colour = "black",
                                     margin = margin(r = 1.5)),
      plot.title      = element_blank(),
      # "topleft" puts the tag in the margin above the panel. A fixed npc
      # position places it over the rotated y-axis title instead.
      plot.tag        = element_text(size = base + 3.5, face = "bold",
                                     hjust = 0, vjust = 1),
      plot.tag.position = "topleft",
      legend.position = "top",
      legend.title    = element_blank(),
      legend.text     = element_text(size = base - 0.5, face = "bold"),
      legend.key.size = unit(3.4, "mm"),
      legend.margin   = margin(0, 0, 1, 0),
      strip.background = element_blank(),
      strip.text      = element_text(size = base, face = "bold"),
      plot.margin     = margin(t = 11, r = 6, b = 4, l = 5)
    )
}
theme_set(theme_thesis())

# ------------------------------------------------------------------- loading --
read_te <- function(path) {
  d <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  # AlphaFold3 writes the top pose of each seed twice, once in its own seed
  # folder and once under collected_models/. Keep one copy.
  if ("collected_copy" %in% names(d)) d <- filter(d, !collected_copy)
  d
}

# Restrict to the nucleophile that closes the macrocycle. Reporting whichever
# nucleophile happens to be nearest measures a different atom in different
# pairings, and gives a substrate with many candidates more chances than one
# with few.
ring_closing <- function(d, residue) {
  out <- filter(d, nuc_residue == residue)
  if (nrow(out) == 0)
    stop("no nucleophile at residue ", residue,
         "; check numbering against the nuc_residue column")
  out
}

add_nac <- function(d) {
  mutate(d, nac = d_ne2 < NE2_CUT & d_attack < ATK_CUT &
           dplyr::between(burgi_dunitz_angle, ANG_LO, ANG_HI))
}

DAT <- FILES %>%
  mutate(raw  = map(file, read_te),
         ring = map2(raw, ring_residue, ring_closing),
         ring = map(ring, add_nac))
get_raw  <- function(k) DAT$raw[[which(DAT$key == k)]]
get_ring <- function(k) DAT$ring[[which(DAT$key == k)]]

# --------------------------------------------------------------- stats utils --
# Rank-biserial correlation, the effect size for a Mann-Whitney test. It is the
# probability that a value drawn from x exceeds one drawn from y, rescaled to
# -1..1. A value of +/-1 means the two distributions do not overlap at all,
# which no p value and no difference in medians conveys.
rank_biserial <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  2 * unname(w$statistic) / (length(x) * length(y)) - 1
}

mw <- function(d, value, group, g1, g2, label) {
  x <- d[[value]][d[[group]] %in% g1]; x <- x[!is.na(x)]
  y <- d[[value]][d[[group]] %in% g2]; y <- y[!is.na(y)]
  # A renamed factor level would otherwise abort the whole run inside
  # wilcox.test. Report the empty group and carry on.
  if (length(x) < 3 || length(y) < 3) {
    warning("Skipped '", label, "': found ", length(x), " and ", length(y),
            " values. Levels present in '", group, "' are: ",
            paste(unique(as.character(d[[group]])), collapse = ", "),
            call. = FALSE)
    return(tibble(test = "Mann-Whitney", comparison = label,
                  statistic = NA_real_, p = NA_real_, effect = NA_real_,
                  effect_type = "not run, empty group",
                  n1 = length(x), n2 = length(y)))
  }
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  tibble(test = "Mann-Whitney", comparison = label,
         statistic = unname(w$statistic), p = w$p.value,
         effect = rank_biserial(x, y), effect_type = "rank-biserial",
         n1 = length(x), n2 = length(y))
}

# The near-attack matrices are mostly zero, so a pale-to-dark ramp leaves the
# informative cells barely distinguishable from the empty ones. Zero is given a
# neutral grey and every non-zero value a saturated viridis hue, with the label
# colour switching so it stays legible on both.
heat_layers <- function(limit = 70) {
  list(geom_tile(colour = "white", linewidth = 1.1),
       geom_text(aes(label = ifelse(pct == 0, "0", sprintf("%.0f", pct)),
                     colour = pct > limit * 0.45),
                 size = 2.8, fontface = "bold", show.legend = FALSE),
       scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")),
       scale_fill_gradientn(
         colours = c("#EDEDED", "#B7E4C7", "#56B4E9", "#0072B2", "#08306B"),
         values  = scales::rescale(c(0, 0.001, 0.25, 0.55, 1)),
         limits  = c(0, limit), name = "Structures Reaching\nNear-Attack (%)",
         guide = guide_colourbar(barwidth = unit(28, "mm"),
                                 barheight = unit(2.6, "mm"),
                                 title.position = "top", title.hjust = 0)),
       theme(axis.line = element_blank(), axis.ticks = element_blank(),
             panel.grid = element_blank()))
}

band_layers <- function(ne2_only = FALSE) {
  l <- list(annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = NE2_CUT,
                     fill = BAND_NEAR, alpha = 0.55))
  if (!ne2_only)
    l <- c(l, list(annotate("rect", xmin = -Inf, xmax = Inf,
                            ymin = NE2_CUT, ymax = 6, fill = BAND_INT, alpha = 0.55)))
  l
}

# stat_boxplot(geom = "errorbar") draws a terminal cap on each whisker. Without
# it ggplot renders the whisker as a bare line, which reads as unfinished next to
# the boxes. Drawn first so the caps sit behind the box.
box_layers <- function(w = 0.62) {
  list(stat_boxplot(geom = "errorbar", width = w * 0.45,
                    linewidth = LW * 0.65, colour = "black",
                    show.legend = FALSE),
       geom_boxplot(outlier.size = 0.45, outlier.alpha = 0.30,
                    outlier.shape = 21, outlier.stroke = 0.25,
                    linewidth = LW * 0.75, width = w, alpha = 0.9,
                    colour = "black", median.linewidth = LW * 1.4),
       scale_y_continuous(expand = expansion(mult = c(0, 0.06))))
}

# Word rasterises an inserted PNG and rescales it, which is why the figures look
# soft on the page. EMF is a vector format Word renders natively at any size, so
# it is written as well when devEMF is available. Insert the .emf into Word and
# keep the .pdf for anything else.
save_fig <- function(plot, name, w, h) {
  ggsave(file.path(OUT, "figures", paste0(name, ".pdf")), plot,
         width = w, height = h, units = "cm", device = cairo_pdf)
  ggsave(file.path(OUT, "figures", paste0(name, ".png")), plot,
         width = w, height = h, units = "cm", dpi = 600)
  if (requireNamespace("devEMF", quietly = TRUE)) {
    devEMF::emf(file.path(OUT, "figures", paste0(name, ".emf")),
                width = w / 2.54, height = h / 2.54, coordDPI = 300,
                emfPlus = TRUE)
    print(plot); grDevices::dev.off()
  }
  message("  wrote figures/", name)
}

STATS <- list()
keep <- function(section, tbl) STATS[[length(STATS) + 1]] <<-
  mutate(tbl, section = section, .before = 1)

# =============================================================================
# 3.1  Confidence is uniformly high
# =============================================================================
message("Section 3.1")

conf_global <- DAT %>%
  transmute(label, d = map(raw, ~ distinct(.x, model_path, .keep_all = TRUE))) %>%
  unnest(d) %>%
  mutate(label = factor(label, levels = FILES$label))

conf_local <- DAT %>%
  transmute(label, d = map(ring, ~ select(.x, nuc_plddt))) %>%
  unnest(d) %>%
  mutate(label = factor(label, levels = FILES$label))

# Category names here are long and the longest single word ("Pristinamycin") is
# already wider than the column, so wrapping cannot prevent collisions. Angling
# them is the workable option; the canvas heights below allow for it.
XLAB <- function(angle = 30)
  theme(axis.text.x = element_text(angle = angle, hjust = 1))

p31a <- ggplot(conf_global, aes(label, iptm, fill = label)) +
  geom_hline(yintercept = 0.7, linetype = "22", linewidth = LW * 0.7,
             colour = "grey30") +
  box_layers() +
  scale_fill_manual(values = PAL_DATASET, guide = "none") +
  labs(x = NULL, y = "Interface Confidence (ipTM)") +
  XLAB(30)

p31b <- ggplot(conf_local, aes(label, nuc_plddt, fill = label)) +
  box_layers() +
  scale_fill_manual(values = PAL_DATASET, guide = "none") +
  labs(x = NULL, y = "Local Confidence at the\nRing-Closing Nucleophile (pLDDT)") +
  XLAB(30)

save_fig(p31a + p31b + plot_annotation(tag_levels = "A"), "Figure_3.1", 18, 10.5)

keep("3.1", conf_global %>% kruskal_test(iptm ~ label) %>%
       transmute(test = "Kruskal-Wallis", comparison = "ipTM across datasets",
                 statistic, p, effect = NA_real_, effect_type = NA_character_,
                 n1 = n, n2 = NA_integer_))
keep("3.1", mw(conf_local, "nuc_plddt", "label",
               "Ramoplanin Family", "Pristinamycin (Native Sequence)",
               "nucleophile pLDDT, ramoplanin family vs pristinamycin"))

conf_summary <- conf_global %>%
  group_by(label) %>%
  summarise(n = n(), ipTM = mean(iptm), pTM = mean(ptm),
            ranking = mean(ranking_score), .groups = "drop") %>%
  left_join(
  conf_local %>% group_by(label) %>%
    summarise(nucleophile_pLDDT = mean(nuc_plddt, na.rm = TRUE), .groups = "drop"),
  by = "label")
write_csv(conf_summary, file.path(OUT, "stats", "Table_3.1.csv"))

# =============================================================================
# 3.2  Ramoplanin family, four substrates on three domains
# =============================================================================
message("Section 3.2")

ENZYME <- c(`218` = "chersinamycin", `220` = "ramoplanin", `225` = "enduracidin")

ramo <- get_ring("ramocase") %>%
  mutate(
    enzyme    = unname(ENZYME[as.character(his_seqid)]),
    substrate = case_when(
      substrate_resname == "chrs" &  str_detect(job, "trial") ~ "Chersinamycin (Natural)",
      substrate_resname == "chrs" & !str_detect(job, "trial") ~ "Synthetic Construct",
      substrate_resname == "endu" ~ "Enduracidin (Natural)",
      substrate_resname == "ramo" ~ "Ramoplanin (Natural)"),
    type    = if_else(substrate == "Synthetic Construct", "synthetic", "natural"),
    cognate = case_when(str_starts(substrate, "Synthetic")     ~ "chersinamycin",
                        str_starts(substrate, "Chersinamycin") ~ "chersinamycin",
                        str_starts(substrate, "Enduracidin")   ~ "enduracidin",
                        TRUE                                   ~ "ramoplanin"),
    pairing = factor(if_else(enzyme == cognate, "Native", "Cross-Domain"),
                     levels = c("Native", "Cross-Domain")),
    substrate = factor(substrate, levels = c(
      "Synthetic Construct", "Chersinamycin (Natural)",
      "Enduracidin (Natural)", "Ramoplanin (Natural)")),
    enzyme = factor(str_to_title(enzyme),
                    levels = c("Chersinamycin", "Enduracidin", "Ramoplanin")))

p32c <- ggplot(ramo, aes(enzyme, d_ne2, fill = pairing)) +
  band_layers(TRUE) + box_layers() +
  facet_wrap(~ substrate, nrow = 1) +
  scale_fill_manual(values = PAL_PAIRING) +
  labs(x = "TE-Domain", y = "Ring-Closing Nucleophile\nto His NE2 (\u00c5)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(size = 7.5))

p32d <- ramo %>%
  group_by(substrate, enzyme) %>%
  summarise(pct = 100 * mean(nac), .groups = "drop") %>%
  mutate(short = fct_rev(factor(str_remove(substrate, " \\(Natural\\)") %>%
                                  str_replace("Synthetic Construct", "Synthetic"),
                                levels = c("Synthetic", "Chersinamycin",
                                           "Enduracidin", "Ramoplanin")))) %>%
  ggplot(aes(enzyme, short, fill = pct)) +
  heat_layers() +
  labs(x = "TE-Domain", y = NULL) +
  coord_fixed(ratio = 0.55) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# patchwork aligns the axis areas of stacked plots, so panel D's long row labels
# push panel C's y-axis title away from its axis. free() releases panel C from
# that alignment; without it the gap cannot be closed by margins alone.
p32c_free <- if ("free" %in% getNamespaceExports("patchwork")) {
  patchwork::free(p32c)
} else {
  message("  patchwork >= 1.2 gives a tighter Figure 3.2 axis; consider upgrading")
  p32c
}
save_fig(p32c_free / p32d + plot_layout(heights = c(1.7, 1)) +
           plot_annotation(tag_levels = list(c("C", "D"))),
         "Figure_3.2", 18, 20)

keep("3.2", ramo %>% mutate(g = paste(substrate, enzyme)) %>%
       kruskal_test(d_ne2 ~ g) %>%
       transmute(test = "Kruskal-Wallis", comparison = "12 pairings",
                 statistic, p, effect = NA_real_, effect_type = NA_character_,
                 n1 = n, n2 = NA_integer_))
keep("3.2", mw(ramo, "d_ne2", "type", "synthetic", "natural",
               "synthetic construct vs three natural peptides"))
keep("3.2", mw(ramo, "d_ne2", "pairing", "Native", "Cross-Domain",
               "native vs cross-domain pairings"))

ramo %>% group_by(substrate, enzyme, pairing) %>%
  summarise(n = n(), median_dNE2 = median(d_ne2), min_dNE2 = min(d_ne2),
            near_attack = sum(nac), .groups = "drop") %>%
  write_csv(file.path(OUT, "stats", "Table_3.2.csv"))

# =============================================================================
# 3.3  Truncation series
# =============================================================================
message("Section 3.3")

LEN <- c(synthpep1 = 18, synthpep2 = 16, synthpep3 = 14,
         synthpep4 = 12, synthpep5 = 10, synthpep6 = 7)

trunc <- get_ring("ramoshort") %>%
  mutate(peptide = factor(substrate_resname, levels = names(LEN)),
         residues = unname(LEN[as.character(substrate_resname)]),
         enzyme  = sub(".*_on_", "", job))

p33c <- ggplot(trunc, aes(peptide, d_ne2)) +
  band_layers(TRUE) +
  stat_boxplot(geom = "errorbar", width = 0.28, linewidth = LW * 0.65,
               colour = "black") +
  geom_boxplot(fill = unname(OKABE["skyblue"]), outlier.size = 0.45,
               outlier.alpha = 0.3, outlier.shape = 21, outlier.stroke = 0.25,
               linewidth = LW * 0.75, width = 0.62, alpha = 0.9,
               colour = "black", median.linewidth = LW * 1.4) +
  stat_summary(fun = median, geom = "line", group = 1,
               colour = unname(OKABE["vermillion"]), linewidth = LW * 1.4) +
  stat_summary(fun = median, geom = "point", group = 1,
               colour = unname(OKABE["vermillion"]), size = 1.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  scale_x_discrete(labels = as.character(LEN)) +
  labs(x = "Peptide Length (Residues)",
       y = "Ring-Closing Nucleophile\nto His NE2 (\u00c5)")

p33d <- trunc %>% filter(d_ne2 < NE2_CUT) %>%
  ggplot(aes(d_attack, burgi_dunitz_angle)) +
  annotate("rect", xmin = 3, xmax = ATK_CUT, ymin = -Inf, ymax = Inf,
           fill = BAND_INT, alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ANG_LO, ymax = ANG_HI,
           fill = BAND_NEAR, alpha = 0.6) +
  geom_hline(yintercept = 105, linetype = "22", linewidth = LW * 0.7,
             colour = "grey30") +
  geom_point(colour = unname(OKABE["vermillion"]), size = 1.5, alpha = 0.75,
             shape = 21, fill = unname(OKABE["vermillion"]), stroke = 0.25) +
  labs(x = "Nucleophile to Electrophilic Carbon (\u00c5)",
       y = "B\u00fcrgi-Dunitz Angle (\u00b0)")

save_fig(p33c + p33d + plot_annotation(tag_levels = list(c("C", "D"))),
         "Figure_3.3", 18, 9.5)

keep("3.3", trunc %>% kruskal_test(d_ne2 ~ peptide) %>%
       transmute(test = "Kruskal-Wallis", comparison = "six peptides",
                 statistic, p, effect = NA_real_, effect_type = NA_character_,
                 n1 = n, n2 = NA_integer_))
sp <- cor.test(trunc$residues, trunc$d_ne2, method = "spearman", exact = FALSE)
keep("3.3", tibble(test = "Spearman", comparison = "residue count vs d_NE2, pooled",
                   statistic = unname(sp$estimate), p = sp$p.value,
                   effect = unname(sp$estimate), effect_type = "rho",
                   n1 = nrow(trunc), n2 = NA_integer_))
# structures from one seed are not independent; repeat on one value per seed
per_seed <- trunc %>% group_by(job, seed) %>%
  summarise(d = min(d_ne2), residues = first(residues), .groups = "drop")
sp2 <- cor.test(per_seed$residues, per_seed$d, method = "spearman", exact = FALSE)
keep("3.3", tibble(test = "Spearman", comparison = "residue count vs d_NE2, per seed",
                   statistic = unname(sp2$estimate), p = sp2$p.value,
                   effect = unname(sp2$estimate), effect_type = "rho",
                   n1 = nrow(per_seed), n2 = NA_integer_))
keep("3.3", trunc %>% dunn_test(d_ne2 ~ peptide, p.adjust.method = "holm") %>%
       transmute(test = "Dunn", comparison = paste(group1, "vs", group2),
                 statistic, p = p.adj, effect = NA_real_,
                 effect_type = NA_character_, n1, n2))

trunc %>% group_by(peptide, residues) %>%
  summarise(n = n(), median_dNE2 = median(d_ne2), min_dNE2 = min(d_ne2),
            near_attack = sum(nac), .groups = "drop") %>%
  write_csv(file.path(OUT, "stats", "Table_3.3.csv"))

# =============================================================================
# 3.4  Streptogramin native substrate
# =============================================================================
message("Section 3.4")

DOM     <- c("pris1", "pris2", "vir1", "vir2", "vir3", "vir4")
DOM_LAB <- c("Pristinamycin 1", "Pristinamycin 2", "Virginiamycin 1",
             "Virginiamycin 2", "Virginiamycin 3", "Virginiamycin 4")
relabel_domain <- function(x) factor(DOM_LAB[match(x, DOM)], levels = DOM_LAB)

pris <- get_ring("pris") %>%
  mutate(domain = relabel_domain(job),
         class  = factor(if_else(str_starts(job, "pris"),
                                 "Pristinamycin", "Virginiamycin"),
                         levels = c("Pristinamycin", "Virginiamycin")))
pris_all <- get_raw("pris") %>% mutate(domain = relabel_domain(job))

p34c <- ggplot(pris, aes(domain, d_ne2, fill = class)) +
  band_layers(TRUE) + box_layers() +
  scale_fill_manual(values = PAL_CLASS) +
  labs(x = "TE-Domain", y = "Ring-Closing Nucleophile\nto His NE2 (\u00c5)") +
  XLAB(30)

p34d <- pris_all %>%
  mutate(nac = d_ne2 < NE2_CUT & d_attack < ATK_CUT &
           dplyr::between(burgi_dunitz_angle, ANG_LO, ANG_HI)) %>%
  mutate(nuc_label = recode(nuc_atom,
           O3 = "O3\nring-closing\nhydroxyl",
           O1 = "O1\nphenolic\n(3HPA cap)",
           N6 = "N6\nmethylamino\n(position 5)",
           .default = nuc_atom)) %>%
  group_by(domain, nuc_label) %>%
  summarise(pct = 100 * mean(nac), .groups = "drop") %>%
  ggplot(aes(domain, nuc_label, fill = pct)) +
  heat_layers() +
  labs(x = "TE-Domain", y = NULL) +
  theme(axis.text.y = element_text(lineheight = 0.85)) +
  XLAB(30)

# The two pristinamycin domains hold the nucleophile at the same distance and
# differ only in approach angle, which is what separates 97 near-attack
# structures from 50. Without this panel that result is prose only.
p34e <- pris %>%
  filter(str_starts(job, "pris")) %>%
  ggplot(aes(burgi_dunitz_angle, fill = domain)) +
  annotate("rect", xmin = ANG_LO, xmax = ANG_HI, ymin = -Inf, ymax = Inf,
           fill = BAND_NEAR, alpha = 0.55) +
  geom_histogram(binwidth = 5, colour = "black", linewidth = LW * 0.5,
                 alpha = 0.85, position = "identity") +
  facet_wrap(~ domain, ncol = 1) +
  scale_fill_manual(values = unname(OKABE[c("blue", "green")]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "B\u00fcrgi-Dunitz Angle (\u00b0)", y = "Structures")

save_fig(p34c + p34d + p34e +
           plot_layout(widths = c(1, 1, 0.85)) +
           plot_annotation(tag_levels = list(c("C", "D", "E"))),
         "Figure_3.4", 24, 10)

keep("3.4", pris %>% kruskal_test(d_ne2 ~ domain) %>%
       transmute(test = "Kruskal-Wallis", comparison = "six domains",
                 statistic, p, effect = NA_real_, effect_type = NA_character_,
                 n1 = n, n2 = NA_integer_))
keep("3.4", mw(pris, "d_ne2", "class", "Pristinamycin", "Virginiamycin",
               "pristinamycin vs virginiamycin domains"))
keep("3.4", pris %>% dunn_test(d_ne2 ~ domain, p.adjust.method = "holm") %>%
       transmute(test = "Dunn", comparison = paste(group1, "vs", group2),
                 statistic, p = p.adj, effect = NA_real_,
                 effect_type = NA_character_, n1, n2))
# Fisher rather than chi-square: four domains have no near-attack structures at
# all, so expected counts fall below five and chi-square is invalid.
ft <- pris %>% count(class, nac) %>%
  pivot_wider(names_from = nac, values_from = n, values_fill = 0) %>%
  column_to_rownames("class") %>% as.matrix() %>% fisher.test()
keep("3.4", tibble(test = "Fisher exact",
                   comparison = "near-attack, pristinamycin vs virginiamycin",
                   statistic = NA_real_, p = ft$p.value, effect = NA_real_,
                   effect_type = NA_character_, n1 = 300L, n2 = 600L))

pris %>% group_by(domain, class) %>%
  summarise(n = n(), his = first(his_seqid), median_dNE2 = median(d_ne2),
            min_dNE2 = min(d_ne2), near_attack = sum(nac), .groups = "drop") %>%
  write_csv(file.path(OUT, "stats", "Table_3.4.csv"))

# =============================================================================
# 3.5  Position-2 variants
# =============================================================================
message("Section 3.5")

VAR <- c(lthr = "L-Thr", ser = "L-Ser", dthr = "D-Thr", dab = "Dab", dap = "Dap")
ana <- get_ring("analogues") %>%
  mutate(variant = factor(unname(VAR[substrate_dir]), levels = unname(VAR)),
         chem = factor(if_else(substrate_dir %in% c("dab", "dap"),
                               "Amine", "Hydroxyl"),
                       levels = c("Hydroxyl", "Amine")))

p35c <- ggplot(ana, aes(variant, d_ne2, fill = chem)) +
  band_layers(TRUE) + box_layers() +
  scale_fill_manual(values = PAL_CHEM) +
  labs(x = "Residue at Position 2",
       y = "Position-2 Nucleophile\nto His NE2 (\u00c5)")

p35d <- ana %>% filter(job == "pris2") %>%
  group_by(variant, chem) %>%
  summarise(contacts = sum(d_ne2 < NE2_CUT), .groups = "drop") %>%
  ggplot(aes(variant, contacts, fill = chem)) +
  geom_col(width = 0.62, colour = "black", linewidth = LW * 0.7) +
  geom_text(aes(label = contacts), vjust = -0.5, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = PAL_CHEM) +
  scale_y_continuous(limits = c(0, 158), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Residue at Position 2",
       y = "Structures on Pristinamycin 2 Within\n3.5 \u00c5 of His NE2 (of 150)")

save_fig(p35c + p35d + plot_annotation(tag_levels = list(c("C", "D"))),
         "Figure_3.5", 18, 9.5)

keep("3.5", ana %>% kruskal_test(d_ne2 ~ variant) %>%
       transmute(test = "Kruskal-Wallis", comparison = "five variants",
                 statistic, p, effect = NA_real_, effect_type = NA_character_,
                 n1 = n, n2 = NA_integer_))
# Each variant is compared against the unmodified analogue, not against the
# other variants, so Dunn's is run against a single control column.
keep("3.5", ana %>% dunn_test(d_ne2 ~ variant, p.adjust.method = "holm") %>%
       filter(group1 == "L-Thr") %>%
       transmute(test = "Dunn vs control", comparison = paste("L-Thr vs", group2),
                 statistic, p = p.adj, effect = NA_real_,
                 effect_type = NA_character_, n1, n2))
for (v in c("L-Ser", "D-Thr", "Dab", "Dap"))
  keep("3.5", mw(ana, "d_ne2", "variant", v, "L-Thr", paste(v, "vs L-Thr")))
ft5 <- ana %>% count(chem, nac) %>%
  pivot_wider(names_from = nac, values_from = n, values_fill = 0) %>%
  column_to_rownames("chem") %>% as.matrix() %>% fisher.test()
keep("3.5", tibble(test = "Fisher exact", comparison = "near-attack, hydroxyl vs amine",
                   statistic = NA_real_, p = ft5$p.value, effect = NA_real_,
                   effect_type = NA_character_, n1 = NA_integer_, n2 = NA_integer_))

ana %>% group_by(variant, chem) %>%
  summarise(n = n(), median_dNE2 = median(d_ne2), min_dNE2 = min(d_ne2),
            NE2_contacts_pris2 = sum(job == "pris2" & d_ne2 < NE2_CUT),
            near_attack = sum(nac), .groups = "drop") %>%
  write_csv(file.path(OUT, "stats", "Table_3.5.csv"))

# =============================================================================
# 3.6  Does confidence identify competent predictions?
# =============================================================================
message("Section 3.6")

# The question can only be asked where competent structures exist. The
# ramoplanin datasets contain none, so the test is restricted to the
# streptogramin native set.
SCORES <- c(nuc_plddt = "Nucleophile\npLDDT", ptm = "pTM",
            ranking_score = "Ranking\nScore", iptm = "ipTM")

auc_of <- function(d, score) {
  if (length(unique(d$nac)) < 2) return(NA_real_)
  suppressMessages(as.numeric(pROC::auc(pROC::roc(d$nac, d[[score]],
                                                  direction = "<", quiet = TRUE))))
}

auc_tbl <- bind_rows(
  map_dfr(names(SCORES), ~ tibble(score = SCORES[[.x]], scope = "Pooled",
                                  auc = auc_of(pris, .x))),
  map_dfr(c("pris1", "pris2"), function(j)
    map_dfr(names(SCORES), ~ tibble(score = SCORES[[.x]],
                                    scope = DOM_LAB[match(j, DOM)],
                                    auc = auc_of(filter(pris, job == j), .x))))
) %>% mutate(score = factor(score, levels = unname(SCORES)),
             scope = factor(scope, levels = c("Pooled", "Pristinamycin 1",
                                              "Pristinamycin 2")))

p36 <- ggplot(auc_tbl, aes(score, auc, fill = scope)) +
  geom_hline(yintercept = 0.5, linetype = "22", linewidth = LW * 0.9,
             colour = "grey20") +
  geom_col(position = position_dodge(0.82), width = 0.74,
           colour = "black", linewidth = LW * 0.7) +
  # Labels sit inside the bar, so they cannot collide with a neighbouring
  # label or with the 0.5 reference line.
  geom_text(aes(label = sprintf("%.2f", auc)),
            position = position_dodge(0.82), vjust = 1.5, size = 2.5,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL_SCOPE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Area Under the ROC Curve")
save_fig(p36, "Figure_3.6", 16, 9)

for (s in names(SCORES)) {
  r <- suppressMessages(pROC::roc(pris$nac, pris[[s]], direction = "<", quiet = TRUE))
  ci <- suppressMessages(pROC::ci.auc(r))
  keep("3.6", tibble(test = "ROC", comparison = paste(SCORES[[s]], "pooled"),
                     statistic = as.numeric(pROC::auc(r)),
                     p = NA_real_, effect = as.numeric(pROC::auc(r)),
                     effect_type = sprintf("AUC (95%% CI %.3f-%.3f)", ci[1], ci[3]),
                     n1 = sum(pris$nac), n2 = sum(!pris$nac)))
  for (j in c("pris1", "pris2")) {
    sub <- filter(pris, job == j)
    keep("3.6", tibble(test = "ROC",
                       comparison = paste(str_replace_all(SCORES[[s]], "\n", " "),
                                          "within", DOM_LAB[match(j, DOM)]),
                       statistic = auc_of(sub, s), p = NA_real_,
                       effect = auc_of(sub, s), effect_type = "AUC",
                       n1 = sum(sub$nac), n2 = sum(!sub$nac)))
  }
}
auc_tbl %>% mutate(score = str_replace_all(score, "\n", " ")) %>%
  write_csv(file.path(OUT, "stats", "Table_3.6.csv"))

# =============================================================================
message("Writing combined statistics")
bind_rows(STATS) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 4))) %>%
  write_csv(file.path(OUT, "stats", "all_statistics.csv"))
message("Done. Figures in ", file.path(OUT, "figures"),
        ", statistics in ", file.path(OUT, "stats"))

# =============================================================================
# Expected values, for checking the script ran correctly
# -----------------------------------------------------------------------------
# 3.1  ipTM means 0.648 / 0.715 / 0.729 / 0.696
#      nucleophile pLDDT 27.3 / 37.8 / 51.9 / 43.0
# 3.2  0 near-attack of 1,800
# 3.3  medians 24.25, 22.33, 18.15, 18.82, 14.52, 10.87
#      Spearman rho 0.864 pooled, 0.849 per seed (n = 903)
#      Kruskal-Wallis H = 2086.9, 0 near-attack of 2,715
# 3.4  medians 2.82, 2.81, 11.99, 11.31, 10.97, 15.26
#      near-attack 50, 97, 0, 0, 0, 0; rank-biserial -1.000; Fisher p = 2.5e-84
# 3.5  medians L-Thr 11.02, L-Ser 11.57, D-Thr 11.42, Dab 12.24, Dap 12.28
#      NE2 contacts on pris2: 129, 78, 55, 0, 0
# 3.6  pooled AUC: nucleophile pLDDT 0.914, pTM 0.744, ranking 0.693, ipTM 0.688
#      within pris1: 0.607, 0.450, 0.487, 0.477
#      within pris2: 0.594, 0.463, 0.499, 0.513
# =============================================================================
