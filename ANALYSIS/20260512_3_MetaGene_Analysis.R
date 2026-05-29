################################################################################
## CLIP Metagene Analysis
## Author: Soon Yi with Antigravity
## Date: May 2026
################################################################################

# === Load Libraries ===
library(dplyr)
library(tidyr)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ggplot2)
library(stringr)

## 1. Basic Setup:
################################################################################
BASE_DIR = "/mnt/2TB_DATA/ANALYSIS_CLIPittyCLIP/"

INPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/2_FILTERED/")

HUR_CoCLIP_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HuR_CoCLIP_PEAK_ANNOTATED_FILTERED.txt"))
HUR_CoCLIP_CTK_TRUNC_FILTERED   = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"))
HUR_CoCLIP_CTK_DEL_FILTERED     = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CTK_DEL_ANNOTATED_FILTERED.txt"))
HUR_CoCLIP_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HUR_CoCLIP_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"))

HUR_eCLIP_K562_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"))
HUR_eCLIP_K562_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HUR_eCLIP_K562_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"))

HNRNPC_BrdU_CLIP2_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_PEAK_ANNOTATED_FILTERED.txt"))
HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED     = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_DEL_ANNOTATED_FILTERED.txt"))
HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_DEL_ANNOTATED_FILTERED.txt"))

HNRNPC_iCLIP_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_PEAK_ANNOTATED_FILTERED.txt"))
HNRNPC_iCLIP_CTK_TRUNC_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_iCLIP_CTK_DEL_FILTERED     = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CTK_DEL_ANNOTATED_FILTERED.txt"))
HNRNPC_iCLIP_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_iCLIP_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"))

HNRNPC_eCLIP_HepG2_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_PEAK_ANNOTATED_FILTERED.txt"))
HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_DEL_ANNOTATED_FILTERED.txt"))

HNRNPC_eCLIP_K562_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"))
HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"))

RBFOX2_BrdU_CLIP_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_PEAK_ANNOTATED_FILTERED.txt"))
RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED   = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"))
RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED     = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_DEL_ANNOTATED_FILTERED.txt"))
RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"))

RBFOX2_eCLIP_HepG2_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_PEAK_ANNOTATED_FILTERED.txt"))
RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_DEL_ANNOTATED_FILTERED.txt"))

RBFOX2_eCLIP_K562_PEAK_FILTERED        = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"))
RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"))
RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED   = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"))

hs_genome = BSgenome.Hsapiens.UCSC.hg38
mm_genome = BSgenome.Mmusculus.UCSC.mm39

################################################################################

## 2. Custom Functions:
################################################################################
motifDensityAroundPeaks = function(peaks_df, genome, motif, window_definition = 100, normalize_by = c("proportion", "raw"), qval = NULL, fixed = TRUE) {
  
  normalize_by = match.arg(normalize_by)
  
  ## Convert RNA motif to DNA:
  dna_motif = gsub("U", "T", motif)
  motif_len = nchar(dna_motif)
  
  ## Optional qvalue filtering:
  if (!is.null(qval)) {
    before = nrow(peaks_df)
    peaks_df = peaks_df %>% filter(qvalue <= qval)
    message("  qvalue <= ", qval, ": ", nrow(peaks_df), " / ", before, " sites kept")
  }
  
  positions    = -window_definition:window_definition
  count_totals = rep(0, length(positions))
  
  ## Build GRanges, center on peak, extend flanks:
  gr = GRanges(
    seqnames = peaks_df$chr,
    ranges   = IRanges(start = peaks_df$start, end = peaks_df$end),
    strand   = peaks_df$strand
  )
  
  ## Center each site and extend by window_definition on each side:
  gr = resize(gr, width = 1, fix = "center")
  gr = resize(gr, width = 2 * window_definition + 1, fix = "center")
  
  ## Keep only standard chromosomes and sites with full window:
  gr = keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr = trim(gr)
  gr = gr[width(gr) == 2 * window_definition + 1]
  
  n_valid = length(gr)
  message("  Sites used after filtering: ", n_valid)
  
  if (n_valid == 0) {
    warning("No valid sites remain after filtering.")
    return(data.frame(midpoint = positions, count = rep(0, length(positions))))
  }
  
  ## Strand-aware getSeq: minus-strand sequences are automatically reverse complemented
  ## so the returned sequence is always 5'->3' RNA orientation:
  seqs = getSeq(genome, gr)
  
  ## Search for dna_motif on each sequence:
  for (i in seq_along(seqs)) {
    matches = matchPattern(dna_motif, seqs[[i]], fixed = fixed)
    if (length(matches) > 0) {
      ## Convert to positions relative to center (center = 0):
      match_starts_rel = start(matches) - window_definition - 1
      valid_pos = match_starts_rel[match_starts_rel >= -window_definition &
                                     match_starts_rel <= window_definition]
      if (length(valid_pos) > 0) {
        count_totals = count_totals + tabulate(valid_pos + window_definition + 1,
                                               nbins = length(positions))
      }
    }
  }
  
  result = data.frame(
    midpoint = positions,
    count    = if (normalize_by == "proportion") count_totals / n_valid else count_totals
  )
  
  return(result)
}

combineDensity = function(..., names) {
  density_list = list(...)
  stopifnot(length(density_list) == length(names))
  
  combined = density_list[[1]][, "midpoint", drop = FALSE]
  colnames(combined) = "position"
  
  for (i in seq_along(density_list)) {
    combined[[names[i]]] = density_list[[i]]$count
  }
  return(combined)
}

plotDensity = function(density_data, columns_list, xaxis_lims = NULL, yaxis_lims = NULL, custom_colors = NULL, custom_linetypes = NULL, custom_linewidths = NULL, densityType = NULL, sampleName = NULL, featureName = NULL, smoothing = NULL) {
  
  plot_data = density_data[, c('position', columns_list)]
  plot_data_long = plot_data %>%
    pivot_longer(cols = {{columns_list}}, names_to = "Data", values_to = "Density")
  plot_data_long$Data = factor(plot_data_long$Data, levels = columns_list)
  
  plot = ggplot(plot_data_long, aes(x = position, y = Density, color = Data)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme_bw() +
    theme(axis.text  = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(y = "Peak Density")
  
  if (!is.null(smoothing)) {
    plot = plot + geom_smooth(span = smoothing, se = FALSE, level = 0.9)
  } else {
    if (!is.null(custom_linetypes) || !is.null(custom_linewidths)) {
      ## Map linetype and/or linewidth per series via aes
      plot_data_long$linetype_val  = custom_linetypes[as.character(plot_data_long$Data)]
      plot_data_long$linewidth_val = custom_linewidths[as.character(plot_data_long$Data)]
      plot = ggplot(plot_data_long, aes(x = position, y = Density, color = Data,
                                        linetype = Data, linewidth = Data)) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
        geom_line() +
        theme_minimal() +
        theme_bw() +
        theme(axis.text  = element_text(size = 14),
              axis.title = element_text(size = 14, face = 'bold')) +
        labs(y = "Peak Density")
      if (!is.null(custom_linetypes))  plot = plot + scale_linetype_manual(values = custom_linetypes)
      if (!is.null(custom_linewidths)) plot = plot + scale_linewidth_manual(values = custom_linewidths)
    } else {
      plot = plot + geom_line(linewidth = 1)
    }
  }
  
  if (!is.null(custom_colors)) {
    plot = plot + scale_color_manual(values = custom_colors)
  }
  if (!is.null(xaxis_lims)) plot = plot + xlim(xaxis_lims)
  if (!is.null(yaxis_lims)) plot = plot + ylim(yaxis_lims)
  
  ## Title/label logic unchanged
  if (is.null(densityType)) {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Density for ', sampleName),
                         x = "Distance from center (nucleotides)")
    } else {
      plot = plot + labs(title = 'Density Plot',
                         x = "Distance from center (nucleotides)")
    }
  } else if (densityType == 'motif_density') {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Motif Density around Peaks for ', sampleName),
                         x = "Distance from peak center (nucleotides)")
    } else {
      plot = plot + labs(title = 'Motif Density around Peaks',
                         x = "Distance from peak center (nucleotides)")
    }
  } else if (densityType == 'feature_metagene') {
    if (!is.null(sampleName)) {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName, ' for ', sampleName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    } else {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    }
  } else if (densityType == 'nucleotide_content') {
    if (!is.null(sampleName)) {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0(featureName, ' for ', sampleName),
                           x = "Distance from center (nucleotides)", y = "Nucleotide Fraction")
      } else {
        plot = plot + labs(title = paste0('Nucleotide Content for ', sampleName),
                           x = "Distance from center (nucleotides)", y = "Nucleotide Fraction")
      }
    } else {
      if (!is.null(featureName)) {
        plot = plot + labs(title = featureName,
                           x = "Distance from center (nucleotides)", y = "Nucleotide Fraction")
      } else {
        plot = plot + labs(title = 'Nucleotide Content',
                           x = "Distance from center (nucleotides)", y = "Nucleotide Fraction")
      }
    }
  }
  
  return(plot)
}
################################################################################

## 3. Motif Density Plot Options
#########################################################################################################################
wd = 100
ylimits = c(0, 0.5)
smoothing = NULL

HUR_Motif    = "UUUUUU"
HNRNPC_Motif = "UUUUUU"
RBFOX2_Motif = "UGCAUG"

plot_names = c("Peak", "CTK Truncation", "CTK Deletion", "CLINK Truncation", "CLINK Deletion")

plot_colors = c(
  "Peak"             = "gray",
  "CTK Truncation"   = "salmon",
  "CTK Deletion"     = "skyblue",
  "CLINK Truncation" = "red",
  "CLINK Deletion"   = "blue"
)

plot_linetypes = c(
  "Peak"             = "solid",
  "CTK Truncation"   = "solid",
  "CTK Deletion"     = "solid",
  "CLINK Truncation" = "solid",
  "CLINK Deletion"   = "solid"
)

plot_linewidths = c(
  "Peak"             = 0.8,
  "CTK Truncation"   = 0.8,
  "CTK Deletion"     = 0.5,
  "CLINK Truncation" = 0.8,
  "CLINK Deletion"   = 0.5
)

eclip_plot_names = c("Peak", "CLINK Truncation", "CLINK Deletion")

eclip_plot_colors = c(
  "Peak"             = "gray",
  "CLINK Truncation" = "red",
  "CLINK Deletion"   = "blue"
)

eclip_plot_linetypes = c(
  "Peak"             = "solid",
  "CLINK Truncation" = "solid",
  "CLINK Deletion"   = "solid"
)

eclip_plot_linewidths = c(
  "Peak"             = 0.8,
  "CLINK Truncation" = 0.8,
  "CLINK Deletion"   = 0.8
)


#########################################################################################################################

## 4. HUR CoCLIP
#########################################################################################################################
DENSITY_HUR_CoCLIP_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HUR_CoCLIP_PEAK_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_CoCLIP_CTK_TRUNC_FILTERED   = motifDensityAroundPeaks(peaks_df = HUR_CoCLIP_CTK_TRUNC_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_CoCLIP_CTK_DEL_FILTERED     = motifDensityAroundPeaks(peaks_df = HUR_CoCLIP_CTK_DEL_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_CoCLIP_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HUR_CoCLIP_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_CoCLIP_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HUR_CoCLIP_CLINK_DEL_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HUR_CoCLIP = combineDensity(
  DENSITY_HUR_CoCLIP_PEAK_FILTERED,
  DENSITY_HUR_CoCLIP_CTK_TRUNC_FILTERED,
  DENSITY_HUR_CoCLIP_CTK_DEL_FILTERED,
  DENSITY_HUR_CoCLIP_CLINK_TRUNC_FILTERED,
  DENSITY_HUR_CoCLIP_CLINK_DEL_FILTERED,
  names = plot_names
)

plotDensity(
  density_data   = DENSITY_HUR_CoCLIP,
  columns_list   = plot_names,
  custom_colors  = plot_colors,
  custom_linetypes = plot_linetypes,
  custom_linewidths = plot_linewidths,
  densityType    = "motif_density",
  sampleName     = "HuR CoCLIP Input",
  yaxis_lims     = c(0, 0.25),
  smoothing      = smoothing
)

#########################################################################################################################

## 5. HUR eCLIP K562
#########################################################################################################################
DENSITY_HUR_eCLIP_K562_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HUR_eCLIP_K562_PEAK_FILTERED,        genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_eCLIP_K562_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HUR_eCLIP_K562_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HUR_eCLIP_K562_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HUR_eCLIP_K562_CLINK_DEL_FILTERED,   genome = hs_genome, motif = HUR_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HUR_eCLIP_K562 = combineDensity(
  DENSITY_HUR_eCLIP_K562_PEAK_FILTERED,
  DENSITY_HUR_eCLIP_K562_CLINK_TRUNC_FILTERED,
  DENSITY_HUR_eCLIP_K562_CLINK_DEL_FILTERED,
  names = eclip_plot_names
)

plotDensity(
  density_data  = DENSITY_HUR_eCLIP_K562,
  columns_list  = eclip_plot_names,
  custom_colors = eclip_plot_colors,
  custom_linetypes = eclip_plot_linetypes,
  custom_linewidths = eclip_plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "HuR eCLIP K562",
  yaxis_lims    = ylimits,
  smoothing     = smoothing
)
#########################################################################################################################

## 6. HNRNPC BrdU CLIP2
#########################################################################################################################
DENSITY_HNRNPC_BrdU_CLIP2_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HNRNPC_BrdU_CLIP2_PEAK_FILTERED,        genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED     = motifDensityAroundPeaks(peaks_df = HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED,     genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HNRNPC_BrdU_CLIP2 = combineDensity(
  DENSITY_HNRNPC_BrdU_CLIP2_PEAK_FILTERED,
  DENSITY_HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED,
  DENSITY_HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED,
  DENSITY_HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED,
  DENSITY_HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED,
  names = plot_names
)

plotDensity(
  density_data  = DENSITY_HNRNPC_BrdU_CLIP2,
  columns_list  = plot_names,
  custom_colors = plot_colors,
  custom_linetypes = plot_linetypes,
  custom_linewidths = plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "HNRNPC Specificity CLIP",
  yaxis_lims    = c(0, 0.5),
  smoothing     = smoothing
)
#########################################################################################################################

## 7. HNRNPC iCLIP
#########################################################################################################################
DENSITY_HNRNPC_iCLIP_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HNRNPC_iCLIP_PEAK_FILTERED,        genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_iCLIP_CTK_TRUNC_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_iCLIP_CTK_TRUNC_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_iCLIP_CTK_DEL_FILTERED     = motifDensityAroundPeaks(peaks_df = HNRNPC_iCLIP_CTK_DEL_FILTERED,     genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_iCLIP_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HNRNPC_iCLIP_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_iCLIP_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_iCLIP_CLINK_DEL_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HNRNPC_iCLIP = combineDensity(
  DENSITY_HNRNPC_iCLIP_PEAK_FILTERED,
  DENSITY_HNRNPC_iCLIP_CTK_TRUNC_FILTERED,
  DENSITY_HNRNPC_iCLIP_CTK_DEL_FILTERED,
  DENSITY_HNRNPC_iCLIP_CLINK_TRUNC_FILTERED,
  DENSITY_HNRNPC_iCLIP_CLINK_DEL_FILTERED,
  names = plot_names
)

plotDensity(
  density_data  = DENSITY_HNRNPC_iCLIP,
  columns_list  = plot_names,
  custom_colors = plot_colors,
  custom_linetypes = plot_linetypes,
  custom_linewidths = plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "HNRNPC iCLIP",
  yaxis_lims    = c(0, 0.5),
  smoothing     = smoothing
)
#########################################################################################################################

## 8. HNRNPC eCLIP HepG2
#########################################################################################################################
DENSITY_HNRNPC_eCLIP_HepG2_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_HepG2_PEAK_FILTERED,        genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HNRNPC_eCLIP_HepG2 = combineDensity(
  DENSITY_HNRNPC_eCLIP_HepG2_PEAK_FILTERED,
  DENSITY_HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED,
  DENSITY_HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED,
  names = eclip_plot_names
)

plotDensity(
  density_data  = DENSITY_HNRNPC_eCLIP_HepG2,
  columns_list  = eclip_plot_names,
  custom_colors = eclip_plot_colors,
  custom_linetypes = eclip_plot_linetypes,
  custom_linewidths = eclip_plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "HNRNPC eCLIP HepG2",
  yaxis_lims    = c(0, 1),
  smoothing     = smoothing
)
#########################################################################################################################

## 9. HNRNPC eCLIP K562
#########################################################################################################################
DENSITY_HNRNPC_eCLIP_K562_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_K562_PEAK_FILTERED,        genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED,   genome = hs_genome, motif = HNRNPC_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_HNRNPC_eCLIP_K562 = combineDensity(
  DENSITY_HNRNPC_eCLIP_K562_PEAK_FILTERED,
  DENSITY_HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED,
  DENSITY_HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED,
  names = eclip_plot_names
)

plotDensity(
  density_data  = DENSITY_HNRNPC_eCLIP_K562,
  columns_list  = eclip_plot_names,
  custom_colors = eclip_plot_colors,
  custom_linetypes = eclip_plot_linetypes,
  custom_linewidths = eclip_plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "HNRNPC eCLIP K562",
  yaxis_lims    = c(0, 1),
  smoothing     = smoothing
)
#########################################################################################################################

## 10. RBFOX2 BrdU CLIP
#########################################################################################################################
DENSITY_RBFOX2_BrdU_CLIP_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = RBFOX2_BrdU_CLIP_PEAK_FILTERED,        genome = mm_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED   = motifDensityAroundPeaks(peaks_df = RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED,   genome = mm_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED     = motifDensityAroundPeaks(peaks_df = RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED,     genome = mm_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED, genome = mm_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED,   genome = mm_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_RBFOX2_BrdU_CLIP = combineDensity(
  DENSITY_RBFOX2_BrdU_CLIP_PEAK_FILTERED,
  DENSITY_RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED,
  DENSITY_RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED,
  DENSITY_RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED,
  DENSITY_RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED,
  names = plot_names
)

plotDensity(
  density_data  = DENSITY_RBFOX2_BrdU_CLIP,
  columns_list  = plot_names,
  custom_colors = plot_colors,
  custom_linetypes = plot_linetypes,
  custom_linewidths = plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "RBFOX2 BrdU CLIP",
  yaxis_lims    = c(0, 0.05),
  smoothing     = smoothing
)
#########################################################################################################################

## 11. RBFOX2 eCLIP HepG2
#########################################################################################################################
DENSITY_RBFOX2_eCLIP_HepG2_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_HepG2_PEAK_FILTERED,        genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED,   genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_RBFOX2_eCLIP_HepG2 = combineDensity(
  DENSITY_RBFOX2_eCLIP_HepG2_PEAK_FILTERED,
  DENSITY_RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED,
  DENSITY_RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED,
  names = eclip_plot_names
)

plotDensity(
  density_data  = DENSITY_RBFOX2_eCLIP_HepG2,
  columns_list  = eclip_plot_names,
  custom_colors = eclip_plot_colors,
  custom_linetypes = eclip_plot_linetypes,
  custom_linewidths = eclip_plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "RBFOX2 eCLIP HepG2",
  yaxis_lims    = c(0, 0.1),
  smoothing     = smoothing
)
#########################################################################################################################

## 12. RBFOX2 eCLIP K562
#########################################################################################################################
DENSITY_RBFOX2_eCLIP_K562_PEAK_FILTERED        = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_K562_PEAK_FILTERED,        genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED, genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")
DENSITY_RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED   = motifDensityAroundPeaks(peaks_df = RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED,   genome = hs_genome, motif = RBFOX2_Motif, window_definition = wd, fixed = TRUE, normalize_by = "proportion")

DENSITY_RBFOX2_eCLIP_K562 = combineDensity(
  DENSITY_RBFOX2_eCLIP_K562_PEAK_FILTERED,
  DENSITY_RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED,
  DENSITY_RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED,
  names = eclip_plot_names
)

plotDensity(
  density_data  = DENSITY_RBFOX2_eCLIP_K562,
  columns_list  = eclip_plot_names,
  custom_colors = eclip_plot_colors,
  custom_linetypes = eclip_plot_linetypes,
  custom_linewidths = eclip_plot_linewidths,
  densityType   = "motif_density",
  sampleName    = "RBFOX2 eCLIP K562",
  yaxis_lims    = c(0, 0.05),
  smoothing     = smoothing
)
#########################################################################################################################

