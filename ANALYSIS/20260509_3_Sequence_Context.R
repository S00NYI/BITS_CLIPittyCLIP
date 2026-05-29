################################################################################
## CLIP Peak Sequence Context Analysis
## Author: Soon Yi with Antigravity
## Date: May 2026
################################################################################

# === Load Libraries ===
library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(ggseqlogo)

## 1. Basic Setup:
################################################################################
BASE_DIR = "/mnt/2TB_DATA/ANALYSIS_CLIPittyCLIP/"

INPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/FILTERED/")


PEAK_MATRIX_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "PEAKS_RBFOX2_BrdU_CLIP_FILTERED", ".txt")
PEAK_MATRIX_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "PEAKS_HuR_CoCLIP_FILTERED", ".txt")
PEAK_MATRIX_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "PEAKS_HNRNPC_iCLIP_FILTERED", ".txt")
PEAK_MATRIX_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "PEAKS_HNRNPC_BrdU_CLIP2_FILTERED", ".txt")

CTK_TRUNC_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "CTK_TRUNC_RBFOX2_BrdU_CLIP_FILTERED", ".txt")
CTK_TRUNC_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "CTK_TRUNC_HuR_CoCLIP_FILTERED", ".txt")
CTK_TRUNC_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "CTK_TRUNC_HNRNPC_iCLIP_FILTERED", ".txt")
CTK_TRUNC_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "CTK_TRUNC_HNRNPC_BrdU_CLIP2_FILTERED", ".txt")

# CTK_DEL_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CTK_DEL_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "HuR_CoCLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CTK_DEL_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "HNRNPC_iCLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CTK_DEL_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2", "_DEL_NORMALIZED_ANNOTATED", ".txt")

CLINK_TRUNC_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "CLINK_TRUNC_RBFOX2_BrdU_CLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "CLINK_TRUNC_HuR_CoCLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "CLINK_TRUNC_HNRNPC_iCLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "CLINK_TRUNC_HNRNPC_BrdU_CLIP2_FILTERED", ".txt")

# CLINK_DEL_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CLINK_DEL_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "HuR_CoCLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CLINK_DEL_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "HNRNPC_iCLIP", "_DEL_NORMALIZED_ANNOTATED", ".txt")
# CLINK_DEL_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2", "_DEL_NORMALIZED_ANNOTATED", ".txt")

OUTPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/SEQ_CONTEXT/")

genome = BSgenome.Hsapiens.UCSC.hg38
genome2 = BSgenome.Mmusculus.UCSC.mm10
################################################################################

## 2. Custom Functions:
################################################################################
# Compute per-position nucleotide frequency around peak centers (vectorized):
nucleotideFrequencyAroundPeaks = function(peaks_df, genome, window = 200, extension = c(0, 0)) {
  
  peaks_df_sense = peaks_df[peaks_df$strand == '+', ]
  peaks_df_antisense = peaks_df[peaks_df$strand == '-', ]
  
  ## Apply strand-aware extension: ext[1] = 5' end, ext[2] = 3' end
  if (nrow(peaks_df_sense) > 0) {
    peaks_df_sense$start = peaks_df_sense$start - extension[1]
    peaks_df_sense$end   = peaks_df_sense$end   + extension[2]
  }
  if (nrow(peaks_df_antisense) > 0) {
    peaks_df_antisense$start = peaks_df_antisense$start - extension[2]
    peaks_df_antisense$end   = peaks_df_antisense$end   + extension[1]
  }
  
  all_seqs = DNAStringSet()
  
  ## Sense strand: extract window around peak center
  if (nrow(peaks_df_sense) > 0) {
    chroms = as.character(peaks_df_sense$chr)
    centers = floor((peaks_df_sense$start + peaks_df_sense$end) / 2)
    
    starts = centers - window
    ends = centers + window - 1
    
    ## Filter valid ranges
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '+')
      seqs = getSeq(genome, gr)
      all_seqs = c(all_seqs, seqs)
    }
  }
  
  ## Antisense strand: extract and reverse complement
  if (nrow(peaks_df_antisense) > 0) {
    chroms = as.character(peaks_df_antisense$chr)
    centers = floor((peaks_df_antisense$start + peaks_df_antisense$end) / 2)
    
    starts = centers - window
    ends = centers + window - 1
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '-')
      seqs = getSeq(genome, gr)
      ## getSeq(strand='-') automatically reverse complements to orient 5'->3' related to RNA
      all_seqs = c(all_seqs, seqs)
    }
  }
  
  if (length(all_seqs) == 0) {
    warning("No valid sequences extracted")
    return(NULL)
  }
  
  ## consensusMatrix: rows = A/C/G/T, columns = positions
  cm = consensusMatrix(all_seqs, as.prob = TRUE)
  
  ## Extract only A, C, G, T rows (ignore N, etc.)
  positions = seq(-window, window - 1)
  
  result = data.frame(
    position = positions,
    A = cm["A", ],
    C = cm["C", ],
    G = cm["G", ],
    U = cm["T", ]   ## T on DNA = U on RNA
  )
  
  return(result)
}

# Plot stacked barplot of nucleotide composition around peaks:
plotNucleotideBar = function(freq_df, sampleName = NULL, window_lim = NULL, nuc_colors = NULL, x_ticks = NULL) {
  
  if (is.null(nuc_colors)) {
    nuc_colors = c("A" = "#4DAF4A", "C" = "#377EB8", "G" = "#FF7F00", "U" = "#E41A1C")
  }
  
  plot_data = freq_df %>%
    pivot_longer(cols = c(A, C, G, U), names_to = "Nucleotide", values_to = "Fraction")
  
  ## Set nucleotide order (bottom to top of stack)
  plot_data$Nucleotide = factor(plot_data$Nucleotide, levels = c("A", "C", "G", "U"))
  
  plot = ggplot(plot_data, aes(x = position, y = Fraction, fill = Nucleotide)) +
    geom_col(position = 'stack', width = 1) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 0.25, color = "gray40", linetype = "dotted", linewidth = 0.3) +
    scale_fill_manual(values = nuc_colors) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold'),
          legend.text = element_text(size = 14),
          panel.grid.minor = element_blank()) +
    labs(x = "Distance from peak center (nt)",
         y = "Nucleotide Fraction")
  
  if (!is.null(sampleName)) {
    plot = plot + ggtitle(paste0("Sequence Context around Peaks: ", sampleName))
  }
  
  if (!is.null(x_ticks) && length(x_ticks) == 2) {
    minor_val = x_ticks[1]
    major_val = x_ticks[2]
    max_pos = max(abs(freq_df$position))
    plot = plot + scale_x_continuous(breaks = seq(-max_pos, max_pos, by = major_val),
                                     minor_breaks = seq(-max_pos, max_pos, by = minor_val))
  }
  
  if (!is.null(window_lim)) {
    plot = plot + xlim(c(-window_lim, window_lim))
  }
  
  return(plot)
}

# Plot sequence logo using ggseqlogo (for narrow window):
plotSeqLogo = function(peaks_df, genome, window = 15, sampleName = NULL, extension = c(0, 0), method = 'bits', x_ticks = NULL, y_lim = NULL, title = NULL, x_label = NULL) {
  
  peaks_df_sense = peaks_df[peaks_df$strand == '+', ]
  peaks_df_antisense = peaks_df[peaks_df$strand == '-', ]
  
  if (nrow(peaks_df_sense) > 0) {
    peaks_df_sense$start = peaks_df_sense$start - extension[1]
    peaks_df_sense$end   = peaks_df_sense$end   + extension[2]
  }
  if (nrow(peaks_df_antisense) > 0) {
    peaks_df_antisense$start = peaks_df_antisense$start - extension[2]
    peaks_df_antisense$end   = peaks_df_antisense$end   + extension[1]
  }
  
  all_seqs = c()
  
  if (nrow(peaks_df_sense) > 0) {
    chroms = as.character(peaks_df_sense$chr)
    centers = floor((peaks_df_sense$start + peaks_df_sense$end) / 2)
    starts = centers - window
    ends = centers + window
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '+')
      seqs = getSeq(genome, gr)
      all_seqs = c(all_seqs, as.character(seqs))
    }
  }
  
  if (nrow(peaks_df_antisense) > 0) {
    chroms = as.character(peaks_df_antisense$chr)
    centers = floor((peaks_df_antisense$start + peaks_df_antisense$end) / 2)
    starts = centers - window
    ends = centers + window
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '-')
      seqs = getSeq(genome, gr)
      all_seqs = c(all_seqs, as.character(seqs))
    }
  }
  
  if (length(all_seqs) == 0) {
    warning("No valid sequences for logo")
    return(NULL)
  }
  
  all_seqs = gsub("T", "U", all_seqs)
  
  ## Resolve title: user override > hardcoded default (with optional sampleName suffix)
  default_prefix = "Sequence Logo around Peak Center"
  plot_title = if (!is.null(title)) title else default_prefix
  if (!is.null(sampleName)) plot_title = paste0(plot_title, ": ", sampleName)
  
  ## Resolve x-axis label: user override > hardcoded default
  default_x_prefix = "Position relative to peak center"
  plot_x_label = if (!is.null(x_label)) x_label else default_x_prefix
  plot_x_label = paste0(plot_x_label, " (±", window, " nt)")
  
  plot = ggseqlogo(all_seqs, method = method, seq_type = 'rna') +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = 'bold'),
          plot.title = element_text(size = 14, face = 'bold'),
          panel.grid.major.y = element_line(color = "grey80"),
          panel.grid.minor.y = element_line(color = "grey92"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle(plot_title) +
    xlab(plot_x_label)
  
  ## Y-axis limits
  if (!is.null(y_lim)) {
    plot = plot + coord_cartesian(ylim = y_lim)
  }
  
  if (!is.null(x_ticks) && length(x_ticks) == 2) {
    minor_val = x_ticks[1]
    major_val = x_ticks[2]
    plot = plot + scale_x_continuous(breaks = seq(1, 2*window + 1, by = major_val),
                                     minor_breaks = seq(1, 2*window + 1, by = minor_val),
                                     labels = seq(-window, window, by = major_val))
  }
  
  return(plot)
}
################################################################################

## 3. Load Pre-filtered Peaks/XL-Sites:
################################################################################
PEAKS_RBFOX2_BrdU_CLIP = fread(PEAK_MATRIX_FILE_RBFOX2_BrdU_CLIP)
PEAKS_HuR_CoCLIP = fread(PEAK_MATRIX_FILE_HuR_CoCLIP)
PEAKS_HNRNPC_iCLIP = fread(PEAK_MATRIX_FILE_HNRNPC_iCLIP)
PEAKS_HNRNPC_BrdU_CLIP2 = fread(PEAK_MATRIX_FILE_HNRNPC_BrdU_CLIP2)

TRUNC_CTK_RBFOX2_BrdU_CLIP = fread(CTK_TRUNC_FILE_RBFOX2_BrdU_CLIP)
TRUNC_CTK_HuR_CoCLIP = fread(CTK_TRUNC_FILE_HuR_CoCLIP)
TRUNC_CTK_HNRNPC_iCLIP = fread(CTK_TRUNC_FILE_HNRNPC_iCLIP)
TRUNC_CTK_HNRNPC_BrdU_CLIP2 = fread(CTK_TRUNC_FILE_HNRNPC_BrdU_CLIP2)

TRUNC_CLINK_RBFOX2_BrdU_CLIP = fread(CLINK_TRUNC_FILE_RBFOX2_BrdU_CLIP)
TRUNC_CLINK_HuR_CoCLIP = fread(CLINK_TRUNC_FILE_HuR_CoCLIP)
TRUNC_CLINK_HNRNPC_iCLIP = fread(CLINK_TRUNC_FILE_HNRNPC_iCLIP)
TRUNC_CLINK_HNRNPC_BrdU_CLIP2 = fread(CLINK_TRUNC_FILE_HNRNPC_BrdU_CLIP2)
# 
# DEL_CTK_RBFOX2_BrdU_CLIP = fread(CTK_DEL_FILE_RBFOX2_BrdU_CLIP)
# DEL_CTK_HuR_CoCLIP = fread(CTK_DEL_FILE_HuR_CoCLIP)
# DEL_CTK_HNRNPC_iCLIP = fread(CTK_DEL_FILE_HNRNPC_iCLIP)
# DEL_CTK_HNRNPC_BrdU_CLIP2 = fread(CTK_DEL_FILE_HNRNPC_BrdU_CLIP2)
# 
# DEL_CLINK_RBFOX2_BrdU_CLIP = fread(CLINK_DEL_FILE_RBFOX2_BrdU_CLIP)
# DEL_CLINK_HuR_CoCLIP = fread(CLINK_DEL_FILE_HuR_CoCLIP)
# DEL_CLINK_HNRNPC_iCLIP = fread(CLINK_DEL_FILE_HNRNPC_iCLIP)
# DEL_CLINK_HNRNPC_BrdU_CLIP2 = fread(CLINK_DEL_FILE_HNRNPC_BrdU_CLIP2)
################################################################################

## 4. Compute Nucleotide Frequencies (Broad View):
################################################################################
broad_window = 100

FREQ_PEAKS_RBFOX2_BrdU_CLIP = nucleotideFrequencyAroundPeaks(PEAKS_RBFOX2_BrdU_CLIP, genome2, broad_window)
FREQ_PEAKS_HuR_CoCLIP = nucleotideFrequencyAroundPeaks(PEAKS_HuR_CoCLIP, genome, broad_window)
FREQ_PEAKS_HNRNPC_iCLIP = nucleotideFrequencyAroundPeaks(PEAKS_HNRNPC_iCLIP, genome, broad_window)
FREQ_PEAKS_HNRNPC_BrdU_CLIP2 = nucleotideFrequencyAroundPeaks(PEAKS_HNRNPC_BrdU_CLIP2, genome, broad_window)

FREQ_TRUNC_CTK_RBFOX2_BrdU_CLIP = nucleotideFrequencyAroundPeaks(TRUNC_CTK_RBFOX2_BrdU_CLIP, genome2, broad_window)
FREQ_TRUNC_CTK_HuR_CoCLIP = nucleotideFrequencyAroundPeaks(TRUNC_CTK_HuR_CoCLIP, genome, broad_window)
FREQ_TRUNC_CTK_HNRNPC_iCLIP = nucleotideFrequencyAroundPeaks(TRUNC_CTK_HNRNPC_iCLIP, genome, broad_window)
FREQ_TRUNC_CTK_HNRNPC_BrdU_CLIP2 = nucleotideFrequencyAroundPeaks(TRUNC_CTK_HNRNPC_BrdU_CLIP2, genome, broad_window)

FREQ_TRUNC_CLINK_RBFOX2_BrdU_CLIP = nucleotideFrequencyAroundPeaks(TRUNC_CLINK_RBFOX2_BrdU_CLIP, genome2, broad_window)
FREQ_TRUNC_CLINK_HuR_CoCLIP = nucleotideFrequencyAroundPeaks(TRUNC_CLINK_HuR_CoCLIP, genome, broad_window)
FREQ_TRUNC_CLINK_HNRNPC_iCLIP = nucleotideFrequencyAroundPeaks(TRUNC_CLINK_HNRNPC_iCLIP, genome, broad_window)
FREQ_TRUNC_CLINK_HNRNPC_BrdU_CLIP2 = nucleotideFrequencyAroundPeaks(TRUNC_CLINK_HNRNPC_BrdU_CLIP2, genome, broad_window)

# FREQ_DEL_CTK_RBFOX2_BrdU_CLIP = nucleotideFrequencyAroundPeaks(DEL_CTK_RBFOX2_BrdU_CLIP, genome2, broad_window)
# FREQ_DEL_CTK_HuR_CoCLIP = nucleotideFrequencyAroundPeaks(DEL_CTK_HuR_CoCLIP, genome, broad_window)
# FREQ_DEL_CTK_HNRNPC_iCLIP = nucleotideFrequencyAroundPeaks(DEL_CTK_HNRNPC_iCLIP, genome, broad_window)
# FREQ_DEL_CTK_HNRNPC_BrdU_CLIP2 = nucleotideFrequencyAroundPeaks(DEL_CTK_HNRNPC_BrdU_CLIP2, genome, broad_window)
# 
# FREQ_DEL_CLINK_RBFOX2_BrdU_CLIP = nucleotideFrequencyAroundPeaks(DEL_CLINK_RBFOX2_BrdU_CLIP, genome2, broad_window)
# FREQ_DEL_CLINK_HuR_CoCLIP = nucleotideFrequencyAroundPeaks(DEL_CLINK_HuR_CoCLIP, genome, broad_window)
# FREQ_DEL_CLINK_HNRNPC_iCLIP = nucleotideFrequencyAroundPeaks(DEL_CLINK_HNRNPC_iCLIP, genome, broad_window)
# FREQ_DEL_CLINK_HNRNPC_BrdU_CLIP2 = nucleotideFrequencyAroundPeaks(DEL_CLINK_HNRNPC_BrdU_CLIP2, genome, broad_window)
################################################################################

## 5. Stacked Barplots (Broad View):
################################################################################
print(plotNucleotideBar(FREQ_PEAKS_RBFOX2_BrdU_CLIP, "RBFOX2_BrdU_CLIP Peak"))
print(plotNucleotideBar(FREQ_PEAKS_HuR_CoCLIP, "HuR_CoCLIP Peak"))
print(plotNucleotideBar(FREQ_PEAKS_HNRNPC_iCLIP, "HNRNPC_iCLIP Peak"))
print(plotNucleotideBar(FREQ_PEAKS_HNRNPC_BrdU_CLIP2, "HNRNPC_BrdU_CLIP2 Peak"))

print(plotNucleotideBar(FREQ_TRUNC_CTK_RBFOX2_BrdU_CLIP, "RBFOX2_BrdU_CLIP CTK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CTK_HuR_CoCLIP, "HuR_CoCLIP CTK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CTK_HNRNPC_iCLIP, "HNRNPC_iCLIP CTK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CTK_HNRNPC_BrdU_CLIP2, "HNRNPC_BrdU_CLIP2 CTK Truncation"))

print(plotNucleotideBar(FREQ_TRUNC_CLINK_RBFOX2_BrdU_CLIP, "RBFOX2_BrdU_CLIP CLINK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CLINK_HuR_CoCLIP, "HuR_CoCLIP CLINK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CLINK_HNRNPC_iCLIP, "HNRNPC_iCLIP CLINK Truncation"))
print(plotNucleotideBar(FREQ_TRUNC_CLINK_HNRNPC_BrdU_CLIP2, "HNRNPC_BrdU_CLIP2 CLINK Truncation"))

# print(plotNucleotideBar(FREQ_DEL_CTK_RBFOX2_BrdU_CLIP, "RBFOX2_BrdU_CLIP CTK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CTK_HuR_CoCLIP, "HuR_CoCLIP CTK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CTK_HNRNPC_iCLIP, "HNRNPC_iCLIP CTK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CTK_HNRNPC_BrdU_CLIP2, "HNRNPC_BrdU_CLIP2 CTK Deletion"))
# 
# print(plotNucleotideBar(FREQ_DEL_CLINK_RBFOX2_BrdU_CLIP, "RBFOX2_BrdU_CLIP CLINK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CLINK_HuR_CoCLIP, "HuR_CoCLIP CLINK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CLINK_HNRNPC_iCLIP, "HNRNPC_iCLIP CLINK Deletion"))
# print(plotNucleotideBar(FREQ_DEL_CLINK_HNRNPC_BrdU_CLIP2, "HNRNPC_BrdU_CLIP2 CLINK Deletion"))

################################################################################

## 6. Sequence Logos (Narrow View):
################################################################################
narrow_window = 50
extension = c(0, 0)
x_ticks = c(5, 10)
# logo_method = 'prob'
logo_method = 'bits'

print(plotSeqLogo(PEAKS_RBFOX2_BrdU_CLIP, genome2, narrow_window, "RBFOX2_BrdU_CLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 1)))
print(plotSeqLogo(PEAKS_HuR_CoCLIP, genome, narrow_window, "HuR_CoCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 1)))
print(plotSeqLogo(PEAKS_HNRNPC_iCLIP, genome, narrow_window, "HNRNPC_iCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 1)))
print(plotSeqLogo(PEAKS_HNRNPC_BrdU_CLIP2, genome, narrow_window, "HNRNPC_BrdU_CLIP2", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 1)))

print(plotSeqLogo(TRUNC_CTK_RBFOX2_BrdU_CLIP, genome2, narrow_window, "RBFOX2_BrdU_CLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.25), title = 'Sequence Logo around CTK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CTK_HuR_CoCLIP, genome, narrow_window, "HuR_CoCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CTK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CTK_HNRNPC_iCLIP, genome, narrow_window, "HNRNPC_iCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CTK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CTK_HNRNPC_BrdU_CLIP2, genome, narrow_window, "HNRNPC_BrdU_CLIP2", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CTK Truncation Site', x_label = 'Position relative to Truncation site'))

print(plotSeqLogo(TRUNC_CLINK_RBFOX2_BrdU_CLIP, genome2, narrow_window, "RBFOX2_BrdU_CLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.25), title = 'Sequence Logo around CLINK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CLINK_HuR_CoCLIP, genome, narrow_window, "HuR_CoCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CLINK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CLINK_HNRNPC_iCLIP, genome, narrow_window, "HNRNPC_iCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CLINK Truncation Site', x_label = 'Position relative to Truncation site'))
print(plotSeqLogo(TRUNC_CLINK_HNRNPC_BrdU_CLIP2, genome, narrow_window, "HNRNPC_BrdU_CLIP2", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 2), title = 'Sequence Logo around CLINK Truncation Site', x_label = 'Position relative to Truncation site'))

# print(plotSeqLogo(DEL_CTK_RBFOX2_BrdU_CLIP, genome2, narrow_window, "RBFOX2_BrdU_CLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CTK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CTK_HuR_CoCLIP, genome, narrow_window, "HuR_CoCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CTK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CTK_HNRNPC_iCLIP, genome, narrow_window, "HNRNPC_iCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CTK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CTK_HNRNPC_BrdU_CLIP2, genome, narrow_window, "HNRNPC_BrdU_CLIP2", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CTK Deletion Site', x_label = 'Position relative to Deletion site'))
# 
# print(plotSeqLogo(DEL_CLINK_RBFOX2_BrdU_CLIP, genome2, narrow_window, "RBFOX2_BrdU_CLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CLINK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CLINK_HuR_CoCLIP, genome, narrow_window, "HuR_CoCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CLINK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CLINK_HNRNPC_iCLIP, genome, narrow_window, "HNRNPC_iCLIP", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CLINK Deletion Site', x_label = 'Position relative to Deletion site'))
# print(plotSeqLogo(DEL_CLINK_HNRNPC_BrdU_CLIP2, genome, narrow_window, "HNRNPC_BrdU_CLIP2", extension = extension, method = logo_method, x_ticks = x_ticks, y_lim = c(0, 0.5), title = 'Sequence Logo around CLINK Deletion Site', x_label = 'Position relative to Deletion site'))

################################################################################


