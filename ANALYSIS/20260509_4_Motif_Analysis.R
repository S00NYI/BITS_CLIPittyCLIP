################################################################################
## CLIP Peak Motif Analysis
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(RBPSpecificity)
library(ggplot2)
library(dplyr)
library(data.table)
library(seqLogo)
library(Biostrings)
library(cowplot)
library(grid)
library(gridExtra)

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

CLINK_TRUNC_FILE_RBFOX2_BrdU_CLIP = paste0(INPUT_DIR, "CLINK_TRUNC_RBFOX2_BrdU_CLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HuR_CoCLIP = paste0(INPUT_DIR, "CLINK_TRUNC_HuR_CoCLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HNRNPC_iCLIP = paste0(INPUT_DIR, "CLINK_TRUNC_HNRNPC_iCLIP_FILTERED", ".txt")
CLINK_TRUNC_FILE_HNRNPC_BrdU_CLIP2 = paste0(INPUT_DIR, "CLINK_TRUNC_HNRNPC_BrdU_CLIP2_FILTERED", ".txt")


OUTPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/MOTIF/")

K_MER = 5

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

################################################################################

## Peak Test 1: ANR, Scramble Off
################################################################################
extensions = c(5, -15)

MOTIF_PEAK_RBFOX2_BrdU_CLIP_ANR = motifEnrichment(PEAKS_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HUR_CoCLIP_ANR = motifEnrichment(PEAKS_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HNRNPC_iCLIP_ANR = motifEnrichment(PEAKS_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ANR= motifEnrichment(PEAKS_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_PEAK_RBFOX2_BrdU_CLIP_ANR$Score); MOTIF_PEAK_RBFOX2_BrdU_CLIP_ANR$MOTIF[MOTIF_PEAK_RBFOX2_BrdU_CLIP_ANR$Score == 1]
1 / median(MOTIF_PEAK_HUR_CoCLIP_ANR$Score); MOTIF_PEAK_HUR_CoCLIP_ANR$MOTIF[MOTIF_PEAK_HUR_CoCLIP_ANR$Score == 1]
1 / median(MOTIF_PEAK_HNRNPC_iCLIP_ANR$Score); MOTIF_PEAK_HNRNPC_iCLIP_ANR$MOTIF[MOTIF_PEAK_HNRNPC_iCLIP_ANR$Score == 1]
1 / median(MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ANR$Score); MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ANR$MOTIF[MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ANR$Score == 1]
################################################################################

## Peak Test 2: ZOOPS, Scramble Off
################################################################################
extensions = c(20, 0)

MOTIF_PEAK_RBFOX2_BrdU_CLIP_ZOOPS = motifEnrichment(PEAKS_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HUR_CoCLIP_ZOOPS = motifEnrichment(PEAKS_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HNRNPC_iCLIP_ZOOPS = motifEnrichment(PEAKS_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ZOOPS= motifEnrichment(PEAKS_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_PEAK_RBFOX2_BrdU_CLIP_ZOOPS$Score); MOTIF_PEAK_RBFOX2_BrdU_CLIP_ZOOPS$MOTIF[MOTIF_PEAK_RBFOX2_BrdU_CLIP_ZOOPS$Score == 1]
1 / median(MOTIF_PEAK_HUR_CoCLIP_ZOOPS$Score); MOTIF_PEAK_HUR_CoCLIP_ZOOPS$MOTIF[MOTIF_PEAK_HUR_CoCLIP_ZOOPS$Score == 1]
1 / median(MOTIF_PEAK_HNRNPC_iCLIP_ZOOPS$Score); MOTIF_PEAK_HNRNPC_iCLIP_ZOOPS$MOTIF[MOTIF_PEAK_HNRNPC_iCLIP_ZOOPS$Score == 1]
1 / median(MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ZOOPS$Score); MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ZOOPS$MOTIF[MOTIF_PEAK_HNRNPC_BrdU_CLIP2_ZOOPS$Score == 1]
################################################################################



## CTK TRUNC Test 1: ANR, Scramble Off
################################################################################
extensions = c(5, -15)

MOTIF_CTK_RBFOX2_BrdU_CLIP_ANR = motifEnrichment(TRUNC_CTK_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HUR_CoCLIP_ANR = motifEnrichment(TRUNC_CTK_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HNRNPC_iCLIP_ANR = motifEnrichment(TRUNC_CTK_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HNRNPC_BrdU_CLIP2_ANR= motifEnrichment(TRUNC_CTK_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_CTK_RBFOX2_BrdU_CLIP_ANR$Score)
1 / median(MOTIF_CTK_HUR_CoCLIP_ANR$Score)
1 / median(MOTIF_CTK_HNRNPC_iCLIP_ANR$Score)
1 / median(MOTIF_CTK_HNRNPC_BrdU_CLIP2_ANR$Score)
################################################################################

## CTK TRUNC Test 2: ZOOPS, Scramble Off
################################################################################
extensions = c(5, -15)

MOTIF_CTK_RBFOX2_BrdU_CLIP_ZOOPS = motifEnrichment(TRUNC_CTK_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HUR_CoCLIP_ZOOPS = motifEnrichment(TRUNC_CTK_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HNRNPC_iCLIP_ZOOPS = motifEnrichment(TRUNC_CTK_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CTK_HNRNPC_BrdU_CLIP2_ZOOPS= motifEnrichment(TRUNC_CTK_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_CTK_RBFOX2_BrdU_CLIP_ZOOPS$Score)
1 / median(MOTIF_CTK_HUR_CoCLIP_ZOOPS$Score)
1 / median(MOTIF_CTK_HNRNPC_iCLIP_ZOOPS$Score)
1 / median(MOTIF_CTK_HNRNPC_BrdU_CLIP2_ZOOPS$Score)
################################################################################




## CLINK TRUNC Test 1: ANR, Scramble Off
################################################################################
extensions = c(5, -15)

MOTIF_CLINK_RBFOX2_BrdU_CLIP_ANR = motifEnrichment(TRUNC_CLINK_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HUR_CoCLIP_ANR = motifEnrichment(TRUNC_CLINK_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HNRNPC_iCLIP_ANR = motifEnrichment(TRUNC_CLINK_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HNRNPC_BrdU_CLIP2_ANR= motifEnrichment(TRUNC_CLINK_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'anr', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_CLINK_RBFOX2_BrdU_CLIP_ANR$Score)
1 / median(MOTIF_CLINK_HUR_CoCLIP_ANR$Score)
1 / median(MOTIF_CLINK_HNRNPC_iCLIP_ANR$Score)
1 / median(MOTIF_CLINK_HNRNPC_BrdU_CLIP2_ANR$Score)
################################################################################

## CLINK TRUNC Test 2: ZOOPS, Scramble Off
################################################################################
extensions = c(5, -15)

MOTIF_CLINK_RBFOX2_BrdU_CLIP_ZOOPS = motifEnrichment(TRUNC_CLINK_RBFOX2_BrdU_CLIP[, .(chr, start, end, name, score, strand)], 'mm39', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HUR_CoCLIP_ZOOPS = motifEnrichment(TRUNC_CLINK_HuR_CoCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HNRNPC_iCLIP_ZOOPS = motifEnrichment(TRUNC_CLINK_HNRNPC_iCLIP[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)
MOTIF_CLINK_HNRNPC_BrdU_CLIP2_ZOOPS= motifEnrichment(TRUNC_CLINK_HNRNPC_BrdU_CLIP2[, .(chr, start, end, name, score, strand)], 'hg38', K = 5, nucleic_acid_type = 'RNA', method = 'zoops', scramble_bkg = F, extension = extensions)

1 / median(MOTIF_CLINK_RBFOX2_BrdU_CLIP_ZOOPS$Score)
1 / median(MOTIF_CLINK_HUR_CoCLIP_ZOOPS$Score)
1 / median(MOTIF_CLINK_HNRNPC_iCLIP_ZOOPS$Score)
1 / median(MOTIF_CLINK_HNRNPC_BrdU_CLIP2_ZOOPS$Score)
################################################################################




## 3. PWM (per condition):
################################################################################
returnPWM = function(data, sd_multipler, plot = TRUE) {
  data = data.frame(data[, c('MOTIF', 'Score')])
  data = data[order(-data[, 'Score']), ]
  rownames(data) = NULL
  
  data = data %>% filter(Score > (mean(data$Score) + sd_multipler*sd(data$Score)))
  
  if (nrow(data) == 0) {
    print("No motifs passed the filter.")
    return(NULL)
  }
  
  motif = DNAStringSet(gsub("U", "T", data$MOTIF))
  pfm_raw = consensusMatrix(motif)
  
  bases = c("A", "C", "G", "T")
  pfm = matrix(0, nrow = 4, ncol = ncol(pfm_raw))
  rownames(pfm) = bases
  colnames(pfm) = colnames(pfm_raw)
  
  common_rows = intersect(bases, rownames(pfm_raw))
  if (length(common_rows) > 0) {
    pfm[common_rows, ] = pfm_raw[common_rows, ]
  }
  
  ppm = prop.table(pfm, margin = 2)
  PWM = makePWM(ppm, alphabet = 'RNA')
  if (plot) seqLogo(PWM, ic.scale = F)
  
  ppm_log_ppm = ppm * log2(ppm)
  ppm_log_ppm[is.nan(ppm_log_ppm)] = 0
  H_l = -colSums(ppm_log_ppm)
  IC_l = 2 - H_l
  
  invisible(list(H_motif = sum(H_l), IC_motif = sum(IC_l), ppm = ppm, PWM = PWM, TopMotif = data$MOTIF[1]))
}

returnPWM(motifs_XL_CatStalling_TRUNC, 2)
returnPWM(motifs_XL_CatInactive_TRUNC, 2)
returnPWM(motifs_XL_CatStalling_DEL, 2)
returnPWM(motifs_XL_CatInactive_DEL, 2)
################################################################################
















