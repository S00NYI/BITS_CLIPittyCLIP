################################################################################
## CLIP Peak Filtering
## Author: Soon Yi with Antigravity
## Date: May 2026
################################################################################

# === Load Libraries ===
library(data.table)
library(tidyverse)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(eulerr)
library(ggExtra)
library(dplyr)
library(patchwork)
library(VennDiagram)
library(RColorBrewer)
library(GenomicRanges)

## 1. Basic Setup:
################################################################################
BASE_DIR = "/mnt/2TB_DATA/ANALYSIS_CLIPittyCLIP/"
INPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/1_ANNOTATED/")
OUTPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/2_FILTERED/")
################################################################################

## 1. Custom Functions
################################################################################
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[[annotation_column]]), row.names = 1)
  if (!is.null(new_column_name)) {
    colnames(temp) = new_column_name
  }
  
  if (!is.null(annotation_to_skip)) {
    temp = temp[rownames(temp) != annotation_to_skip, , drop = FALSE]
  }
  
  if (!is.null(fraction)) {
    temp = temp / sum(temp)
  }
  
  return(temp)
}
fillAnnotation = function(annotation_counts, annotation_list) {
  temp = data.frame(Sample = numeric(length(annotation_list)))
  rownames(temp) = annotation_list
  temp2 = merge(temp, annotation_counts, by = "row.names", all = TRUE)
  temp2[is.na(temp2)] = 0
  temp2$Sample = NULL
  rownames(temp2) = temp2$Row.names
  temp2 = temp2[annotation_list, -1, drop = FALSE]
  return(temp2)
}
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, color_map = NULL, counts = NULL, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list),
                aes(fill = Annotation, y = Freq, x = Source)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title  = element_text(hjust = 0.5),
      axis.text   = element_text(size = 14),
      axis.title  = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)
    )
  
  if (!is.null(color_map)) {
    plot = plot + scale_fill_manual(values = color_map)
  } else {
    plot = plot + scale_fill_brewer(palette = "Set3")
  }
  
  if (!is.null(counts)) {
    counts_df = data.frame(
      Source = names(counts),
      n      = unname(counts),
      stringsAsFactors = FALSE
    )
    counts_df$Source = factor(counts_df$Source, levels = sample_list)
    plot = plot + geom_text(
      data        = counts_df,
      mapping     = aes(x = Source, y = 1.02, label = n),
      inherit.aes = FALSE,
      size        = 4.5,
      vjust       = 0
    ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),
        limits = c(0, 1.12)
      )
  }
  
  if (!is.null(y_lim)) {
    plot = plot + ylim(y_lim)
  }
  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by = y_tick),
                                     limits = c(0, y_lim[2]))
  }
  
  return(plot)
}
################################################################################


## 2. Load annotated Peak/CITS/CIMS matrix:
################################################################################
HUR_CoCLIP_PEAK         = fread(paste0(INPUT_DIR, "HuR_CoCLIP_PEAK_ANNOTATED.txt"))
HUR_CoCLIP_CTK_TRUNC    = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CTK_TRUNC_ANNOTATED.txt"))
HUR_CoCLIP_CTK_TRUNC    = HUR_CoCLIP_CTK_TRUNC %>% mutate(pval = as.numeric(sub(".*\\[P=([^]]+)\\]$", "\\1", name)))
HUR_CoCLIP_CTK_DEL      = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CTK_DEL_ANNOTATED.txt"))
HUR_CoCLIP_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CLINK_TRUNC_ANNOTATED.txt"))
HUR_CoCLIP_CLINK_DEL    = fread(paste0(INPUT_DIR, "HuR_CoCLIP_CLINK_DEL_ANNOTATED.txt"))

HUR_eCLIP_K562_PEAK         = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_PEAK_ANNOTATED.txt"))
HUR_eCLIP_K562_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"))
HUR_eCLIP_K562_CLINK_DEL    = fread(paste0(INPUT_DIR, "HuR_eCLIP_K562_CLINK_DEL_ANNOTATED.txt"))

HNRNPC_BrdU_CLIP2_PEAK        = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_PEAK_ANNOTATED.txt"))
HNRNPC_BrdU_CLIP2_CTK_TRUNC   = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_TRUNC_ANNOTATED.txt"))
HNRNPC_BrdU_CLIP2_CTK_TRUNC   = HNRNPC_BrdU_CLIP2_CTK_TRUNC %>% mutate(pval = as.numeric(sub(".*\\[P=([^]]+)\\]$", "\\1", name)))
HNRNPC_BrdU_CLIP2_CTK_DEL     = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_DEL_ANNOTATED.txt"))
HNRNPC_BrdU_CLIP2_CLINK_TRUNC = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_TRUNC_ANNOTATED.txt"))
HNRNPC_BrdU_CLIP2_CLINK_DEL   = fread(paste0(INPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_DEL_ANNOTATED.txt"))

HNRNPC_iCLIP_PEAK         = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_PEAK_ANNOTATED.txt"))
HNRNPC_iCLIP_CTK_TRUNC    = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CTK_TRUNC_ANNOTATED.txt"))
HNRNPC_iCLIP_CTK_TRUNC    = HNRNPC_iCLIP_CTK_TRUNC %>% mutate(pval = as.numeric(sub(".*\\[P=([^]]+)\\]$", "\\1", name)))
HNRNPC_iCLIP_CTK_DEL      = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CTK_DEL_ANNOTATED.txt"))
HNRNPC_iCLIP_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CLINK_TRUNC_ANNOTATED.txt"))
HNRNPC_iCLIP_CLINK_DEL    = fread(paste0(INPUT_DIR, "HNRNPC_iCLIP_CLINK_DEL_ANNOTATED.txt"))

HNRNPC_eCLIP_HepG2_PEAK         = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_PEAK_ANNOTATED.txt"))
HNRNPC_eCLIP_HepG2_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED.txt"))
HNRNPC_eCLIP_HepG2_CLINK_DEL    = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_DEL_ANNOTATED.txt"))

HNRNPC_eCLIP_K562_PEAK         = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_PEAK_ANNOTATED.txt"))
HNRNPC_eCLIP_K562_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"))
HNRNPC_eCLIP_K562_CLINK_DEL    = fread(paste0(INPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_DEL_ANNOTATED.txt"))

RBFOX2_BrdU_CLIP_PEAK         = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_PEAK_ANNOTATED.txt"))
RBFOX2_BrdU_CLIP_CTK_TRUNC    = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_TRUNC_ANNOTATED.txt"))
RBFOX2_BrdU_CLIP_CTK_TRUNC    = RBFOX2_BrdU_CLIP_CTK_TRUNC %>% mutate(pval = as.numeric(sub(".*\\[P=([^]]+)\\]$", "\\1", name)))
RBFOX2_BrdU_CLIP_CTK_DEL      = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_DEL_ANNOTATED.txt"))
RBFOX2_BrdU_CLIP_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_TRUNC_ANNOTATED.txt"))
RBFOX2_BrdU_CLIP_CLINK_DEL    = fread(paste0(INPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_DEL_ANNOTATED.txt"))

RBFOX2_eCLIP_HepG2_PEAK         = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_PEAK_ANNOTATED.txt"))
RBFOX2_eCLIP_HepG2_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED.txt"))
RBFOX2_eCLIP_HepG2_CLINK_DEL    = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_DEL_ANNOTATED.txt"))

RBFOX2_eCLIP_K562_PEAK         = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_PEAK_ANNOTATED.txt"))
RBFOX2_eCLIP_K562_CLINK_TRUNC  = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"))
RBFOX2_eCLIP_K562_CLINK_DEL    = fread(paste0(INPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_DEL_ANNOTATED.txt"))

################################################################################


## 3. Filter annotated Peak/CITS/CIMS matrix:
################################################################################
CTK_sig = 0.001
CLINK_sig = 0.01

HUR_CoCLIP_PEAK_FILTERED             = HUR_CoCLIP_PEAK %>% filter(NormedTC_HuR_WT >= quantile(NormedTC_HuR_WT)[3])
HUR_CoCLIP_CTK_TRUNC_FILTERED        = HUR_CoCLIP_CTK_TRUNC %>% filter(score > median(score))  %>% filter(pval <= CTK_sig)
HUR_CoCLIP_CTK_DEL_FILTERED          = HUR_CoCLIP_CTK_DEL %>% filter(score > median(score))  %>% filter(FDR <= CTK_sig)
HUR_CoCLIP_CLINK_TRUNC_FILTERED      = HUR_CoCLIP_CLINK_TRUNC %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)
HUR_CoCLIP_CLINK_DEL_FILTERED        = HUR_CoCLIP_CLINK_DEL %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)

HUR_eCLIP_K562_PEAK_FILTERED         = HUR_eCLIP_K562_PEAK %>% filter(NormedTC_K562_ELAVL1 >= mean(NormedTC_K562_ELAVL1)) %>% filter(BC_K562_ELAVL1 == 2) 
HUR_eCLIP_K562_CLINK_TRUNC_FILTERED  = HUR_eCLIP_K562_CLINK_TRUNC %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)
HUR_eCLIP_K562_CLINK_DEL_FILTERED    = HUR_eCLIP_K562_CLINK_DEL %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)

HNRNPC_BrdU_CLIP2_PEAK_FILTERED          = HNRNPC_BrdU_CLIP2_PEAK %>% filter(NormedTC_HNRNPC_WT >= quantile(NormedTC_HNRNPC_WT)[3])
HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED     = HNRNPC_BrdU_CLIP2_CTK_TRUNC %>% filter(score > median(score))  %>% filter(pval <= CTK_sig)
HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED       = HNRNPC_BrdU_CLIP2_CTK_DEL %>% filter(score > median(score))  %>% filter(FDR <= CTK_sig)
HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED   = HNRNPC_BrdU_CLIP2_CLINK_TRUNC %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)
HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED     = HNRNPC_BrdU_CLIP2_CLINK_DEL %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)

HNRNPC_iCLIP_PEAK_FILTERED               = HNRNPC_iCLIP_PEAK %>% filter(NormedTC_HNRNPC_WT >= quantile(NormedTC_HNRNPC_WT)[3])
HNRNPC_iCLIP_CTK_TRUNC_FILTERED          = HNRNPC_iCLIP_CTK_TRUNC %>% filter(pval <= CTK_sig)
HNRNPC_iCLIP_CTK_DEL_FILTERED            = HNRNPC_iCLIP_CTK_DEL %>% filter(FDR <= CTK_sig)
HNRNPC_iCLIP_CLINK_TRUNC_FILTERED        = HNRNPC_iCLIP_CLINK_TRUNC %>% filter(qvalue <= CLINK_sig)
HNRNPC_iCLIP_CLINK_DEL_FILTERED          = HNRNPC_iCLIP_CLINK_DEL %>% filter(qvalue <= CLINK_sig)

HNRNPC_eCLIP_HepG2_PEAK_FILTERED         = HNRNPC_eCLIP_HepG2_PEAK %>% filter(NormedTC_HepG2_HNRNPC >= mean(NormedTC_HepG2_HNRNPC)) %>% filter(BC_HepG2_HNRNPC == 2) 
HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED  = HNRNPC_eCLIP_HepG2_CLINK_TRUNC %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)
HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED    = HNRNPC_eCLIP_HepG2_CLINK_DEL %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)

HNRNPC_eCLIP_K562_PEAK_FILTERED          = HNRNPC_eCLIP_K562_PEAK %>% filter(NormedTC_K562_HNRNPC >= mean(NormedTC_K562_HNRNPC)) %>% filter(BC_K562_HNRNPC == 2)
HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED   = HNRNPC_eCLIP_K562_CLINK_TRUNC %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)
HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED     = HNRNPC_eCLIP_K562_CLINK_DEL %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)

RBFOX2_BrdU_CLIP_PEAK_FILTERED           = RBFOX2_BrdU_CLIP_PEAK %>% filter(NormedTC_RBFOX2 >= quantile(NormedTC_RBFOX2)[3])
RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED      = RBFOX2_BrdU_CLIP_CTK_TRUNC %>% filter(score > median(score))  %>% filter(pval <= CTK_sig)
RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED        = RBFOX2_BrdU_CLIP_CTK_DEL %>% filter(score > median(score))  %>% filter(FDR <= CTK_sig)
RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED    = RBFOX2_BrdU_CLIP_CLINK_TRUNC %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)
RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED      = RBFOX2_BrdU_CLIP_CLINK_DEL %>% filter(score > median(score))  %>% filter(qvalue <= CLINK_sig)

RBFOX2_eCLIP_HepG2_PEAK_FILTERED         = RBFOX2_eCLIP_HepG2_PEAK %>% filter(NormedTC_HepG2_RBFOX2 >= mean(NormedTC_HepG2_RBFOX2)) %>% filter(BC_HepG2_RBFOX2 == 2)
RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED  = RBFOX2_eCLIP_HepG2_CLINK_TRUNC %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)
RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED    = RBFOX2_eCLIP_HepG2_CLINK_DEL %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)

RBFOX2_eCLIP_K562_PEAK_FILTERED          = RBFOX2_eCLIP_K562_PEAK %>% filter(NormedTC_K562_RBFOX2 >= mean(NormedTC_K562_RBFOX2)) %>% filter(BC_K562_RBFOX2 == 2)
RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED   = RBFOX2_eCLIP_K562_CLINK_TRUNC %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)
RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED     = RBFOX2_eCLIP_K562_CLINK_DEL %>% filter(score > mean(score) + 2*sd(score))  %>% filter(qvalue <= CLINK_sig)


fwrite(HUR_CoCLIP_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HuR_CoCLIP_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HUR_CoCLIP_CTK_TRUNC_FILTERED,   paste0(OUTPUT_DIR, "HuR_CoCLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"),   sep = "\t")
fwrite(HUR_CoCLIP_CTK_DEL_FILTERED,     paste0(OUTPUT_DIR, "HuR_CoCLIP_CTK_DEL_ANNOTATED_FILTERED.txt"),     sep = "\t")
fwrite(HUR_CoCLIP_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HuR_CoCLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HUR_CoCLIP_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HuR_CoCLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(HUR_eCLIP_K562_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HuR_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HUR_eCLIP_K562_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HuR_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HUR_eCLIP_K562_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HuR_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(HNRNPC_BrdU_CLIP2_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_TRUNC_ANNOTATED_FILTERED.txt"),   sep = "\t")
fwrite(HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED,     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_DEL_ANNOTATED_FILTERED.txt"),     sep = "\t")
fwrite(HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(HNRNPC_iCLIP_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HNRNPC_iCLIP_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HNRNPC_iCLIP_CTK_TRUNC_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"),   sep = "\t")
fwrite(HNRNPC_iCLIP_CTK_DEL_FILTERED,     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CTK_DEL_ANNOTATED_FILTERED.txt"),     sep = "\t")
fwrite(HNRNPC_iCLIP_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HNRNPC_iCLIP_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(HNRNPC_eCLIP_HepG2_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(HNRNPC_eCLIP_K562_PEAK_FILTERED,        paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(RBFOX2_BrdU_CLIP_PEAK_FILTERED,        paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED,   paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_TRUNC_ANNOTATED_FILTERED.txt"),   sep = "\t")
fwrite(RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED,     paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_DEL_ANNOTATED_FILTERED.txt"),     sep = "\t")
fwrite(RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(RBFOX2_eCLIP_HepG2_PEAK_FILTERED,        paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

fwrite(RBFOX2_eCLIP_K562_PEAK_FILTERED,        paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_PEAK_ANNOTATED_FILTERED.txt"),        sep = "\t")
fwrite(RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED, paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_TRUNC_ANNOTATED_FILTERED.txt"), sep = "\t")
fwrite(RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED,   paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_DEL_ANNOTATED_FILTERED.txt"),   sep = "\t")

################################################################################

## 4_0 Annotation list
################################################################################
Genic_Annotation_List = c(
  "3'UTR", "3'UTR-distal", "5'UTR", "CDS",
  "First intron", "Proximal intronic", "Deep intronic",
  "Retained intron:coding", "Other exon"
)

NonGenic_Annotation_List = c(
  "lncRNA", "miRNA", "YRNA",
  "rRNA*", "snRNA*", "snoRNA*", "tRNA*",
  "Alu", "Other SINE", "LINE", "LTR/ERV", "Satellite/LC",
  "Intergenic"
)

color_map = c(
  ## Genic
  "3'UTR"                  = "#d9a528",
  "3'UTR-distal"           = "#f6d062",
  "5'UTR"                  = "#f57f20",
  "CDS"                    = "#e21f26",
  "First intron"           = "#076936",
  "Proximal intronic"      = "#2da260",
  "Deep intronic"          = "#76c376",
  "Retained intron:coding" = "#c8e4bf",
  "Other exon"             = "#bcbec0",
  ## Non-genic
  "lncRNA"                 = "#6b3f98",
  "miRNA"                  = "#b25a29",
  "YRNA"                   = "#f57f20",
  "rRNA*"                  = "#0b539d",
  "snRNA*"                 = "#3183be",
  "snoRNA*"                = "#6baed6",
  "tRNA*"                  = "#bfd7e7",
  "Alu"                    = "#c61d8a",
  "Other SINE"             = "#f069a0",
  "LINE"                   = "#f59eb5",
  "LTR/ERV"                = "#fce0ec",
  "Satellite/LC"           = "#e6e7e8",
  "Intergenic"             = "#6d6e71"
)
################################################################################

## 4_1 HuR CoCLIP Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HUR_CoCLIP_PEAK        = countAnnotation(HUR_CoCLIP_PEAK_FILTERED, "annotation", "HUR_CoCLIP_PEAK")
DIST_HUR_CoCLIP_CTK_TRUNC   = countAnnotation(HUR_CoCLIP_CTK_TRUNC_FILTERED, "annotation", "HUR_CoCLIP_CTK_TRUNC")
DIST_HUR_CoCLIP_CTK_DEL     = countAnnotation(HUR_CoCLIP_CTK_DEL_FILTERED, "annotation", "HUR_CoCLIP_CTK_DEL")
DIST_HUR_CoCLIP_CLINK_TRUNC = countAnnotation(HUR_CoCLIP_CLINK_TRUNC_FILTERED, "annotation", "HUR_CoCLIP_CLINK_TRUNC")
DIST_HUR_CoCLIP_CLINK_DEL   = countAnnotation(HUR_CoCLIP_CLINK_DEL_FILTERED, "annotation", "HUR_CoCLIP_CLINK_DEL")

DIST_HUR_CoCLIP_PEAK        = fillAnnotation(DIST_HUR_CoCLIP_PEAK, rownames(DIST_HUR_CoCLIP_CTK_TRUNC))
DIST_HUR_CoCLIP_CTK_DEL     = fillAnnotation(DIST_HUR_CoCLIP_CTK_DEL, rownames(DIST_HUR_CoCLIP_CTK_TRUNC))
DIST_HUR_CoCLIP_CLINK_TRUNC = fillAnnotation(DIST_HUR_CoCLIP_CLINK_TRUNC, rownames(DIST_HUR_CoCLIP_CTK_TRUNC))
DIST_HUR_CoCLIP_CLINK_DEL   = fillAnnotation(DIST_HUR_CoCLIP_CLINK_DEL, rownames(DIST_HUR_CoCLIP_CTK_TRUNC))

CLIP_List = c("HUR_CoCLIP_PEAK", "HUR_CoCLIP_CTK_TRUNC", "HUR_CoCLIP_CLINK_TRUNC", "HUR_CoCLIP_CTK_DEL", "HUR_CoCLIP_CLINK_DEL")

DIST_HUR_CoCLIP_COMBINED = cbind(
  DIST_HUR_CoCLIP_PEAK,
  DIST_HUR_CoCLIP_CTK_TRUNC,
  DIST_HUR_CoCLIP_CLINK_TRUNC,
  DIST_HUR_CoCLIP_CTK_DEL,
  DIST_HUR_CoCLIP_CLINK_DEL
)
DIST_HUR_CoCLIP_COMBINED$Annotation = rownames(DIST_HUR_CoCLIP_COMBINED)

DIST_HUR_CoCLIP_COMBINED = DIST_HUR_CoCLIP_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HUR_CoCLIP_COMBINED$Source = factor(DIST_HUR_CoCLIP_COMBINED$Source, levels = CLIP_List)

DIST_HUR_CoCLIP_COMBINED_GENIC = DIST_HUR_CoCLIP_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HUR_CoCLIP_COMBINED_GENIC = tapply(DIST_HUR_CoCLIP_COMBINED_GENIC$Freq, DIST_HUR_CoCLIP_COMBINED_GENIC$Source, sum)

DIST_HUR_CoCLIP_COMBINED_NONGENIC = DIST_HUR_CoCLIP_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_DIST_HUR_CoCLIP_COMBINED_NONGENIC = tapply(DIST_HUR_CoCLIP_COMBINED_NONGENIC$Freq, DIST_HUR_CoCLIP_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HUR_CoCLIP_COMBINED_GENIC = plotStackedBar(
  DIST_HUR_CoCLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HuR CoCLIP Input Genic Annotation Distributions")

PLOT_DIST_HUR_CoCLIP_COMBINED_NONGENIC = plotStackedBar(
  DIST_HUR_CoCLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HuR CoCLIP Input Non-Genic Annotation Distributions")

# print(PLOT_DIST_HUR_CoCLIP_COMBINED_GENIC)
# print(PLOT_DIST_HUR_CoCLIP_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HUR_CoCLIP_COMBINED_GENIC = DIST_HUR_CoCLIP_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HUR_CoCLIP_COMBINED_NONGENIC = DIST_HUR_CoCLIP_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HUR_CoCLIP_COMBINED_GENIC = plotStackedBar(
  PROP_HUR_CoCLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HUR_CoCLIP_COMBINED_GENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HuR CoCLIP Input Genic Annotation Proportions")

PLOT_PROP_HUR_CoCLIP_COMBINED_NONGENIC = plotStackedBar(
  PROP_HUR_CoCLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_DIST_HUR_CoCLIP_COMBINED_NONGENIC,  
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HuR CoCLIP Input Non-Genic Annotation Proportions")

# print(PLOT_PROP_HUR_CoCLIP_COMBINED_GENIC)
# print(PLOT_PROP_HUR_CoCLIP_COMBINED_NONGENIC)
################################################################################

## 4_2 HuR eCLIP K562 Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HUR_eCLIP_K562_PEAK          = countAnnotation(HUR_eCLIP_K562_PEAK_FILTERED, "annotation", "HUR_eCLIP_K562_PEAK")
DIST_HUR_eCLIP_K562_CLINK_TRUNC   = countAnnotation(HUR_eCLIP_K562_CLINK_TRUNC_FILTERED, "annotation", "HUR_eCLIP_K562_CLINK_TRUNC")
DIST_HUR_eCLIP_K562_CLINK_DEL     = countAnnotation(HUR_eCLIP_K562_CLINK_DEL_FILTERED, "annotation", "HUR_eCLIP_K562_CLINK_DEL")

CLIP_List = c("HUR_eCLIP_K562_PEAK", "HUR_eCLIP_K562_CLINK_TRUNC", "HUR_eCLIP_K562_CLINK_DEL")

DIST_HUR_eCLIP_K562_COMBINED = cbind(
  DIST_HUR_eCLIP_K562_PEAK,
  DIST_HUR_eCLIP_K562_CLINK_TRUNC,
  DIST_HUR_eCLIP_K562_CLINK_DEL
)
DIST_HUR_eCLIP_K562_COMBINED$Annotation = rownames(DIST_HUR_eCLIP_K562_COMBINED)

DIST_HUR_eCLIP_K562_COMBINED = DIST_HUR_eCLIP_K562_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HUR_eCLIP_K562_COMBINED$Source = factor(DIST_HUR_eCLIP_K562_COMBINED$Source, levels = CLIP_List)

DIST_HUR_eCLIP_K562_COMBINED_GENIC = DIST_HUR_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HUR_eCLIP_K562_COMBINED_GENIC = tapply(DIST_HUR_eCLIP_K562_COMBINED_GENIC$Freq, DIST_HUR_eCLIP_K562_COMBINED_GENIC$Source, sum)

DIST_HUR_eCLIP_K562_COMBINED_NONGENIC = DIST_HUR_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_HUR_eCLIP_K562_COMBINED_NONGENIC = tapply(DIST_HUR_eCLIP_K562_COMBINED_NONGENIC$Freq, DIST_HUR_eCLIP_K562_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HUR_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  DIST_HUR_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HuR eCLIP K562 Genic Annotation Distributions")

PLOT_DIST_HUR_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  DIST_HUR_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HuR eCLIP K562 Non-Genic Annotation Distributions")

# print(PLOT_DIST_HUR_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_DIST_HUR_eCLIP_K562_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HUR_eCLIP_K562_COMBINED_GENIC = DIST_HUR_eCLIP_K562_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HUR_eCLIP_K562_COMBINED_NONGENIC = DIST_HUR_eCLIP_K562_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HUR_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  PROP_HUR_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HUR_eCLIP_K562_COMBINED_GENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HuR eCLIP K562 Genic Annotation Proportions")

PLOT_PROP_HUR_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  PROP_HUR_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HUR_eCLIP_K562_COMBINED_NONGENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HuR eCLIP K562 Non-Genic Annotation Proportions")

# print(PLOT_PROP_HUR_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_PROP_HUR_eCLIP_K562_COMBINED_NONGENIC)
################################################################################

## 4_3 HNRNPC Specificity CLIP Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HNRNPC_BrdU_CLIP2_PEAK        = countAnnotation(HNRNPC_BrdU_CLIP2_PEAK_FILTERED, "annotation", "HNRNPC_BrdU_CLIP2_PEAK")
DIST_HNRNPC_BrdU_CLIP2_CTK_TRUNC   = countAnnotation(HNRNPC_BrdU_CLIP2_CTK_TRUNC_FILTERED, "annotation", "HNRNPC_BrdU_CLIP2_CTK_TRUNC")
DIST_HNRNPC_BrdU_CLIP2_CTK_DEL     = countAnnotation(HNRNPC_BrdU_CLIP2_CTK_DEL_FILTERED, "annotation", "HNRNPC_BrdU_CLIP2_CTK_DEL")
DIST_HNRNPC_BrdU_CLIP2_CLINK_TRUNC = countAnnotation(HNRNPC_BrdU_CLIP2_CLINK_TRUNC_FILTERED, "annotation", "HNRNPC_BrdU_CLIP2_CLINK_TRUNC")
DIST_HNRNPC_BrdU_CLIP2_CLINK_DEL   = countAnnotation(HNRNPC_BrdU_CLIP2_CLINK_DEL_FILTERED, "annotation", "HNRNPC_BrdU_CLIP2_CLINK_DEL")

DIST_HNRNPC_BrdU_CLIP2_CLINK_DEL        = fillAnnotation(DIST_HNRNPC_BrdU_CLIP2_CLINK_DEL, rownames(DIST_HNRNPC_BrdU_CLIP2_PEAK))

CLIP_List = c("HNRNPC_BrdU_CLIP2_PEAK", "HNRNPC_BrdU_CLIP2_CTK_TRUNC", "HNRNPC_BrdU_CLIP2_CLINK_TRUNC", "HNRNPC_BrdU_CLIP2_CTK_DEL", "HNRNPC_BrdU_CLIP2_CLINK_DEL")

DIST_HNRNPC_BrdU_CLIP2_COMBINED = cbind(
  DIST_HNRNPC_BrdU_CLIP2_PEAK,
  DIST_HNRNPC_BrdU_CLIP2_CTK_TRUNC,
  DIST_HNRNPC_BrdU_CLIP2_CLINK_TRUNC,
  DIST_HNRNPC_BrdU_CLIP2_CTK_DEL,
  DIST_HNRNPC_BrdU_CLIP2_CLINK_DEL
)
DIST_HNRNPC_BrdU_CLIP2_COMBINED$Annotation = rownames(DIST_HNRNPC_BrdU_CLIP2_COMBINED)

DIST_HNRNPC_BrdU_CLIP2_COMBINED = DIST_HNRNPC_BrdU_CLIP2_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HNRNPC_BrdU_CLIP2_COMBINED$Source = factor(DIST_HNRNPC_BrdU_CLIP2_COMBINED$Source, levels = CLIP_List)

DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC = DIST_HNRNPC_BrdU_CLIP2_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HNRNPC_BrdU_CLIP2_COMBINED_GENIC = tapply(DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC$Freq, DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC$Source, sum)

DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC = DIST_HNRNPC_BrdU_CLIP2_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC = tapply(DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC$Freq, DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC = plotStackedBar(
  DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC Specificity CLIP Genic Annotation Distributions")

PLOT_DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC = plotStackedBar(
  DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC Specificity CLIP Non-Genic Annotation Distributions")

# print(PLOT_DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC)
# print(PLOT_DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HNRNPC_BrdU_CLIP2_COMBINED_GENIC = DIST_HNRNPC_BrdU_CLIP2_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC = DIST_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_GENIC = plotStackedBar(
  PROP_HNRNPC_BrdU_CLIP2_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_BrdU_CLIP2_COMBINED_GENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC Specificity CLIP Genic Annotation Proportions")

PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC = plotStackedBar(
  PROP_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC Specificity CLIP Non-Genic Annotation Proportions")

# print(PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_GENIC)
# print(PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC)
################################################################################

## 4_4 HNRNPC iCLIP Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HNRNPC_iCLIP_PEAK        = countAnnotation(HNRNPC_iCLIP_PEAK_FILTERED, "annotation", "HNRNPC_iCLIP_PEAK")
DIST_HNRNPC_iCLIP_CTK_TRUNC   = countAnnotation(HNRNPC_iCLIP_CTK_TRUNC_FILTERED, "annotation", "HNRNPC_iCLIP_CTK_TRUNC")
DIST_HNRNPC_iCLIP_CTK_DEL     = countAnnotation(HNRNPC_iCLIP_CTK_DEL_FILTERED, "annotation", "HNRNPC_iCLIP_CTK_DEL")
DIST_HNRNPC_iCLIP_CLINK_TRUNC = countAnnotation(HNRNPC_iCLIP_CLINK_TRUNC_FILTERED, "annotation", "HNRNPC_iCLIP_CLINK_TRUNC")
DIST_HNRNPC_iCLIP_CLINK_DEL   = countAnnotation(HNRNPC_iCLIP_CLINK_DEL_FILTERED, "annotation", "HNRNPC_iCLIP_CLINK_DEL")

DIST_HNRNPC_iCLIP_CTK_TRUNC   = fillAnnotation(DIST_HNRNPC_iCLIP_CTK_TRUNC, rownames(DIST_HNRNPC_iCLIP_PEAK))
DIST_HNRNPC_iCLIP_CTK_DEL     = fillAnnotation(DIST_HNRNPC_iCLIP_CTK_DEL, rownames(DIST_HNRNPC_iCLIP_PEAK))
DIST_HNRNPC_iCLIP_CLINK_TRUNC = fillAnnotation(DIST_HNRNPC_iCLIP_CLINK_TRUNC, rownames(DIST_HNRNPC_iCLIP_PEAK))
DIST_HNRNPC_iCLIP_CLINK_DEL   = fillAnnotation(DIST_HNRNPC_iCLIP_CLINK_DEL, rownames(DIST_HNRNPC_iCLIP_PEAK))

CLIP_List = c("HNRNPC_iCLIP_PEAK", "HNRNPC_iCLIP_CTK_TRUNC", "HNRNPC_iCLIP_CLINK_TRUNC", "HNRNPC_iCLIP_CTK_DEL", "HNRNPC_iCLIP_CLINK_DEL")

DIST_HNRNPC_iCLIP_COMBINED = cbind(
  DIST_HNRNPC_iCLIP_PEAK,
  DIST_HNRNPC_iCLIP_CTK_TRUNC,
  DIST_HNRNPC_iCLIP_CLINK_TRUNC,
  DIST_HNRNPC_iCLIP_CTK_DEL,
  DIST_HNRNPC_iCLIP_CLINK_DEL
)
DIST_HNRNPC_iCLIP_COMBINED$Annotation = rownames(DIST_HNRNPC_iCLIP_COMBINED)

DIST_HNRNPC_iCLIP_COMBINED = DIST_HNRNPC_iCLIP_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HNRNPC_iCLIP_COMBINED$Source = factor(DIST_HNRNPC_iCLIP_COMBINED$Source, levels = CLIP_List)

DIST_HNRNPC_iCLIP_COMBINED_GENIC = DIST_HNRNPC_iCLIP_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HNRNPC_iCLIP_COMBINED_GENIC = tapply(DIST_HNRNPC_iCLIP_COMBINED_GENIC$Freq, DIST_HNRNPC_iCLIP_COMBINED_GENIC$Source, sum)

DIST_HNRNPC_iCLIP_COMBINED_NONGENIC = DIST_HNRNPC_iCLIP_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_HNRNPC_iCLIP_COMBINED_NONGENIC = tapply(DIST_HNRNPC_iCLIP_COMBINED_NONGENIC$Freq, DIST_HNRNPC_iCLIP_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HNRNPC_iCLIP_COMBINED_GENIC = plotStackedBar(
  DIST_HNRNPC_iCLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC iCLIP Genic Annotation Distributions")

PLOT_DIST_HNRNPC_iCLIP_COMBINED_NONGENIC = plotStackedBar(
  DIST_HNRNPC_iCLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC iCLIP Non-Genic Annotation Distributions")

# print(PLOT_DIST_HNRNPC_iCLIP_COMBINED_GENIC)
# print(PLOT_DIST_HNRNPC_iCLIP_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HNRNPC_iCLIP_COMBINED_GENIC = DIST_HNRNPC_iCLIP_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HNRNPC_iCLIP_COMBINED_NONGENIC = DIST_HNRNPC_iCLIP_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HNRNPC_iCLIP_COMBINED_GENIC = plotStackedBar(
  PROP_HNRNPC_iCLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_iCLIP_COMBINED_GENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC iCLIP Genic Annotation Proportions")

PLOT_PROP_HNRNPC_iCLIP_COMBINED_NONGENIC = plotStackedBar(
  PROP_HNRNPC_iCLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_iCLIP_COMBINED_NONGENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "HNRNPC iCLIP Non-Genic Annotation Proportions")

# print(PLOT_PROP_HNRNPC_iCLIP_COMBINED_GENIC)
# print(PLOT_PROP_HNRNPC_iCLIP_COMBINED_NONGENIC)
################################################################################

## 4_5 HNRNPC eCLIP HepG2 Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HNRNPC_eCLIP_HepG2_PEAK          = countAnnotation(HNRNPC_eCLIP_HepG2_PEAK_FILTERED, "annotation", "HNRNPC_eCLIP_HepG2_PEAK")
DIST_HNRNPC_eCLIP_HepG2_CLINK_TRUNC   = countAnnotation(HNRNPC_eCLIP_HepG2_CLINK_TRUNC_FILTERED, "annotation", "HNRNPC_eCLIP_HepG2_CLINK_TRUNC")
DIST_HNRNPC_eCLIP_HepG2_CLINK_DEL     = countAnnotation(HNRNPC_eCLIP_HepG2_CLINK_DEL_FILTERED, "annotation", "HNRNPC_eCLIP_HepG2_CLINK_DEL")

DIST_HNRNPC_eCLIP_HepG2_CLINK_TRUNC      = fillAnnotation(DIST_HNRNPC_eCLIP_HepG2_CLINK_TRUNC, rownames(DIST_HNRNPC_eCLIP_HepG2_PEAK))
DIST_HNRNPC_eCLIP_HepG2_CLINK_DEL = fillAnnotation(DIST_HNRNPC_eCLIP_HepG2_CLINK_DEL, rownames(DIST_HNRNPC_eCLIP_HepG2_PEAK))

CLIP_List = c("HNRNPC_eCLIP_HepG2_PEAK", "HNRNPC_eCLIP_HepG2_CLINK_TRUNC", "HNRNPC_eCLIP_HepG2_CLINK_DEL")

DIST_HNRNPC_eCLIP_HepG2_COMBINED = cbind(
  DIST_HNRNPC_eCLIP_HepG2_PEAK,
  DIST_HNRNPC_eCLIP_HepG2_CLINK_TRUNC,
  DIST_HNRNPC_eCLIP_HepG2_CLINK_DEL
)
DIST_HNRNPC_eCLIP_HepG2_COMBINED$Annotation = rownames(DIST_HNRNPC_eCLIP_HepG2_COMBINED)

DIST_HNRNPC_eCLIP_HepG2_COMBINED = DIST_HNRNPC_eCLIP_HepG2_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HNRNPC_eCLIP_HepG2_COMBINED$Source = factor(DIST_HNRNPC_eCLIP_HepG2_COMBINED$Source, levels = CLIP_List)

DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC = DIST_HNRNPC_eCLIP_HepG2_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HNRNPC_eCLIP_HepG2_COMBINED_GENIC = tapply(DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC$Freq, DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC$Source, sum)

DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC = DIST_HNRNPC_eCLIP_HepG2_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC = tapply(DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC$Freq, DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC = plotStackedBar(
  DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP HepG2 Genic Annotation Distributions")

PLOT_DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC = plotStackedBar(
  DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP HepG2 Non-Genic Annotation Distributions")

# print(PLOT_DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC)
# print(PLOT_DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HNRNPC_eCLIP_HepG2_COMBINED_GENIC= DIST_HNRNPC_eCLIP_HepG2_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC = DIST_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_GENIC = plotStackedBar(
  PROP_HNRNPC_eCLIP_HepG2_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"), counts = COUNTS_HNRNPC_eCLIP_HepG2_COMBINED_GENIC,
  "HNRNPC eCLIP HepG2 Genic Annotation Proportions")

PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC = plotStackedBar(
  PROP_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,  counts = COUNTS_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP HepG2 Non-Genic Annotation Proportions")

# print(PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_GENIC)
# print(PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC)
################################################################################

## 4_6 HNRNPC eCLIP K562 Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_HNRNPC_eCLIP_K562_PEAK          = countAnnotation(HNRNPC_eCLIP_K562_PEAK_FILTERED, "annotation", "HNRNPC_eCLIP_K562_PEAK")
DIST_HNRNPC_eCLIP_K562_CLINK_TRUNC   = countAnnotation(HNRNPC_eCLIP_K562_CLINK_TRUNC_FILTERED, "annotation", "HNRNPC_eCLIP_K562_CLINK_TRUNC")
DIST_HNRNPC_eCLIP_K562_CLINK_DEL     = countAnnotation(HNRNPC_eCLIP_K562_CLINK_DEL_FILTERED, "annotation", "HNRNPC_eCLIP_K562_CLINK_DEL")

DIST_HNRNPC_eCLIP_K562_CLINK_DEL = fillAnnotation(DIST_HNRNPC_eCLIP_K562_CLINK_DEL, rownames(DIST_HNRNPC_eCLIP_K562_CLINK_TRUNC))

CLIP_List = c("HNRNPC_eCLIP_K562_PEAK", "HNRNPC_eCLIP_K562_CLINK_TRUNC", "HNRNPC_eCLIP_K562_CLINK_DEL")

DIST_HNRNPC_eCLIP_K562_COMBINED = cbind(
  DIST_HNRNPC_eCLIP_K562_PEAK,
  DIST_HNRNPC_eCLIP_K562_CLINK_TRUNC,
  DIST_HNRNPC_eCLIP_K562_CLINK_DEL
)
DIST_HNRNPC_eCLIP_K562_COMBINED$Annotation = rownames(DIST_HNRNPC_eCLIP_K562_COMBINED)

DIST_HNRNPC_eCLIP_K562_COMBINED = DIST_HNRNPC_eCLIP_K562_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_HNRNPC_eCLIP_K562_COMBINED$Source = factor(DIST_HNRNPC_eCLIP_K562_COMBINED$Source, levels = CLIP_List)

DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC = DIST_HNRNPC_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_HNRNPC_eCLIP_K562_COMBINED_GENIC = tapply(DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC$Freq, DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC$Source, sum)

DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC = DIST_HNRNPC_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_HNRNPC_eCLIP_K562_COMBINED_NONGENIC = tapply(DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC$Freq, DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map, 
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP K562 Genic Annotation Distributions")

PLOT_DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP K562 Non-Genic Annotation Distributions")

# print(PLOT_DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_HNRNPC_eCLIP_K562_COMBINED_GENIC= DIST_HNRNPC_eCLIP_K562_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_HNRNPC_eCLIP_K562_COMBINED_NONGENIC = DIST_HNRNPC_eCLIP_K562_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  PROP_HNRNPC_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_eCLIP_K562_COMBINED_GENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP K562 Genic Annotation Proportions")

PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  PROP_HNRNPC_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_HNRNPC_eCLIP_K562_COMBINED_NONGENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "HNRNPC eCLIP K562 Non-Genic Annotation Proportions")

# print(PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_NONGENIC)
################################################################################

## 4_7 RBFOX2 BrdU-CLIP Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_RBFOX2_BrdU_CLIP_PEAK        = countAnnotation(RBFOX2_BrdU_CLIP_PEAK_FILTERED, "annotation", "RBFOX2_BrdU_CLIP_PEAK")
DIST_RBFOX2_BrdU_CLIP_CTK_TRUNC   = countAnnotation(RBFOX2_BrdU_CLIP_CTK_TRUNC_FILTERED, "annotation", "RBFOX2_BrdU_CLIP_CTK_TRUNC")
DIST_RBFOX2_BrdU_CLIP_CTK_DEL     = countAnnotation(RBFOX2_BrdU_CLIP_CTK_DEL_FILTERED, "annotation", "RBFOX2_BrdU_CLIP_CTK_DEL")
DIST_RBFOX2_BrdU_CLIP_CLINK_TRUNC = countAnnotation(RBFOX2_BrdU_CLIP_CLINK_TRUNC_FILTERED, "annotation", "RBFOX2_BrdU_CLIP_CLINK_TRUNC")
DIST_RBFOX2_BrdU_CLIP_CLINK_DEL   = countAnnotation(RBFOX2_BrdU_CLIP_CLINK_DEL_FILTERED, "annotation", "RBFOX2_BrdU_CLIP_CLINK_DEL")

DIST_RBFOX2_BrdU_CLIP_CTK_DEL     = fillAnnotation(DIST_RBFOX2_BrdU_CLIP_CTK_DEL, rownames(DIST_RBFOX2_BrdU_CLIP_PEAK))
DIST_RBFOX2_BrdU_CLIP_CLINK_DEL   = fillAnnotation(DIST_RBFOX2_BrdU_CLIP_CLINK_DEL, rownames(DIST_RBFOX2_BrdU_CLIP_PEAK))

CLIP_List = c("RBFOX2_BrdU_CLIP_PEAK", "RBFOX2_BrdU_CLIP_CTK_TRUNC", "RBFOX2_BrdU_CLIP_CLINK_TRUNC", "RBFOX2_BrdU_CLIP_CTK_DEL", "RBFOX2_BrdU_CLIP_CLINK_DEL")

DIST_RBFOX2_BrdU_CLIP_COMBINED = cbind(
  DIST_RBFOX2_BrdU_CLIP_PEAK,
  DIST_RBFOX2_BrdU_CLIP_CTK_TRUNC,
  DIST_RBFOX2_BrdU_CLIP_CLINK_TRUNC,
  DIST_RBFOX2_BrdU_CLIP_CTK_DEL,
  DIST_RBFOX2_BrdU_CLIP_CLINK_DEL
)
DIST_RBFOX2_BrdU_CLIP_COMBINED$Annotation = rownames(DIST_RBFOX2_BrdU_CLIP_COMBINED)

DIST_RBFOX2_BrdU_CLIP_COMBINED = DIST_RBFOX2_BrdU_CLIP_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_RBFOX2_BrdU_CLIP_COMBINED$Source = factor(DIST_RBFOX2_BrdU_CLIP_COMBINED$Source, levels = CLIP_List)

DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC = DIST_RBFOX2_BrdU_CLIP_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_RBFOX2_BrdU_CLIP_COMBINED_GENIC = tapply(DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC$Freq, DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC$Source, sum)

DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC = DIST_RBFOX2_BrdU_CLIP_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC = tapply(DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC$Freq, DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC = plotStackedBar(
  DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "RBFOX2 BrdU CLIP Genic Annotation Distributions")

PLOT_DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC = plotStackedBar(
  DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "RBFOX2 BrdU CLIP Non-Genic Annotation Distributions")

# print(PLOT_DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC)
# print(PLOT_DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_RBFOX2_BrdU_CLIP_COMBINED_GENIC = DIST_RBFOX2_BrdU_CLIP_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC = DIST_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_GENIC = plotStackedBar(
  PROP_RBFOX2_BrdU_CLIP_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_RBFOX2_BrdU_CLIP_COMBINED_GENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "RBFOX2 BrdU CLIP Genic Annotation Proportions")

PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC = plotStackedBar(
  PROP_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC,
  c("PEAK","CTK\nTRUNC", "CLINK\nTRUNC", "CTK\nDEL", "CLINK\nDEL"),
  "RBFOX2 BrdU CLIP Non-Genic Annotation Proportions")

# print(PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_GENIC)
# print(PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC)
################################################################################

## 4_8 RBFOX2 eCLIP HepG2 Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_RBFOX2_eCLIP_HepG2_PEAK          = countAnnotation(RBFOX2_eCLIP_HepG2_PEAK_FILTERED, "annotation", "RBFOX2_eCLIP_HepG2_PEAK")
DIST_RBFOX2_eCLIP_HepG2_CLINK_TRUNC   = countAnnotation(RBFOX2_eCLIP_HepG2_CLINK_TRUNC_FILTERED, "annotation", "RBFOX2_eCLIP_HepG2_CLINK_TRUNC")
DIST_RBFOX2_eCLIP_HepG2_CLINK_DEL     = countAnnotation(RBFOX2_eCLIP_HepG2_CLINK_DEL_FILTERED, "annotation", "RBFOX2_eCLIP_HepG2_CLINK_DEL")

DIST_RBFOX2_eCLIP_HepG2_CLINK_TRUNC   = fillAnnotation(DIST_RBFOX2_eCLIP_HepG2_CLINK_TRUNC, rownames(DIST_RBFOX2_eCLIP_HepG2_PEAK))
DIST_RBFOX2_eCLIP_HepG2_CLINK_DEL     = fillAnnotation(DIST_RBFOX2_eCLIP_HepG2_CLINK_DEL, rownames(DIST_RBFOX2_eCLIP_HepG2_PEAK))

CLIP_List = c("RBFOX2_eCLIP_HepG2_PEAK", "RBFOX2_eCLIP_HepG2_CLINK_TRUNC", "RBFOX2_eCLIP_HepG2_CLINK_DEL")

DIST_RBFOX2_eCLIP_HepG2_COMBINED = cbind(
  DIST_RBFOX2_eCLIP_HepG2_PEAK,
  DIST_RBFOX2_eCLIP_HepG2_CLINK_TRUNC,
  DIST_RBFOX2_eCLIP_HepG2_CLINK_DEL
)
DIST_RBFOX2_eCLIP_HepG2_COMBINED$Annotation = rownames(DIST_RBFOX2_eCLIP_HepG2_COMBINED)

DIST_RBFOX2_eCLIP_HepG2_COMBINED = DIST_RBFOX2_eCLIP_HepG2_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_RBFOX2_eCLIP_HepG2_COMBINED$Source = factor(DIST_RBFOX2_eCLIP_HepG2_COMBINED$Source, levels = CLIP_List)

DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC = DIST_RBFOX2_eCLIP_HepG2_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_RBFOX2_eCLIP_HepG2_COMBINED_GENIC = tapply(DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC$Freq, DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC$Source, sum)

DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC = DIST_RBFOX2_eCLIP_HepG2_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC = tapply(DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC$Freq, DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC = plotStackedBar(
  DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "RBFOX2 eCLIP HepG2 Genic Annotation Distributions")

PLOT_DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC = plotStackedBar(
  DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "RBFOX2 eCLIP HepG2 Non-Genic Annotation Distributions")

# print(PLOT_DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC)
# print(PLOT_DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_RBFOX2_eCLIP_HepG2_COMBINED_GENIC= DIST_RBFOX2_eCLIP_HepG2_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC = DIST_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_GENIC = plotStackedBar(
  PROP_RBFOX2_eCLIP_HepG2_COMBINED_GENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"), counts = COUNTS_RBFOX2_eCLIP_HepG2_COMBINED_GENIC,
  "RBFOX2 eCLIP HepG2 Genic Annotation Proportions")

PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC = plotStackedBar(
  PROP_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"), counts = COUNTS_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC,
  "RBFOX2 eCLIP HepG2 Non-Genic Annotation Proportions")

# print(PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_GENIC)
# print(PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC)
################################################################################

## 4_9 RBFOX2 eCLIP K562 Visualization
################################################################################
## Peak Counts Distribution Stacked Bar Graph
DIST_RBFOX2_eCLIP_K562_PEAK          = countAnnotation(RBFOX2_eCLIP_K562_PEAK_FILTERED, "annotation", "RBFOX2_eCLIP_K562_PEAK")
DIST_RBFOX2_eCLIP_K562_CLINK_TRUNC   = countAnnotation(RBFOX2_eCLIP_K562_CLINK_TRUNC_FILTERED, "annotation", "RBFOX2_eCLIP_K562_CLINK_TRUNC")
DIST_RBFOX2_eCLIP_K562_CLINK_DEL     = countAnnotation(RBFOX2_eCLIP_K562_CLINK_DEL_FILTERED, "annotation", "RBFOX2_eCLIP_K562_CLINK_DEL")

DIST_RBFOX2_eCLIP_K562_CLINK_DEL     = fillAnnotation(DIST_RBFOX2_eCLIP_K562_CLINK_DEL, rownames(DIST_RBFOX2_eCLIP_K562_PEAK))

CLIP_List = c("RBFOX2_eCLIP_K562_PEAK", "RBFOX2_eCLIP_K562_CLINK_TRUNC", "RBFOX2_eCLIP_K562_CLINK_DEL")

DIST_RBFOX2_eCLIP_K562_COMBINED = cbind(
  DIST_RBFOX2_eCLIP_K562_PEAK,
  DIST_RBFOX2_eCLIP_K562_CLINK_TRUNC,
  DIST_RBFOX2_eCLIP_K562_CLINK_DEL
)
DIST_RBFOX2_eCLIP_K562_COMBINED$Annotation = rownames(DIST_RBFOX2_eCLIP_K562_COMBINED)

DIST_RBFOX2_eCLIP_K562_COMBINED = DIST_RBFOX2_eCLIP_K562_COMBINED %>%
  gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>%
  dplyr::select(Source, Freq, Annotation)
DIST_RBFOX2_eCLIP_K562_COMBINED$Source = factor(DIST_RBFOX2_eCLIP_K562_COMBINED$Source, levels = CLIP_List)

DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC = DIST_RBFOX2_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% Genic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = Genic_Annotation_List))

COUNTS_RBFOX2_eCLIP_K562_COMBINED_GENIC = tapply(DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC$Freq, DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC$Source, sum)

DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC = DIST_RBFOX2_eCLIP_K562_COMBINED %>% 
  filter(Annotation %in% NonGenic_Annotation_List) %>%
  mutate(Annotation = factor(Annotation, levels = NonGenic_Annotation_List))

COUNTS_RBFOX2_eCLIP_K562_COMBINED_NONGENIC = tapply(DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC$Freq, DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC$Source, sum)

PLOT_DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "RBFOX2 eCLIP K562 Genic Annotation Distributions")

PLOT_DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "RBFOX2 eCLIP K562 Non-Genic Annotation Distributions")

# print(PLOT_DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC)

## Add proportional plot for all samples:
PROP_RBFOX2_eCLIP_K562_COMBINED_GENIC = DIST_RBFOX2_eCLIP_K562_COMBINED_GENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PROP_RBFOX2_eCLIP_K562_COMBINED_NONGENIC = DIST_RBFOX2_eCLIP_K562_COMBINED_NONGENIC %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  ungroup()

PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_GENIC = plotStackedBar(
  PROP_RBFOX2_eCLIP_K562_COMBINED_GENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_RBFOX2_eCLIP_K562_COMBINED_GENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"), 
  "RBFOX2 eCLIP K562 Genic Annotation Proportions")

PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_NONGENIC = plotStackedBar(
  PROP_RBFOX2_eCLIP_K562_COMBINED_NONGENIC,
  CLIP_List, color_map = color_map, counts = COUNTS_RBFOX2_eCLIP_K562_COMBINED_NONGENIC,
  c("PEAK", "CLINK\nTRUNC", "CLINK\nDEL"),
  "RBFOX2 eCLIP K562 Non-Genic Annotation Proportions")

# print(PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_GENIC)
# print(PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_NONGENIC)
################################################################################





################################################################################
PLOT_PROP_HUR_CoCLIP_COMBINED_GENIC
PLOT_PROP_HUR_eCLIP_K562_COMBINED_GENIC

PLOT_PROP_HUR_CoCLIP_COMBINED_NONGENIC
PLOT_PROP_HUR_eCLIP_K562_COMBINED_NONGENIC




PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_GENIC
PLOT_PROP_HNRNPC_iCLIP_COMBINED_GENIC
PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_GENIC
PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_GENIC

PLOT_PROP_HNRNPC_BrdU_CLIP2_COMBINED_NONGENIC
PLOT_PROP_HNRNPC_iCLIP_COMBINED_NONGENIC
PLOT_PROP_HNRNPC_eCLIP_HepG2_COMBINED_NONGENIC
PLOT_PROP_HNRNPC_eCLIP_K562_COMBINED_NONGENIC




PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_GENIC
PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_GENIC
PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_GENIC

PLOT_PROP_RBFOX2_BrdU_CLIP_COMBINED_NONGENIC
PLOT_PROP_RBFOX2_eCLIP_HepG2_COMBINED_NONGENIC
PLOT_PROP_RBFOX2_eCLIP_K562_COMBINED_NONGENIC
################################################################################





