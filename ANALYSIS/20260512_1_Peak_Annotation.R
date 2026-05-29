################################################################################
## CLIP Peak Normalization, Biological Complexity, and Annotation
## Author: Soon Yi with Antigravity and claude
## Date: May 2026
## v2: Aligned with JL annotation scheme (first-wins engine, 25-level
##     priority, intron subcategories, TE splitting, rescue step)
################################################################################

# === Load Libraries ===
library(data.table)
library(tibble)
library(tidyverse)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)
library(AnnotationDbi)
library(GenomeInfoDb)

## 1. Basic Setup:
################################################################################
BASE_DIR   = "/mnt/2TB_DATA/ANALYSIS_CLIPittyCLIP/"
INPUT_DIR  = paste0(BASE_DIR, "ANALYSIS/")
OUTPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/1_ANNOTATED/")
ANNO_DIR   = "/mnt/2TB_DATA/Annotations/"

# GTF_FILE   = paste0(ANNO_DIR, "Human/gencode.v49.primary_assembly.basic.annotation.gtf")
# RMSK_FILE  = paste0(ANNO_DIR, "Human/Rmsk_hg38_2026Jan.txt")
# TXDB_CACHE = paste0(ANNO_DIR, "Human/gencode.v49.txdb.sqlite")  # separate cache — don't reuse the mouse one

GTF_FILE    = paste0(ANNO_DIR, "Mouse/gencode.vM38.primary_assembly.annotation.gtf")
RMSK_FILE   = paste0(ANNO_DIR, "Mouse/Rmsk_mm39_2026May.txt")
TXDB_CACHE  = paste0(ANNO_DIR, "Mouse/gencode.vM38.txdb.sqlite")  # built once, loaded thereafter

PROXIMAL_WINDOW = 300L   # nt from splice site defining "Proximal intronic"
DISTAL_UTR_EXT  = 10000L # nt past most-distal 3'UTR end for "3'UTR-distal"
RESCUE_WINDOW   = 10000L # nt around any transcript body for deep-intergenic rescue
################################################################################

## 2. Load Data:
################################################################################

## ── TxDb (build once; reload from cache on subsequent runs) ──────────────────
if (file.exists(TXDB_CACHE)) {
  message("Loading TxDb from cache: ", TXDB_CACHE)
  txdb = loadDb(TXDB_CACHE)
} else {
  message("Building TxDb from GTF (~1 min)...")
  txdb = makeTxDbFromGFF(GTF_FILE, format = "gtf")
  saveDb(txdb, TXDB_CACHE)
  message("TxDb cached to: ", TXDB_CACHE)
}

## ── GTF (needed for transcript biotype subsetting) ───────────────────────────
message("Importing GTF...")
gtf_all = import(GTF_FILE)

## ── RepeatMasker ─────────────────────────────────────────────────────────────
repMask_Main = read.delim(RMSK_FILE)
repMask_Main = repMask_Main[, c("genoName","genoStart","genoEnd","strand","repName","repClass","repFamily")]
colnames(repMask_Main) = c("chr","start","end","strand","name","repClass","repFamily")
## Fix rmsk strand encoding: "C" (complement) → "-"; anything else non-standard → "*"
repMask_Main$strand = ifelse(repMask_Main$strand == "C",  "-", repMask_Main$strand)
repMask_Main$strand = ifelse(repMask_Main$strand %in% c("+","-"), repMask_Main$strand, "*")

## ── BiomaRt ───────────────────────────────────────────────────────────────────
# mart.hs = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                   host    = "https://may2025.archive.ensembl.org",
#                   dataset = "hsapiens_gene_ensembl")

mart.hs = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   host    = "https://may2025.archive.ensembl.org",
                   dataset = "mmusculus_gene_ensembl")
################################################################################

## 3. Build Annotation Features:
################################################################################

## ── Group A: Genic features ──────────────────────────────────────────────────
message("Building Group A: Genic features...")

feat_3utr = unlist(threeUTRsByTranscript(txdb))
feat_5utr = unlist(fiveUTRsByTranscript(txdb))
feat_cds  = unlist(cdsBy(txdb, by = "tx"))
feat_exon = unlist(exonsBy(txdb, by = "gene"))

## 3'UTR-distal: extend DISTAL_UTR_EXT nt past the most-distal annotated 3'UTR
## end, anchored per gene (not per isoform) to avoid overlapping longer-isoform
## annotated 3'UTRs.
tx2gene_map = AnnotationDbi::select(
  txdb,
  keys    = keys(txdb, keytype = "TXID"),
  columns = "GENEID",
  keytype = "TXID"
)
tx2gene_map = data.frame(
  TXID   = as.integer(tx2gene_map$TXID),
  GENEID = as.character(tx2gene_map$GENEID),
  stringsAsFactors = FALSE
)

utr3_names = names(feat_3utr)   # save names before stripping

utr3_df = as.data.frame(unname(feat_3utr))[, c("seqnames", "start", "end", "strand")]
utr3_df$seqnames = as.character(utr3_df$seqnames)
utr3_df$strand   = as.character(utr3_df$strand)
utr3_df$TXID     = as.integer(sub("\\.[0-9]+$", "", utr3_names))
utr3_df = utr3_df %>%
  left_join(tx2gene_map, by = "TXID") %>%
  filter(!is.na(GENEID)) %>%
  mutate(GENEID = as.character(GENEID),
         start  = as.integer(start),
         end    = as.integer(end))

distal_plus  = utr3_df %>% filter(strand == "+") %>%
  group_by(GENEID) %>% dplyr::slice(which.max(.data$end)) %>% ungroup()

distal_minus = utr3_df %>% filter(strand == "-") %>%
  group_by(GENEID) %>% dplyr::slice(which.min(.data$start)) %>% ungroup()

feat_3utr_distal = suppressWarnings(c(
  GRanges(seqnames = distal_plus$seqnames,  strand = "+",
          ranges = IRanges(start = distal_plus$end + 1L,
                           end   = distal_plus$end + DISTAL_UTR_EXT)),
  GRanges(seqnames = distal_minus$seqnames, strand = "-",
          ranges = IRanges(start = pmax(1L, distal_minus$start - DISTAL_UTR_EXT),
                           end   = distal_minus$start - 1L))
))
rm(tx2gene_map, utr3_df, distal_plus, distal_minus)

## TSS-proximal: 1 kb upstream / 100 bp downstream of TSS per gene
feat_tss_proximal = suppressWarnings(
  promoters(genes(txdb), upstream = 1000L, downstream = 100L)
)

## Retained introns (split by host gene biotype)
ri_all         = gtf_all[!is.na(gtf_all$transcript_biotype) &
                           gtf_all$transcript_biotype == "retained_intron" &
                           gtf_all$type == "transcript"]
feat_ri_coding = ri_all[ri_all$gene_type == "protein_coding"]
feat_ri_noncod = ri_all[ri_all$gene_type != "protein_coding"]

## All introns (reduced across all transcripts)
message("  Extracting all introns (~30 s)...")
all_introns_flat    = unlist(intronsByTranscript(txdb))
all_introns_reduced = suppressWarnings(
  keepStandardChromosomes(reduce(all_introns_flat), pruning.mode = "coarse")
)

## First intron: for each transcript, the intron immediately after exon 1
## (minimum genomic start on +, maximum start on -).
message("  Extracting first introns...")
intron_names = names(all_introns_flat)   # save before stripping

intron_df = as.data.frame(unname(all_introns_flat))[, c("seqnames", "start", "end", "strand")]
intron_df$seqnames = as.character(intron_df$seqnames)
intron_df$strand   = as.character(intron_df$strand)
intron_df$tx_id    = sub("\\.[0-9]+$", "", intron_names)

first_df = bind_rows(
  intron_df %>% filter(strand == "+")           %>% group_by(tx_id) %>% dplyr::slice(which.min(.data$start)) %>% ungroup(),
  intron_df %>% filter(strand == "-")           %>% group_by(tx_id) %>% dplyr::slice(which.max(.data$start)) %>% ungroup(),
  intron_df %>% filter(!strand %in% c("+","-")) %>% group_by(tx_id) %>% dplyr::slice(1L) %>% ungroup()
)
feat_first_intron = suppressWarnings(keepStandardChromosomes(
  reduce(GRanges(seqnames = first_df$seqnames,
                 ranges   = IRanges(start = first_df$start, end = first_df$end),
                 strand   = first_df$strand)),
  pruning.mode = "coarse"
))
rm(intron_df, first_df)

## Non-first introns → proximal (±PROXIMAL_WINDOW nt of splice site) or deep
feat_non_first_intron = suppressWarnings(
  setdiff(all_introns_reduced, feat_first_intron)
)

make_proximal_windows = function(intron_gr, window = PROXIMAL_WINDOW) {
  if (length(intron_gr) == 0) return(GRanges())
  w5 = GRanges(seqnames = seqnames(intron_gr), strand = strand(intron_gr),
               ranges = IRanges(start = start(intron_gr),
                                end   = pmin(start(intron_gr) + window - 1L, end(intron_gr))))
  w3 = GRanges(seqnames = seqnames(intron_gr), strand = strand(intron_gr),
               ranges = IRanges(start = pmax(end(intron_gr) - window + 1L, start(intron_gr)),
                                end   = end(intron_gr)))
  suppressWarnings(reduce(c(w5, w3)))
}

feat_proximal_intronic = make_proximal_windows(feat_non_first_intron)
feat_deep_intronic     = suppressWarnings(
  setdiff(feat_non_first_intron, feat_proximal_intronic)
)

## Transcript bodies extended ±RESCUE_WINDOW (for deep-intergenic rescue step)
tx_gr  = transcripts(txdb)
tx_ext = suppressWarnings(GRanges(
  seqnames = seqnames(tx_gr), strand = "*",
  ranges   = IRanges(start = pmax(1L, start(tx_gr) - RESCUE_WINDOW),
                     end   = end(tx_gr) + RESCUE_WINDOW)
))

## Gene bodies (for biomaRt gene-name lookup downstream)
genes_gr  = genes(txdb)
genes_bed = data.frame(
  chr    = as.character(seqnames(genes_gr)),
  start  = start(genes_gr) - 1L,
  end    = end(genes_gr),
  name   = unlist(genes_gr$gene_id),
  score  = 0,
  strand = as.character(strand(genes_gr))
)

message("Group A done.")

## ── Group B: Discrete ncRNA genes ────────────────────────────────────────────
## Called BEFORE intronic categories so that lncRNA/miRNA gene bodies are not
## mis-classified as "Deep intronic" (lncRNA loci are largely intronic).
message("Building Group B: Discrete ncRNA genes...")

make_ncrna_gr = function(biotype_val) {
  idx = !is.na(gtf_all$transcript_biotype) &
    gtf_all$transcript_biotype == biotype_val &
    gtf_all$type == "transcript"
  gtf_all[idx]
}

feat_mirna  = make_ncrna_gr("miRNA")
feat_lncrna = make_ncrna_gr("lncRNA")

## YRNA: misc_RNA entries whose gene name matches RNY[0-9P]
feat_yrna = gtf_all[
  gtf_all$type == "transcript" &
    !is.na(gtf_all$gene_name) &
    grepl("^RNY[0-9P]", gtf_all$gene_name, ignore.case = TRUE)
]

## Other ncRNA: misc_RNA (excluding YRNA) + scaRNA
## scaRNA is absorbed here per Joe's scheme (no separate category)
feat_misc = gtf_all[
  gtf_all$type == "transcript" &
    !is.na(gtf_all$transcript_biotype) &
    gtf_all$transcript_biotype == "misc_RNA" &
    !(!is.na(gtf_all$gene_name) & grepl("^RNY[0-9P]", gtf_all$gene_name, ignore.case = TRUE))
]
feat_scarna = make_ncrna_gr("scaRNA")
feat_other_ncrna = suppressWarnings(c(feat_misc, feat_scarna))
rm(feat_misc, feat_scarna)

message("Group B done.")

## ── Group C: Structural RNAs ─────────────────────────────────────────────────
## (* = multi-mapping concern with unique-mapped reads)
message("Building Group C: Structural RNAs...")
feat_rrna     = make_ncrna_gr("rRNA")
feat_snrna    = make_ncrna_gr("snRNA")
feat_snorna   = make_ncrna_gr("snoRNA")
feat_trna_gtf = make_ncrna_gr("tRNA")
message("Group C done.")

## ── Group D: Transposable elements (RepeatMasker) ────────────────────────────
## Called BEFORE broad intronic categories so that TE-overlapping intronic peaks
## are labelled as the TE, not "Deep intronic".
## NOTE: ignore.strand = TRUE used for TEs (see annotation engine below) because
##   rmsk strand annotation is imperfect (many old/degenerate copies are "*").
message("Building Group D: TEs from RepeatMasker...")

make_rmsk_gr = function(class_filter = NULL, family_filter = NULL, negate_family = FALSE) {
  sub_df = repMask_Main
  if (!is.null(class_filter))  sub_df = sub_df[sub_df$repClass %in% class_filter, ]
  if (!is.null(family_filter)) {
    if (negate_family) sub_df = sub_df[!(sub_df$repFamily %in% family_filter), ]
    else               sub_df = sub_df[   sub_df$repFamily %in%  family_filter,  ]
  }
  if (nrow(sub_df) == 0) return(GRanges())
  ## BED coords are 0-based: add 1 to start for GRanges (1-based)
  GRanges(seqnames  = sub_df$chr,
          ranges    = IRanges(start = sub_df$start + 1L, end = sub_df$end),
          strand    = sub_df$strand,
          name      = sub_df$name,
          repClass  = sub_df$repClass,
          repFamily = sub_df$repFamily)
}

feat_alu        = make_rmsk_gr("SINE", "Alu")
feat_other_sine = make_rmsk_gr("SINE", "Alu", negate_family = TRUE)
feat_line       = make_rmsk_gr("LINE")
feat_ltr        = make_rmsk_gr("LTR")
feat_satellite  = make_rmsk_gr(c("Satellite","Low_complexity","Simple_repeat"))
feat_trna_rmsk  = make_rmsk_gr("tRNA")

## Combined tRNA: GTF entries + RepeatMasker tRNA-derived SINEs
feat_trna = suppressWarnings(
  if (length(feat_trna_gtf) > 0 && length(feat_trna_rmsk) > 0)
    c(feat_trna_gtf, feat_trna_rmsk)
  else if (length(feat_trna_gtf) > 0) feat_trna_gtf
  else feat_trna_rmsk
)

message("Group D done.")
message("All annotation features built.")
################################################################################

## 4. Job Table:
################################################################################
## Add or remove rows to run any combination of files in a single pass.
## Output directories are created automatically if they don't exist.

# jobs = data.frame(
#   input = c(
#     "/mnt/2TB_DATA/ANALYSIS_HUR/coCLIP/20260512_HuR_CLIP_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/coCLIP/20260512_HuR_CLIP_CTK/5_CTK_Analysis/HuR_WT/CITS/HuR_WT_CITS.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/coCLIP/20260512_HuR_CLIP_CTK/5_CTK_Analysis/HuR_WT/CIMS/HuR_WT_CIMS_del.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/coCLIP/20260512_HuR_CLIP_CLINK/5_Clink/GROUP_HuR_WT/HuR_WT_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/coCLIP/20260512_HuR_CLIP_CLINK/5_Clink/GROUP_HuR_WT/HuR_WT_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/SPECIFICITY/20260512_HNRNPC_CLIP_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/SPECIFICITY/20260512_HNRNPC_CLIP_CTK/5_CTK_Analysis/HNRNPC_WT/CITS/HNRNPC_WT_CITS.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/SPECIFICITY/20260512_HNRNPC_CLIP_CTK/5_CTK_Analysis/HNRNPC_WT/CIMS/HNRNPC_WT_CIMS_del.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/SPECIFICITY/20260512_HNRNPC_CLIP_CLINK/5_Clink/GROUP_HNRNPC_WT/HNRNPC_WT_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/SPECIFICITY/20260512_HNRNPC_CLIP_CLINK/5_Clink/GROUP_HNRNPC_WT/HNRNPC_WT_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/iCLIP/20260512_HNRNPC_iCLIP_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/iCLIP/20260512_HNRNPC_iCLIP_CTK/5_CTK_Analysis/HNRNPC_WT/CITS/HNRNPC_WT_CITS.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/iCLIP/20260512_HNRNPC_iCLIP_CTK/5_CTK_Analysis/HNRNPC_WT/CIMS/HNRNPC_WT_CIMS_del.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/iCLIP/20260512_HNRNPC_iCLIP_CLINK/5_Clink/GROUP_HNRNPC_WT/HNRNPC_WT_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/iCLIP/20260512_HNRNPC_iCLIP_CLINK/5_Clink/GROUP_HNRNPC_WT/HNRNPC_WT_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_HEPG2_STAR_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_HEPG2_STAR_CLINK/5_Clink/GROUP_HepG2_HNRNPC/HepG2_HNRNPC_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_HEPG2_STAR_CLINK/5_Clink/GROUP_HepG2_HNRNPC/HepG2_HNRNPC_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_K562_STAR_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_HNRNPC/K562_HNRNPC_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HNRNPC/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_HNRNPC/K562_HNRNPC_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_HUR/eCLIP/20260511_K562_STAR_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_ELAVL1/K562_ELAVL1_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_HUR/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_ELAVL1/K562_ELAVL1_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_HEPG2_STAR_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_HEPG2_STAR_CLINK/5_Clink/GROUP_HepG2_RBFOX2/HepG2_RBFOX2_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_HEPG2_STAR_CLINK/5_Clink/GROUP_HepG2_RBFOX2/HepG2_RBFOX2_deletions.bed",
#     
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_K562_STAR_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_RBFOX2/K562_RBFOX2_truncations.bed",
#     "/mnt/2TB_DATA/ANALYSIS_RBFOX2/eCLIP/20260511_K562_STAR_CLINK/5_Clink/GROUP_K562_RBFOX2/K562_RBFOX2_deletions.bed"
#   ),
#   output = c(
#     paste0(OUTPUT_DIR, "HuR_CoCLIP_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_CoCLIP_CTK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_CoCLIP_CTK_DEL_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_CoCLIP_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_CoCLIP_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CTK_DEL_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_BrdU_CLIP2_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CTK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CTK_DEL_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_iCLIP_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_HepG2_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HNRNPC_eCLIP_K562_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "HuR_eCLIP_K562_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "HuR_eCLIP_K562_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_HepG2_CLINK_DEL_ANNOTATED.txt"),
#     
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_PEAK_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_TRUNC_ANNOTATED.txt"),
#     paste0(OUTPUT_DIR, "RBFOX2_eCLIP_K562_CLINK_DEL_ANNOTATED.txt")
#     
#   ),
#   stringsAsFactors = FALSE
# )

jobs = data.frame(
  input = c(
    "/mnt/2TB_DATA/ANALYSIS_RBFOX2/CZ/20260512_RBFOX2_BrdU_CLIP_CLINK/4_PEAKS/COMBINED_PEAKS/COMBINED_PEAK_MATRIX.txt",
    "/mnt/2TB_DATA/ANALYSIS_RBFOX2/CZ/20260512_RBFOX2_BrdU_CLIP_CTK/5_CTK_Analysis/RBFOX2/CITS/RBFOX2_CITS.txt",
    "/mnt/2TB_DATA/ANALYSIS_RBFOX2/CZ/20260512_RBFOX2_BrdU_CLIP_CTK/5_CTK_Analysis/RBFOX2/CIMS/RBFOX2_CIMS_del.txt",
    "/mnt/2TB_DATA/ANALYSIS_RBFOX2/CZ/20260512_RBFOX2_BrdU_CLIP_CLINK/5_Clink/GROUP_RBFOX2/RBFOX2_truncations.bed",
    "/mnt/2TB_DATA/ANALYSIS_RBFOX2/CZ/20260512_RBFOX2_BrdU_CLIP_CLINK/5_Clink/GROUP_RBFOX2/RBFOX2_deletions.bed"
  ),
  output = c(
    paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_PEAK_ANNOTATED.txt"),
    paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_TRUNC_ANNOTATED.txt"),
    paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CTK_DEL_ANNOTATED.txt"),
    paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_TRUNC_ANNOTATED.txt"),
    paste0(OUTPUT_DIR, "RBFOX2_BrdU_CLIP_CLINK_DEL_ANNOTATED.txt")
  ),
  stringsAsFactors = FALSE
)


################################################################################

## 5. Annotation Engine — constants defined once outside the loop:
################################################################################
##
## Priority (first-wins; a peak assigned to level N is invisible to N+1 onward):
##
##  1. 3'UTR                  ← wins over everything
##  2. 3'UTR-distal
##  3. 5'UTR
##  4. CDS
##
##  5. miRNA                  ← ncRNA genes before introns (so lncRNA intronic
##  6. lncRNA                   peaks aren't called "Deep intronic")
##  7. YRNA
##  8. rRNA*
##  9. snRNA*
## 10. snoRNA*
## 11. tRNA*
## 12. Other ncRNA
##
## 13. Alu                    ← TEs before broad intron categories (so an Alu
## 14. Other SINE               inside an intron is called "Alu", not "Intronic")
## 15. LINE
## 16. LTR/ERV
## 17. Satellite/LC
##
## 18. Retained intron:coding ← specific biotype before general intron classes
## 19. Retained intron:ncRNA
## 20. First intron
## 21. Proximal intronic
## 22. Deep intronic
## 23. Other exon
## 24. TSS-proximal
##
## 25. Rescue: Deep intronic  ← unassigned but within RESCUE_WINDOW of any transcript
##     Rescue: Intergenic     ← truly gene-distal
##
## Strand handling:
##   ignore.strand = FALSE for genic / ncRNA features (CLIP peaks are strand-specific)
##   ignore.strand = TRUE  for TE features (rmsk strand annotation is imperfect;
##                                          many old/degenerate copies annotated as "*")

## assign_hits() is defined once here and reused across all loop iterations.
## It writes into `labels` and `assigned` via <<-, which resolve to whatever
## iteration's environment they were initialized in (top-level for loop = global).
assign_hits = function(feat_gr, label, ignore.strand = FALSE) {
  if (length(feat_gr) == 0) return(invisible(NULL))
  hits = suppressWarnings(
    findOverlaps(peaksGR, feat_gr, minoverlap = 1L, ignore.strand = ignore.strand)
  )
  idx = unique(queryHits(hits))
  idx = idx[!assigned[idx]]
  if (length(idx) == 0) return(invisible(NULL))
  labels[idx]   <<- label
  assigned[idx] <<- TRUE
}

ANNOT_LEVELS = c(
  ## Group A — Genic
  "3'UTR", "3'UTR-distal", "5'UTR", "CDS",
  "Retained intron:coding", "Retained intron:ncRNA",
  "First intron", "Proximal intronic", "Deep intronic",
  "Other exon", "TSS-proximal",
  ## Group B — Discrete ncRNA
  "miRNA", "lncRNA", "YRNA", "Other ncRNA",
  ## Group C — Structural RNAs
  "rRNA*", "snRNA*", "snoRNA*", "tRNA*",
  ## Group D — TEs
  "Alu", "Other SINE", "LINE", "LTR/ERV", "Satellite/LC",
  ## Intergenic
  "Intergenic"
)
################################################################################

## 6. Loop: Annotate each file in the job table:
################################################################################
for (i in seq_len(nrow(jobs))) {
  
  PEAK_MATRIX_FILE = jobs$input[i]
  OUTPUT_FILE      = jobs$output[i]
  message(sprintf("[%d/%d] Processing: %s", i, nrow(jobs), basename(PEAK_MATRIX_FILE)))
  
  ## Create output directory if it doesn't exist
  dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
  
  ## ── Load peaks ─────────────────────────────────────────────────────────────
  peaksMatrix = fread(PEAK_MATRIX_FILE)
  colnames(peaksMatrix)[1:6] = c("chr", "start", "end", "name", "score", "strand")
  
  ## BED coords are 0-based: add 1 to start for GRanges (1-based)
  peaksGR = GRanges(
    seqnames = peaksMatrix$chr,
    ranges   = IRanges(start = peaksMatrix$start + 1L, end = peaksMatrix$end),
    strand   = peaksMatrix$strand
  )
  
  ## ── First-wins annotation ──────────────────────────────────────────────────
  ## Reset per-iteration state; assign_hits() writes into these via <<-
  labels   = rep(NA_character_, length(peaksGR))
  assigned = rep(FALSE,         length(peaksGR))
  
  ## Group A: high-confidence mRNA features
  assign_hits(feat_3utr,        "3'UTR")
  assign_hits(feat_3utr_distal, "3'UTR-distal")
  assign_hits(feat_5utr,        "5'UTR")
  assign_hits(feat_cds,         "CDS")
  
  ## Group B: discrete ncRNA genes (before introns)
  assign_hits(feat_mirna,       "miRNA")
  assign_hits(feat_lncrna,      "lncRNA")
  assign_hits(feat_yrna,        "YRNA")
  
  ## Group C: structural RNAs
  assign_hits(feat_rrna,        "rRNA*")
  assign_hits(feat_snrna,       "snRNA*")
  assign_hits(feat_snorna,      "snoRNA*")
  assign_hits(feat_trna,        "tRNA*")
  assign_hits(feat_other_ncrna, "Other ncRNA")
  
  ## Group D: TEs (ignore.strand — rmsk strand annotation is imperfect)
  assign_hits(feat_alu,         "Alu",          ignore.strand = TRUE)
  assign_hits(feat_other_sine,  "Other SINE",   ignore.strand = TRUE)
  assign_hits(feat_line,        "LINE",         ignore.strand = TRUE)
  assign_hits(feat_ltr,         "LTR/ERV",      ignore.strand = TRUE)
  assign_hits(feat_satellite,   "Satellite/LC", ignore.strand = TRUE)
  
  ## Group A continued: intron subcategories
  assign_hits(feat_ri_coding,         "Retained intron:coding")
  assign_hits(feat_ri_noncod,         "Retained intron:ncRNA")
  assign_hits(feat_first_intron,      "First intron")
  assign_hits(feat_proximal_intronic, "Proximal intronic")
  assign_hits(feat_deep_intronic,     "Deep intronic")
  assign_hits(feat_exon,              "Other exon")
  assign_hits(feat_tss_proximal,      "TSS-proximal")
  
  ## Rescue: unassigned peaks
  remaining = which(!assigned)
  if (length(remaining) > 0) {
    near_tx = suppressWarnings(
      countOverlaps(peaksGR[remaining], tx_ext, ignore.strand = TRUE)
    ) > 0
    labels[remaining[ near_tx]] = "Deep intronic"
    labels[remaining[!near_tx]] = "Intergenic"
  }
  
  peaksMatrix$annotation = labels
  
  ## ── Gene lookup ────────────────────────────────────────────────────────────
  genes_hits       = as.data.frame(findOverlaps(query = peaksGR, subject = genes_gr, minoverlap = 1L, select = "first"))
  peaksMatrix$gene = genes_bed[genes_hits[, 1], "name"]
  
  peaksMatrix$gene_base = sub("\\.[0-9]+$", "", peaksMatrix$gene)
  
  gene_names = getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = na.omit(unique(peaksMatrix$gene_base)),
    mart       = mart.hs
  )
  
  peaksMatrix = peaksMatrix %>%
    left_join(gene_names, by = c("gene_base" = "ensembl_gene_id"), relationship = "many-to-many") %>%
    dplyr::select(-gene_base)
  
  ## ── Grouped annotation & factor ordering ───────────────────────────────────
  peaksMatrix = peaksMatrix %>%
    mutate(
      annotation = factor(annotation, levels = ANNOT_LEVELS),
      grouped_annotation = case_when(
        annotation %in% c("3'UTR", "3'UTR-distal", "5'UTR", "CDS",
                          "Retained intron:coding", "Retained intron:ncRNA",
                          "First intron", "Proximal intronic", "Deep intronic",
                          "Other exon", "TSS-proximal")                              ~ "Genic",
        annotation %in% c("miRNA", "lncRNA", "YRNA", "Other ncRNA",
                          "rRNA*", "snRNA*", "snoRNA*", "tRNA*")                     ~ "ncRNA",
        annotation %in% c("Alu", "Other SINE", "LINE", "LTR/ERV", "Satellite/LC")   ~ "TE",
        annotation == "Intergenic"                                                    ~ "Intergenic",
        TRUE                                                                          ~ NA_character_
      )
    )
  
  ## ── Save ───────────────────────────────────────────────────────────────────
  fwrite(peaksMatrix, OUTPUT_FILE, sep = "\t")
  message("  → Written to: ", OUTPUT_FILE)
  
}

message("All ", nrow(jobs), " files processed.")
################################################################################