################################################################################
## CLIP Peak Normalization, Biological Complexity, and Annotation
## Author: Soon Yi
## Date: May 2026
################################################################################

# === Load Libraries ===
library(data.table)
library(tibble)
library(tidyverse)
library(tidyr)
library(readr)
library(stringr)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)

## 1. Basic Setup:
################################################################################
BASE_DIR = "/mnt/2TB_DATA/ANALYSIS_CLIPittyCLIP/"
INPUT_DIR = paste0(BASE_DIR, "ANALYSIS/")
OUTPUT_DIR = paste0(BASE_DIR, "ANALYSIS_OUTPUT/")
ANNO_DIR = "/mnt/2TB_DATA/Annotations/"

# GTF_FILE = paste0(ANNO_DIR, "Human/gencode.v49.primary_assembly.basic.annotation.gtf")
# RMSK_FILE = paste0(ANNO_DIR, "Human/Rmsk_hg38_2026Jan.txt")

GTF_FILE = paste0(ANNO_DIR, "Mouse/gencode.vM38.primary_assembly.annotation.gtf")
RMSK_FILE = paste0(ANNO_DIR, "Mouse/Rmsk_mm39_2026May.txt")

rmchr = function(gr) {
  seqlevels(gr) = gsub("^chr", "", seqlevels(gr))
  return(gr)
}
################################################################################

## 2. Load Data:
################################################################################
gtf = import(GTF_FILE)
gtf_df = as.data.frame(gtf)

## Use the same GTF to make txdb
txdb = makeTxDbFromGFF(GTF_FILE, format="gtf")

repMask_Main = read.delim(RMSK_FILE)
repMask_Main = repMask_Main[, c('genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass', 'repFamily')]
colnames(repMask_Main) = c('chr', 'start', 'end', 'strand', 'name', 'repClass', 'repFamily')

# mart.hs = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                   host = "https://may2025.archive.ensembl.org",
#                   dataset = "hsapiens_gene_ensembl")

mart.hs = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  host = "https://may2025.archive.ensembl.org",
                  dataset = "mmusculus_gene_ensembl")
################################################################################

## Subset the data_frame by gene/transcript types:
################################################################################
protCoding = subset(gtf_df, gene_type == "protein_coding")
miR = subset(gtf_df, transcript_type == "miRNA" & type == "transcript")
lncRNA = subset(gtf_df, transcript_type == "lncRNA" & type == "transcript")
rRNA = subset(gtf_df, transcript_type == "rRNA" & type == "transcript")
snoRNA = subset(gtf_df, transcript_type == "snoRNA" & type == "transcript")
scaRNA = subset(gtf_df, transcript_type == "scaRNA" & type == "transcript")
snRNA = subset(gtf_df, transcript_type == "snRNA" & type == "transcript")
miscRNA = subset(gtf_df, transcript_type == "misc_RNA" & type == "transcript")
Prot_retained_int = subset(gtf_df, transcript_type == "retained_intron" & gene_type == "protein_coding" & type == "transcript")
nc_retained_int = subset(gtf_df, transcript_type == "retained_intron" & gene_type != "protein_coding" & type == "transcript")

## Make Granges for GTF derived categories above:
miR.gr = GRanges(seqnames=miR$seqnames, ranges=IRanges(start=miR$start, end=miR$end, names=miR$gene_name), strand=miR$strand)
lncRNA.gr = GRanges(seqnames=lncRNA$seqnames, ranges=IRanges(start=lncRNA$start, end=lncRNA$end, names=lncRNA$gene_name), strand=lncRNA$strand)
rRNA.gr = GRanges(seqnames=rRNA$seqnames, ranges=IRanges(start=rRNA$start, end=rRNA$end, names=rRNA$gene_name), strand=rRNA$strand)
snoRNA.gr = GRanges(seqnames=snoRNA$seqnames, ranges=IRanges(start=snoRNA$start, end=snoRNA$end, names=snoRNA$gene_name), strand=snoRNA$strand)
scaRNA.gr = GRanges(seqnames=scaRNA$seqnames, ranges=IRanges(start=scaRNA$start, end=scaRNA$end, names=scaRNA$gene_name), strand=scaRNA$strand)
snRNA.gr = GRanges(seqnames=snRNA$seqnames, ranges=IRanges(start=snRNA$start, end=snRNA$end, names=snRNA$gene_name), strand=snRNA$strand)
miscRNA.gr = GRanges(seqnames=miscRNA$seqnames, ranges=IRanges(start=miscRNA$start, end=miscRNA$end, names=miscRNA$gene_name), strand=miscRNA$strand)
Prot_retained_int.gr = GRanges(seqnames=Prot_retained_int$seqnames, ranges=IRanges(start=Prot_retained_int$start, end=Prot_retained_int$end, names=Prot_retained_int$gene_name), strand=Prot_retained_int$strand)
nc_retained_int.gr = GRanges(seqnames=nc_retained_int$seqnames, ranges=IRanges(start=nc_retained_int$start, end=nc_retained_int$end, names=nc_retained_int$gene_name), strand=nc_retained_int$strand)


## Extract features of interest from hg38 TxDb:
exons = unique(unlist(exonsBy(txdb, "tx", use.names=T)))    
CDS = unique(unlist(cdsBy(txdb, "tx", use.names=T)))
introns = unique(unlist(intronsByTranscript(txdb, use.names=T)))
fiveUTRs = unique(unlist(fiveUTRsByTranscript(txdb,use.names=T)))
threeUTRs = unique(unlist(threeUTRsByTranscript(txdb,use.names=T)))

## RepeatMasker features:
tRNA.bed = subset(repMask_Main, repClass == 'tRNA')
LINE.bed = subset(repMask_Main, repClass == 'LINE')
LC_SR.bed = subset(repMask_Main, (repClass == 'Low_complexity' | repClass == 'Simple_repeat'))
LTR.bed = subset(repMask_Main, repClass == 'LTR')
Satellite.bed = subset(repMask_Main, repClass == 'Satellite')
SINE.bed = subset(repMask_Main, repClass == 'SINE')

tRNA.gr = GRanges(seqnames=tRNA.bed$chr, ranges=IRanges(start=tRNA.bed$start, end=tRNA.bed$end, names=tRNA.bed$name), strand=tRNA.bed$strand)
LINE.gr = GRanges(seqnames=LINE.bed$chr, ranges=IRanges(start=LINE.bed$start, end=LINE.bed$end, names=LINE.bed$name), strand=LINE.bed$strand)
LC_SR.gr = GRanges(seqnames=LC_SR.bed$chr, ranges=IRanges(start=LC_SR.bed$start, end=LC_SR.bed$end, names=LC_SR.bed$name), strand=LC_SR.bed$strand)
LTR.gr = GRanges(seqnames=LTR.bed$chr, ranges=IRanges(start=LTR.bed$start, end=LTR.bed$end, names=LTR.bed$name), strand=LTR.bed$strand)
Satellite.gr = GRanges(seqnames=Satellite.bed$chr, ranges=IRanges(start=Satellite.bed$start, end=Satellite.bed$end, names=Satellite.bed$name), strand=Satellite.bed$strand)
SINE.gr = GRanges(seqnames=SINE.bed$chr, ranges=IRanges(start=SINE.bed$start, end=SINE.bed$end, names=SINE.bed$name), strand=SINE.bed$strand)

## Extract gene regions for annotation of genes:
genes = genes(txdb)
genes.bed = data.frame(chr=seqnames(genes), start=start(genes)-1, end=end(genes), name=unlist(genes$gene_id), score=0, strand=strand(genes))
################################################################################

## 3. Annotation - Feature Setup:
################################################################################
##
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "PEAK/HuR_CoCLIP_PEAK", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "PEAK/HuR_CoCLIP_PEAK", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "PEAK/HNRNPC_iCLIP_PEAK", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "PEAK/HNRNPC_iCLIP_PEAK", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "PEAK/HNRNPC_BrdU_CLIP2_PEAK", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "PEAK/HNRNPC_BrdU_CLIP2_PEAK", "_NORMALIZED_ANNOTATED.txt")


##
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/DEL/HuR_CoCLIP_DEL", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/DEL/HuR_CoCLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/DEL/HNRNPC_iCLIP_DEL", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/DEL/HNRNPC_iCLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/DEL/HNRNPC_BrdU_CLIP2_DEL", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/DEL/HNRNPC_BrdU_CLIP2_DEL", "_NORMALIZED_ANNOTATED.txt")


##
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/TRUNC/HuR_CoCLIP_TRUNC", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/TRUNC/HuR_CoCLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/TRUNC/HNRNPC_iCLIP_TRUNC", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/TRUNC/HNRNPC_iCLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/TRUNC/HNRNPC_BrdU_CLIP2_TRUNC", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/TRUNC/HNRNPC_BrdU_CLIP2_TRUNC", "_NORMALIZED_ANNOTATED.txt")


## USE MOUSE GENOME FOR RBFOX2:
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "PEAK/RBFOX2_BrdU_CLIP_PEAK", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "PEAK/RBFOX2_BrdU_CLIP_PEAK", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/DEL/RBFOX2_BrdU_CLIP_DEL", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/DEL/RBFOX2_BrdU_CLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/TRUNC/RBFOX2_BrdU_CLIP_TRUNC", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/TRUNC/RBFOX2_BrdU_CLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")



##
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/DEL/HuR_CoCLIP_DEL", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/DEL/HuR_CoCLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/DEL/HNRNPC_iCLIP_DEL", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/DEL/HNRNPC_iCLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/DEL/HNRNPC_BrdU_CLIP2_DEL", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/DEL/HNRNPC_BrdU_CLIP2_DEL", "_NORMALIZED_ANNOTATED.txt")


##
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/TRUNC/HuR_CoCLIP_TRUNC", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/TRUNC/HuR_CoCLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/TRUNC/HNRNPC_iCLIP_TRUNC", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/TRUNC/HNRNPC_iCLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/TRUNC/HNRNPC_BrdU_CLIP2_TRUNC", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/TRUNC/HNRNPC_BrdU_CLIP2_TRUNC", "_NORMALIZED_ANNOTATED.txt")


## USE MOUSE GENOME FOR RBFOX2:
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/DEL/RBFOX2_BrdU_CLIP_DEL", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/DEL/RBFOX2_BrdU_CLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/TRUNC/RBFOX2_BrdU_CLIP_TRUNC", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/TRUNC/RBFOX2_BrdU_CLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")

## USE MOUSE GENOME FOR RBFOX2 REGULAR CLIP:
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "PEAK/RBFOX2_Reg_CLIP_PEAK", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "PEAK/RBFOX2_Reg_CLIP_PEAK", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/DEL/RBFOX2_Reg_CLIP_DEL", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/DEL/RBFOX2_Reg_CLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CTK/TRUNC/RBFOX2_Reg_CLIP_TRUNC", ".txt"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CTK/TRUNC/RBFOX2_Reg_CLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")
# PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/DEL/RBFOX2_Reg_CLIP_DEL", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/DEL/RBFOX2_Reg_CLIP_DEL", "_NORMALIZED_ANNOTATED.txt")
PEAK_MATRIX_FILE = paste0(INPUT_DIR, "CLINK/TRUNC/RBFOX2_Reg_CLIP_TRUNC", ".bed"); OUTPUT_FILE = paste0(OUTPUT_DIR, "CLINK/TRUNC/RBFOX2_Reg_CLIP_TRUNC", "_NORMALIZED_ANNOTATED.txt")


peaksMatrix = fread(PEAK_MATRIX_FILE)
colnames(peaksMatrix)[1:6] = c("chr", "start", "end", "name", "score", "strand")
################################################################################

## 4. Annotation - Overlap Detection:
################################################################################
## Create GRanges object for peaks
peaksGR = GRanges(seqnames=peaksMatrix$chr, ranges=IRanges(start=peaksMatrix$start, end=peaksMatrix$end, names=peaksMatrix$name), strand=peaksMatrix$strand)

## Overlap with features:
fiveUTRs.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=fiveUTRs,  minoverlap=1, select="first"))
threeUTRs.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=threeUTRs,  minoverlap=1, select="first"))
CDS.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=CDS,  minoverlap=1, select="first"))
introns.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=introns,  minoverlap=1, select="first"))
genes.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=genes,  minoverlap=1, select="first"))
tRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=tRNA.gr,  minoverlap=1, select="first"))
LINE.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LINE.gr,  minoverlap=1, select="first"))
LC_SR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LC_SR.gr,  minoverlap=1, select="first"))
LTR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LTR.gr,  minoverlap=1, select="first"))
Satellite.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=Satellite.gr,  minoverlap=1, select="first"))
SINE.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=SINE.gr,  minoverlap=1, select="first"))
miR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=miR.gr,  minoverlap=1, select="first"))
lncRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=lncRNA.gr,  minoverlap=1, select="first"))
rRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=rRNA.gr,  minoverlap=1, select="first"))
snoRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=snoRNA.gr,  minoverlap=1, select="first"))
scaRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=scaRNA.gr,  minoverlap=1, select="first"))
snRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=snRNA.gr,  minoverlap=1, select="first"))
miscRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=miscRNA.gr,  minoverlap=1, select="first"))
Prot_retained_int.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=Prot_retained_int.gr,  minoverlap=1, select="first"))
nc_retained_int.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=nc_retained_int.gr,  minoverlap=1, select="first"))

## Assign feature data to matrix:

peaksMatrix$fiveUTRs = ifelse(is.na(fiveUTRs.peaks[,1]), NA, "5'UTR")
peaksMatrix$threeUTRs = ifelse(is.na(threeUTRs.peaks[,1]), NA, "3'UTR")
peaksMatrix$CDS = ifelse(is.na(CDS.peaks[,1]), NA, "CDS")
peaksMatrix$introns = ifelse(is.na(introns.peaks[,1]), NA, "intron")
peaksMatrix$gene = genes.bed[genes.peaks[,1], "name"]

peaksMatrix$tRNA1 = ifelse(is.na(tRNA.peaks[,1]), NA, "tRNA")
peaksMatrix$LINE1 = ifelse(is.na(LINE.peaks[,1]), NA, "TE")
peaksMatrix$LTR1 = ifelse(is.na(LTR.peaks[,1]), NA, "Other")
peaksMatrix$LC_SR1 = ifelse(is.na(LC_SR.peaks[,1]), NA, "Other")
peaksMatrix$Satellite1 = ifelse(is.na(Satellite.peaks[,1]), NA, "Other")
peaksMatrix$SINE1 = ifelse(is.na(SINE.peaks[,1]), NA, "TE")
peaksMatrix$miR1 = ifelse(is.na(miR.peaks[,1]), NA, "miRNA")
peaksMatrix$lncRNA1 = ifelse(is.na(lncRNA.peaks[,1]), NA, "lncRNA")
peaksMatrix$rRNA1 = ifelse(is.na(rRNA.peaks[,1]), NA, "rRNA")
peaksMatrix$snoRNA1 = ifelse(is.na(snoRNA.peaks[,1]), NA, "snoRNA")
peaksMatrix$scaRNA1 = ifelse(is.na(scaRNA.peaks[,1]), NA, "scaRNA")
peaksMatrix$snRNA1 = ifelse(is.na(snRNA.peaks[,1]), NA, "snRNA")
peaksMatrix$miscRNA1 = ifelse(is.na(miscRNA.peaks[,1]), NA, "Other")
peaksMatrix$Prot_retained_int1 = ifelse(is.na(Prot_retained_int.peaks[,1]), NA, "CDS_Retained_intron")
peaksMatrix$nc_retained_int1 = ifelse(is.na(nc_retained_int.peaks[,1]), NA, "ncRNA_Retained_intron")
################################################################################

## 5. Annotation - Consolidation:
################################################################################
# Columns to collapse
annoCols = c("fiveUTRs", "threeUTRs", "CDS", "introns", "tRNA1", "LTR1", 
            "LINE1", "SINE1", "Satellite1", "LC_SR1", "miR1", "lncRNA1", 
            "rRNA1", "snoRNA1", "scaRNA1", "snRNA1", "miscRNA1", 
            "Prot_retained_int1", "nc_retained_int1")

# Collapse all overlapping features into a single column using vectorized unite
peaksMatrix = peaksMatrix %>%
  unite("annotation", all_of(annoCols), sep = "|", na.rm = TRUE, remove = FALSE)

# Priority list
priority_list = c("5'UTR", "CDS", "3'UTR", "intron",
                  "miRNA", "lncRNA",  "rRNA", "snoRNA",  "scaRNA",  "snRNA",  "miscRNA",  "tRNA",  "TE",  "Other",  
                  "CDS_Retained_intron", "ncRNA_Retained_intron", 
                  "unannotated")

# Extract the highest priority term for finalized_annotation
peaksMatrix = peaksMatrix %>%
  mutate(finalized_annotation = sapply(strsplit(annotation, "\\|"), function(terms) {
    for (term in priority_list) {
      if (term %in% terms) {
        return(term)
      }
    }
    return(NA)
  }))

# Map terms to grouped_annotation
peaksMatrix = peaksMatrix %>%
  mutate(grouped_annotation = ifelse(finalized_annotation %in% c("5'UTR", "CDS", "3'UTR", "intron", "Other", "unannotated"), finalized_annotation,
                                     ifelse(finalized_annotation %in% c("miRNA", "lncRNA", "rRNA", "scaRNA", "snRNA", "miscRNA", "tRNA", "snoRNA"), "ncRNA",
                                            ifelse(finalized_annotation %in% c("CDS_Retained_intron", "ncRNA_Retained_intron"), "retained_intron", "unannotated"))),
         annotation_count = sapply(strsplit(annotation, "\\|"), length))

peaksMatrix[peaksMatrix == '' | peaksMatrix == 'NA'] = NA

# Cleanup intermediate annotation columns
drops = c('fiveUTRs', 'threeUTRs', 'CDS', 'introns',
          'tRNA1', 'LINE1', 'LTR1', 'LC_SR1', 'Satellite1', 'SINE1',
          'miR1', 'lncRNA1', 'rRNA1', 'snoRNA1', 'scaRNA1', 'snRNA1', 'miscRNA1',
          'Prot_retained_int1', 'nc_retained_int1')
setDT(peaksMatrix); peaksMatrix[, (drops) := NULL]

## Use biomaRt for gene mapping:
gene_names = getBM(attributes = c("ensembl_gene_id_version", "external_gene_name"),
                   filters = "ensembl_gene_id_version",
                   values = peaksMatrix$gene,
                   mart = mart.hs)

peaksMatrix = peaksMatrix %>% left_join(gene_names, by = c("gene" = "ensembl_gene_id_version"), relationship = "many-to-many")
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(is.na(finalized_annotation), 'unannotated', finalized_annotation))
################################################################################

## 6. Save Results:
################################################################################
fwrite(peaksMatrix, OUTPUT_FILE, sep = "\t")
################################################################################
