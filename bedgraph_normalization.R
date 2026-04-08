## CoCLIP Analysis: 
## BedGraph Normalization by Mapped Read Depth
## Written by Soon Yi
## Last Edit: 2024-01-06
## Use this after running normalized_bedgraph.sh

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

baseDir = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/normalize_BEDGRAPH/'
bgDir = paste0(baseDir, 'normalized_combined_bedgraph')
setwd(bgDir)
bgFiles = list.files(bgDir)

for (bgFile in bgFiles) {
  out_bg = str_replace(bgFile, '.bedgraph', '.normalized.bedGraph')
  
  expType = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][1]
  
  compartment = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][2]
  if (compartment == 'NLS') {compartment = 'Nuc'} else if (compartment == 'NES') {compartment = 'Cyto'} else if (compartment == 'G3BP') {compartment = 'SG'}
  
  condition = str_split(str_split(bgFile, "\\.")[[1]][1], "_")[[1]][3]
  if (condition == 'Ars') {condition = 'Stress'} else {condition = 'Mock'}
  
  strand = str_split(bgFile, "\\.")[[1]][4]
  if (strand == 'pos') {strand = '(+)'} else {strand = '(-)'}
  
  bg = read.delim(bgFile, header = F)
  bg = bg[2:nrow(bg), ]
  bg = bg %>% mutate(final = (rowSums(bg[, colnames(bg)[4:length(colnames(bg))]]) / length(4:length(colnames(bg)))))
  # bg = bg %>% filter(final != 0)
  bg = bg[, c('V1', 'V2', 'V3', 'final')]
  
  bg_header = paste0('track type=bedGraph name=', expType, '_', compartment, '_', condition, '_', strand)
  
  writeLines(bg_header, paste0(baseDir, '/finished_bedgraph/', out_bg))
  write.table(bg, paste0(baseDir, '/finished_bedgraph/', out_bg), row.names = F, col.names = F, quote = F, append = T)
}

