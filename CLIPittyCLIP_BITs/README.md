# BITs
Collection of custom <ins>**B**</ins>io<ins>**I**</ins>nformatics <ins>**T**</ins>ools and <ins>**S**</ins>cripts.

Last Update: 2026-04-08

----
## CLIPittyCLIP BITs
- Scripts associated with CLIPittyClip but not necessarily with CoCLIP Analysis.

#### AllowedBarcodeMismatches.py:
  - Python script to check how many mismatches are allowed for barcode demultiplexing.

#### peak_grooming.py:
  - Python script to organize peak matrix from CLIPittyClip pipeline.
  - Renames and reorganizes columns to make the matrix more readable.
  - Calculates biological complexities and total peak heights.
  - Normalizes raw tag counts with total mapped reads (from separate file).

#### bedgraph_normalization.R and .sh
  - Command line and R script for making mapped depth normalized bedGraphs.
  - Requires table containing mapped depth information.
    
#### fastx_barcode_splitter_custom.pl
  - Modified version of fastx_barcode_splitter.pl for demultiplexing.
  - This script was used to generate GEO submitted demultiplexed fastq.gz files.
  - Usage example (for generating demultiplexed fastq.gz files):
      - cat Specificity_Pool.fastq | perl fastx_barcode_splitter_custom.pl --bcfile Specificity_Pool_BC.txt --prefix ~/Specificity/Sequencing/GEO/demux/demux_ --suffix ".fastq" --bol --mismatches 2
      - gzip *.fastq


----
## CoCLIP Analysis
CoCLIP analysis scripts have been moved to standalone repo:
- https://github.com/S00NYI/BITS_CoCLIP
