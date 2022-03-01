# Load Packages
suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(tidyverse)
})

#' Exon boundary coordinates
#'
#' \code{exonCoords} creates a dataframe containing exon start/stop coordinates for each transcript of a gene
#' 
#' @param txdb a TxDb object created from a reference genome (.gtf or .gff) using GenomicFeatures::makeTxDbFromGFF()
#' @param gene gene name corresponding to chosen reference (e.g., "ZDHHC8")
exonCoords <- function(txdb, gene) {
  
  # Get chromosome identifier and gene coordinates from txdb
  genes_txdb <- genes(txdb, single.strand.genes.only = FALSE)
  gene_info <- data.frame(genes_txdb) %>%
    filter(group_name == gene & str_starts(seqnames, "NC"))
  
  chrID <- as.character(gene_info[,3])
  coords <- c(gene_info[,4], gene_info[,5])
  
  # Get dataframe with start/end coords for all mRNA transcripts of the gene 
  gene_tx <- data.frame(transcripts(txdb)) %>%
    filter(seqnames == chrID) %>%
    filter((start >= coords[1]) & (end <= coords[2])) %>%
    filter(str_starts(tx_name, "NM") |
             str_starts(tx_name, "XM")) %>% 
    dplyr::select(tx_name) %>%
    list(.$tx_name) %>%
    .[[2]]
  
  # Get dataframe with start/end coords for all exons in all transcripts
  gr <- exonsBy(txdb, by = "tx", use.names = TRUE)[gene_tx]
  df <- data.frame(gr) %>%
    rename(txID = group_name,
           chrID = seqnames) %>% 
    select(-c(1,8,9))
  
  # Reorder so exons ordered as they appear from left-to-right (for GRanges)
  df <- df %>%
    group_by(txID) %>% 
    arrange(start, .by_group = TRUE) %>% 
    mutate(exon_num = row_number())
  
  # Get number of confirmed transcripts vs. all transcripts
  gene_tx <- df$txID %>%
    unique(.) %>%
    as.data.frame(.) %>%
    filter(str_starts(., pattern = "NM")) %>%
    .$.
  
  confirmed_transcripts <- length(gene_tx)
  
  all_tx <- df$txID %>%
    unique(.) %>%
    length(.)
  
  # Select all confirmed transcripts plus a subset of predicted transcripts
  if (confirmed_transcripts < 2 & all_tx > 1) {
    gene_tx <- df$txID %>%
      unique(.)
    if (length(gene_tx) > 8) {
      gene_tx <- gene_tx %>% 
        .[1:8]
    }
  }
  
  # Subset dataframe to only include transcripts selected above
  df <- df %>%
    filter(txID %in% gene_tx)
  
  return(df)
}

#' Junctions GRanges
#' 
#' \code{junctionsGRanges} creates a GRanges object of junction ranges for use in sashimi plot filter
#'
#' @param exoncoords the dataframe produced using exonCoords()
junctionsGRanges <- function(exoncoords) {
  
  chrID <- as.character(exoncoords[1,2])
  
  rows_count <- nrow(exoncoords)
  boundaries_start <- vector()
  boundaries_end <- vector()
  
  # List all unique pairings of start and stop coordinates for junctions
  i = 1
  while (i < rows_count) {
    if (exoncoords[i, 5] != exoncoords[i + 1, 5]) {
      i = i + 1
    }
    start <- exoncoords[i, 2][[1]] + 1
    end <- exoncoords[i + 1, 1][[1]] - 1
    
    check <- match(start, boundaries_start)
    if (!is.na(check)) {
      if (boundaries_end[check] != end) {
        boundaries_start <- append(boundaries_start, start)
        boundaries_end <- append(boundaries_end, end)
      }
    }
    if (is.na(check)) {
      boundaries_start <- append(boundaries_start, start)
      boundaries_end <- append(boundaries_end, end)
    }
    
    i = i + 1
  }
  
  range <- GRanges(chrID, IRanges(start = boundaries_start, end = boundaries_end))
  
  return(range)
}

#' Coverage and junctions plot
#' 
#' \code{coverageJunctionPlot} creates a plot of read coverage and junctions per provided bam, aligned to isoform models
#'
#' @param txdb a TxDb object created from a reference genome (.gtf or .gff) using GenomicFeatures::makeTxDbFromGFF()
#' @param gene gene name corresponding to chosen reference (e.g., "ZDHHC8")
#' @param bams a vector of maximum length 4 where each item is the absolute path to an indexed .bam file
#' @param bamtitles a vector of maximum length 4 where each item is a title for provided bams (must match order of bams vector)
coverageJunctionPlot <- function(txdb, gene, bams, bamtitles) {
  
  # Create genome axis track
  genomeAxis <-
    GenomeAxisTrack(
      name = "MyAxis",
      col = "lightsteelblue4",
      fontcolor = "lightsteelblue4",
      add35 = TRUE,
      add53 = TRUE
    )
  
  # Get exon coordinates
  exoncoords <- exonCoords(txdb, gene)
  chrID <- as.character(exoncoords[1,2])
  rows_count <- nrow(exoncoords)
  
  # Extract list of transcripts to plot
  gene_tx <- exoncoords$transcript %>%
    unique(.)
  
  tx_num <- length(gene_tx)
  
  # Get junction coordinates (for sashimiFilter)
  if (nrow(exoncoords) > 1) {
    junctions <- junctionsGRanges(exoncoords, chrID)
    plottypes <- c('coverage', 'sashimi')
  } else {
    junctions <- NULL
    plottypes <- 'coverage'
  }
  
  # Prepare plot for each BAM
  i = 1
  bamplots <- vector()
  
  for (bam in bams) {
    if (!is.null(bam)) {
      # Create plot for coverage and sashimi
      bamplot <- AlignmentsTrack(
        bam,
        name = bamtitles[i],
        cex = 2,
        background.title = "white",
        sashimiFilter = junctions,
        sashimiFilterTolerance = 2L,
        col.axis = "lightsteelblue4",
        col.title = "lightsteelblue4",
        type = plottypes
      )
      bamplots <- append(bamplots, bamplot)
      
      # Update index
      i = i + 1
    }
  }
  
  # Create gene models
  gr <- exonsBy(txdb, by = "tx", use.names = TRUE)[gene_tx]
  gr <- unlist(gr)
  elementMetadata(gr)$transcript <- names(gr)
  gene_models <-
    Gviz::GeneRegionTrack(
      gr,
      showId = TRUE,
      options(ucscChromosomeNames = FALSE),
      just.group = "above",
      transcriptAnnotation = "transcript",
      name = "Gene Model",
      background.title = "white",
      col.axis = "lightsteelblue4",
      col.title = "lightsteelblue4",
      fill = "darkgrey",
      fontcolor.group = "lightsteelblue4",
      col.line = "lightsteelblue4"
    )
  
  
  # Indicate tracks to plot
  tracks <- c(genomeAxis)
  for (plot in bamplots) {
    tracks <- append(tracks, plot)
  }
  tracks <- append(tracks, gene_models)
  
  # Select track sizes
  tracksizes <- c(1)
  for (bam in bams) {
    if (!is.null(bam)) {
      tracksizes <- append(tracksizes, 3)
    }
  }
  if (tx_num == 8) {
    tracksizes <- append(tracksizes, 5)
  } else if (tx_num %in% 5:7) {
    tracksizes <- append(tracksizes, 4)
  } else if (tx_num %in% 2:4) {
    tracksizes <- append(tracksizes, 3)
  } else {
    tracksizes <- append(tracksizes, 2)
  }
  
  # Put the figure together
  options(ucscChromosomeNames = FALSE)
  fig <- plotTracks(
    trackList = tracks,
    main = gene,
    showId = TRUE,
    transcriptAnnotation = "transcript",
    chromosome = chrID,
    sizes = tracksizes,
    from = coords[1],
    to = coords[2],
    fill = "lightsteelblue",
    col.sashimi = "lightsteelblue4",
    col.main = "lightsteelblue4"
  )
  
  return(fig)
}