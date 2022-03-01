# Load Packages
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(tidyverse)
})

#' Index all BAMs in a directory
#'
#' \code{indexBAMs} creates a .bai (index) file for all .bam files in a directory
#' 
#' @param bamdir absolute path to the directory where .bam files are stored
indexBAMs <- function(bamdir) {
  
  # Create list of file paths for the BAMs
  bamfiles <- list.files(bamdir, pattern="*.bam")
  bampaths <- paste0(bamdir, "/", bamfiles)
  
  # Ignoring any previously created index files
  bampaths_nobai <- str_subset(bampaths, ".bai", negate=TRUE)
  baifiles <- bampaths %>% 
    str_subset(".bai") %>% 
    str_remove(".bai")
  
  # Index each BAM file
  i=0
  for (file in bampaths_nobai) {
    if (!(file %in% baifiles)) {
      indexBam(file)
    }
    print(i)
    i=i+1
  }
}

#' Make subclass BAMs
#' 
#' \code{makeSubclassBAMs} creates subclass-level BAM files by merging BAMs in subclass-specific folders
#' 
#' This function assumes there is a directory containing STAR outputs divided by subclass, so \code{bamdir}
#' could be ".../VISp" and it should contain folders like ".../VISp/pvalb_outputs" that were populated using STAR
#' 
#' @param bamdir absolute path to directory containing folders per subclass with STAR outputs
#' @param outputdir absolute path to desired location for merged BAMs
#' @param subclasses a vector of subclass names, should match folder names in bamdir minus the "_outputs" suffix
#' @param txdb a TxDb object created from a reference genome (.gtf or .gff) using GenomicFeatures::makeTxDbFromGFF()
#' @param gene gene name corresponding to chosen reference (i.e., "Zdhhc8" for GRCm39)
#' @param index whether BAMs in the directory need to be indexed
makeSubclassBAMs <- function(bamdir, outputdir, subclasses, txdb, gene, index=FALSE) {
  
  # Get chromosome identifier and gene coordinates from txdb
  genes_txdb <- genes(txdb, single.strand.genes.only = FALSE)
  gene_info <- data.frame(genes_txdb) %>%
    filter(group_name == gene & str_starts(seqnames, "NC"))
  
  chrID <- as.character(gene_info[,3])
  coords <- c(gene_info[,4], gene_info[,5])
  
  # Set up lists of subclass directories and output names
  subclasspaths <- c()
  outputpaths <- c()
  
  for (subclass in subclasses) {
    subclasspaths <- append(subclasspaths, paste0(bamdir, "/", subclass,
                                                  "_outputs/STAR_results/coord_bams"))
    outputpaths <- append(outputpaths, paste0(outputdir, "/", subclass, ".bam"))
  }
  
  for (i in 1:length(subclasspaths)) {
    # Create list of file paths for the BAMs
    bamfiles <- list.files(subclasspaths[i])
    bampaths <- paste0(subclasspaths[i], "/", bamfiles)
    
    # Ignoring any previously created index files
    bampaths_nobai <- str_subset(bampaths, ".bai", negate=TRUE)
    set.seed(12345) # For reproducibility
    bampaths_nobai <- bampaths_nobai %>% sample(size=100)
    
    # If index=TRUE, index each BAM file
    if (index) {
      i=0
      for (file in bampaths_nobai) {
        indexBam(file)
        print(i)
        i=i+1
      }
    }
    
    # Merge all BAMs into one which only contains gene of interest
    mergeBam(files=bampaths_nobai, destination=outputpaths[i], 
             overwrite=TRUE, region = GRanges(chrID, IRanges(coords[1], coords[2])), 
             indexDestination=TRUE)
  }
}

#' Reads to isoform proportions
#'
#' \code{getReads} creates a dataframe of cell type-specific isoform proportions based on read coverage
#' 
#' @param bams a vector of absolute paths to indexed .bam files, assumed to be merged subclass BAMs
#' @param subclasses a vector of subclass names matching the order in \code{bams}
getReads <- function(bams, subclasses) {
  
  # Prepare full reads dataframe
  df <- data.frame(matrix(ncol=2, nrow=0))
  colnames(df) <- c("EV", "OV")
  
  # Set read ranges
  full_range <- GRanges(chrID, IRanges(20145229, 20148007))
  ov_specific_range <- GRanges(chrID, IRanges(20145229, 20147021))
  
  i=1
  for (bam in bams) {
    # Create BamFile object
    bf <- BamFile(bam, asMates=TRUE)
    
    # Get reads per range
    full <- data.frame(union=assays(summarizeOverlaps(full_range, bf, singleEnd=FALSE, fragments=TRUE, 
                                                      param=ScanBamParam(which=full_range)))$counts)
    ov_specific <- data.frame(union=assays(summarizeOverlaps(ov_specific_range, bf, singleEnd=FALSE, 
                                                             fragments=TRUE, 
                                                             param=ScanBamParam(which=ov_specific_range)))$counts)
    
    # Compute reads as percentage of whole final exon
    read_percent_ov <- 100*(rowSums(ov_specific))/(rowSums(full))
    read_percent_ev <- 100-read_percent_ov
    
    rowname <- subclasses[i]
    i=i+1
    
    newrow <- data.frame(subclass=rowname, EV=read_percent_ev, OV=read_percent_ov)
    df <- rbind(df, newrow)
  }
  
  return(df)
}

#' Junctions to isoform proportions
#' 
#' \code{getJunctions} creates a dataframe of cell type-specific isoform proportions based on isoform-specific junctions
#' 
#' @param bams a vector of absolute paths to indexed .bam files, assumed to be merged subclass BAMs
#' @param subclasses a vector of subclass names matching the order in \code{bams}
getJunctions <- function(bams, subclasses) {
  
  # Prepare full junction dataframe
  df <- data.frame(matrix(ncol=2, nrow=0))
  colnames(df) <- c("EV", "OV")
  
  i=1
  for (bam in bams) {
    # Create GAlignments object
    galign <- readGAlignmentPairs(file=bam, index=bam, strandMode=1)
    
    # Get junction info
    junctions <- summarizeJunctions(galign)
    
    # Create dataframe from junction info
    junction_df <- data.frame(c(as.data.frame(junctions@ranges), as.data.frame(junctions$score)))
    
    # Filter dataframe for junctions of interest
    junction_df <- junction_df %>% filter(start==20143757) %>% filter(end==20145228 | end==20147020)
    
    ov <- junction_df %>% filter(end==20145228) %>% .[1,4]
    ev <- junction_df %>% filter(end==20147020) %>% .[1,4]
    
    if (is.na(ov)){
      ov<-0
    }
    if (is.na(ev)){
      ev<-0
    }
    
    junction_percent_ov <- 100*ov/(ev+ov)
    
    if (is.na(junction_percent_ov)){
      junction_percent_ov<-0
    }
    
    junction_percent_ev <- 100*ev/(ev+ov)
    
    if (is.na(junction_percent_ev)){
      junction_percent_ev<-0
    }
    
    rowname <- subclasses[i]
    i=i+1
    
    newrow <- data.frame(subclass=rowname, EV=junction_percent_ev, OV=junction_percent_ov)
    df <- rbind(df, newrow)
  }
  
  return(df)
}