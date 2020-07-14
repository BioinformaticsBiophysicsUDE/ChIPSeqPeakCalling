#domain.merger.function
require(data.table)
require(GenomicRanges)
require(parallel)
require(foreach)
library(doParallel)

options(scipen=999)


#' Title
#'
#' @param analysis.bed 
#' @param set 
#'
#' @return
#' @export
#'
#' @examples
domainMergerSeq <- function(analysis.bed, set){
  options(scipen=999)
  if(set == "broad"){
    settings = c(bins = 800, window = 5600, treshold= 4000)
  }
  if(set == "narrow"){
    settings = c(bins = 200, window = 1000, treshold= 800)
  }
  # number of bins in window
  factor <-  settings[2]/settings[1]
  
  # load analysis. bed file
  analysis.df <- as.data.frame(fread(analysis.bed))
  names(analysis.df)<- c("chrom", "start", "end", "id")
  
  # convert data frame to GRanges
  #bed files start 100 end 200 means 100 - 199, substract 1 from end for conversion to GRanges
  analysis.df$end <- analysis.df$end - 1
  analysis.Gr <- makeGRangesFromDataFrame(analysis.df, keep.extra.columns =T)
  
  # name for merged domain bed file
  sample <-  unlist(strsplit(analysis.bed, split = "/"))
  sample <- unlist(strsplit(sample[length(sample)], split = ".bed"))[1]
  bed.file <- paste(sample, ".w", settings[2], ".t", settings[3], ".edgeFilter.bed", sep ="")
  
  ## split GRanges by chromosome
  Gr.list <- split(analysis.Gr, seqnames(analysis.Gr))
  chrom.names <- names(Gr.list)
  
  ### loop over chromosomes
  #k = 1:length(Gr.list)
  final.bed <- foreach(k = 1:length(Gr.list), .combine = "rbind")%do%{
    print(k)
    Gr.chrom <- Gr.list[[k]]
    ranges.chrom <- ranges(Gr.chrom)
    
    ## not the whole chromosomes are covered with domains
    # first position with peak on chromosome
    domain.start <- start(ranges.chrom[1])
    # last position with peak on chromosome
    domain.end <- end(ranges.chrom[length(ranges.chrom)])
    # range on chromosome covered with peaks
    domain.range <- domain.end - domain.start
    # number of sliding windows required to cover the complete range
    window.count <- round(domain.range / settings[1]) - factor +1
    # convert ranges to coverage, 1 peak, 0 no peak for each position
    coverage.chrom <- coverage(ranges.chrom)
    
    new.ranges= c()
    tmp <- foreach(i = 1:window.count)%do%{
      tmp.start <- domain.start + (i -1 )*settings[1]
      tmp.end <- tmp.start + settings[2] -1
      # coverage of window
      tmp.window <- window(coverage.chrom, start=tmp.start, end = tmp.end)
      #+ check if coverage is over specified window
      if(sum(tmp.window) >= settings[3]){
        ## check if bins at edges are empty
        if(tmp.window@values[1] != 0 &&  tmp.window@values[length(tmp.window@values)] != 0){
          ## if windows with empty edges are included, will broaden peaks, if its a broad peak,
          # will cought by next sliding window
          new.ranges <- append(new.ranges, IRanges(tmp.start, tmp.end))
        }
      }
    }
    
    if(length(new.ranges) > 0){
      # overlapping ranges are merged
      new.ranges.reduced <- reduce(new.ranges)
      ### union with original data -> get small isolated peaks
      ranges.union <- union(new.ranges.reduced, ranges.chrom)
      ranges.union.df <- as.data.frame(ranges.union)
      new.bed <- data.frame(chrom = chrom.names[k], start= ranges.union.df$start, end = ranges.union.df$end +1)
    }else{
      ranges.df <- as.data.frame(ranges.chrom)
      new.bed <- data.frame(chrom = chrom.names[k], start= ranges.df$start, end = ranges.df$end +1)
    }
    write.table(new.bed, bed.file, quote = F, sep = "\t", row.names = F, col.names = F, append = T)
    new.bed
  }
  return(final.bed)
}
