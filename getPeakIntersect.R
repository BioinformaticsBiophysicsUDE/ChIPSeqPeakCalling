# for merged domains
getIntersectVectNoID <- function(analysis.files, save.name){
  require(data.table)
  require(GenomicRanges)
  options(scipen=999)
  
  nb.files <- length(analysis.files)
  if(nb.files < 2){
    return()
  }
  if(nb.files >= 2){
    # load file1
    file1.analysis.df <- as.data.frame(fread(analysis.files[[1]]))
    names(file1.analysis.df)<- c("chrom", "start", "end")
    #bed files start 100 end 200 means 100 - 199, substract 1 from end for conversion to GRanges
    file1.analysis.df$end <- file1.analysis.df$end - 1
    file1.analysis.Gr <- makeGRangesFromDataFrame(file1.analysis.df, keep.extra.columns =T)
    # load file2
    file2.analysis.df <- as.data.frame(fread(analysis.files[[2]]))
    names(file2.analysis.df)<- c("chrom", "start", "end")
    file2.analysis.df$end <- file2.analysis.df$end - 1
    file2.analysis.Gr <- makeGRangesFromDataFrame(file2.analysis.df, keep.extra.columns =T)
    ## get intersect
    overlap <- intersect(file1.analysis.Gr, file2.analysis.Gr)
  }
  if(nb.files == 3){
    
    file1.analysis.Gr <- overlap
    # load file2
    file2.analysis.df <- as.data.frame(fread(analysis.files[[3]]))
    names(file2.analysis.df)<- c("chrom", "start", "end")
    file2.analysis.df$end <- file2.analysis.df$end - 1
    file2.analysis.Gr <- makeGRangesFromDataFrame(file2.analysis.df, keep.extra.columns =T)
    ## get intersect
    overlap <- intersect(file1.analysis.Gr, file2.analysis.Gr)
  }
  
  overlap.df <- as.data.frame(overlap)
  overlap.df$seqnames <- as.character(overlap.df$seqnames)
  ## save as bed.file
  save.name1<- paste(save.name, ".bed", sep = "")
  write.table(cbind(overlap.df$seqnames, overlap.df$start, overlap.df$end +1), sep = "\t", file = save.name1, quote = F, 
              row.names = F, col.names = F)
  
  options(scipen=0)
  return(overlap)
}
