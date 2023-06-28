#!/usr/local/bin/Rscript --vanilla

#paste0 function
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

#default variables
#bin sizes: 30 and 50kb
binSizesToUse=c(30,50)

#retrieve arguments
args<-commandArgs(TRUE)

#help
help <- function(){
  cat("\nAnalyze shallow whole genome sequencing (sWGS) BAM files using QDNAseq\n")
  cat("\nUsage: analyzeWGS_QDNAseq.R \n")
  cat("-i : input BAM file\n")
  cat("-l : SLX library ID\n")
  cat("-p : sample prefix\n")
  q()
}

#save values of arguments
if(length(args)==0 || !is.na(charmatch("-help",args))){
  help()
  } else {
    for(ii in 1:length(args)){
      if(grepl("^-",args[ii]) && args[ii] != "-"){
        if(ii+1<=length(args) && (!grepl("^-",args[ii+1]) || args[ii+1]=="-")){
          assign(gsub("-","",args[ii]),args[ii+1])
          } else {
            assign(gsub("-","",args[ii]),1)
          }
        }
    }
  }

#check if BAM file exist
if(exists("i")){
  if(file.exists(i)){
  inputBAM=i
  } else {
    cat("Input file does not exist\n")
    q()
    }
  } else {
    cat("No input specified\n")
    q()
  }

#check if user supplied SLX ID
if(exists("l")){
  seqlibname=l
  cat(paste("sample SLX ID: ",seqlibname,"\n"))
  } else {
    cat(paste("sample SLX ID required\n"))
    q()
  }

#check if user supplied sample prefix
if(exists("p")){
  sampleID=p
  cat(paste("sample prefix: ",sampleID,"\n"))
  } else {
    cat(paste("sample prefix required\n"))
    q()
  }

#load library
library(QDNAseq)
library(methods)	#by default, library methods are not loaded in Rscript, although it is in interactive R


#iterate by bin_size
for(bin_size in binSizesToUse)
{

  ############################
  #download pre-created bins
  ############################
  bins <- getBinAnnotations(binSize=bin_size)

  ############################
  #     load bam file
  ############################
  readCounts <- binReadCounts(bins,bamfiles=inputBAM)

  ############################
  #plot raw read counts per bin
  ############################
  png(file=paste0("/XXXX/QDNA_results/",seqlibname,"/RawReadsPerBin/",seqlibname,"_",sampleID,"_",bin_size,"kb.png"))
  plot(readCounts,logTransform=FALSE,ylim=c(-50, 400))
  dev.off()

  ############################
  #apply filters
  ############################
  readCountsFiltered <- applyFilters(readCounts,chromosomes='Y')

  ############################
  #plot isobarPlot
  ############################
  png(file=paste0("/XXXX/QDNA_results/",seqlibname,"/IsoBarPlot/",seqlibname,"_",sampleID,"_",bin_size,"kb.png"))
  isobarPlot(readCountsFiltered)
  dev.off()

  ###############################################
  #estimate for GC content and mappability
  ###############################################
  readCountsFiltered <- estimateCorrection(readCountsFiltered)

  ############################
  #plot noise plot
  ############################
  png(file=paste0("/XXXX/QDNA_results/",seqlibname,"/NoisePlot/",seqlibname,"_",sampleID,"_",bin_size,"kb.png"))
  noisePlot(readCountsFiltered)
  dev.off()

  ###############################################
  #corect for GC content and mappability
  ###############################################
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

  ###############################################
  #   plot copy numbers and export
  ###############################################
  png(file=paste0("/XXXX/QDNA_results/",seqlibname,"/CopyNumbers/",seqlibname,"_",sampleID,"_",bin_size,"kb.png"))
  plot(copyNumbersSmooth)
  dev.off()

  ###############################################
  #     segment, call CNA and plot
  ###############################################
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth,transformFun="sqrt",smoothBy=1L)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

  png(file=paste0("/XXXX/QDNA_results/",seqlibname,"/CopyNumbersSegmented/",seqlibname,"_",sampleID,"_",bin_size,"kb.png"))
  plot(copyNumbersSegmented)
  dev.off()

  ###############################################
  #call aberrations and export
  ###############################################
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  exportBins(copyNumbersCalled, file=paste0("/XXXX/sWGS/QDNA_results/",seqlibname,"/copyNumbersCalled/igv/",seqlibname,"_",sampleID,"_",bin_size,"kb.igv"), format="igv")
  exportBins(copyNumbersCalled, file=paste0("/XXXX/sWGS/QDNA_results/",seqlibname,"/copyNumbersCalled/bed/",seqlibname,"_",sampleID,"_",bin_size,"kb.bed"), format="bed")

}


cat("DONE!!!\n")




#######
