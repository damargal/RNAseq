#' Get count of fragments assgined to each gene
#' 
#' Get raw feature counts assgined to each gene 
#' in group of bam files using featureCounts (Rsubread)
#' 
#' 
#'    RNA-seq workflow - gene-level exploratory analysis and differential expression
#'    https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#summarizing-an-rna-seq-experiment-as-a-count-matrix
#'  
#' required params:
#' @param sampleFile A vector with the path to the file with a table (tab-separated) sample information 
#' @param filenames A vector with the paths to bam files
#' @param gtf A verctor with the path to the annotation file
#' 
#' @return counts numeric matrix normalized by library size and feature length.
#'

gene_counts <- function(sampleTableFile='',
                        filenames='',
                        gtf='',
                        allowMultiOverlap=FALSE,
                        countMultiMappingReads=FALSE,
                        fraction=FALSE,
                        isPairedEnd=TRUE,
                        minMQS=0,
                        primaryOnly=FALSE,
                        strandSpecific='stranded',
                        nthreads='1',
                        reportReads=FALSE,
                        tmpDir=".", 
                        requireBothEndsMapped=FALSE,
                        checkFragLength=FALSE,
                        minFragLength=50,
                        maxFragLength=600,
                        ignoreDup=FALSE){
  require("Rsamtools")
  require("Rsubread")
  
  #load sample table
  sampleTable <- read.csv(sampleTableFile,row.names=1, sep="\t")
  
  #check if files exists
  file.exists(filenames)
  #load bamfiles (Rsamtools)
  bamfiles <- BamFileList(filenames, yieldSize=2000000)
  #check chromosome names
  seqinfo(bamfiles[1])

  fc <- featureCounts(files=filenames, 
                      annot.ext=gtf, 
                      isGTFAnnotationFile=TRUE,
                      allowMultiOverlap=allowMultiOverlap,
                      countMultiMappingReads=countMultiMappingReads,
                      fraction=fraction,
                      isPairedEnd=isPairedEnd,
                      minMQS=minMQS,
                      primaryOnly=primaryOnly,
                      strandSpecific=strandSpecific,
                      nthreads=nthreads,
                      reportReads=reportReads,
                      tmpDir=tmpDir,
                      requireBothEndsMapped=requireBothEndsMapped,
                      checkFragLength=checkFragLength,
                      minFragLength=minFragLength,
                      maxFragLength=maxFragLength,
                      ignoreDup=ignoreDup)
  
  colnames(fc$counts) <- row.names(sampleTable)
  
  return(fc)
  
}
