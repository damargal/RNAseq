#' Get count of fragments assgined to each gene
#' 
#' Get raw feature counts assgined to each gene
#' in group of bam files using featureCounts (Rsubread)
#' 
#'    RNA-seq workflow - gene-level exploratory analysis and differential expression
#'    https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#summarizing-an-rna-seq-experiment-as-a-count-matrix
#'    
#' @param filenames A vector with the paths to bam files
#' @param sampleFile A vector with the path to the file with a table (tab-separated) sample information 
#' @param gtffile A verctor with the path to the annotation file
#' @return counts numeric matrix normalized by library size and feature length.
#'

gene_counts <- function(samples=sampleFile,bams=bamfiles,gtf=gtffile){
  require("Rsamtools")
  require("Rsubread")
  
  #load sample table
  sampleTable <- read.csv(samplefile,row.names=1, sep="\t")
  
  #check if files exists
  file.exists(filenames)
  #load bamfiles (Rsamtools)
  bamfiles <- BamFileList(filenames, yieldSize=2000000)
  #check chromosome names
  seqinfo(bamfiles[1])
  
  #create a matrix with read/fragment counts assigned to each gene
  #to include more options for counting reads, such as:
  # library specificity
  # to handgle multiple mapping reads
  # ...
  fc <- featureCounts(files=filenames, 
                      annot.ext=gtffile, 
                      isGTFAnnotationFile=TRUE,
                      isPairedEnd=TRUE)
  
  colnames(fc$counts) <- row.names(sampleTable)
  
  return(fc)
  
}