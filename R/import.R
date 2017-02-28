##### load a .seg file into a GRanges #####

#' CNV.import
#' 
#' @description import segments where !is.na(pval) from a .seg file into GRanges
#'
#' @param seg       the .seg filename
#' @param genome    genome against which segments were annotated ("hg19")
#' @param name      .seg file column to use as $name metadata ("ID")
#' @param score     .seg file column to use as $score metadata ("seg.median")
#'
#' @return \code{GRanges} class.
#' @details This makes the output of CNV.write easier to merge and annotate.
#' @examples
#' 
#' os <- system.file("extdata", "osteo.hg19.seg", package="conumee")
#' osteo <- CNV.import(os, genome="hg19")
#' 
#' nbl <- system.file("extdata", "neuro.hg19.seg", package="conumee")
#' neuro <- CNV.import(nbl, genome="hg19") 
#' 
#' merged <- sort(c(osteo, neuro))
#' show(merged)
#'
#' # library(Homo.sapiens)
#' # library(VariantAnnotation)
#' # annotated <- locateVariants(merged, Homo.sapiens)
#' # show(annotated)
#'
#' @author Tim Triche, Jr.\email{tim.triche@@gmail.com}
#' @export
CNV.import <- function(seg, genome="hg19", name="ID", score="seg.median") { 
  segs <- subset(read.table(seg, header=TRUE), !is.na(pval))
  names(segs) <- base::sub(name, "name", names(segs))
  names(segs) <- base::sub(score, "score", names(segs))
  firstTwo <- c("score", "name")
  lastTwo <- c("num.mark", "pval")
  segs <- segs[, c(firstTwo, setdiff(names(segs),c(firstTwo,lastTwo)), lastTwo)]
  if (nrow(segs) > 0) {
    gr <- sort(makeGRangesFromDataFrame(segs, keep=TRUE))
    names(gr) <- gr$name
    return(gr)
  } else { 
    return(GRanges())
  }
}
