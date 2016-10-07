##### ANNOTATION methods #####

#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICmanifest
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom rtracklayer import
NULL

#' CNV.create_anno
#' @description Create annotations for CNV analysis.
#' @param bin_minprobes numeric. Minimum number of probes per bin. Bins are interatively merged with neighboring bin until minimum number is reached.
#' @param bin_minsize numeric. Minimum size of a bin.
#' @param bin_maxsize numeric. Maximum size of a bin. Merged bins that are larger are filtered out.
#' @param array_type character. One of \code{450k}, \code{EPIC}, or \code{overlap}. Defaults to \code{450k}.
#' @param chrXY logical. Should chromosome X and Y be included in the analysis?
#' @param exclude_regions GRanges object or path to bed file containing genomic regions to be excluded.
#' @param detail_regions GRanges object or path to bed file containing genomic regions to be examined in detail.
#' @return \code{CNV.anno} object.
#' @details This function collects all annotations required for CNV analysis using Illumina 450k or EPIC arrays. The output \code{CNV.anno} object is not editable. Rerun \code{CNV.create_anno} to change parameters.
#' @examples
#' # create annotation object
#' anno <- CNV.create_anno()
#' anno
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
CNV.create_anno <- function(bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06, 
    array_type = "450k", chrXY = FALSE, exclude_regions = NULL, detail_regions = NULL) {
    object <- new("CNV.anno")
    object@date <- date()
    
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    object@args <- as.list(sapply(unique(names(c(a1, a2))), function(an) if (is.element(an, 
        names(a2))) 
        a2[[an]] else a1[[an]], simplify = FALSE))
    
    if (is.null(array_type)) {
      array_type <- "450k"
    }
    if (!is.element(array_type, c("450k", "EPIC", "overlap"))) {
      stop("array_type must be on of 450k, EPIC, or overlap")
    }
    
    if (chrXY) {
        object@genome <- data.frame(chr = paste("chr", c(1:22, "X", "Y"), 
            sep = ""), stringsAsFactors = FALSE)
    } else {
        object@genome <- data.frame(chr = paste("chr", 1:22, sep = ""), 
            stringsAsFactors = FALSE)
    }
    rownames(object@genome) <- object@genome$chr
    
    message("using genome annotations from UCSC")
    tbl.chromInfo <- tbl_ucsc$chromInfo[match(object@genome$chr, tbl_ucsc$chromInfo$chrom), 
        "size"]
    object@genome$size <- tbl.chromInfo
    
    tbl.gap <- tbl_ucsc$gap[is.element(tbl_ucsc$gap$chrom, object@genome$chr), 
        ]
    object@gap <- sort(GRanges(as.vector(tbl.gap$chrom), IRanges(tbl.gap$chromStart + 
        1, tbl.gap$chromEnd), seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
    
    tbl.cytoBand <- tbl_ucsc$cytoBand[is.element(tbl_ucsc$cytoBand$chrom, 
        object@genome$chr), ]
    # find the gap that is overlapping with the end of the last p-band, use
    # center of that gap for indicating centromers in the genome plots
    pq <- sapply(split(tbl.cytoBand$chromEnd[grepl("p", tbl.cytoBand$name)], 
        as.vector(tbl.cytoBand$chrom[grepl("p", tbl.cytoBand$name)])), 
        max)
    object@genome$pq <- start(resize(subsetByOverlaps(object@gap, GRanges(names(pq), 
        IRanges(pq, pq))), 1, fix = "center"))
    
    probes450k <- probesEPIC <- GRanges()
    if (is.element(array_type, c("450k", "overlap"))) {
      message("getting 450k annotations")
      probes450k <- sort(minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
    }
    if (is.element(array_type, c("EPIC", "overlap"))) {
      message("getting EPIC annotations")
      probesEPIC <- sort(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
    }
    if (array_type == "overlap") {
      probes <- sort(subsetByOverlaps(probes450k, probesEPIC))
    } else {
      probes <- unique(sort(c(probes450k, probesEPIC)))
    }
    
    # CpG probes only
    probes <- probes[substr(names(probes), 1, 2) == "cg" & is.element(as.vector(seqnames(probes)), 
        object@genome$chr)]
    object@probes <- sort(GRanges(as.vector(seqnames(probes)), ranges(probes), 
        seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
    message(" - ", length(object@probes), " probes used")
    
    if (!is.null(exclude_regions)) {
      message("importing regions to exclude from analysis")
      if (class(exclude_regions) == "GRanges") {
        object@exclude <- GRanges(as.vector(seqnames(exclude_regions)), 
                                  ranges(exclude_regions), seqinfo = Seqinfo(object@genome$chr, 
                                                                             object@genome$size))
        values(object@exclude) <- values(exclude_regions)
        object@exclude <- sort(object@exclude)
      } else {
        object@exclude <- sort(rtracklayer::import(exclude_regions, 
                                                   seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
      }
    } else {
      object@exclude <- GRanges(seqinfo = Seqinfo(object@genome$chr, 
                                                  object@genome$size))
    }
    
    if (!is.null(detail_regions)) {
      message("importing regions for detailed analysis")
      if (class(detail_regions) == "GRanges") {
        object@detail <- GRanges(as.vector(seqnames(detail_regions)), 
                                 ranges(detail_regions), seqinfo = Seqinfo(object@genome$chr, 
                                                                           object@genome$size))
        if (any(grepl("name", names(values(detail_regions))))) {
          values(object@detail)$name <- values(detail_regions)[[grep("name", 
                                                                     names(values(detail_regions)))[1]]]
        }
        if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
          values(object@detail)$thick <- values(detail_regions)[[grep("IRanges", 
                                                                      sapply(values(detail_regions), class))[1]]]
        }
        object@detail <- sort(object@detail)
      } else {
        object@detail <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(object@genome$chr, 
                                                                                    object@genome$size)))
      }
      if (!is.element("name", names(values(object@detail)))) {
        stop("detailed region bed file must contain name column.")
      }
      if (!all(table(values(object@detail)$name) == 1)) {
        stop("detailed region names must be unique.")
      }
    } else {
      object@detail <- GRanges(seqinfo = Seqinfo(object@genome$chr, object@genome$size))
    }
    if (!is.element("thick", names(values(object@detail)))) {
      values(object@detail)$thick <- resize(ranges(object@detail), fix = "center", 
                                            1e+06)
    }
    
    message("creating bins")
    anno.tile <- CNV.create_bins(hg19.anno = object@genome, bin_minsize = bin_minsize, 
        hg19.gap = object@gap, hg19.exclude = object@exclude)
    message(" - ", length(anno.tile), " bins created")
    
    message("merging bins")
    object@bins <- CNV.merge_bins(hg19.anno = object@genome, hg19.tile = anno.tile, 
        bin_minprobes = bin_minprobes, hg19.probes = object@probes, bin_maxsize = bin_maxsize)
    message(" - ", length(object@bins), " bins remaining")
    
    return(object)
}

#' CNV.create_bins
#' @description Split genome into bins of defined size.
#' @param hg19.anno foo
#' @param bin_minsize foo
#' @param hg19.gap foo
#' @param hg19.exclude foo
#' @return \code{GRanges} object.
CNV.create_bins <- function(hg19.anno, bin_minsize = 50000, hg19.gap, hg19.exclude) {
    hg19.tile <- sort(tileGenome(Seqinfo(hg19.anno$chr, hg19.anno$size), 
        tilewidth = bin_minsize, cut.last.tile.in.chrom = TRUE))
    # setdiff for gaps (on every second window to avoid merging)
    hg19.tile <- sort(c(setdiff(hg19.tile[seq(1, length(hg19.tile), 2)], 
        hg19.gap), setdiff(hg19.tile[seq(2, length(hg19.tile), 2)], hg19.gap)))
    # setdiff for exluded regions
    hg19.tile <- sort(c(setdiff(hg19.tile[seq(1, length(hg19.tile), 2)], 
        hg19.exclude), setdiff(hg19.tile[seq(2, length(hg19.tile), 2)], 
        hg19.exclude)))
    return(hg19.tile)
}

#' CNV.merge_bins
#' @description Merge bins containing less than the defined number probes with neighboring bin containing fewer probes.
#' @param hg19.anno foo
#' @param hg19.tile foo
#' @param bin_minprobes foo
#' @param hg19.probes foo
#' @param bin_maxsize foo
#' @param verbose foo
#' @return \code{GRanges} object.
CNV.merge_bins <- function(hg19.anno, hg19.tile, bin_minprobes = 20, hg19.probes, 
    bin_maxsize = 5e+06, verbose = FALSE) {
    values(hg19.tile)$probes <- countOverlaps(hg19.tile, hg19.probes)
    hg19.tile.df <- as.data.frame(hg19.tile)[, c("seqnames", "start", "end", 
        "probes")]
    hg19.tile.df$seqnames <- as.vector(hg19.tile.df$seqnames)  # not factor
    
    hg19.tile.df.bin <- do.call(rbind, lapply(split(hg19.tile.df, hg19.tile.df$seqnames), 
        function(mdf) {
            while (min(mdf$probes) < bin_minprobes) {
                mw <- which(mdf$probes == min(mdf$probes))[1]
                mwn <- NA
                mwns <- Inf
                # left
                if (is.element(mdf[mw, "start"] - 1, mdf[, "end"])) {
                  mwn <- mw - 1
                  mwns <- mdf[mw - 1, "probes"]
                  # }
                }
                # right
                if (is.element(mdf[mw, "end"] + 1, mdf[, "start"])) {
                  if (mdf[mw + 1, "probes"] < mwns) {
                    mwn <- mw + 1
                    mwns <- mdf[mw + 1, "probes"]
                  }
                }
                if (is.na(mwn)) {
                  if (verbose) 
                    message(paste(mdf[mw, 1:3], collapse = "-"), " has only ", 
                      mdf[mw, "probes"], " probes and cannot be merged - remove")
                  mdf <- mdf[-mw, ]
                } else {
                  # merge
                  mdf[mwn, "start"] <- min(mdf[c(mwn, mw), "start"])
                  mdf[mwn, "end"] <- max(mdf[c(mwn, mw), "end"])
                  mdf[mwn, "probes"] <- sum(mdf[c(mwn, mw), "probes"])
                  mdf <- mdf[-mw, ]
                }
            }
            return(mdf)
        }))
    
    hg19.tile.bin <- sort(GRanges(hg19.tile.df.bin$seqnames, IRanges(hg19.tile.df.bin$start, 
        hg19.tile.df.bin$end), seqinfo = seqinfo(hg19.tile)))
    hg19.tile.bin <- hg19.tile.bin[width(hg19.tile.bin) <= bin_maxsize]
    
    values(hg19.tile.bin)$probes <- countOverlaps(hg19.tile.bin, hg19.probes)
    values(hg19.tile.bin)$midpoint <- as.integer(start(hg19.tile.bin) + 
        (end(hg19.tile.bin) - start(hg19.tile.bin))/2)
    # values(hg19.tile.bin)$offset <-
    # hg19.anno[as.vector(seqnames(hg19.tile.bin)),
    # 'offset']+values(hg19.tile.bin)$midpoint
    
    names(hg19.tile.bin) <- paste(as.vector(seqnames(hg19.tile.bin)), formatC(unlist(lapply(table(seqnames(hg19.tile.bin)), 
        seq)), width = nchar(max(table(seqnames(hg19.tile.bin)))), format = "d", 
        flag = "0"), sep = "-")
    return(hg19.tile.bin)
} 
