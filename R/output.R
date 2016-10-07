##### OUTPUT methods #####

.cumsum0 <- function(x, left = TRUE, right = FALSE, n = NULL) {
    xx <- c(0, cumsum(as.numeric(x)))
    if (!left) 
        xx <- xx[-1]
    if (!right) 
        xx <- head(xx, -1)
    names(xx) <- n
    xx
}

#' CNV.genomeplot
#' @description Create CNV plot for the whole genome or chromosomes.
#' @param object \code{CNV.analysis} object.
#' @param chr character vector. Which chromomsomes to plot. Defaults to \code{'all'}.
#' @param chrX logical. Plot values for chrX? Defaults to \code{TRUE}. Set \code{CNV.create_anno(chrXY = FALSE)} if chrX and Y should not be included at all.
#' @param chrY logical. Plot values for chrY? Defaults to \code{TRUE}.
#' @param centromere logical. Show dashed lines at centromeres? Defaults to \code{TRUE}.
#' @param detail logical. If available, include labels of detail regions? Defaults to \code{TRUE}.
#' @param main character. Title of the plot. Defaults to sample name.
#' @param ylim numeric vector. The y limits of the plot. Defaults to \code{c(-1.25, 1.25)}.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param cols character vector. Colors to use for plotting intensity levels of bins. Centered around 0. Defaults to \code{c('red', 'red', 'lightgrey', 'green', 'green')}.
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return \code{NULL}.
#' @details This method provides the functionality for generating CNV plots for the whole genome or defined chromosomes. Bins are shown as dots, segments are shown as lines. See parameters for more information.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#' 
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.genomeplot", function(object, ...) {
    standardGeneric("CNV.genomeplot")
})

#' @rdname CNV.genomeplot
setMethod("CNV.genomeplot", signature(object = "CNV.analysis"), function(object, 
    chr = "all", chrX = TRUE, chrY = TRUE, centromere = TRUE, detail = TRUE, 
    main = NULL, ylim = c(-1.25, 1.25), set_par = TRUE, cols = c("red", 
        "red", "lightgrey", "green", "green")) {
    # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
    if (length(object@bin) == 0) 
        stop("bin unavailable, run CNV.bin")
    # if(length(object@detail) == 0) stop('bin unavailable, run
    # CNV.detail')
    if (length(object@seg) == 0) 
        stop("bin unavailable, run CNV.seg")
    
    if (set_par) {
        mfrow_original <- par()$mfrow
        mar_original <- par()$mar
        oma_original <- par()$oma
        par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
    }
    
    if (is.null(main)) 
        main <- object@name
    if (chr[1] == "all") {
        chr <- object@anno@genome$chr
    } else {
        chr <- intersect(chr, object@anno@genome$chr)
    }
    chr.cumsum0 <- .cumsum0(object@anno@genome[chr, "size"], n = chr)
    if (!chrX & is.element("chrX", names(chr.cumsum0))) 
        chr.cumsum0["chrX"] <- NA
    if (!chrY & is.element("chrY", names(chr.cumsum0))) 
        chr.cumsum0["chrY"] <- NA
    
    plot(NA, xlim = c(0, sum(as.numeric(object@anno@genome[chr, "size"])) - 
        0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA, 
        ylab = NA, main = main)
    abline(v = .cumsum0(object@anno@genome[chr, "size"], right = TRUE), 
        col = "grey")
    if (centromere) 
        abline(v = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr, 
            "pq"], col = "grey", lty = 2)
    axis(1, at = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr, 
        "size"]/2, labels = object@anno@genome[chr, "chr"], las = 2)
    if (all(ylim == c(-1.25, 1.25))) {
        axis(2, at = round(seq(-1.2, 1.2, 0.4), 1), las = 2)
    } else {
        axis(2, las = 2)
    }
    
    # ratio
    bin.ratio <- object@bin$ratio - object@bin$shift
    bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
    bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
    bin.ratio.cols <- apply(colorRamp(cols)((bin.ratio + max(abs(ylim)))/(2 * 
        max(abs(ylim)))), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    
    lines(chr.cumsum0[as.vector(seqnames(object@anno@bins))] + values(object@anno@bins)$midpoint, 
        bin.ratio, type = "p", pch = 16, cex = 0.75, col = bin.ratio.cols)
    
    for (i in seq(length(object@seg$summary$seg.median))) {
        lines(c(object@seg$summary$loc.start[i] + chr.cumsum0[object@seg$summary$chrom[i]], 
            object@seg$summary$loc.end[i] + chr.cumsum0[object@seg$summary$chrom[i]]), 
            rep(min(ylim[2], max(ylim[1], object@seg$summary$seg.median[i])), 
                2) - object@bin$shift, col = "darkblue", lwd = 2)
    }
    
    # detail
    if (detail & length(object@detail) > 0) {
        detail.ratio <- object@detail$ratio - object@bin$shift
        detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
        detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
        detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) | 
            detail.ratio < -0.85
        
        lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
              + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
            detail.ratio, type = "p", pch = 16, col = "darkblue")
        text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
             + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
            ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", 
                values(object@anno@detail)$name, sep = ""), adj = c(0, 
                0.5), srt = 90, col = "darkblue")
        text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
             + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
            ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, 
                "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "darkblue")
    }
    
    if (set_par) 
        par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

#' CNV.detailplot
#' @description Create CNV plot for detail region.
#' @param object \code{CNV.analysis} object.
#' @param name character. Name of detail region to plot.
#' @param yaxt character. Include y-axis? \code{'l'}: left, \code{'r'}: right, \code{'n'}: no. Defaults to \code{'l'}.
#' @param ylim numeric vector. The y limits of the plot. Defaults to \code{c(-1.25, 1.25)}.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param cols character vector. Colors to use for plotting intensity levels of bins. Centered around 0. Defaults to \code{c('red', 'red', 'lightgrey', 'green', 'green')}.
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return \code{NULL}.
#' @details This method provides the functionality for generating detail regions CNV plots. Probes are shown as dots, bins are shown as lines. See parameters for more information.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#' 
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detailplot", function(object, ...) {
    standardGeneric("CNV.detailplot")
})

#' @rdname CNV.detailplot
setMethod("CNV.detailplot", signature(object = "CNV.analysis"), function(object, 
    name, yaxt = "l", ylim = c(-1.25, 1.25), set_par = TRUE, cols = c("red", 
        "red", "lightgrey", "green", "green")) {
    if (!is.element(name, values(object@anno@detail)$name)) 
        stop("detail_name not in list of detail regions.")
    
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    if (length(object@bin) == 0) 
        stop("bin unavailable, run CNV.bin")
    if (length(object@detail) == 0) 
        stop("bin unavailable, run CNV.detail")
    # if(length(object@seg) == 0) stop('bin unavailable, run CNV.seg')
    
    if (set_par) {
        mfrow_original <- par()$mfrow
        mar_original <- par()$mar
        oma_original <- par()$oma
        par(mfrow = c(1, 1), mar = c(8, 4, 4, 4), oma = c(0, 0, 0, 0))
    }
    
    detail.gene <- object@anno@detail[match(name, values(object@anno@detail)$name)]
    detail.region <- detail.gene
    ranges(detail.region) <- values(detail.gene)$thick
    
    plot(NA, xlim = c(start(detail.region), end(detail.region)), ylim = ylim, 
        xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = values(detail.gene)$name)
    axis(1, at = mean(c(start(detail.region), end(detail.region))), labels = as.vector(seqnames(detail.region)), 
        tick = 0, las = 1)
    axis(1, at = start(detail.region), labels = format(start(detail.region), 
        big.mark = ",", scientific = FALSE), las = 2, padj = 1)
    axis(1, at = end(detail.region), labels = format(end(detail.region), 
        big.mark = ",", scientific = FALSE), las = 2, padj = 0)
    if (yaxt != "n") 
        if (all(ylim == c(-1.25, 1.25))) {
            axis(ifelse(yaxt == "r", 4, 2), at = round(seq(-1.2, 1.2, 0.4), 
                1), las = 2)
        } else {
            axis(ifelse(yaxt == "r", 4, 2), las = 2)
        }
    axis(3, at = c(start(detail.gene), end(detail.gene)), labels = NA)
    
    detail.bins <- names(object@bin$ratio)[as.matrix(findOverlaps(detail.region, 
        object@anno@bins, maxgap = width(detail.region)))[, 2]]
    detail.probes <- names(object@anno@probes)[as.matrix(findOverlaps(detail.region, 
        object@anno@probes, maxgap = width(detail.region)))[, 2]]
    
    detail.ratio <- object@fit$ratio[detail.probes] - object@bin$shift
    detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
    detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
    detail.ratio.cols <- apply(colorRamp(cols)((detail.ratio + max(abs(ylim)))/(2 * 
        max(abs(ylim)))), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    names(detail.ratio.cols) <- names(detail.ratio)
    lines(start(object@anno@probes[detail.probes]), detail.ratio[detail.probes], 
        type = "p", pch = 4, cex = 0.75, col = detail.ratio.cols[detail.probes])
    
    anno.bins.detail <- object@anno@bins[detail.bins]
    anno.bins.ratio <- object@bin$ratio[detail.bins] - object@bin$shift
    anno.bins.ratio[anno.bins.ratio > ylim[2]] <- ylim[2]
    anno.bins.ratio[anno.bins.ratio < ylim[1]] <- ylim[1]
    lines(as.vector(rbind(rep(start(anno.bins.detail), each = 2), rep(end(anno.bins.detail), 
        each = 2))), as.vector(rbind(NA, anno.bins.ratio, anno.bins.ratio, 
        NA)), col = "darkblue", lwd = 2)
    
    if (set_par) 
        par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

#' CNV.detailplot_wrap
#' @description Create CNV plot for all detail regions.
#' @param object \code{CNV.analysis} object.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param main character. Title of the plot. Defaults to sample name.
#' @param ... Additional paramters supplied to \code{CNV.detailplot}.
#' @return \code{NULL}.
#' @details This method is a wrapper of the \code{CNV.detailplot} method to plot all detail regions.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#' 
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detailplot_wrap", function(object, ...) {
    standardGeneric("CNV.detailplot_wrap")
})

#' @rdname CNV.detailplot_wrap
setMethod("CNV.detailplot_wrap", signature(object = "CNV.analysis"), function(object, 
    set_par = TRUE, main = NULL, ...) {
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    if (length(object@bin) == 0) 
        stop("bin unavailable, run CNV.bin")
    if (length(object@detail) == 0) 
        stop("bin unavailable, run CNV.detail")
    # if(length(object@seg) == 0) stop('bin unavailable, run CNV.seg')
    
    if (set_par) {
        mfrow_original <- par()$mfrow
        mar_original <- par()$mar
        oma_original <- par()$oma
        par(mfrow = c(1, length(object@anno@detail) + 2), mar = c(8, 0, 
            4, 0), oma = c(0, 0, 4, 0))
    }
    
    frame()
    for (i in seq(length(object@anno@detail))) {
        if (i == 1) {
            CNV.detailplot(object, name = values(object@anno@detail)$name[i], 
                yaxt = "l", set_par = FALSE, ...)
        } else if (i == length(object@anno@detail)) {
            CNV.detailplot(object, name = values(object@anno@detail)$name[i], 
                yaxt = "r", set_par = FALSE, ...)
        } else {
            CNV.detailplot(object, name = values(object@anno@detail)$name[i], 
                yaxt = "n", set_par = FALSE, ...)
        }
    }
    frame()
    
    if (is.null(main)) 
        main <- object@name
    title(main, outer = TRUE)
    
    if (set_par) 
        par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})


#' CNV.write
#' @description Output CNV analysis results as table.
#' @param object \code{CNV.analysis} object.
#' @param file Path where output file should be written to. Defaults to \code{NULL}: No file is written, table is returned as data.frame object.
#' @param what character. This should be (an unambiguous abbreviation of) one of \code{'probes'}, \code{'bins'}, \code{'detail'} or \code{'segments'}. Defaults to \code{'segments'}.
#' @param ... Additional parameters (\code{CNV.write} generic, currently not used).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#' 
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @return if parameter \code{file} is not supplied, the table is returned as a \code{data.frame} object.
#' @export
setGeneric("CNV.write", function(object, ...) {
    standardGeneric("CNV.write")
})

#' @rdname CNV.write
setMethod("CNV.write", signature(object = "CNV.analysis"), function(object, 
    file = NULL, what = "segments") {
    w <- pmatch(what, c("probes", "bins", "detail", "segments"))
    if (w == 1) {
        if (length(object@fit) == 0) 
            stop("fit unavailable, run CNV.fit")
        if (!is.null(file)) 
            if (!grepl(".igv$", file)) 
                warning("filename does not end in .igv")
        x <- data.frame(Chromosome = as.vector(seqnames(object@anno@probes)), 
            Start = start(object@anno@probes) - 1, End = end(object@anno@probes), 
            Feature = names(object@anno@probes), Value = round(object@fit$ratio - 
                object@bin$shift, 3), row.names = NULL)
        colnames(x) <- sub("Value", object@name, colnames(x))
    } else if (w == 2) {
        if (length(object@bin) == 0) 
            stop("bin unavailable, run CNV.bin")
        if (!is.null(file)) 
            if (!grepl(".igv$", file)) 
                warning("filename does not end in .igv")
        x <- data.frame(Chromosome = as.vector(seqnames(object@anno@bins)), 
            Start = start(object@anno@bins), End = end(object@anno@bins), 
            Feature = names(object@anno@bins), Value = round(object@bin$ratio - 
                object@bin$shift, 3), row.names = NULL)
        colnames(x) <- sub("Value", object@name, colnames(x))
    } else if (w == 3) {
        if (length(object@detail) == 0) 
            stop("detail unavailable, run CNV.bin")
        if (!is.null(file)) 
            if (!grepl(".txt$", file)) 
                warning("filename does not end in .txt")
        x <- data.frame(chr = as.vector(seqnames(object@anno@detail)), 
            start = start(object@anno@detail), end = end(object@anno@detail), 
            name = names(object@detail$probes), sample = object@name, probes = object@detail$probes, 
            value = round(object@detail$ratio - object@bin$shift, 3), row.names = NULL)
    } else if (w == 4) {
        if (length(object@seg) == 0) 
            stop("seg unavailable, run CNV.bin")
        if (!is.null(file)) 
            if (!grepl(".seg$", file)) 
                warning("filename does not end in .seg")
        # seg format, last numeric column is used in igv
        x <- data.frame(ID = object@name, chrom = object@seg$summary$chrom, 
            loc.start = object@seg$summary$loc.start, loc.end = object@seg$summary$loc.end, 
            num.mark = object@seg$summary$num.mark, bstat = object@seg$p$bstat, 
            pval = object@seg$p$pval, seg.mean = round(object@seg$summary$seg.mean - 
                object@bin$shift, 3), seg.median = round(object@seg$summary$seg.median - 
                object@bin$shift, 3), row.names = NULL)
    } else {
        stop("value for what is ambigious.")
    }
    if (is.null(file)) {
        return(x)
    } else {
        write.table(x, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
    }
}) 
