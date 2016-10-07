##### PROCESSING methods #####

#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (single sample).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param name character. Optional parameter to set query sample name.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"), 
    function(query, ref, anno, name = NULL, intercept = TRUE) {
        if (ncol(query@intensity) == 0) 
            stop("query intensities unavailable, run CNV.load")
        if (ncol(ref@intensity) == 0) 
            stop("reference set intensities unavailable, run CNV.load")
        
        if (ncol(query@intensity) != 1) 
            stop("query contains more than one sample.")
        if (ncol(ref@intensity) == 1) 
            warning("reference set contains only a single sample. use more samples for better results.")
        
        p <- names(anno@probes)  # ordered by location
        if (!all(is.element(p, rownames(query@intensity)))) 
            stop("query intensities not given for all probes.")
        if (!all(is.element(p, rownames(ref@intensity)))) 
            stop("reference set intensities not given for all probes.")
        
        object <- new("CNV.analysis")
        object@date <- date()
        object@fit$args <- list(intercept = intercept)
        
        if (!is.null(name)) {
            names(object) <- name
        } else {
            names(object) <- colnames(query@intensity)
        }
        object@anno <- anno
        
        r <- cor(query@intensity[p, ], ref@intensity[p, ])[1, ] < 0.99
        if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
        if (intercept) {
            ref.fit <- lm(y ~ ., data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, r]))
        } else {
            ref.fit <- lm(y ~ . - 1, data = data.frame(y = query@intensity[p, 
                1], X = ref@intensity[p, r]))
        }
        object@fit$coef <- ref.fit$coefficients
        
        ref.predict <- predict(ref.fit)
        ref.predict[ref.predict < 1] <- 1
        
        object@fit$ratio <- log2(query@intensity[p, 1]/ref.predict[p])
        object@fit$noise <- sqrt(mean((object@fit$ratio[-1] - object@fit$ratio[-length(object@fit$ratio)])^2, 
            na.rm = TRUE))
        
        return(object)
    })


#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin is calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
    standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    
    o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
    o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]], 
        probe = names(object@anno@probes)[o1[, "subjectHits"]], stringsAsFactors = FALSE)
    
    object@bin$ratio <- sapply(split(object@fit$ratio[o2[, "probe"]], o2[, 
        "bin"]), median, na.rm = TRUE)[names(object@anno@bins)]
    object@bin$shift <- optim(0, function(s) median(abs(object@bin$ratio - 
        s), na.rm = TRUE), method = "Brent", lower = -100, upper = 100)$par
    
    return(object)
})


#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
    standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
    if (length(object@fit) == 0) 
        stop("fit unavailable, run CNV.fit")
    # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')
    
    if (length(object@anno@detail) == 0) {
        message("no detail regions provided, define using CNV.create_anno")
    } else {
        d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
        d2 <- data.frame(detail = values(object@anno@detail)$name[d1[, 
            "queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]), 
            stringsAsFactors = FALSE)
        
        object@detail$ratio <- sapply(split(object@fit$ratio[d2[, "probe"]], 
            d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
        object@detail$probes <- table(d2[, 1])[values(object@anno@detail)$name]
    }
    return(object)
})


#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package).
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis. 
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#' 
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#' 
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#' 
#' # general information
#' x
#' show(x)
#' 
#' # coefficients of linear regression
#' coef(x)
#' 
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object, 
    alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo", 
    undo.SD = 2.2, verbose = 0, ...) {
    # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
    if (length(object@bin) == 0) 
        stop("bin unavailable, run CNV.bin")
    # if(length(object@detail) == 0) stop('bin unavailable, run
    # CNV.detail')
    
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))), 
        c("object", "verbose")), function(an) if (is.element(an, names(a2))) 
        a2[[an]] else a1[[an]], simplify = FALSE))
    
    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[names(object@anno@bins)], 
        chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint, 
        data.type = "logratio", sampleid = "sampleid")
    x2 <- DNAcopy::segment(x = x1, verbose = verbose, min.width = min.width, 
        nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, 
        ...)
    object@seg$summary <- DNAcopy::segments.summary(x2)
    object@seg$summary$chrom <- as.vector(object@seg$summary$chrom)  # DNAcopy will factor chrom names. is there another way? 
    object@seg$p <- DNAcopy::segments.p(x2)
    object@seg$p$chrom <- as.vector(object@seg$p$chrom)
    
    return(object)
}) 
