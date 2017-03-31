##### COERCIONS #####

setAs("CNV.analysis", "GRanges", function(from) {
    res <- sort(subset(as(slot(from, "seg")[["p"]], "GRanges"), !is.na(pval)))
    res$score <- res$seg.mean
    return(res)
})
