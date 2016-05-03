##### included DATA objects #####

#' tcgaBRCA.sentrix2name
#' @name tcgaBRCA.sentrix2name
#' @description Named vector for Sentrix ID to TCGA ID conversion of breast cancer example data (see README).
#' @details Based on \code{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/jhu-usc.edu/humanmethylation450/methylation/jhu-usc.edu_BRCA.HumanMethylation450.aux.1.8.0/BRCA.mappings.csv}.
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
NULL

#' tbl_ucsc
#' @name tbl_ucsc
#' @description UCSC tables required for creating annotation object.
#' @details Imported using \code{rtracklayer::browserSession('UCSC')}: \code{chromInfo}, \code{gap}, \code{cytoBand}.
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
NULL
# mySession <- rtracklayer::browserSession('UCSC')
# rtracklayer::genome(mySession) <- 'hg19' tbl_ucsc <- list()
# tbl_ucsc$chromInfo <-
# rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, table =
# 'chromInfo')) tbl_ucsc$gap <-
# rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track =
# 'gap', table = 'gap')) tbl_ucsc$cytoBand <-
# rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track =
# 'cytoBand', table = 'cytoBand')) devtools::use_data(tbl_ucsc,
# internal = TRUE)

#' exclude_regions
#' @name exclude_regions
#' @description Example of genomic regions to exclude (e.g. known polymorphic regions).
#' @details Imported using \code{rtracklayer}. Raw data stored in \code{inst/extdata/exclude_regions.bed}.
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
NULL
# exclude_regions <- import('inst/extdata/exclude_regions.bed')
# save(exclude_regions, file = 'data/exclude_regions.rda')

#' detail_regions
#' @name detail_regions
#' @description Example of genomic regions to be analyzed in detail (e.g. candidate oncogenes/TSGs).
#' @details Imported using \code{rtracklayer}. Raw data stored in \code{inst/extdata/detail_regions.bed}.
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
NULL
# detail_regions.gene <- import('inst/extdata/detail_regions.bed')
# detail_regions.promoter <- flank(detail_regions.gene, width = 2000)
# detail_regions <- punion(detail_regions, detail_regions.flank)
# values(detail_regions)$name <- values(detail_regions.gene)$name
# values(detail_regions)$thick <- resize(ranges(detail_regions), 1E6,
# fix = 'center') score(detail_regions) <-
# countOverlaps(detail_regions, anno@probes)
# values(detail_regions)$probes_gene <-
# countOverlaps(detail_regions.gene, anno@probes)
# values(detail_regions)$probes_promoter <-
# countOverlaps(detail_regions.promoter, anno@probes)
# save(detail_regions, file = 'data/detail_regions.rda') 
