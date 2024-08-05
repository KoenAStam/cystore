#' Phenotype labelling
#'
#' @description
#' returns a list with three elements `$phenotype_table`, `$marker_bins` and
#' `$marker_quantiles`.
#'
#' * phenotype_table, a `matrix`. Rows represent the clusterID and columns
#'   represent the markers. Each entry is either 1 (positive), -1 (negative) or 0 (either)
#'
#' * marker_bins, a `list` with two elements, the $background and $clusters specific
#'   bins
#'
#' * marker_quantiles, a `list` with n (number of markers) elements. Each element
#'   is a matrix with 41 rows (quantiles of 2.5%) and p (number of clusters). Thus each
#'   entry is a specific quantile of a specific cluster for the given marker.
#'
#'
#' @param clusterID character vector with identifier for each cell to a cluster
#' @param cellExpr data.frame with numerical cell expression. Rows represents
#' the cells that correspond with the clusterID. Columns represent the markers
#' that will be used for phenotyping.
#'
#'
#' @return a list, check description
#'
#' @details I'll explain the algorithm later
#'
#' @keywords phenotype, markers, labels
#'
#'
#'
#' @export
getPhenotype <- function(clusterID, cellExpr){
    markers <- colnames(cellExpr)
    clusters <- unique(clusterID)

    # for background quantile we don't need more than 100.000 cells
    sampleCells <- sample(nrow(cellExpr), 100000, replace=T)

    # Containers to save data
    phenotype <- matrix(0, nrow=length(clusters), ncol=length(markers),
                        dimnames=list(clusters, markers))
    marker_quantiles <- setNames(vector("list", length=length(markers)), markers)
    marker_bins <- setNames(vector("list", length=length(markers)), markers)


    # Quantile per Marker
    for(i in seq_along(markers)){
        marker <- markers[i]
        message("\rProcessing: ",
                marker, " (", i,"/", length(markers), ")", sep="")

        # split data per cluster
        sv <- split(cellExpr[,marker], clusterID)


        # Binning of expression
        bins <- setNames(vector("list", length=2), c("background", "clusters"))
        markerExpr <- cellExpr[sampleCells, marker]
        binBreaks <- seq(0, ceiling(max(markerExpr)*10), length.out=30)/10
        cutExpr <- cut(markerExpr, breaks=binBreaks, include.lowest = TRUE)
        bins$background <- as.data.frame.table(table(cutExpr) / length(sampleCells))
        bs <- lapply(sv,
                     function(x){
                         cx <- cut(x, breaks=binBreaks, include.lowest=TRUE)
                         table(cx) / length(x)
                     })
        bins$clusters <- do.call(cbind, bs)
        marker_bins[[markers[i]]] <- bins


        # Calculating quantiles
        qs <- lapply(sv, quantile, probs=seq(0, 1, 0.025), names=F)
        qm <- do.call(cbind, qs)
        row.names(qm) <- paste0(seq(0,1,0.025) * 100, "%")
        marker_quantiles[[markers[i]]] <- qm


        # Cluster quantiles in pos and negative
        hclust2 <- cutree(hclust(dist(t(qm))), 2)
        med1 <- mean(qm["50%", hclust2 == 1])
        med2 <- mean(qm["50%", hclust2 == 2])
        posCluster <- which.max(c(med1, med2))
        negCluster <- which.min(c(med1, med2))
        lowCut <- mean(qm["5%", hclust2 == posCluster])
        highCut <- mean(qm["95%", hclust2 == negCluster])


        # set gates
        if(lowCut >= highCut){
            gates <- data.frame(pos = hclust2 == posCluster &
                                    qm["50%",] > lowCut,
                                neg = hclust2 == negCluster &
                                    qm["50%",] < highCut)
        } else if (lowCut <= highCut){
            gates <- data.frame(pos = hclust2 == posCluster |
                                    qm["50%",] > lowCut,
                                neg = hclust2 == negCluster |
                                    qm["50%",] < highCut)
        }

        # set labels
        pos <- numeric(ncol(qm))
        pos <- pos + gates[row.names(phenotype), 1]
        pos <- pos + -gates[row.names(phenotype), 2]
        phenotype[,marker] <- pos
    }

    return(list(phenotype_table = phenotype,
                marker_bins = marker_bins,
                marker_quantiles = marker_quantiles))
}
