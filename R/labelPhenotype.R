#' Phenotype labeller
#'
#' @description Function that labels the phenotype table generated by \code{\link[CyStore]{getPhenotype}}
#'
#' @param labels character
#'
#' @return list
#'
#' @keywords labels, markers, phenotype, scheme
#' @export
labelPhenotype <- function(phenotype_table){

  label_list <- CyStore::label_list # this should be added as parameter later

  # global parameters
  lb <- setNames(c("dim", "+"), c(1, 2))
  markers <- colnames(phenotype_table)

  # Setup labels with cluster_long
  labels <- setNames(vector("list", length=nrow(phenotype_table)),
                     row.names(phenotype_table))
  for(i in row.names(phenotype_table)){
    x <- phenotype_table[i,]
    labels[[i]]$cluster_long <- paste0(markers[(x + 1) %in% names(lb)], lb[x + 1])
  }
  long_clusters <- lapply(labels, unlist, use.names=F) # save copy


  for(cl in label_list){

    # Positive/negative markers
    sp <- apply(phenotype_table[,markers %in% cl$pos, drop=F] > 0, 1, all)
    sn <- apply(phenotype_table[,markers %in% cl$neg, drop=F] < 0, 1, all)

    # Not positive/negative markers
    snp <- apply(phenotype_table[,markers %in% cl$not_pos, drop=F] != 1, 1, all)
    snn <- apply(phenotype_table[,markers %in% cl$not_neg, drop=F] != -1, 1, all)

    scores <- apply(cbind(sp, sn, snp, snn), 1, all)

    if(!any(scores)){
      message(paste0(cl$label, " not found in phenotype table"))
    }

    for(i in names(scores)){
      if(scores[i]){
        if(!cl$label %in% labels[[i]]$subset){
          labels[[i]]$subset <- c(labels[[i]]$subset, cl$label)
          labels[[i]]$lineage <- c(labels[[i]]$lineage, cl$lineage)
        }
      }
    }

    clusters_in_subset <- names(scores)[scores]
    reduce_markers <- Reduce(intersect, long_clusters[clusters_in_subset])

    for(i in clusters_in_subset){
      middleMarkers <- setdiff(long_clusters[[i]], c(reduce_markers, "CD45"))

      labels[[i]]$cluster_middle <- paste(
        labels[[i]]$subset,
        paste(middleMarkers, collapse=" "),
        sep = " ")
      labels[[i]]$cluster_short <- paste(
        labels[[i]]$subset,
        paste(middleMarkers[grep("[+^]", middleMarkers)], collapse=" "),
        sep = " ")

    }

  }
  labels
}
