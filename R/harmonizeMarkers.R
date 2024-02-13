#' Harmonize cellular marker labels
#'
#' @description Harmonization of marker labels. Each element of a given character
#' vector will be evaluated and assigned a common label.
#'
#' @param markers character
#'
#'
#'
#' @return character
#'
#' @keywords panel, markers, labels
#' @export
harmonizeMarkers <- function(markers){
    HGNC <- CyStore::HGNC
    LUNC <- CyStore::LUNC

    new_markers <- character(length(markers))

    for(i in seq_along(markers)){
        message("Matching ", markers[i])

        ## First we match to CD symbol database
        m <- match(markers[i], HGNC$cd)

        # No perfect match, check CD* in label
        if(is.na(m)){
            if(grepl("CD[[:alnum:]]{1,7}", markers[i])){
                r <- regmatches(markers[i], regexpr("CD[[:alnum:]]{1,7}", markers[i]))
                m <- match(r, HGNC$cd)
            }
        }

        if(!is.na(m)){
            new_markers[i] <- HGNC$cd[m]
            message("-> ", new_markers[i], "\n", appendLF = FALSE)
            next
        }


        ## Second, match to HGNC symbol
        message("-> HGNC (symbol) ", appendLF=FALSE)
        m2 <- match(toupper(markers[i]), toupper(HGNC$symbol))

        if(!is.na(m2)){
            # add CD label if able else HGNC symbol
            if(HGNC$cd[m2] != ""){
                new_markers[i] <- HGNC$cd[m2]
            } else {
                new_markers[i] <- HGNC$symbol[m2]
            }
            message("-> ", new_markers[i], "\n", appendLF = FALSE)
            next
        }

        ## Lastly, match to department personal database
        message("-> Department ",appendLF = FALSE)
        m3 <- match(toupper(markers[i]), toupper(LUNC$fcs_parameters_desc))

        if(!is.na(m3)){
            new_markers[i] <- LUNC$symbol[m3]
            message("-> ", new_markers[i], "\n", appendLF=FALSE)
            next
        }

        # Still no match throw warning
        message("-> !No match!")
    }

    return(new_markers)
}
