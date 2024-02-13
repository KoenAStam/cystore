#' Draft panel from CyTOF template
#'
#' @description Function that takes a template file used for the CyTOF. It will
#' convert this template to a data frame containing all relevant information.
#' It will also add a column with harmonized markers see
#' \link[CyStore]{harmonizeMarkers.R} for more information
#'
#' @param template character, location of file
#'
#' @importFrom XML xmlToDataFrame xmlParse
#'
#' @return data.frame
#'
#' @keywords panel, markers, labels, template, cytof
#' @export
templateToPanel <- function(template){

  if(!file.exists(template)){
    stop("File not found, check the spelling of the file name.")
  }

  tem <- xmlToDataFrame(xmlParse(template))
  panel <- tem[!is.na(tem$Mass), c("Mass", "Element", "Target")]

  # re-order panel
  panel$Mass <- as.numeric(panel$Mass)
  panel <- panel[order(panel$Mass),]

  # Label empty channels and beads
  panel$empty_channel <- as.numeric(panel$Target == "")
  panel$beads <- as.numeric(grepl("bead", tolower(panel$Target)))

  # Add channel metal
  panel$channel_metal <- paste0(panel$Element, panel$Mass, "Di")

  # Harmonize labels
  panel$channel_label <- harmonizeMarkers(panel$Target)

  return(panel)
}
