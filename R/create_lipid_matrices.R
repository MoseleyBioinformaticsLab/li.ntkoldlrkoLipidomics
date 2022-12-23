##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
create_lipid_matrices <- function() {
  
  fatty_acid_sat = rbind(
    data.frame(saturation = "SAT", values = c(
      paste0(seq(12, 24, by = 2), ":0"))),
    data.frame(saturation = "ODD", values = 
                 c("15:0", "17:0")),
    data.frame(saturation = "MUFA", values = 
                 paste0(seq(14, 24, by = 2), ":1")),
    data.frame(saturation = "PUFA", values = 
                 c("18:2", "20:2", "20:3", "20:4", "22:4", "18:3", "18:4", "20:5", "22:5", "22:6")))
  
  fatty_acid_sat$saturation = factor(fatty_acid_sat$saturation, levels = c("SAT", "ODD", "MUFA", "PUFA"), ordered = TRUE)
  
  fatty_acid_classes = rbind(
    data.frame(
      type = "Membrane", 
      lipid_class = c("PC", "PE", "PI", "LPC", "LPE")),
    data.frame(
      type = "Neutral", 
      lipid_class = c("CE", "TAG", "DAG", "MAG")),
    data.frame(
      type = "Ceramide",
      lipid_class = c("CER", "DCER", "HCER", "LCER", "SM")
    ))
  
  fatty_acid_matrix = create_matrix(fatty_acid_classes$lipid_class, fatty_acid_sat$values)
  
  fatty_acid_classes$type = factor(fatty_acid_classes$type, levels = c("Membrane", "Neutral", "Ceramide"), ordered = TRUE)
  
  sphingolipids_sat = rbind(
    data.frame(saturation = "SAT", values = c(
      paste0(seq(14, 26, by = 2), ":0"))),
    data.frame(saturation = "MUFA", values = c(
      paste0(seq(18, 26, by = 2), ":1")))
  )
  
  sphingolipids_sat$saturation = factor(sphingolipids_sat$saturation, levels = c("SAT", "MUFA"), ordered = TRUE)
  
  sphingolipids_classes = data.frame(type = "",
                                     lipid_class = c("DCER", "CER", "HCER", "LCER", "SM"))
  sphingolipids_matrix = create_matrix(sphingolipids_classes$lipid_class, sphingolipids_sat$values)
  
  complex_lipid_classes = rbind(
    data.frame(
      type = "Polar", 
      lipid_class = c("PC", "PE", "PI", "LPC", "LPE")),
    data.frame(
      type = "Neutral", 
      lipid_class = c("CE", "TAG", "DAG", "MAG")))
  
  complex_lipid_sat = fatty_acid_sat
  
  complex_lipid_matrix = create_matrix(complex_lipid_classes$lipid_class, complex_lipid_sat$values)
  
  complex_lipid_classes$type = factor(complex_lipid_classes$type, level = c("Polar", "Neutral"), ordered = TRUE)
  
  list(fatty_acids = list(
    classes = fatty_acid_classes,
    saturation = fatty_acid_sat,
    matrix = fatty_acid_matrix),
    sphingolipids = list(
      classes = sphingolipids_classes,
      saturation = sphingolipids_sat,
      matrix = sphingolipids_matrix),
    complex_lipids = list(
      classes = complex_lipid_classes,
      saturation = complex_lipid_sat,
      matrix = complex_lipid_matrix))
}

create_matrix = function(classes, sat) {
  lipid_matrix = matrix(NA, nrow = length(classes),
                        ncol = length(sat))
  rownames(lipid_matrix) = classes
  colnames(lipid_matrix) = sat
  lipid_matrix
  
}
