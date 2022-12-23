##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
create_graphs <- function() {
  
  fa_graph = readr::read_delim(
    'from  | to 
  #------|----
  Diet & Lipogenesis | 12:0
  12:0 | 14:0
  Diet & Lipogenesis | 16:0
  Diet & Lipogenesis | 18:0
  14:0 | 16:0
  14:0 | 14:1
  14:1 | OMEGA-5
  16:0 | 16:1
  16:1 | 18:1
  18:1 | OMEGA-7
  16:0 | 18:0
  18:0 | 18:1
  18:1 | 20:1
  20:1 | 22:1
  22:1 | 24:1
  24:1 | OMEGA-9
  18:0 | 20:0
  20:0 | 22:0
  22:0 | 24:0
  24:0 | Saturated
  Exogenous | 18:3
  Exogenous | 18:2
  18:3 | 18:4
  18:4 | 20:4
  20:4 | 20:5
  20:5 | 22:5
  22:5 | 22:6
  22:6 | OMEGA-3
  18:2 | 18:3
  18:3 | 20:3
  20:3 | 20:4
  20:4 | 22:4
  22:4 | 22:5
  22:5 | OMEGA-6
  18:2 | 20:2
', delim = '|', trim_ws = TRUE, comment = "#")
  
  fa_graph2 = tidygraph::as_tbl_graph(fa_graph)
  
  sp_graph = readr::read_delim(
    'from  | to 
  #------|----
  Serine + Palmitoyl-CoA | 3SA
  3SA | SA
  SA  | SA1P
  SA  | DCER
  DCER | SA
  DCER | CER
  CER | SM
  CER | SP
  CER | HCER
  SM | CER
  SP | CER
  SP | S1P
  S1P | SP
  S1P | Hydrolysis Pathway
  Hydrolysis Pathway | S1P
  HCER | CER
  HCER | LCER
  LCER | HCER
  LCER | GM
  GM | De Novo Pathway
', delim = '|', trim_ws = TRUE, comment = "#")
  sp_graph2 = tidygraph::as_tbl_graph(sp_graph)
  
  storage_graph = readr::read_delim(
    'from  | to 
  #------|----
  TAG | DAG
  DAG | TAG
  DAG | MAG
  MAG | DAG
  CE  | Membrane CHOL
  Membrane CHOL | CE
  Cyto DAG | TAG
', delim = '|', trim_ws = TRUE, comment = "#")
  
  storage_graph2 = tidygraph::as_tbl_graph(storage_graph)
  
  complex_graph = readr::read_delim(
    'from  | to 
  #------|----
  G3P | LPA
  DHAP | LPA
  LPA | PA
  PA | DAG
  DAG | PC
  PC | LPC
  LPC | PC
  DAG | PE
  PE | LPE
  LPE | PE
  PE | PC
  PE | PS
  PS | LPS
  LPS | PS
  PA | CDP-DAG
  CDP-DAG | PS
  CDP-DAG | PI
  PI | LPI
  LPI | PI
  CDP-DAG | PG
  PG | LPG
  LPG | PG
  PG | CL
  CL | LCL
  LCL | CL
  CDP-DAG | CL
', delim = '|', trim_ws = TRUE, comment = "#")
  complex_graph2 = tidygraph::as_tbl_graph(complex_graph)
  
  list(fatty_acids = fa_graph2,
       sphingolipids = sp_graph2,
       complex_lipid_storage = storage_graph2,
       complex_lipid = complex_graph2)
  
}
