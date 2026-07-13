args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 3:4)
version <- args[1]
input <- readRDS(args[2])
if (version == "legacy") {
  .libPaths(c(args[4], .libPaths()))
  library(scDNS)
  create <- scDNS::CreatScDNSobject
  step1 <- scDNS::scDNS_1_CalDivs
  step2 <- scDNS::scDNS_2_creatNEAModel_v2
  step3 <- scDNS::scDNS_3_GeneZscore_v2
  step4 <- scDNS::scDNS_4_scContribution
} else {
  library(spDNS)
  create <- spDNS::CreatScDNSobject
  step1 <- spDNS::scDNS_1_CalDivs
  step2 <- spDNS::scDNS_2_creatNEAModel_v2
  step3 <- spDNS::scDNS_3_GeneZscore_v2
  step4 <- spDNS::scDNS_4_scContribution
}
set.seed(123)
obj <- create(input$counts, input$data, input$network, input$group,
              k = 5, n.grid = 12, n.coarse = 6, parallel.sz = 1,
              n.randNet = 5000)
obj <- step1(obj)
obj <- step2(obj, n.randNet = 5000, repTime = 1)
obj <- step3(obj)
obj <- step4(obj, sigGene = input$genes, q.th = 1)
saveRDS(list(Network = obj@Network, Zscore = obj@Zscore,
             scZscore = obj@scZscore), args[3])
