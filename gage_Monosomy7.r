#!/app/easybuild/software/R/3.5.1-foss-2016b-fh1/bin/Rscript

#Jenny Smith 
#2/20/19

library(methods)

setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2019.02.19_Monosomy7_DEGs')
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/GAGE_GSEA_Function.r")


# mono7vsAML <- readRDS("TARGET_AML_RBD_Monosomy7_vs_OtherAMLs.RDS")
mono7vsNBM <- readRDS("TARGET_AML_RBD_Monosomy7_vs_NBM.RDS")

filename <- "TARGET_AML_1031_Monosomy7_vs_OtherAML_GAGE"
filename2 <- "TARGET_AML_1031_Monosomy7_vs_NBM_GAGE"

C2.KEGG <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.cp.kegg.v6.0.symbols.RDS")



print("starting1")

# GSA <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML,
#               method="voom",
#               type="expn",
#               geneset=C2.KEGG)
#   saveRDS(GSA,file=paste0(filename, "_expn_C2.KEGG.RDS"))
#   rm(GSA)
#   gc()
  
GSA <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM,
                            method="voom",
                            type="expn",
                            geneset=C2.KEGG)
  saveRDS(GSA,file=paste0(filename2, "_expn_C2.KEGG.RDS"))
  rm(GSA)
  gc()

print("done1")


print("starting2")


C2.All <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.all.v6.0.symbols.RDS")

# GSA.C2.All <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML, 
#                      method="voom",
#                      type="expn",
#                      geneset=C2.All)
# saveRDS(GSA.C2.All,file=paste0(filename, "_expn_C2.All.RDS"))
# rm(GSA.C2.All)
# gc()

GSA.C2.All <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM, 
                                 method="voom",
                                 type="expn",
                                 geneset=C2.All)
saveRDS(GSA.C2.All,file=paste0(filename2, "_expn_C2.All.RDS"))
rm(GSA.C2.All)
gc()

print("done2")


print("starting3")

# GSA.KEGG <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML,
#                    method="voom",
#                    type="expn",
#                    geneset=NULL)
#   saveRDS(GSA.KEGG,file=paste0(filename, "_expn_HSA.KEGG.RDS"))
#   rm(GSA.KEGG)
#   gc()

GSA.KEGG <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM,
                                 method="voom",
                                 type="expn",
                                 geneset=NULL)
  saveRDS(GSA.KEGG,file=paste0(filename2, "_expn_HSA.KEGG.RDS"))
  rm(GSA.KEGG)
  gc() 
print("done3")

C5 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c5_list_SetSize_GE.50_LE.300_v6.1.symbols.RDS")

print("starting4")
# GSA.GO.BioProcess <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML,
#                             method="voom",
#                             type="expn", 
#                             geneset=C5[["c5.bp"]])
#   saveRDS(GSA.GO.BioProcess, file=paste0(filename, "_expn_C5.BioProcess_SetSize50to300.RDS"))
#   rm(GSA.GO.BioProcess)
#   gc()

GSA.GO.BioProcess <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM,
                                        method="voom",
                                        type="expn", 
                                        geneset=C5[["c5.bp"]])
saveRDS(GSA.GO.BioProcess, file=paste0(filename2, "_expn_C5.BioProcess_SetSize50to300.RDS"))
rm(GSA.GO.BioProcess)
gc()


print("done4")



print("starting5")
# GSA.GO.CellComp <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML, 
#                           method="voom",
#                           type="expn", 
#                           geneset=C5[["c5.cc"]])
#   saveRDS(GSA.GO.CellComp, file=paste0(filename, "_expn_C5.CellComp_SetSize50to300.RDS"))
#   rm(GSA.GO.CellComp)
#   gc()

GSA.GO.CellComp <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM, 
                                      method="voom",
                                      type="expn", 
                                      geneset=C5[["c5.cc"]])
  saveRDS(GSA.GO.CellComp, file=paste0(filename2, "_expn_C5.CellComp_SetSize50to300.RDS"))
  rm(GSA.GO.CellComp)
  gc()
print("done5")    



print("starting6")
# GSA.GO.MolFunc <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsAML, 
#                          method="voom",
#                          type="expn",
#                          geneset=C5[["c5.mf"]])
#   saveRDS(GSA.GO.MolFunc, file=paste0(filename, "_expn_C5.MolFunc_SetSize50to300.RDS"))
#   rm(GSA.GO.MolFunc)
#   gc()

GSA.GO.MolFunc <- gage_from_pipeline(twoGroups_DEGs.res=mono7vsNBM, 
                                     method="voom",
                                     type="expn",
                                     geneset=C5[["c5.mf"]])
  saveRDS(GSA.GO.MolFunc, file=paste0(filename2, "_expn_C5.MolFunc_SetSize50to300.RDS"))
  rm(GSA.GO.MolFunc)
  gc()
print("done6")










