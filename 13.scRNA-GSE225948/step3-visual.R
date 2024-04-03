rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
sce.all.int = readRDS('2-harmony/sce.all_int.rds')
sp='mouse'
colnames(sce.all.int@meta.data) 
table(sce.all.int$RNA_snn_res.0.8)
load('./phe_by_paper.Rdata')
sce.all = sce.all.int
rownames(phe) = colnames(sce.all)
sce.all@meta.data = phe
sel.clust = "celltype"
sce.all <- SetIdent(sce.all, value = sel.clust)
table(sce.all@active.ident) 
DimPlot(sce.all)
colnames(sce.all@meta.data)
DimPlot(sce.all,group.by = 'parent')
 
cg = stringr::str_to_title(c('ALDH16A1','GGCX','SLC33A1'))
cg = cg[cg %in% rownames(sce.all)]
cg
ct = as.data.frame(sce.all@assays$RNA$counts[cg,])
ct = as.data.frame(t(ct ))
head(ct)
table(ct$Aldh16a1 > 0 ,sce.all$tissue)
table(ct$Aldh16a1 > 0 ,sce.all$sex)
table(ct$Aldh16a1 > 0 ,sce.all$age)
kp = sce.all$tissue=='brain'
table(ct$Aldh16a1[kp] > 0 ,sce.all$treatment[kp])
kp = sce.all$tissue !='brain'
table(ct$Aldh16a1[kp] > 0 ,sce.all$treatment[kp])

table(sce.all$tissue,sce.all$celltype)

colnames(sce.all@meta.data)
df = cbind(ct,sce.all@meta.data[16:21])
save(df,file = 'df.Rdata')

colnames(sce.all@meta.data)
for (g in cg) {
  # g=cg[1]
  library(patchwork)
  # p = VlnPlot(sce.all, g, pt.size = 0) + 
  #   VlnPlot(sce.all,g,    group.by = "parent",pt.size = 0) 
  # print(p)
  # ggsave(paste0('all-choose-',g,'-VlnPlot.pdf'),width = 10,height = 8) 
  # 
  # FeaturePlot(sce.all, g ,order = T,raster = T)
  # FeaturePlot(sce.all, g )
  FeaturePlot(sce.all, g,split.by = "treatment")
  ggsave(paste0('FeaturePlot-',g,'-treatment.pdf'),width = 9,height = 6) 
  FeaturePlot(sce.all, g,split.by = "sex")
  ggsave(paste0('FeaturePlot-',g,'-sex.pdf'),width = 7,height = 6) 
  FeaturePlot(sce.all, g,split.by = "age")
  ggsave(paste0('FeaturePlot-',g,'-age.pdf'),width = 12,height = 6) 
  FeaturePlot(sce.all, g,split.by = "tissue")
  ggsave(paste0('FeaturePlot-',g,'-tissue.pdf'),width = 7,height = 6) 
}

genes_to_check = cg
pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce.all, cg,) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap.pdf",width = 16,height = 15)

pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce.all, cg,order = T,raster = T) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap2.pdf",width = 16,height = 15)


# pbmc_small <- BuildClusterTree(object = sce.all.int)
# plot(Tool(object = pbmc_small, slot = 'BuildClusterTree'))
# plot(pbmc_small@tools$BuildClusterTree)