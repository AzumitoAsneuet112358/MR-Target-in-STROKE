rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(ggsci) 
library(patchwork) 
library(ggsci)
library(ggpubr)
library(RColorBrewer) 

getwd()
dir.create('paper-figures/')
setwd('paper-figures/')


## 1. 载入，并且理解数据 ---- 
sce.all.int = readRDS('../2-harmony/sce.all_int.rds')
sp='mouse' 
load('../phe_by_paper.Rdata')
sce.all = sce.all.int
rownames(phe) = colnames(sce.all)
sce.all@meta.data = phe
sel.clust = "celltype"
sce.all <- SetIdent(sce.all, value = sel.clust)
table(sce.all@active.ident) 
DimPlot(sce.all) 

table(Idents(sce.all))  
unique(sce.all$celltype)
length(unique(sce.all$celltype))
sce.all=sce.all[,sce.all$celltype != '13']

ord = c( "cycle" ,  'endo' 
          , 'Microglial',"Mac" ,"Mono","NK" ,'Tcells','Bcells')
all(ord %in% unique(sce.all$celltype))
sce.all$celltype = factor(sce.all$celltype ,levels = ord)
table(sce.all$celltype)
DimPlot(sce.all) 

## 2. 检查内置分群 ----
sce =sce.all
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),paletteer_d("awtools::a_palette"),paletteer_d("awtools::mpalette"))
p2_umap=DimPlot(sce, reduction = "umap",
                # split.by = "orig.ident",
                group.by = "celltype",label = T,label.size = 0,cols = color)
p2_umap
ggsave(p2_umap,filename = "umap_by_celltype_without_label.pdf",
       width = 10,height = 6,units = "cm")


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  
# ord = c( "cycle" ,  'endo' 
#          , 'Microglial',"Mac" ,"Mono","NK" ,'Tcells','Bcells')
selected_genes <-  c(    'MKI67' , 'TOP2A',  
                      'PECAM1', 'VWF',  ## endo 
                     'CD68', 'CD163', 'CD14',  'C1QA',  'C1QB', # mac  # myeloids  
                     'S100A9', 'S100A8', 'MMP19',# monocyte
                     "GNLY","GZMA","NKG7","GZMK", # NK
                     'CD3D', 'CD3E', 'CD4','CD8A',  # Tcells 
                     'CD19', 'CD79A', 'MS4A1' , # Bcells 
                     'IGHG1', 'MZB1', 'SDC1'
                   
) 
selected_genes <-  unique(str_to_title(selected_genes))
p <- DotPlot(sce, features = selected_genes,
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="check_marker_by_celltype.pdf")


## 3. 各式各样的umap  ----
table(sce$tissue)
table(sce$sex)
table(sce$age)
table(sce$treatment)


plot1 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,repel = T)+  
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
plot2 <- DimPlot(sce, reduction = "umap", group.by = "celltype",label = F,cols = color,repel = T) +
  theme(legend.key.size = unit(0.02, 'cm')) + theme(plot.title = element_text(size=12),
                                                    axis.title = element_text(size=10),
                                                    axis.text = element_text(size=10),
                                                    legend.text = element_text(size=10))
plot3 <- DimPlot(sce, reduction = "umap", group.by = "tissue",label = T,cols = sample(color,6),repel = T) +
  theme(plot.title = element_text(size=12),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10))

table(sce$tissue)
table(sce$sex)
table(sce$age)
table(sce$treatment)

p_list = lapply(c('tissue','sex','age','treatment'), function(x){
  DimPlot(sce, reduction = "umap", group.by = x ,
          label = T,cols = sample(color,6),repel = T) +
    theme(plot.title = element_text(size=12),
          axis.title = element_text(size=10),
          axis.text = element_text(size=10),
          legend.text = element_text(size=10))
})
p_list[[1]]
fig1 <- plot1|plot2|plot3
fig1

fig1 = (p_list[[1]]+p_list[[2]])/(p_list[[3]]+p_list[[4]])
ggsave(fig1,filename = "fig1.pdf",units = "cm",width = 22,height = 22)

sce$group = sce$tissue
plot4 <- DimPlot(sce, reduction = "umap", group.by = "celltype",split.by = "group",
                 label = F,cols = color) + theme(plot.title = element_text(size=12),
                                                 axis.title = element_text(size=10),
                                                 axis.text = element_text(size=10))

plot4
ggsave(plot4,filename = "fig2.pdf",units = "cm",width = 22,height = 10)


fig2 <- (plot1|plot2|plot3)/ plot4 + plot_layout(heights = c(1:2))
fig2
ggsave(fig2,filename = "fig1-2.pdf",units = "cm",width = 22,height = 15)


plot5 <- DimPlot(sce, reduction = "umap", group.by = "celltype",
                 split.by = "orig.ident",ncol = 3,label = F,cols = color) 
plot5
ggsave(plot5,filename = "fig2-2.pdf",units = "cm",width = 22,height = 40)



## 4. 已知的基因  ----
pro = 'fig3'
library(dplyr)  
library(paletteer)
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))

cg = c( 'PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
        'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
        'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
        'IFNG', 'CCL4', 'CCL3' ,
        'PRF1' , 'NKG7',  'KLRB1','NCR1', # NK 
        'CD19', 'CD79A', 'MS4A1' , 
        'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
        'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
        'IGHG1', 'MZB1', 'SDC1',
        'CD68', 'CD163', 'CD14', 
        'MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7', 
        'MAF','APOE','FOLR2','RELB','BST2','BATF3',
        'TPSAB1' , 'TPSB2',  # mast cells,
        'RCVRN','FPR1' , 'ITGAM' ,
        'C1QA',  'C1QB',  # mac
        'S100A9', 'S100A8', 'MMP19',# monocyte
        'FCGR3A','JCHAIN',
        'LAMP3', 'IDO1','IDO2',## DC3 
        'CD1E','CD1C','Cd209a','Clec9a', # DC2
        'FGF7','MME', 'ACTA2', ## fibo 
        'DCN', 'LUM',  'GSN' , ##  fibo 
        'PDGFRB', 'PDGFRA',
        'CSPG4','GJB2', 'RGS5','ITGA7',
        'MKI67' , 'TOP2A', 
        'PECAM1', 'VWF',  ## endo 
        'EPCAM' , 'KRT19',  
        'AGER','SFTPA1','SCGB1A1','KRT17','TPPP3',
        'KRT4','KRT14','KRT8','KRT18',
        'PROM1', 'ALDH1A1' )
cg
library(stringr)
genes_to_check = unique(str_to_title(cg))

# 基本上代替热图的小提琴图
table(Idents(sce.all))
p_all_markers=DotPlot(sce.all, 
                      group.by  = "celltype",
                      features = genes_to_check,
                      scale = T,assay='RNA' ) + coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggplot2::ggsave(paste0(pro,'_all_genes_to_check-DotPlot.pdf'),
                height = 15,width = 6)

 
p1 <- VlnPlot(sce.all,features =  selected_genes,
              group.by  = "celltype",
              flip = T,stack = T,cols = color)
p1 + NoLegend()
ggplot2::ggsave(paste0(pro,'_genes_to_check-VlnPlot-heatmap.pdf'),width = 5 )

sce.Scale <- ScaleData(subset(sce.all,downsample=100),features =  selected_genes  )  
DoHeatmap(sce.Scale,
          features =  selected_genes,
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(paste0(pro,filename = "_genes_to_check-pheatmap.pdf") )


## 4. COSG的基因  ---- 
table(Idents(sce))
table(Idents(sce.all))

if(T){
  sce= sce.all
  pro = 'cosg_celltype_'
  library(COSG)
  marker_cosg <- cosg(
    sce,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  # width <-0.006*dim(sce)[2];width
  # height <- 0.25*length(top_10)+4.5;height
  
  width <- 15+0.5*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_10);height
  
  # levels(  Idents(sce) )  = 0:length(levels(  Idents(sce) ))
  DoHeatmap( subset(sce,downsample=100), top_10 , 
             size=3)
  
  ggsave(filename=paste0(pro,'DoHeatmap_check_top10_markers_by_clusters.pdf') ,
         # limitsize = FALSE,
         units = "cm",width=width,height=height)
  width <- 8+0.6*length(unique(Idents(sce)));width
  height <- 8+0.2*length(top_10);height
  DotPlot(sce, features = top_10 ,
          assay='RNA'  )  + coord_flip() +FontSize(y.text = 4)
  ggsave(paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'),
         units = "cm",width=width,height=height)
  
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
  width <- 15+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  
  DoHeatmap( subset(sce,downsample=100), top_3 ,
             size=3)
  ggsave(filename=paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf') ,
         units = "cm",width=width,height=height)
  
  width <- 8+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  DotPlot(sce, features = top_3 ,
          assay='RNA'  )  + coord_flip()
  ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),width=width,height=height)
  
  
}
 

getwd()
pro = 'cosg_celltype_'
load(file = paste0(pro,'_marker_cosg.Rdata'))

top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10))) 
sce.Scale <- ScaleData(subset(sce.all,downsample=100),features =  top_10  )  
table(sce.Scale$celltype)
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
 
table(Idents(sce.Scale)) 

ll =  as.list(as.data.frame(apply(marker_cosg$names,2,head,10)))
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
ll
DoHeatmap(sce.Scale,
          features = unlist(ll),
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(filename = "fig4-top10-marker_pheatmap.pdf",units = "cm",width = 25,height = 25)



### 4.1 . 针对top基因绘制样本分类的小提琴图 ---- 
ll
sce$group = sce$tissue
pdf('fig4-compare_sub_top10_by_VlnPlot.pdf')
lapply(names(ll), function(x){
  cg = ll[[x]]
  VlnPlot(sce,features = cg,group.by = "celltype",split.by = "group",flip = T,stack = T,
          split.plot =T )+scale_fill_d3()
  
})
dev.off()


### 4.2 . 针对top基因绘制FeaturePlot ---- 
pdf('fig4-FeaturePlot-top10.pdf',width = 8,height = 8)
lapply(names(ll), function(x){
  # x=1
  cg = ll[[x]]
  FeaturePlot(sce,features = cg ,ncol = 3,order = T,raster = T ) & 
    theme(plot.title = element_text(size=12),
          axis.title = element_text(size=10),
          axis.text = element_text(size=10)) & 
    scale_colour_continuous(low = 'lightgrey',high = 'red') 
  
})
dev.off()


### 4.3 .. 针对top基因绘制小提琴图 ----
p1 <- VlnPlot(sce ,features =  unique(unlist(ll)),flip = T,
              stack = T,
              group.by = "celltype",
              split.by = "celltype")+
  scale_color_npg()

p1
ggsave(p1,filename = "fig4_Vlnplot_top10.pdf",width = 8,height = 25)




## 5 看比例  ---- 

load(file = '../phe_by_paper.Rdata')
head(phe)
phe$group=phe$tissue
table(phe$group,phe$orig.ident)


### 5.1 批量组合不同变量（分组，样品）看细胞类型比例  ---- 
cal_table = function(x,y,prefix ){
  # x = phe$orig.ident
  # y = phe$celltype
  library(sur)
  library(reshape2)
  tbl =  table(x,y)
  pdf(paste0(prefix,'-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( tbl )
  dev.off() 
  df = dcast(as.data.frame(tbl),x~y)
  head(df)
  write.csv(  df ,paste0(prefix,'-table.csv'))
  
  # ptb = round(sur::percent.table(x,y),2)
  ptb = round(100*tbl/rowSums(df[,-1]),2)
  
  pdf(paste0(prefix,'-percent-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( ptb )
  dev.off()
  write.csv(  dcast(as.data.frame(ptb),x~y) ,paste0(prefix,'-percent-table.csv')) 
  
}
cal_table(phe$orig.ident,phe$celltype,prefix = 'celltype-vs-orig.ident')
cal_table(phe$group,phe$celltype,prefix = 'celltype-vs-group')

cal_table(phe$group,phe$seurat_clusters,prefix = 'seurat_clusters-vs-group')

cal_table(phe$group,phe$RNA_snn_res.0.8,prefix = 'RNA_snn_res.0.8-vs-group')


### 5.2 不同分组的细胞类型比例折线图  ---- 
x='celltype';y='group' 
plot_data <- data.frame(table(phe[, y ],
                              phe[, x ]))
head(plot_data)
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
colnames(plot_data) <- c("group","cluster","Freq","Total","Percentage") 
head(plot_data)
ggplot(plot_data,aes(group,Percentage,group= cluster,color =cluster))+geom_point(size=4)+
  geom_line(position = position_dodge(0.1),cex=2)+theme_bw()+theme_test(base_size = 30) 

ggsave("celltype-vs-group-percent-plot.pdf")

### 5.3 不同分组的细胞类型比例箱线图  ---- 
x = phe$orig.ident
y = phe$celltype
tbl =  table(x,y)
df = dcast(as.data.frame(tbl),x~y)
ptb = round(100*tbl/rowSums(df[,-1]),2)
df = as.data.frame(ptb)
colnames(df)=c('orig.ident','celltype','Percentage')
head(df)
gpinfo=unique(phe[,c('orig.ident','group')]);gpinfo
df=merge(df,gpinfo,by='orig.ident')
head(df)

p <- ggplot(df,aes(x=group,y=Percentage))+
  geom_boxplot(outlier.alpha=0, aes(fill=celltype))+facet_grid(~celltype,scales = "free")+
  geom_jitter(aes(x = group, y = Percentage,color=celltype))+ 
  scale_color_npg() +
  scale_fill_npg() +
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold",angle = 45, hjust=1, vjust=0.5,size = 14), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill = "steelblue3", linetype = "solid"),
        strip.text = element_text(face = "bold", color = "black",hjust = 0.5, size = 12,),
        plot.title=element_text(face="bold.italic",size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(face = "bold",size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p   
pro='Percentage-of-each--celltype-in-orig.ident-by-group'
ggsave(paste0(pro,"_boxplot.pdf"),width = 12,height = 6)




