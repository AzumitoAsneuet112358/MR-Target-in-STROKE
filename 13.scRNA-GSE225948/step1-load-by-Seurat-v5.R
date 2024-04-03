###
### Create: Jianming Zeng
### Date:  2023-12-31  
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2023-12-31   First version 
### 


rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
getwd()

###### step1: 导入数据 ######   
# 付费环节 800 元人民币
# 参考：https://mp.weixin.qq.com/s/tw7lygmGDAbpzMTx57VvFw

dir='GSE225948_RAW' 
samples=list.files( dir,pattern = '_counts.csv.gz' )
samples  
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  ct=ct[,-1]
  sce=CreateSeuratObject(counts =  ct ,
                         #project =  gsub('_counts.csv.gz','',strsplit(pro,'_')[[1]][2]),
                         project =  gsub('_counts.csv.gz','',pro),
                         # min.cells = 5,
                         # min.features = 300,
                         )
  
  return(sce)
})
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  gsub('_counts.csv.gz','',samples) ) 
names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )
 
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident) 
sce.all

# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
sp='mouse'
getwd()

###### step3: harmony整合多个单细胞样品 ######
 
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../')
  
} 
###### step4:  看标记基因库 ######
# 原则上分辨率是需要自己肉眼判断，取决于个人经验
# 为了省力，我们直接看  0.1和0.8即可

table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.8) 

getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.5')
setwd('check-by-0.5')
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.8')
setwd('check-by-0.8')
sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

last_markers_to_check


###### step5: 确定单细胞亚群生物学名字 ######
# 一般来说，为了节省工作量，我们选择0.1的分辨率进行命名
# 因为命名这个步骤是纯人工 操作
# 除非0.1确实分群太粗狂了，我们就选择0.8 

## 如果是皮肤
if(F){
  
  # keratinocyte：SFN
  #melanocyte1 "MLANA" 
  selected_genes=c('MUC4', 'PI3', 'SIX3', # nose
                   'SCGB1A1', 'TFF3',   'KRT1',   'KRT10',    
                   'KRT5', 'TP63',   'DLK2',
                   'MKI67', 'TOP2A', 'CDC20' ,
                   'DMD','PUS7','PGAP1','IFI44L',"MLANA" ,'SFN',
                   'MUC5AC' , 'MUC5B',# secretory cells
                   'FOXJ1', 'TPPP3',  'SNTN',  # multiciliated cells 
                   'DEUP1', 'FOXN4',  'CDC20B' # deuterosomal cells
  )
  p1 <- DotPlot(sce.all.int, features = unique(str_to_upper(selected_genes)),
                assay='RNA' ,group.by = 'RNA_snn_res.0.1'  )  + coord_flip() # +th
  
  p1
  
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.8,sce.all.int$orig.ident))
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  
}
## 如果是大脑
if(F){
  
  astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3") 
  endothelial = c("CLDN5", "ABCB1", "EBF1") 
  excitatory = c("CAMK2A", "CBLN2", "LDB2") 
  inhibitory = c("GAD1", "LHFPL3", "PCDH15") 
  microglia = c("C3", "LRMDA", "DOCK8") 
  oligodendrocytes = c("MBP", "PLP1", "ST18")
  # 下面的 OPC是 上面的 oligodendrocytes 的前体细胞 
  OPC='Tnr,Igsf21,Neu4,Gpr17'
  Ependymal='Cfap126,Fam183b,Tmem212,pifo,Tekt1,Dnah12'
  pericyte=c(  'DCN', 'LUM',  'GSN' ,'FGF7','MME', 'ACTA2','RGS5')
  
  gene_list = list(
    Astro = astrocytes,
    Endo = endothelial,
    Excit = excitatory,
    Inhib = inhibitory,
    Mic = microglia,
    Oligo = oligodendrocytes,
    OPC= str_to_upper(trimws(strsplit(OPC,',')[[1]])),
    Ependymal= str_to_upper(trimws(strsplit(Ependymal,',')[[1]])) ,
    peri = pericyte
  )
  gene_list = lapply(gene_list , str_to_title)
  
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.1') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('brain-0.1.pdf',width=12)
  p2 = DotPlot( sce.all.int,  features = gene_list, 
                group.by = 'RNA_snn_res.0.5') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p2    
  ggsave('brain-0.5.pdf',width=12)
  
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.8,sce.all.int$orig.ident))
  gplots::balloonplot(table(sce.all.int$RNA_snn_res.0.1,sce.all.int$orig.ident))
  
  Astrocytes = c("Slc1a2", "Slc1a3","Gpc5", "Prex2","Atp1a2") 
  Epithelial = c("Kitl","Thsd4","Gbp3","Rps8","Tshz2") 
  Oligodendrosytes = c("Mbp","St18","Rnf220","Plp1","Slc24a2") 
  Ependymal = c( "Dnah6","Dnah12","Cfap44","Spag17","Spag16") 
  Oligodendrosyte_progenitor = c("Lhfpl3","Tnr","Dscam","Xylt1","Pcdh15") 
  Contaminating_neurons = c("Meg3","Snhg11","Syt1","Cntnap2","Nrg3") 
  Immune_cells = c( "Ctsd","Hexb","Inpp5d","Ccr5","Cx3cr1")
  Choroid_plexus_epithelial_cells = c("Htr2c","Ttr","Otx2os1","Enpp2","Trpm3")
  Meningeal_cells=c("Fbxl7","Eya2","Bnc2","Cped1","Tmtc1")
  Endothelial =c("Flt1","Slco1a4","Mecom","Adgrl4","Cldn5")
  
  gene_list = list(Astrocytes,Epithelial,Oligodendrosytes,Ependymal,Oligodendrosyte_progenitor,
                   Contaminating_neurons,Immune_cells,Choroid_plexus_epithelial_cells,Meningeal_cells,Endothelial)
 
  names(gene_list)=trimws( strsplit('Astrocytes,Epithelial,Oligodendrosytes,Ependymal,Oligodendrosyte_progenitor,
                   Contaminating_neurons,Immune_cells,Choroid_plexus_epithelial_cells,Meningeal_cells,Endothelial',',')[[1]])
  p_all_markers=DotPlot(sce.all.int, 
                        group.by = 'RNA_snn_res.0.5',
                        features = gene_list,
                        scale = T,assay='RNA' )+  
    theme(axis.text.x=element_text(angle=45,hjust = 1))
  p_all_markers
  ggsave('check_paper_markers_RNA_snn_res.0.5.pdf',
         height = 8,width = 12)
  
}

source('scRNA_scripts/lib.R')
sce.all.int = readRDS('2-harmony/sce.all_int.rds')
sp='mouse'
colnames(sce.all.int@meta.data) 
table(sce.all.int$RNA_snn_res.0.8)
# pbmc_small <- BuildClusterTree(object = sce.all.int)
# plot(Tool(object = pbmc_small, slot = 'BuildClusterTree'))
# plot(pbmc_small@tools$BuildClusterTree)

# 付费环节 800 元人民币
# 如果是手动给各个单细胞亚群命名
if(F){
  sce.all.int
  celltype=data.frame(ClusterID=0:13 ,
                      celltype= 0:13 ) 
  #定义细胞亚群        
  celltype[celltype$ClusterID %in% c( 1 ,7 ),2]='Microglial'  
  celltype[celltype$ClusterID %in% c( 12,2 ),2]='Mono'   
  celltype[celltype$ClusterID %in% c( 3,4,8,11),2]='Mac'   
  celltype[celltype$ClusterID %in% c(5),2]='Tcells' 
  celltype[celltype$ClusterID %in% c( 6 ),2]='endo' 
  celltype[celltype$ClusterID %in% c( 0 ),2]='Bcells'
  celltype[celltype$ClusterID %in% c( 9 ),2]='NK'
  celltype[celltype$ClusterID %in% c( 10 ),2]='cycle'
  
  head(celltype)
  celltype
  table(celltype$celltype)
  sce.all.int@meta.data$celltype = "NA"
  
  for(i in 1:nrow(celltype)){
    sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  Idents(sce.all.int)=sce.all.int$celltype
  

}

# 如果前面成功的给各个细胞亚群命名了
# 就可以运行下面的代码
if("celltype" %in% colnames(sce.all.int@meta.data ) ){
  
  sel.clust = "celltype"
  sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
  table(sce.all.int@active.ident) 
  
  dir.create('check-by-celltype')
  setwd('check-by-celltype')
  source('../scRNA_scripts/check-all-markers.R')
  setwd('../') 
  getwd()
  phe=sce.all.int@meta.data
  save(phe,file = 'phe.Rdata')

}

pdf('celltype-vs-orig.ident.pdf',width = 12)
gplots::balloonplot(table(sce.all.int$celltype,
                          sce.all.int$orig.ident))
dev.off()
 

###### step6: 单细胞亚群比例差异  ######
# 付费环节 800元人民币
###### step7: 单细胞亚群表达量差异分析  ######
# 付费环节 800 元人民币

load(file =  'check-by-0.5/qc-_marker_cosg.Rdata')
head(marker_cosg)
## Top10 genes
library(dplyr)  
cat(paste0('cluster',0:ncol(marker_cosg$names),':',
           unlist(apply(marker_cosg$names,2,function(x){
             paste(head(x),collapse=',')
           })),'\n'))
a=data.table::fread('ACT_Annotation results_top1-by-0.5.txt')
 