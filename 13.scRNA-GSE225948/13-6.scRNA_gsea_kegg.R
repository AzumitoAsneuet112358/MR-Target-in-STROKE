rm(list = ls()) 
options(stringsAsFactors = F)
getOption('timeout')
options(timeout=10000)
options(scipen = 20)
library(enrichplot)
library(ggplot2)
a<-read.csv("comp1_kegg.gsea.csv",row.names = 1)
class(a)
load("gsea_kk.Rdata")
b<-kk@result
class(b)

# up pathways
e<-a[c("hsa05010","hsa04142","hsa01240","hsa05012","hsa05016","hsa05014"),]
f<-b[c("hsa05010","hsa04142","hsa01240","hsa05012","hsa05016","hsa05014"),]
identical(e$enrichmentScore,f$enrichmentScore)
identical(e$NES,f$NES)
kk@result[c("hsa05010","hsa04142","hsa01240","hsa05012","hsa05016","hsa05014"),]<-e
up<-gseaplot2(kk, 
              geneSetID = c("hsa05010","hsa04142","hsa01240","hsa05012","hsa05016","hsa05014"),
              pvalue_table = TRUE)
up
ggsave("figures/gsea_up_01.pdf",
       w=15,
       h=15,
       dpi = 600)

# down pathways
g<-a[c("hsa04727","hsa04728","hsa04724","hsa04360","hsa04151") ,]
h<-b[c("hsa04727","hsa04728","hsa04724","hsa04360","hsa04151") ,]
kk@result[c("hsa04727","hsa04728","hsa04724","hsa04360","hsa04151") ,]<-g
down<-gseaplot2(kk, 
                geneSetID = c("hsa04727","hsa04728","hsa04724","hsa04360","hsa04151"),
                pvalue_table = TRUE)

down
ggsave("figures/gsea_down_01.pdf",
       w=15,
       h=15,
       dpi = 600)

