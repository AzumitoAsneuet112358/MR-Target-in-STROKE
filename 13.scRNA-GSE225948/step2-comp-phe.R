dir='GSE225948_RAW/'
samples=list.files( dir ,pattern = 'metadata')
samples 
ctList = lapply(samples,function(pro){ 
  # pro=samples[7] 
  print(pro)  
  ct <- data.table::fread( file.path(dir,pro),
                           data.table = F)
  ct[1:4,1:4] 
  return(ct)
}) 
do.call(rbind,lapply(ctList, dim))
phe2 = do.call(rbind,ctList)
table(phe2$parent)
rownames(phe2)=phe2$V1

load(file = 'phe.Rdata')
table(phe$celltype)

head(rownames(phe))
head(rownames(phe2))
tmp = stringr::str_split( gsub('aged_','',rownames(phe) ),'[-_]',simplify = T)
rownames(phe)=tmp[,4]
rownames(phe)=gsub('[.][0-9]','',rownames(phe))

ids = intersect(rownames(phe),rownames(phe2))
dim(phe);
dim(phe2);
length(ids)
g1 = phe[ids,'celltype']
g2 = phe2[ids,'parent']
 
gplots::balloonplot(
  table(g1 ,g2 )
)

tmp = phe2[match(rownames(phe),rownames(phe2)),]
colnames(tmp)
phe = cbind(phe,tmp[,c("tissue","sex","age" ,"treatment" ,"parent" )])
save(phe,file = 'phe_by_paper.Rdata')


