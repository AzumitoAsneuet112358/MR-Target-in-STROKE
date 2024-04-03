library(ggplot2)

library(dplyr)

# 数据读入
data <- read.table("data.txt",header=T,sep="\t")

# 自定义颜色
color <- rep("#ae4531", 12)
color[which(data$Score < 0)] <- "#2f73bb"
data$color<-color
# 添加坐标轴标签
up <- data %>% filter(Score>0)
down<- data %>% filter(Score<0)

# 作图
ggplot(data=data,aes(x=reorder(Pathway,Score),Score,fill=color))+

geom_bar(stat = "identity")+

  theme_classic()+

  ylim(-12,12)+

  coord_flip()+

  scale_fill_manual(values = c("#2f73bb","#ae4531"))+

  ggtitle("Basal-like Subtype \n FDR < 0.0001")+

  ylab("Normalized Enrichment Score")+

  #geom_hline(yintercept=0,linetype=1)+

  geom_segment(aes(y = 0, yend = 0,x = 0, xend = 13))+  #geom_segment比genom_line更为灵活，因为可以定义起始位置

  theme(

    legend.position = "none", #消除图例

    plot.title = element_text(hjust = 0.5), #标题居中

    axis.line.y = element_blank(), #删除Y轴

    axis.title.y = element_blank(), #删除Y轴标题

    axis.ticks.y = element_blank(), #删除Y轴刻度

    axis.text.y = element_blank()  #删除Y轴标签

      )+

  geom_text(data=up,aes(x = Pathway, y = 0, label = Pathway), hjust=1,size = 4)+  #添加score是正的label

  geom_text(data=down,aes(x = Pathway, y = 0, label = Pathway), hjust=0,size = 4)+ #添加score是负的label

  geom_segment(aes(y = -1, yend = -10,x = 12.8, xend = 12.8),arrow = arrow(length = unit(0.2, "cm"),type="closed"),size=0.5)+

  geom_segment(aes(y = 1, yend = 10,x = 12.8, xend = 12.8),arrow = arrow(length = unit(0.2, "cm"),type="closed"),size=0.5)+   #添加了2个箭头

  annotate("text", x = 12.8 , y = -12,label = "CP",colour="black",size=5)+  #添加CP和MEP文字注释

  annotate("text", x = 12.8 , y = 12,label = "MEP",colour="black",size=5)+

  theme(text = element_text(size = 15))

ggsave()