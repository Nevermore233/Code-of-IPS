rm(list=ls())
options(stringsAsFactors = F)
gc()

#GSE113585
##### 预处理 #####
# ensembl to symbol
gse113585 = read.csv(file = "GSE113585_log2FPKM.csv", header = T, row.names = 1)
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

x = matrix(NA, nrow = nrow(gse113585), ncol = 1)

for (i in 1:nrow(gse113585)) {
  x[i,1] = unlist(strsplit(rownames(gse113585)[i], "[.]"))[1] 
  
}

genes <-x[,1]
a = genes[!duplicated(genes)]
gse113585 = gse113585[!duplicated(genes), ]
rownames(gse113585) = a

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes, mart= mart)
b = G_list[!duplicated(G_list$hgnc_symbol),]
gse113585 = gse113585[b$ensembl_gene_id, ]
rownames(gse113585) = b$hgnc_symbol

log(8,2)
2^log(8,2)
c = 2^gse113585

load(file = "gp18.Rdata")
source(file = "method.R")
library(ggplot2)
library("RColorBrewer")

data = generate_testset2(c, gp18)
data2 = data[,-16]
for (i in 1:nrow(data2)) {
  data2[i, which(data2[i,] == 0)] = 2
}
for (i in 1:nrow(data2)) {
  data2[i, which(data2[i,] == 1)] = 4
}

write.csv(data2, file = "leidatu.csv")
data2 = read.csv(file = "leidatu.csv", header = T, row.names = 1)
for(i in 1:nrow(data2)){
  data2[i,] = as.numeric(unlist(data2[i,]))
}
data2 = as.matrix(data2)

# 画图
#雷达坐标系的定义

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{  theta <- match.arg(theta, c("x", "y"))
r <- if (theta == "x") 
  "y"
else "x"
ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
        direction = sign(direction),
        is_linear = function(coord) TRUE)}
  
# 三个颜色
# brewer.pal(7,"Set1")[1] 1-5
# brewer.pal(7,"Set1")[2] 6-9
# brewer.pal(7,"Set1")[5] 10-16
# brewer.pal(7,"Set1")[7] 17-24

for (i in 17:24) {
  print(i)
  label_data<-data.frame(car=c(1:17),
                         id=c(1:17) ,
                         value=data2[i,])
  
  AddRow<-c(NA,nrow(label_data)+1,label_data[1,ncol(label_data)])
  mydata<-rbind(label_data,AddRow)
  
  myAngle<- 360- 360 * (label_data$id-1) /nrow(label_data)
  pic = ggplot() + 
    geom_polygon(data=mydata,aes(x=id, y=value),color = brewer.pal(7,"Set1")[7], fill=brewer.pal(7,"Set1")[7],alpha=0.06)+
    geom_point(data=mydata,aes(x=id, y=value),size=5,shape=21,color = 'black', fill=brewer.pal(7,"Set1")[7])+
    # coord_polar() + #实现为图3-8-1(c) 的圆形雷达图
    coord_radar()+
    scale_x_continuous(breaks =label_data$id,labels=label_data$car)+
    ylim(0,5)+
    theme_light()+
    theme(axis.text.x=element_text(size = 15,colour="black",angle = myAngle),
          panel.border = element_blank(),# 去除外框
          legend.title=element_blank(),#去除图例
          axis.text.y = element_blank()#去除Y坐标
          )+
    xlab(NULL)+
    ylab(NULL)
  
  ggsave(pic, file = paste("leidatu", i, ".pdf", sep = ""), width=12, height=10) 
}
