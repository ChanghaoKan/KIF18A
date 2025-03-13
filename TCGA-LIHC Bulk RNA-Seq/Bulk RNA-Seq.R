####准备工作####
#设置工作目录
setwd("xena")
#快捷键
#Ctrl+s 保存
#Ctrl+z 返回
#Ctrl+f 搜索
#alt+-  <- 
#ctrl + shift + m %>%           

#更新r
install.packages("installr")
library(installr)
updateR()
#安装包与加载包
#install.packages("tidyverse")#""中间加入要安装包的名称
library(tidyverse)#每次重新打开R都要library一下——运行包

#对照tegf肿瘤缩写表下载对应肿瘤数据库
#xena官网 https://xenabrowser.net/datapages/

#文件的读取
#读取TSV文件仅替换单引号中间内容
counts1 = read.table(file = 'TCGA-LIHC.htseq_counts.tsv', sep = '\t',
                     header = TRUE) 
#将文件counts第一行作为的rownames（行名） #Alt - 为 <- 
rownames(counts1) <- counts1[,1] 
#将counts1变为不带第一行的counts1 #-1不要第一行
counts1 = counts1[,-1]
#substr函数——substr(s, first, last）first和last为截取的起始位置
#table函数——table #colnames列名，14，16位为01A 01B 02A 02B 11A
table(substr(colnames(counts1),14,16))
#c("01A","11A")——集合包含01A和11A——01A肿瘤，11A正常
#%in%符号用于判断是否属于
counts1 <- counts1[,substr(colnames(counts1),14,16)%in% c("01A","11A")]
#展示counts1中列名14到16位  
table(substr(colnames(counts1),14,16))
#保留行名前15位（前15为基因名）
rownames(counts1) <- substr(rownames(counts1),1,15)
#ceiling——取整函数，根据xena中表格数据为log2x-1逆推x
counts <- ceiling(2^(counts1)-1)
#文件的输出
#输出为文本
write.table(counts,"counts_LIHC.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#输出为表格
write.csv(counts, file = "counts_LIHC.csv")








####差异分析####
#读取文本文件
counts <- read.table("counts_LIHC.txt",sep = "\t",row.names = 1,check.names = F,
                     stringsAsFactors = F,header = T)
#加载基因注释文件
Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,
                      stringsAsFactors = F,header = T,row.names = 1)
#美元符号$代表提取列
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交集 intersect函数 <- 取交集
comgene <- intersect(rownames(counts),rownames(Ginfo)) 
counts <- counts[comgene,]
#class()判断数据类型
class(counts)   #data <- 数据
class(comgene)  #character <-字符串 
#将两个表格数据顺序内容进行比较，相同则输出true，便于更改基因名称  
Ginfo <- Ginfo[comgene,]#（x,y）x <- 行，y <- 列  #按照comgene的顺序提取Ginfo
a <- rownames(counts)
b <- rownames(Ginfo)
identical(a,b)
#新增Gene Symbol
counts$Gene <- as.character(Ginfo$genename)   #存储为字符串as.character
#去重复!duplicated
counts <- counts[!duplicated(counts$Gene),] 
#将行名变为Gene Symbol  
rownames(counts) <- counts$Gene   
#计算行列各有多少
ncol(Ginfo)#——col——列
nrow(Ginfo)#——row——行
#去除最后一列  
counts <- counts[,-ncol(counts)]   
#保存
write.table(counts, file = "LIHC_counts_mRNA_all.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#读取counts文件
counts <- read.table("LIHC_counts_mRNA_all.txt")
#保存癌症患者的counts
tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <- counts[,tumor]
write.table(counts_01A, file = "LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#差异分析
library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
#分组
counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == 
                                            "01A","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")
#跑差异分析的代码
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
#确认肿瘤与正常组织对比
resultsNames(dds)
#输出结果  
res <- results(dds)
#一定要保存  
save(res,file = "LIHC_DEG.rda")#一定要保存！
#将所有的基因的差异分析进行保存
exprSet <- as.data.frame(res)

write.csv(exprSet, file = "LIHC_deseq2_all.csv")
#进行差异分析
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)#根据自己需要,abs绝对值在几倍

#DEG:differentially expressed genes 差异表达基因
#读取差异基因文件
#设置工作目录！！！
setwd("xena")
#读取包
library(tidyverse)
library(DESeq2)
#直接在文件夹双击——报错？

#通过load函数读取rda文件  
load("LIHC_DEG.rda")
#读取后筛选差异基因
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)#根据自己需要,abs绝对值在几倍
#输出为文本
write.table(res_deseq2,"LIHC_abs1.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#输出为表格
write.csv(res_deseq2, file = "LIHC_abs1.csv")




#加载包
BiocManager::install("clusterProfiler") 
BiocManager::install("enrichplot")  #画图
BiocManager::install("org.Hs.eg.db") #基因注释
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
#读取差异分析结果
df = read.table("LIHC_abs2.txt",header = T) #读入txt
# df = read.csv("gene_diff.csv",header = T) #读入csv
#将行名作为SYMBOL一列导入矩阵  
df <- df %>% rownames_to_column("SYMBOL")
#提取symbol与logfc
df1 <- df %>% select(1, 3)
#进行转换
df_id<-bitr(df1$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
#合并数据集id与SYMBOL
df_all<-merge(df1,df_id,by="SYMBOL",all=F)#使用merge合并
#按照log2FoldChange降序排序
df_all_sort <- df_all[order(df_all$log2FoldChange, decreasing = T),]
#按顺序提取log2FoldChange这列，并加上对应的ENTREZID 
gene_fc = df_all_sort$log2FoldChange
names(gene_fc) <- df_all_sort$ENTREZID 
#KEGG Pathway
KEGG <- gseKEGG(gene_fc, organism = "hsa")
#按照enrichment score从高到低排序  
sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]
#保存结果 
write.table(sortKEGG,"gsea_sortKEGG.txt") 
#作图
gseaplot2(
  x, #gseaResult object，即GSEA结果
  core_enrichment,#富集的ID编号
  title = "", #标题
  color = "green",#GSEA线条颜色
  base_size = 11,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = FALSE, #是否添加 pvalue table
  ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
)
gseaplot2(KEGG, "hsa04657", color = "firebrick", rel_heights=c(1, .2, .6))
paths <- c("hsa04657", "hsa04976", "hsa00980")#选取你需要展示的通路ID
gseaplot2(KEGG,paths, pvalue_table = TRUE)

#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
#保存文件  
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_ICA1.rda")#以ICA1为差异作图
#单个图绘制
library(enrichplot)
gseaplot2(KEGG,1,color="red")
gseaplot2(KEGG,3,color="red",pvalue_table = T)




####GSEA####
#加载包
library(ggplot2)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), 
       function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
if(!require(GSVA))BiocManager::install('GSVA')
library(GSVA)
library(GSEABase)
#读取差异分析文件 <- LIHC_abs2.txt  
deg=read.table("LIHC_abs1.txt",sep="\t",header=T,check.names=F)
#设置分组值
logFC_t=0.5#存疑，但无较大影响
deg <-DEG 
#单设一列分为高低组  
deg$g=ifelse(deg$pvalue>0.05,'stable',
             ifelse( deg$log2FoldChange > logFC_t,'UP',
                     ifelse( deg$log2FoldChange < -logFC_t,'DOWN','stable') )
)
#查看高低组各有多少  
table(deg$g)
#提取第一列并命名为symbol
deg$symbol=deg[,1]
deg$symbol=rownames(deg)
#进行转换得到数据集id
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
#提取原有SYMBOL列表  
DEG=deg
#使用merge函数合并数据集id与SYMBOL  
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
#按照log2FoldChange降序排序
data_all_sort <- DEG %>% 
  arrange(desc(log2FoldChange))
#把foldchange按照从大到小提取出来
geneList = data_all_sort$log2FoldChange
#给上面提取的foldchange加上对应上ENTREZID
names(geneList) <- data_all_sort$ENTREZID
#查看genelist  
head(geneList)
#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 2,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
#run！！！
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
#af为最后gsea富集结果
af=as.data.frame(kk2@result)
#输出并保存  
write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)
#排序后分别取GSEA结果的前5个和后5个
num=10
#保存为pdf  
pdf(paste0("2.","down_GSEA.pdf"),width = 8,height = 8)
#绘图 <- 下调
gseaplot2(kk2, geneSetID = rownames(kk2@result)
          [head(order(kk2@result$enrichmentScore),num)])
#关闭图片  
dev.off()``
#保存为pdf
pdf(paste0("2.","up_GSEA.pdf"),width = 8,height = 8)
#绘图 <- 上调
gseaplot2(kk2, geneSetID = rownames(kk2@result)
          [tail(order(kk2@result$enrichmentScore),num)])
#关闭图片  
dev.off()
#排序后取前5个和后5个一起展示
num=5
#保存为pdf
pdf(paste0("2.","all_GSEA.pdf"),width = 10,height = 10)
#绘图
gseaplot2(kk2, geneSetID = rownames(kk2@result)
          [c(head(order(kk2@result$enrichmentScore),num),
             tail(order(kk2@result$enrichmentScore),num))])
#关闭图片  
dev.off()
#单独展示,自行保存
gseaplot2(kk2,
          title = "name",  #设置标题
          "hsa04721", #绘制hsa04658通路的结果，通路名称与编号对应
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, 
          pvalue_table = T) # 显示p值
#山脊图，展示10个，自行保存
library(stringr)
kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
ridgeplot(kk2,showCategory = 20)




####整理fpkm文件####
#与counts几乎相同，fpkm不需进行log转换
#counts数据用来差异分析，fpkm用来后续分析
#加载包
library(tidyverse)
#设置工作目录
setwd("xena")
#读取tsv文件
fpkm1 = read.table(file = 'TCGA-LIHC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
#提取第一行为行名
rownames(fpkm1) <- fpkm1[,1]
#去除第一列
fpkm1 = fpkm1[,-1]
#展示每一列14~16 即01A等数量
table(substr(colnames(fpkm1),14,16))
#提取01A，11A
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
#展示01A，11A数量  
table(substr(colnames(fpkm1),14,16))
#只要行名1~15位
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1
#加载基因注释文件
Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,
                      stringsAsFactors = F,header = T,row.names = 1)
#只要编码RNA
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] 
#取行名交集
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
#取fpkm中行名为comgene所含的
fpkm <- fpkm[comgene,]
#取Ginfo中行名为comgene所含的
Ginfo <- Ginfo[comgene,]
#新增Gene Symbol  
fpkm$Gene <- as.character(Ginfo$genename)
#去重复
fpkm <- fpkm[!duplicated(fpkm$Gene),]
#将行名变为Gene Symbol  
rownames(fpkm) <- fpkm$Gene   
#去除最后一列Gene
fpkm <- fpkm[,-ncol(fpkm)]   
#保存所以患者的fpkm文件
write.table(fpkm, file = "LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#保存癌症患者的fpkm文件
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)




####作图前准备####
#设置工作目录  
setwd("xena")
#加载包  
library(tidyverse)
#读取文件  
fpkm_01A <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,
                       check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,
                       check.names = F,stringsAsFactors = F,header = T)
#读取差异分析
load("LIHC_DEG.rda")
#进行差异分析(abs自己设置)
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
#提取目标基因
#法一（基因量较少）
gene <- c("CDC25C","CENPF","BIRC5","UBE2C","CDK1","BUB1","PLK1","EZH2",
          "PBK","KIF18A","MAD2L1","EPHA2","ETV4","UCHL1","EGF","SERPINE1",
          "HGF","SPTA1","MMP3")
#法二（基因量较多）
#读取想要基因集
gene60 = read.csv("gene60.csv",header = T)
#改行名并取交集
rownames(gene60) <- gene60$GENE
gene_60 <- intersect(rownames(gene60),rownames(fpkm_01A))
#得到想要基因
a <- fpkm_01A[gene,]
b <- fpkm_11A[gene,]
#t转换—行列转换——变成数据矩阵
a <- t(a)
b <- t(b)
class(a)
#变成数据框
a <- as.data.frame(a)
b <- as.data.frame(b)
#变成数据框——运用传导符  %>%  cltrl+shift+M 
#a <- a %>% t() %>% as.data.frame()
#b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "gene_19_01A.csv")
write.csv(b, file = "gene_19_11A.csv")
#Graphpad














####ROC####
#读取生存信息tsv文件
setwd("ROC")
library(tidyverse)
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
#gsub函数替换surv中sample的—换为.
surv$sample <- gsub("-",".",surv$sample)
#提取surv中sample作为行名
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#保存整理好的生存信息
write.table(surv, file = "survival.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,
                   check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#提取上次cox作图的10个基因
gene <- c("ICA1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()
#整合表达谱与生存信息
exp_sur <- cbind(exp10,surv)
write.table(exp_sur, file = "exp_sur.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#准备R包
library(ROCR)
library(rms)
#构建ROC预测模型
#提取exp_sur中ICA1基因和生存状态
ROC1 <- prediction(exp_sur$ICA1,exp_sur$OS) 
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR值
AUC <- performance(ROC1,"auc")   #计算曲线下面积(AUC)值
#看data AUC中 y vaule中值
AUC<- 0.5604839 #改 根据结果对AUC进行赋值 <- 0.6 0.7以上最好

#绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)   #绘制对角线
dev.off()




####timeROC####
#设置工作目录
setwd("timeROC")
#加载包
install.packages("timeROC")
install.packages("survival")
library(timeROC)
library(survival)
library(tidyverse)
#数据的整理与载入
exp_sur <- read.table("exp_sur.txt", header=T,sep="\t", 
                      check.names=F, row.names=1)
#将生存时间转为年  
exp_sur$OS.time <- exp_sur$OS.time/365
#提取其中肿瘤患者数据
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",]
#保存文件  
write.table(exp_sur_01A, file = "exp_sur_01A.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#构建ROC曲线函数
ROC3 <- timeROC(T=exp_sur_01A$OS.time,   #结局时间
                delta=exp_sur_01A$OS,   #结局指标
                marker=exp_sur_01A$ICA1,   #预测变量——基因ICA1
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
#查看模型变量信息  
ROC3   

#绘制ROC曲线
plot(ROC3,
     time=1, col="red")   #time是时间点，col是线条颜色
plot(ROC3,
     time=3, col="green", add=TRUE)   #add指是否添加在上一张图中
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   #添加标签信息

dev.off()




####TCGA差异分析热图与火山图####
#设置工作目录  
setwd("xena")
#加载包
library(tidyverse)
#读取fpkm全部基因组文件 
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,
                  check.names = F,stringsAsFactors = F,header = T)
#通过load函数读取rda文件  
load("LIHC_DEG.rda") 
#将其中p<0.05内容作为数据形式命名为deg  
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
#设置logFC截断值为1  
logFC_cutoff <- 1
#分组 > 1，上调；< 1，下调  
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
#新增change列且如果为K1附带上调，K2附带下调，若均不是为稳定
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
#展示change列有多少组内容
table(DEG$change)
#与特定基因取交集
#读取想要基因集
gene_anoikis = read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv",header = F)
#改行名并取交集
rownames(gene_anoikis) <- gene_anoikis$V1
gene_anoikis <- intersect(rownames(gene_anoikis),rownames(DEG))
#得到想要基因的cox
DEG <- DEG[gene_anoikis,]
###差异分析热图（适合基因量较少）
#读取包
library(pheatmap)
#提取deg中上调和下调部分 
cg = rownames(DEG)[DEG$change !="NOT"]
#将exp中上调和下调部分存储为exp_diff  
exp_diff <- exp[cg,]
#分组
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),
                  levels = c("N","T"))
group_list <- sort(group_list, decreasing = FALSE)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,#表达谱
         annotation_col=annotation_col,#分组量
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)
dev.off()
write.csv(exp_diff, file = "exp_diff.csv")


#带基因注释的火山图
#加载包
library(readxl)
library(ggrepel)
# 读取差异基因数据集
exprSet <- read.csv("LIHC_deseq2_all.csv")
# 设置阈值，整理数据
# 阈值不同，结果不同
cutoff_padj = 0.05
cutoff_logFC = 1
exprSet$Sig <- ifelse(exprSet$padj < cutoff_padj & 
                        abs(exprSet$log2FoldChange) >= cutoff_logFC, 
                      ifelse(exprSet$log2FoldChange > cutoff_logFC ,'Up','Down'),'no-DEGs')
#展示change列有多少组内容
table(exprSet$Sig)

ggplot(exprSet, aes(x = log2FoldChange, y = -log10(padj), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-16, 16)) + 
  # 辅助线
  geom_vline(xintercept=c(-cutoff_logFC,cutoff_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cutoff_padj),
             lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2FoldChange", y="-log10 (p.adjust)") +
  # 主题
  theme_bw() +
  # 标题
  ggtitle("p.adjust vs log2FoldChange") +
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  ) +  
  # 给点标上基因名
  geom_text_repel(
    # 可以设置跟上面不同的阈值，用数值替换即可
    data = subset(exprSet, exprSet$padj < cutoff_padj & abs(exprSet$log2FoldChange) >= cutoff_logFC),
    aes(label = Gene.Symbol), size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )



###差异基因火山图
library(ggpubr)
library(ggthemes)
#将其设为负对数（原因未知，数据对就行，猜测为p本身身为小数取负对数便于分析）
DEG$logP <- -log10(DEG$padj)
#绘图，x轴为logFC，y轴为logp  
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()
#增加基因上下调信息
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()
#添加分界线
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")#分界线更改需要跟随截断值
dev.off()
#添加基因标签信息
#新加一列label 
DEG$Label = ""
#对差异基因的p值进行从小到大的排序  
DEG <- DEG[order(DEG$padj), ]   
#将行名提取到基因中  
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes
#作图
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

dev.off()







####根据分组做生存分析####
###根据某一个基因的高低表达分组，并看生存有无差异  
#设置工作目录  
setwd("survival")
group <- read.csv("生存大于30分组结果.csv")
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
#gsub函数替换surv中sample的—换为.
surv$sample <- gsub("-",".",surv$sample)
#将sample换为行名
rownames(surv) <- surv$sample
#选取生存时间大于30天
surv1 = surv$OS.time >= 30;table(surv1)
surv2 = !(is.na(surv$OS.time)|is.na(surv$OS));table(surv2)
surv = surv[surv1&surv2,]
#输出为表格
write.csv(surv, file = "surv_筛选后大于30天数据.csv")
#将患者名变为行名
rownames(group) <- group$X
#提取group和surv中相同行名
com<- intersect(rownames(group),rownames(surv)) 
#筛选出group中行名与surv中相同的
group <- group[com,]
#筛选出surv中行名与group中相同的
surv <- surv[com,]
#判断顺序是否相同
identical(rownames(group),rownames(surv))
#surv中新增group
surv$group <- as.character(group$Group)
#将生存时间转换为月  
surv$OS.time <- surv$OS.time/30
#输出为表格
write.csv(surv, file = "surv_筛选后大于30天含分组数据.csv")
#判断数据形式  
class(surv$group)
#展示分组情况  
table(surv$group)
#加载包
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
#计算p值 <- 看右侧  
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
#拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
#绘制生存曲线
#方法1
#设置颜色，坐标
plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity(%)"
)
#添加标签
legend("topright",
       title = "Group",
       c("Low", "High"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))
#添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",
                                                               round(pValue, 3))))
text(25, 0.2, p.lab)
dev.off()


#方法2 <- 好看
#加载包
library(survminer)
#添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",
                                                               round(pValue, 3))))
#作图  
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "Set1", # 配色采用jco（黄蓝）
           legend.labs = c("C1", "C2"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = T, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
ydev.off()




####根据基因高低组做生存分析####
# 加载包
library(survival)
library(survminer)

# 可修改参数 
gene <- "KIF18A"              # 只需在此处修改基因名称
output_title <- "Survival Analysis"  # 自定义图表标题

# 数据读取与预处理
surv <- read.csv("surv大于30天的counts01A数据(带预后信息)调整后.csv")

# 转换为生存时间（月）并分组
surv$OS.time <- surv$OS.time / 30
surv$group <- ifelse(
  surv[[gene]] < median(surv[[gene]]),  # 动态引用基因列
  "Low", 
  "High"
)
surv$group <- factor(surv$group, levels = c("High", "Low"))

# 生存分析
fit <- survfit(
  Surv(OS.time, OS) ~ group, 
  data = surv
)

# 绘制生存曲线
ggsurvplot(
  fit,
  data = surv,
  title = output_title,        # 使用自定义标题
  pval = TRUE,                 
  pval.coord = c(30, 0.1),     # 调整P值显示位置
  conf.int = TRUE,             
  risk.table = TRUE,           
  palette = "Set1",             
  legend.labs = c("High", "Low"),
  legend.title = gene,         # 图例标题显示基因名称
  xlab = "Time (Months)",      
  ylab = "Survival Probability (%)",
  xlim = c(0, 120),            
  break.time.by = 20,          
  surv.median.line = "hv",     
  risk.table.y.text = FALSE,   
  ggtheme = theme_minimal()    # 更简洁的主题
)





####特定基因差异分析####
#设置工作目录  
setwd("ICA1_deg")
#加载包
library(DESeq2)
library(tidyverse)
#读取数据
#counts用来差异分析 <- 只能用于差异分析
counts_01A <- read.table("LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,
                         check.names = F,stringsAsFactors = F,header = T)
#fpkm用来分组
exp <- read.table("LIHC_fpkm_mRNA_01A.txt", sep = "\t",row.names = 1,
                  check.names = F,header = T)
#取两者交集  
com <- intersect(colnames(counts_01A),colnames(exp))
#提取，依据逻辑顺序使两组顺序相同  
exp <- exp[,com]
counts_01A <- counts_01A[,com]
#判断两个列名是否相同  
identical(colnames(counts_01A),colnames(exp))
#每次运行只改这个基因名  
gene <- "ICA1"
#提取exp中gene所对行，并进行数值取中位数  
med=median(as.numeric(exp[gene,]))
#根据表达量相较于中位数分组
conditions=data.frame(sample=colnames(exp),group=factor
                      (ifelse(exp[gene,]>med,"high","low"),
                        levels = c("low","high"))) %>% 
  column_to_rownames("sample")
#差异分析 <- 时间较久
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
#确定比较两组
resultsNames(dds)
res <- results(dds)
#保存  
save(res,file="res_deseq2_ICA1.Rda")
#加载
load("res_deseq2_ICA1.Rda")
#挑选差异基因  
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
#输出为文本
write.table(res_deseq2,"LIHC_res_deseq2_ICA1_abs1.txt",sep = "\t",row.names = T,
            col.names = NA,quote = F)
#输出为表格
write.csv(res_deseq2, file = "LIHC_res_deseq2_ICA1_abs1.csv")




####特定基因差异分析GO####
#设置工作目录  
setwd("ICA1_fuji")
#加载包
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
#读取差异分析结果
load("res_deseq2_ICA1.rda")
#自选差异  
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.6, padj < 0.05)
#将行名转换成列，并命名为基因名
DEG <- DEG %>% rownames_to_column("Gene")
#将基因名称转换为数字形式
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
#将其加入DEG最后一行
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,#对应基因数字
                OrgDb = org.Hs.eg.db, 
                ont = "all",#三组分：分子功能，细胞组分，生物过程
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
#读取go分析结果  
ego_res <- ego@result
#保存go分析结果  
save(ego,ego_res,file = "GO_ICA1_DEG_ABS_0.6.Rdata")
#保存go分析整体结果
write.csv(ego_res,file="egO_ICA1_res_0.6.csv")
#可视化
#柱状图
barplot(ego, showCategory = 20,color = "pvalue")
#气泡图
dotplot(ego, showCategory = 20)
#分类展示 <- 三组分：分子功能，细胞组分，生物过程
#柱状图 
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
#气泡图  
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')




####特定基因差异分析KEGG####
#设置工作目录  
setwd("ICA1_KEGG")
#读取包
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
#读取差异分析结果 
load("res_deseq2_ICA1.Rda")
#自选差异比例  
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.6, padj < 0.05)
#将行名变成列给其命名为GENE 
DEG <- DEG %>% rownames_to_column("Gene")
#根据org.Hs.eg.db包将GENE转换为ENTREZID形式
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
#将其导入到deg中  
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#KEGG分析
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
#提取kk结果  
kk_res <- kk@result
#保存  
save(kk,kk_res,file = "KEGG_ICA1_DEG_ABS_0.6.Rdata")
#保存go分析整体结果
write.csv(kk_res,file="KEGG_ICA1_res_0.6.csv")
#加载
load("KEGG_ICA1_DEG_ABS_0.6.Rdata")
#柱状图
barplot(kk, showCategory = 20,color = "pvalue")
#气泡图
dotplot(kk, showCategory = 20)

dev.off()




####特定基因差异分析GSEA####
#内容同kegg
#设置工作目录  
setwd("GSEA_ICA1")
#读取包
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
#读取差异分析结果 
load("res_deseq2_ICA1.Rda")
#自选差异比例  
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.6, padj < 0.05)
#将行名变成列给其命名为GENE 
DEG <- DEG %>% rownames_to_column("Gene")


































#根据org.Hs.eg.db包将GENE转换为ENTREZID形式
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
#将其导入到deg中  
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
#设置文件夹
msigdb_GMTs <- "msigdb_v7.0_GMTs"
#设置所用文件   #c2.all.v7.0.entrez.gmt 或 c5.all.v7.0.entrez.gmt
msigdb <- "c5.all.v7.0.entrez.gmt"   
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
#无脑run  
geneList = DEG[,3]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
#降序排列  
geneList = sort(geneList, decreasing = TRUE)
#设置随机数种子  
set.seed(1)
#GSEA分析  004
KEGG<-GSEA(geneList,TERM2GENE = kegmt)
#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
#保存文件  
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
write.csv(KEGG_result_df, file = "GSEA_MSigDb_C5_result.csv")
save(KEGG,KEGG_result_df,file = "GSEA_deg_ICA1.rda")#以ICA1为差异作图
#单个图绘制
library(enrichplot)
gseaplot2(KEGG,1,color="red")
gseaplot2(KEGG,3,color="red",pvalue_table = T)

#汇总结果
gseaplot2(KEGG, geneSetID = c(1), subplots = 1:3)
gseaplot2(KEGG, geneSetID = c(1,3), subplots = 1:3)
gseaplot2(KEGG, geneSetID = 1:3, subplots = 1)
gseaplot2(KEGG, geneSetID = 10:26, subplots = 1:3)

dev.off()




##### ConsensusClusterPlus聚类 #####
#设置工作目录  
setwd("2_cluster")
#加载包
library(tidyverse)
library(ConsensusClusterPlus)
surv <- read.csv("surv_筛选后大于30天数据.csv")
#读取fpkm01A数据  
exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,
                  check.names = F,stringsAsFactors = F,header = T)
#将患者设为行名
rownames(surv) <- surv$X
#进行t转换
surv <- t(surv)
#取表达谱和生存的患者相同交集
com<- intersect(colnames(exp),colnames(surv)) 
#要表达谱中共有的患者
exp <- exp[,com]
#要生存中共有的患者
surv <- surv[,com]
#判断是否相同
identical(colnames(exp),colnames(surv))
#输出文件
write.csv(exp,file = "surv大于30天的fpkm01A数据.csv")
write.csv(surv,file = "surv大于30天的有表达的生存数据(源于fpkm).csv")

d=as.matrix(exp)
#挑选所选基因  
gene0 = read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv",header=F) 
gene=as.matrix(gene0)
d <- d[gene,]
write.csv(d,file = "ConsensusClusterPlus聚类原始数据 surv大于30天.csv")
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:41],] #2是基因个数，和所选基因相关 
d = sweep(d,1, apply(d,1,median,na.rm=T))
title=("JULEI") ##文件夹输出图片的位置
set.seed(1) #我发现设不设置种子都一样
results = ConsensusClusterPlus(d,maxK=9,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance=
                                 "pearson",seed=1,plot="pdf")
#验证可以不跑  
results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]
#画另一组图片
icl = calcICL(results,title=title,plot="pdf")
#其中数字2由consensus中Delta area中突变值分 <- 一般分为两组
group<-results[[2]][["consensusClass"]]
group<-as.data.frame(group)
group$group <- factor(group$group,levels=c(1,2))
#保存group文件
save(group,file = "group_41.Rda")
exp_gene <- exp[gene,]
#绘制ConsensusClusterPlus后的热图
library(pheatmap)
group <- group %>% 
  rownames_to_column("sample")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(exp_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% 
  column_to_rownames("sample") %>% t(.)
pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,
         fontsize_col=3)
pheatmap(c,annotation = annotation,
         annotation_colors = list(group = c("1" ="#01468b","2"= "#ee0000")),
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,cluster_row = F,
         fontsize_col=3)
dev.off()




####挑选共有的counts数据进行聚类间差异分析####
surv <- read.csv("surv_筛选后大于30天数据.csv")
#读取fpkm01A数据  
exp <- read.table("LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,
                  check.names = F,stringsAsFactors = F,header = T)
#将患者设为行名
rownames(surv) <- surv$X
#进行t转换
surv <- t(surv)
#取表达谱和生存的患者相同交集
com<- intersect(colnames(exp),colnames(surv)) 
#要表达谱中共有的患者
exp <- exp[,com]
#要生存中共有的患者
surv <- surv[,com]
#判断是否相同
identical(colnames(exp),colnames(surv))
#输出文件
write.csv(exp,file = "surv大于30天的counts01A数据.csv")
write.csv(surv,file = "surv大于30天的有表达的生存数据(源于counts).csv")
#读取预后筛选后的表达谱数据
exp <- read.csv("surv大于30天的counts01A数据.csv")
#将exp患基因作为行名
rownames(exp) <- exp$X
#去除第一行
exp <- exp[,-1]
#读取分组数据
group <- read.csv("生存大于30分组结果.csv")
#将group患者名作为行名
rownames(group) <- group$X
#进行t转换
group <- t(group)
identical(colnames(exp),colnames(group))
#取两者交集  
com <- intersect(colnames(exp),colnames(group))
#展示01A和11A有多少  
table(substr(com,14,16))
#提取交集部分
exp <- exp[,com]
group <- group[,com]
identical(colnames(exp),colnames(group))
#保存待整理数据
write.csv(exp,"surv大于30天的counts01A数据_待合并.csv")
write.csv(group,"surv大于30天的group数据_待合并.csv")



####计算患者免疫评分（免疫浸润情况）与肿瘤纯度#####
#加载包
library(utils) #这个包应该不用下载，自带的
#下载estimate包
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
expr <- read.table("surv大于30天的fpkm01A数据.txt",sep = "\t",row.names = 1,check.names = F,
                   stringsAsFactors = F,header = T)
#计算免疫评分
filterCommonGenes(input.f = "surv大于30天的fpkm01A数据.txt",   #输入文件名
                  output.f = "surv大于30天的fpkm01A数据.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("surv大于30天的fpkm01A数据.gct",   #刚才的输出文件名
              "surv大于30天的fpkm01A数据_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台
#输出每个样品的打分

#对免疫评分进行整理  
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))
rownames(result) <- colnames(expr)
#保存文件  
write.table(result, file = "surv大于30天的fpkm01A数据_estimate_score.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F) # 保存并覆盖得分
#读取免疫评分结果
result <- read.table("surv大于30天的fpkm01A数据_estimate_score.txt",sep = "\t",
                     row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#读取患者分组信息文件
patient_groups <- read.csv("生存大于30分组结果.csv", header = TRUE)

# 读取患者分组信息文件
TCGA_TME <- result
#将分组的行名换为患者名字
rownames(patient_groups) <- patient_groups$X
# 数据预处理
#判断患者分组与患者的ciber是否匹配对应
identical(rownames(patient_groups),rownames(TCGA_TME))
#提取患者的分组信息加入ciber结果
TCGA_TME$group <- patient_groups$Group
#将行名新增为样本列
TCGA_TME$sample <- row.names(TCGA_TME)
TME_TME_new <- melt(TCGA_TME)
colnames(TME_TME_new) <- c("Group", "Sample", "variable", "value")
# 设置主题
mytheme <- 
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5),
    axis.title = element_text(size = 12, color = "black"), 
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 12)
  )
# 出图
box_TME_new <- ggplot(TME_TME_new, aes(x = Group, y = value)) + 
  labs(y = "ESTIMATE", x = NULL, title = NULL) +  
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, size = 0.4,
               outlier.alpha = 0.5, outlier.shape = 16, outlier.size = 2) + 
  theme_bw() + mytheme + 
  scale_fill_manual(values = c("#EF6262","#1D5B79")) + 
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ variable, scales = "free", ncol = 4) + 
  stat_compare_means(aes(group = Group),
                     label = "p.format",
                     method = "wilcox.test",
                     size = 3.5,
                     hide.ns = TRUE)
#出图
box_TME_new


####cibersort分析22种免疫细胞####
#设置工作目录
setwd("cibersort")   
#安装包  install.packages('e1071')
#安装包  install.packages('parallel')
#install.packages("BiocManager")
BiocManager::install("preprocessCore")
# 读取包
library(e1071)
library(parallel)
library(preprocessCore)
library(ggplot2)
# 脚本勿动！！！！！！！！！！！！！！！！！！！！
source("CIBERSORT.R")
# 加载免疫细胞注释
sig_matrix <- "LM22.txt"
# 肿瘤患者表达谱 <- 通过表达谱中免疫细胞特殊表达值算比例
mixture_file <- 'surv大于30天的fpkm01A数据.txt'
# 运算比例 <- 时间较久
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm = 100, QN = TRUE)
# 保存中间文件
save(res_cibersort, file = "res_cibersort_surv大于30天的fpkm01A.Rdata")
# 后三列无用保存前22列
res_cibersort <- res_cibersort[, 1:22]
# 去除丰度全为0的细胞
ciber.res <- res_cibersort[, colSums(res_cibersort) > 0]
# 保存文件
write.csv(ciber.res, "res_cibersort_surv大于30天的fpkm01A 丰度大于0.csv")
#读取文件
ciber.res <- read.csv("res_cibersort_surv大于30天的fpkm01A 丰度大于0.csv", row.names = 1 , check.names = FALSE)
#读取患者分组信息文件
patient_groups <- read.csv("生存大于30分组结果.csv", header = TRUE)

####总体样本直方图####
# 加载所需的R包
library(tidyverse)  # 用于数据处理和可视化
library(RColorBrewer)  # 用于生成颜色调色板
library(dplyr)  # 用于数据处理
# 创建颜色调色板
mypalette <- colorRampPalette(brewer.pal(8, "Set1"))
#将ciber.res转换为数据框，并添加"Sample"列作为行名
dat <- ciber.res %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
#将数据从宽格式转换为长格式，将"Cell_type"列作为键，"Proportion"列作为值
  gather(key = Cell_type, value = Proportion, -Sample)
#根据"Cell_type"列的中位数对"Cell_type"进行排序
a <- dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)
#将"Cell_type"列转换为因子，并按照排序后的顺序设置水平
dat$Cell_type <- factor(dat$Cell_type, levels = a)
#绘制箱线图
ggplot(dat, aes(Cell_type, Proportion, fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21, color = "black",) + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))






####基因表达与ciber的相关性####
#设置工作目录
setwd("cor")
#安装包install.packages("corrplot")
#加载包
library(corrplot)
library(tidyverse)
#读取表达谱数据
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,
                   check.names = F,stringsAsFactors = F,header = T)
#提取上次cox所得十个基因
gene <- c("ICA1")
#取表达谱中cox所得十个基因数据  
exp <- expr[gene,]
#将其t转换，并变为数据框形式
exp <- exp %>% t() %>% as.data.frame()
# 读取cibersort文件
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,
                    check.names = F,stringsAsFactors = F,header = T)
#消除无用列 
ciber <- ciber[,1:22]
#判断两者行名是否相同
identical(rownames(ciber),rownames(exp))
#判断两者是否为数值形式
class(exp$ICA1)
class(ciber$`B cells naive`)
#计算相关性
cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
#将行名变为基因名  
rownames(cor)<-colnames(exp)
#计算p值
cor_res <- cor.mtest(cor,conf.level = 0.95)#置信区间
corrplot(cor,
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.6,#图中数字大小
         tl.cex = 0.7,#图例上字母大小
         cl.align = "l")
dev.off()



####4大免疫细胞####
#加载所需包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
#读取数据
ciber.res <- read.csv("res_cibersort_surv大于30天的fpkm01A 丰度大于0.csv", row.names = 1, check.names = FALSE)
TCGA_TME_four <- as.data.frame(ciber.res[, 1:20])
immCell_four_type <- read.table("Cibersort_four_types.txt", header = TRUE, row.names = NULL, sep = "\t")
colnames(TCGA_TME_four) <- immCell_four_type$Immune.cells
# 读取患者分组信息文件
patient_groups <- read.csv("生存大于30分组结果.csv", header = TRUE)
#将分组的行名换为患者名字
rownames(patient_groups) <- patient_groups$X
# 数据预处理
#判断患者分组与患者的ciber是否匹配对应
identical(rownames(patient_groups),rownames(TCGA_TME_four))
#提取患者的分组信息加入ciber结果
TCGA_TME_four$group <- patient_groups$Group
#将行名新增为样本列
TCGA_TME_four$sample <- row.names(TCGA_TME_four)
TME_four_new <- melt(TCGA_TME_four)
colnames(TME_four_new) <- c("Group", "Sample", "Immune.cells", "Composition")
TCGA_TME_four_new2 <- left_join(TME_four_new, immCell_four_type, by = "Immune.cells") %>% 
  group_by(Sample, Group, Types) %>%
  summarize(Sum = sum(Composition))
# 设置主题
mytheme <- 
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5),
    axis.title = element_text(size = 12, color = "black"), 
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 12)
  )
# 出图
box_four_immtypes <- ggplot(TCGA_TME_four_new2, aes(x = Group, y = Sum)) + 
  labs(y = "Cell composition", x = NULL, title = NULL) +  
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, size = 0.4,
               outlier.alpha = 0.5, outlier.shape = 16, outlier.size = 2) + 
  theme_bw() + mytheme + 
  scale_fill_manual(values = c("#EF6262","#1D5B79")) + 
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ Types, scales = "free", ncol = 4) + 
  stat_compare_means(aes(group = Group),
                     label = "p.format",
                     method = "wilcox.test",
                     size = 3.5,
                     hide.ns = TRUE)
#出图
box_four_immtypes




####cluster分组CIBERSORT画箱式图 <- 免疫细胞比例差异####
#加载包  
library(tidyverse)
#取22种免疫细胞
ciber.res <- read.csv("res_cibersort_surv大于30天的fpkm01A 丰度大于0.csv", row.names = 1 , check.names = FALSE)
#提取分组信息
group <- read.csv("生存大于30分组结果.csv")
rownames(group) <- group$X
#判断其行名顺序是否一致  
identical(rownames(ciber.res),rownames(group))
class(group$Group)
#把group中group放入a中  
ciber.res$group <- group$Group
#将a的行名变为一列  
ciber.res <- ciber.res %>% rownames_to_column("sample")
#保存处理后数据
write.csv(ciber.res,"cluster分组CIBERSORT画箱式图原始数据(可用于热图和箱式图绘制).csv")
library(ggsci)
library(tidyr)
library(ggpubr)
#绘图无脑跑
p <- gather(ciber.res,key=CIBERSORT,value = Proportion,-c(group,sample))
# 自定义分组颜色
my_colors <- c("C1" = "#EF6262", "C2" = "#1D5B79")

# 在ggboxplot()函数中使用scale_fill_manual()来设置自定义颜色
ggboxplot(p, x = "CIBERSORT", y = "Proportion",
          fill = "group",outlier.size = 0.8 ,outlier.alpha = 0.4,outlier.alpha = 0.8,outlier.shape = 16) +
  scale_fill_manual(values = my_colors) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns"))) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
          )






####ssGSEA <- 免疫细胞表达量差异####
#设置工作目录  
setwd("ssGSEA_AY")
#读取包  
BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)
###无脑run
#准备细胞marker
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
#将cellMarker文件列名的第2个修改为celltype  
colnames(cellMarker)[2] <- "celltype"
#将cellMarker文件以celltype为分组拆分成list数据格式
type <- split(cellMarker,cellMarker$celltype)
#处理data.tables列表通常比使用group by参数按组对单个data.table进行操作要慢得多
#将list中每个celltype中的基因进行合并 
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
#保存中间文件
save(cellMarker,file = "cellMarker_ssGSEA.Rdata")

#表达量矩阵的准备
#行是基因，列是样本
#读取表达文件 <- 更换样本记得改名字  
expr <- data.table::fread("surv大于30天的fpkm01A数据.txt",data.table = F)   
rownames(expr) <- expr[,1]   #将第一列作为行名
expr <- expr[,-1]   #去除第一列
expr <- as.matrix(expr)   #将expr转换为矩阵格式

#使用ssGSEA量化免疫浸润
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
#作为数据形式，并进行t转换  
a <- gsva_data %>% t() %>% as.data.frame()
#读取分组信息  
load("group_AY.Rda")
#判断行名顺序  
identical(rownames(a),rownames(group))
#提取group分组信息并作为a中group一列  
a$group <- group$Group
a <- a %>% rownames_to_column("sample")
#保存  
write.table(a,"surv大于30天的fpkm01A的ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#读取包  
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))
#作图  
ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()




####COX回归分析####
#设置工作目录
setwd("cox")
#安装包install.packages("survival")
#安装包install.packages("forestplot")
#加载包  
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
#gsub函数替换surv中sample的—换为.
surv$sample <- gsub("-",".",surv$sample)
#将sample换为行名
rownames(surv) <- surv$sample
#删掉无用列
surv <- surv[,-1]
surv <- surv[,-2]
#选取生存时间大于30天
surv1 = surv$OS.time >= 30;table(surv1)
surv2 = !(is.na(surv$OS.time)|is.na(surv$OS));table(surv2)
surv = surv[surv1&surv2,]
#读取转录组表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,
                   check.names = F,stringsAsFactors = F,header = T)
#取两者交集  
comgene <- intersect(colnames(expr),rownames(surv))
#展示01A和11A有多少  
table(substr(comgene,14,16))
#提取交集部分
expr <- expr[,comgene]
surv <- surv[comgene,]
#表达数据整理完毕
#读取tcga差异分析结果
load("LIHC_DEG.rda")
#进行差异分析(abs自己设置)  
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
#整合（行名为患者编号，列名为基因差异）
#提取expr中行名为res_deseq2中包含的内容并进行横竖转换，保存为数据形式
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
#cbind按列将两个数据框合为一个 <- col
#rbind按行将两个数据框合为一个 <- row
surv.expr <- cbind(surv,deg_expr)
#保存生存和表达患者数据
write.csv(surv.expr, file = "surv30_expr1_Padj005.csv")
#Cox分析
#hr大于1：基因表达越高，患者存活时间越低，危险因素 <- 暴露因素为阳性事件发生的促进因素
#hr小于1：基因表达越高，患者存活时间越长，保护因素 <- 暴露因素为阳性事件发生的阻碍因素
#hr等于1：基因表达水平与患者生存事件无相关性       <- 暴露因素与阳性事件发生无关
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- 
    rbind.data.frame(Coxoutput,
                     data.frame(gene = g,
                                HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                lower = as.numeric(coxSummary$conf.int[,3][1]),
                                upper = as.numeric(coxSummary$conf.int[,4][1]),
                                stringsAsFactors = F),
                     stringsAsFactors = F)
}

#保存全部基因cox结果
write.table(Coxoutput, file = "cox results all.txt",sep = "\t",row.names = F,
            col.names = T,
            quote = F)
#筛选top基因
setwd("cox")
#读取全部基因coxoutput
Coxoutput <- read.table("cox results all.txt",sep = "\t",row.names = 1,
                        check.names = F,
                        stringsAsFactors = F,header = T)
#读取想要基因集
gene = read.csv("GENE126新（2023.3.16）.csv",header = T)
#改行名并取交集
rownames(gene) <- gene$GENE
coxgene <- intersect(rownames(gene),rownames(Coxoutput))
#得到想要基因的cox
Coxoutput <- Coxoutput[coxgene,]
#输出目的基因cox回归结果
write.table(Coxoutput, file = "cox results_LIHC abs1_Anoikis0.4 all.txt",sep = "\t",
            row.names = T,
            col.names = T,
            quote = T)
#读取全部目的基因cox回归结果
Coxoutput <- read.table("cox results_LIHC abs1_Anoikis0.4 all.txt",sep = "\t",
                        row.names = 1,
                        check.names = F,
                        stringsAsFactors = F,header = T)
#p值截取点0.05
pcutoff <- 0.05
HR <- 1.1
#取出p值小于阈值的基因
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] 
#取出HR小于阈值的基因
topgene <- topgene[which(topgene$HR > HR),] 
#将行名变为其中一列
topgene$Gene <- rownames(topgene)
#保存最终cox基因结果
write.csv(topgene, file = "cox results_LIHC abs1_Anoikis0.4 p0.05.csv")

##################################对cox results_LIHC abs1_Anoikis0.4 p0.05文件进行HR指数的筛选
#加载包
library(ggplot2)
library(dplyr)
# 读取数据
cox_results <- read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后.csv")
# 提取相关信息
hazard_ratio <- cox_results$HR
lower_ci <- cox_results$lower
upper_ci <- cox_results$upper
p_value <- cox_results$pvalue
gene <- cox_results$Gene
# 创建数据框
df <- data.frame(gene = gene, hazard_ratio = as.numeric(hazard_ratio), 
                 lower_ci = as.numeric(lower_ci), upper_ci = as.numeric(upper_ci), p_value = p_value)
# 定义函数根据p值返回对应的星号
get_pvalue_label <- function(p_value) {
  ifelse(p_value < 0.001, "***", 
         ifelse(p_value >= 0.001 & p_value < 0.01, "**", 
                                        ifelse(p_value >= 0.01 & p_value < 0.05, "*", "")))
}
# 添加星号标签列
df <- df %>% mutate(pvalue_label = get_pvalue_label(p_value))
# 绘制森林图
forest_plot <- ggplot(df, aes(x = hazard_ratio, y = gene)) +
  geom_point(aes(color = hazard_ratio), shape = 15, size = 4) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Hazard Ratio", y = "Gene", title = "Forest Plot") +
  scale_color_gradient(low = "blue", high = "red", name = "HR") +
  geom_text(aes(label = pvalue_label), vjust = -0.2, hjust = -0.2, size = 5, color = "black") +
  theme_minimal()
# 显示森林图
print(forest_plot)




####LASSO COX 回归模型的建立####
#构建内部训练集和内部验证集
#读取病人表达谱和生存信息数据
  data = read.table("ConsensusClusterPlus聚类原始数据 surv大于30天 含生存状态 用于预后模型构建.txt",header = T)
#对data进行t转换，要求行名为患者，列名为基因
  data <- t(data)
#保存转置后的数据
  write.csv(data,"surv大于30天的counts01A数据(带预后信息)调整后.csv")
#加载包
  library(survival)
  library(survminer)
  library(glmnet)
  library(survivalROC)
  library(caret)
  library(coin)
#读取调整后病人表达谱和生存信息数据
  data = read.csv("surv大于30天的counts01A数据(带预后信息)调整后.csv",header = T)
#在data前加一排样本编号，便于切割数据
  sample <- 1:nrow(data)
#合并样本和样本编号
  data <- cbind(sample,data)
#设计随机数种子
  set.seed(287)
#7:3比例分割数据
  ind <- createDataPartition(data$sample,p = 0.7,list = F)
#以随机数对data数据分为训练集和验证集
  train <- data[ind,]
  test <- data[-ind,]
#保存训练集和验证集
  write.csv(train,file = "train_TCGA_seed287.csv")
  write.csv(test,file = "test_TCGA_seed287.csv")




#进行内部训练集训练
#读取训练集和内部验证集
  train <- read.csv("train_TCGA_seed287.csv")
  test <- read.csv("test_TCGA_seed287.csv")
#以数据形式将训练集数据读取
  train <- as.data.frame(train)
  test <- as.data.frame(test)
#基因表达从第几列开始就开始的列数，这里是从第6列开始（不同读取方式可能不同要注意）
  x<-as.matrix(train[,6:dim(train)[2]])
#这里注意生存的列名，如果研究的RFS，根据列名替换掉
  y<-data.matrix(survival::Surv(train$OS.time,train$OS))
#进行lasso回归
  cv.fit<-cv.glmnet(x,y,type.measure = "deviance",family = "cox")
#绘制LASSO交叉验证曲线
  plot(cv.fit)
#进行LASSO权重分析
  cv.fit1<-glmnet(x,y,type.measure = "deviance",family = "cox")
#进行LASSO权重分析曲线
  plot(cv.fit1,label = TRUE)
#计算每个样本的风险评分
  risk.score <- predict(cv.fit$glmnet.fit,
                        newx = x,
                        s=cv.fit$lambda.min, 
                        type="link")
  train$risk.score<-risk.score
#输出含风险因子且带预后的基因表达谱
  write.csv(train,file = "含风险因子评分且带预后的训练集表达谱.csv")
#输出每个基因的风险评分
  lasso.coef <- predict(cv.fit, s = cv.fit$lambda.min, type = "coefficients")
#展示每个基因的风险评分
  lasso.coef
#将lasso获得的每个基因的风险评分进行保存
  lasso <- as.data.frame(as.matrix(lasso.coef))
  write.csv(lasso,file = "lasso基因权重.csv")




#ROC验证模型
#1年生存率,也可改为1年生存率 365*1
  ROC1<- survivalROC(Stime=train$OS.time, 
                    status=train$OS, 
                    marker = train$risk.score, 
                    predict.time =365*1, #1年生存率,也可改为1年生存率 365*1
                    method="KM")
#3年生存率,也可改为3年生存率 365*3
  ROC3<- survivalROC(Stime=train$OS.time, 
                     status=train$OS, 
                     marker = train$risk.score, 
                     predict.time =365*3, #3年生存率,也可改为3年生存率 365*3
                   method="KM")
#5年生存率,也可改为5年生存率 365*5
  ROC5<- survivalROC(Stime=train$OS.time, 
                     status=train$OS, 
                     marker = train$risk.score, 
                     predict.time =365*5, #5年生存率,也可改为5年生存率 365*5
                     method="KM")
#查看1，3，5年对应的AUC值
  ROC1$AUC;ROC3$AUC;ROC5$AUC
#绘制1年、3年和5年的ROC曲线
  plot(ROC1$FP, ROC1$TP, xlim = c(0, 1), ylim = c(0, 1), col = "#377eb8", type = "l",
       xlab = "False Positive Rate", ylab = "True Positive Rate", 
       main = "Time-dependent ROC Curve (LIHC Train)", lwd = 2, 
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)
#添加3年和5年的ROC曲线
  lines(x = ROC3$FP, y = ROC3$TP, lwd = 2, type = "l", col = "#00ba38", lty = 2)
  lines(x = ROC5$FP, y = ROC5$TP, lwd = 2, type = "l", col = "#fb8d3d", lty = 3)
#添加图例
  legend("bottomright", bty = "n",
         col = c("#377eb8", "#00ba38", "#fb8d3d"),
         legend = c(paste("AUC (1 year) =", round(ROC1$AUC, 2)), 
                    paste("AUC (3 year) =", round(ROC3$AUC, 2)), 
                    paste("AUC (5 year) =", round(ROC5$AUC, 2))),
         lty = c(1, 2, 3), lwd = 2, cex = 1, text.col = "black")
#添加对角线
  abline(0, 1, col = "gray", lty = 2)

  
#将风险因子排序并分高低风险组
  train$risk.score <- ifelse(train$risk.score>median(train$risk.score),"High Risk Score","Low Risk Score")
#将生存时间转换为月  
  train$OS.time <- train$OS.time/30
  fit <- survfit(Surv(OS.time,OS) ~ train$risk.score, data = train)
#作图
  ggsurvplot(fit,
             data = train,
             pval = T,
             conf.int = F, # 显示置信区间
             risk.table = TRUE, # 显示风险表
             risk.table.col = "strata",
             palette = "Set1", # 配色采用jco（黄蓝）
             legend.labs = c("High Risk Score","Low Risk Score"), # 图例
             size = 1,
             xlim = c(0,120), # x轴长度，一般为0-10年
             break.time.by = 20, # x轴步长为20个月
             legend.title = "Riskscore",
             surv.median.line = "hv", # 是否添加垂直和水平的中位生存
             ylab = "Survival probability (%)", # 修改y轴标签
             xlab = "Time (Months)", # 修改x轴标签
             ncensor.plot = T, # 显示删失图块
             ncensor.plot.height = 0.25,
             risk.table.y.text = FALSE)

#进行风险因子关联图绘制
  library(ggrisk)
  library(survival)
  library(survminer)
  library(rms)
#提取用于风险测试列
  train <- train[,c("OS.time","OS","KIF18A","SPP1","PLK1","SLC2A2","EZH2","MMP3")]
#进行风险关联图绘制
  fit <- cph(Surv(OS.time,OS) ~ KIF18A + SPP1 + PLK1 + SLC2A2 + EZH2 + MMP3, train)
#绘图
  ggrisk(fit, cutoff.value = "median",
         cutoff.show = FALSE,#不显示cutoff值
         color.A=c(low='#1D5B79',high='#EF6262'),#A图中点的颜色
         color.B=c(code.0='#1D5B79',code.1='#EF6262'), #B图中点的颜色
         color.C = c(low = "#1D5B79", median = "#ffffff", high = "#EF6262"), #C图中表达量的颜色
         )
  
#进行单因素cox回归
  cox <- coxph(Surv(OS.time,OS)~risk.score,train)
  summary(cox)
#记录结果
#                            exp(coef) exp(-coef) lower .95 upper .95
#risk.scoreLow Risk Score    0.3867      2.586    0.2442    0.6125

#Score (logrank) test = 17.62  on 1 df,   p=0.00003

#低风险组死亡风险是高风险组的0.3867倍，P为0.00003




#进行内部验证集验证
#加载lasso获得风险评分公式
  test$risk.score <- 0.173030087*test$KIF18A+
                     0.063695246*test$SPP1+
                     0.014878703*test$PLK1-
                     0.037960223*test$SLC2A2+
                     0.001804856*test$EZH2+
                     0.121272810*test$MMP3
#输出含风险因子且带预后的基因表达谱
  write.csv(test,file = "含风险因子评分且带预后的验证集表达谱.csv")
####仙桃学术绘风险因子关联图
#1年生存率,也可改为3年生存率 365*3
  ROC1<- survivalROC(Stime=test$OS.time, 
                    status=test$OS, 
                    marker = test$risk.score, 
                    predict.time =365*1, #1年生存率,也可改为1年生存率 365*1
                    method="KM")
#3年生存率,也可改为3年生存率 365*3
  ROC3<- survivalROC(Stime=test$OS.time, 
                     status=test$OS, 
                     marker = test$risk.score, 
                     predict.time =365*3, #3年生存率,也可改为3年生存率 365*3
                     method="KM")
#5年生存率,也可改为5年生存率 365*5
  ROC5<- survivalROC(Stime=test$OS.time, 
                     status=test$OS, 
                     marker = test$risk.score, 
                     predict.time =365*5, #5年生存率,也可改为5年生存率 365*5
                     method="KM")
#查看1，3，5年对应的AUC值
  ROC1$AUC;ROC3$AUC;ROC5$AUC
#绘制1年、3年和5年的ROC曲线
  plot(ROC1$FP, ROC1$TP, xlim = c(0, 1), ylim = c(0, 1), col = "#377eb8", type = "l",
       xlab = "False Positive Rate", ylab = "True Positive Rate", 
       main = "Time-dependent ROC Curve (LIHC Test)", lwd = 2, 
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)
#添加3年和5年的ROC曲线
  lines(x = ROC3$FP, y = ROC3$TP, lwd = 2, type = "l", col = "#00ba38", lty = 2)
  lines(x = ROC5$FP, y = ROC5$TP, lwd = 2, type = "l", col = "#fb8d3d", lty = 3)
#添加图例
  legend("bottomright", bty = "n",
         col = c("#377eb8", "#00ba38", "#fb8d3d"),
         legend = c(paste("AUC (1 year) =", round(ROC1$AUC, 2)), 
                    paste("AUC (3 year) =", round(ROC3$AUC, 2)), 
                    paste("AUC (5 year) =", round(ROC5$AUC, 2))),
         lty = c(1, 2, 3), lwd = 2, cex = 1, text.col = "black")
#添加对角线
  abline(0, 1, col = "gray", lty = 2)


#根据风险得分将验证集分高低组
  test$risk.score <- ifelse(test$risk.score>median(test$risk.score),"High Risk Score","Low Risk Score")
#将生存时间转换为月  
  test$OS.time <- test$OS.time/30
  fit <- survfit(Surv(OS.time,OS) ~ test$risk.score, data = test)
#作图
  ggsurvplot(fit,
             data = test,
             pval = T,
             conf.int = F, # 显示置信区间
             risk.table = TRUE, # 显示风险表
             risk.table.col = "strata",
             palette = "Set1", # 配色采用jco（黄蓝）
             legend.labs = c("High Risk Score","Low Risk Score"), # 图例
             size = 1,
             xlim = c(0,120), # x轴长度，一般为0-10年
             break.time.by = 20, # x轴步长为20个月
             legend.title = "Riskscore",
             surv.median.line = "hv", # 是否添加垂直和水平的中位生存
             ylab = "Survival probability (%)", # 修改y轴标签
             xlab = "Time (Months)", # 修改x轴标签
             ncensor.plot = T, # 显示删失图块
             ncensor.plot.height = 0.25,
             risk.table.y.text = FALSE)

#进行cox回归
  cox <- coxph(Surv(OS.time,OS)~risk.score,test)
  summary(cox)

#                          exp(coef) exp(-coef) lower .95 upper .95
#risk.scoreLow Risk Score    0.4578      2.184    0.2356    0.8898
#Score (logrank) test = 5.57  on 1 df,   p=0.02




#对整体数据进行分析
#加入风险评分
  data$risk.score <- 0.173030087*data$KIF18A+
                     0.063695246*data$SPP1+
                     0.014878703*data$PLK1-
                     0.037960223*data$SLC2A2+
                     0.001804856*data$EZH2+
                     0.121272810*data$MMP3
#输出含风险因子且带预后的基因表达谱
  write.csv(data,file = "含风险因子评分且带预后的TCGA所有样本表达谱.csv")
  data <- read.csv("含风险因子评分且带预后的TCGA所有样本表达谱.csv")
#1年生存率,也可改为3年生存率 365*3
  ROC1<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$risk.score, 
                     predict.time =365*1, #1年生存率,也可改为1年生存率 365*1
                     method="KM")
#3年生存率,也可改为3年生存率 365*3
  ROC3<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$risk.score, 
                     predict.time =365*3, #3年生存率,也可改为3年生存率 365*3
                     method="KM")
#5年生存率,也可改为5年生存率 365*5
  ROC5<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$risk.score, 
                     predict.time =365*5, #5年生存率,也可改为5年生存率 365*5
                     method="KM")
#查看1，3，5年对应的AUC值
  ROC1$AUC;ROC3$AUC;ROC5$AUC
#绘制1年、3年和5年的ROC曲线
  plot(ROC1$FP, ROC1$TP, xlim = c(0, 1), ylim = c(0, 1), col = "#377eb8", type = "l",
       xlab = "False Positive Rate", ylab = "True Positive Rate", 
       main = "Time-dependent ROC Curve (LIHC)", lwd = 2, 
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)
#添加3年和5年的ROC曲线
  lines(x = ROC3$FP, y = ROC3$TP, lwd = 2, type = "l", col = "#00ba38", lty = 2)
  lines(x = ROC5$FP, y = ROC5$TP, lwd = 2, type = "l", col = "#fb8d3d", lty = 3)
#添加图例
  legend("bottomright", bty = "n",
         col = c("#377eb8", "#00ba38", "#fb8d3d"),
         legend = c(paste("AUC (1 year) =", round(ROC1$AUC, 2)), 
                    paste("AUC (3 year) =", round(ROC3$AUC, 2)), 
                    paste("AUC (5 year) =", round(ROC5$AUC, 2))),
         lty = c(1, 2, 3), lwd = 2, cex = 1, text.col = "black")
#添加对角线
  abline(0, 1, col = "gray", lty = 2)

#根据风险因子分高低组
  data$risk.score <- ifelse(data$risk.score>median(data$risk.score),"High Risk Score","Low Risk Score")
#保存风险因子高低组的样本数据，用于绘制桑基图
  write.csv(data,file = "含风险因子评分且带预后的TCGA所有样本表达谱，用于绘制桑基图.csv")
#将生存时间转换为月  
  data$OS.time <- data$OS.time/30
  fit <- survfit(Surv(OS.time,OS) ~ data$risk.score, data = data)
#作图
  ggsurvplot(fit,
             data = data,
             pval = T,
             conf.int = F, # 显示置信区间
             risk.table = TRUE, # 显示风险表
             risk.table.col = "strata",
             palette = "Set1", # 配色采用jco（黄蓝）
             legend.labs = c("High Risk Score","Low Risk Score"), # 图例
             size = 1,
             xlim = c(0,120), # x轴长度，一般为0-10年
             break.time.by = 20, # x轴步长为20个月
             legend.title = "Riskscore",
             surv.median.line = "hv", # 是否添加垂直和水平的中位生存
             ylab = "Survival probability (%)", # 修改y轴标签
             xlab = "Time (Months)", # 修改x轴标签
             ncensor.plot = T, # 显示删失图块
             ncensor.plot.height = 0.25,
             risk.table.y.text = FALSE)
  

  cox <- coxph(Surv(OS.time,OS)~risk.score,data)
  summary(cox)




#使用外部验证集
#加载必要的R包
  library(GEOquery)          #用于从GEO数据库获取数据
  library(limma)             #用于统计分析和数据标准化
  library(hthgu133a.db)      #用于将探针ID转换为基因符号
  library(hgu133a2.db)

#下载并加载GEO数据集GSE14520
gset <- getGEO('GSE14520',
               destdir=".", 
               AnnotGPL=FALSE, 
               getGPL=FALSE)

#获取第一个平台的表达矩阵
  expr_matrix_1 <- exprs(gset[[1]])    #使用exprs函数获取样本表达矩阵
#获取第一个平台的临床信息
  pd1 <- pData(gset[[1]])              #使用pData函数获取样本临床信息
#获取第二个平台的表达矩阵
  expr_matrix_2 <- exprs(gset[[2]])    #使用exprs函数获取样本表达矩阵
#获取第二个平台的临床信息
  pd2 <- pData(gset[[2]])              #使用pData函数获取样本临床信息

#转换第一个平台的探针ID到基因符号
#获取探针ID到基因符号的映射表
  probe_to_gene_1 <- toTable(hthgu133aSYMBOL) 
#过滤掉没有基因符号的行
  probe_to_gene_1 <- probe_to_gene_1[probe_to_gene_1$symbol != '',] 
#仅保留存在于表达矩阵中的探针ID
  probe_to_gene_1 <- probe_to_gene_1[probe_to_gene_1$probe_id %in% rownames(expr_matrix_1),]
#更新表达矩阵，仅保留有基因符号的探针
  expr_matrix_1 <- expr_matrix_1[probe_to_gene_1$probe_id,]
#将表达矩阵的行名更新为基因符号
  rownames(expr_matrix_1) <- probe_to_gene_1$symbol

#转换第二个平台的探针ID到基因符号
#获取探针ID到基因符号的映射表
  probe_to_gene_2 <- toTable(hgu133a2SYMBOL)
#过滤掉没有基因符号的行
  probe_to_gene_2 <- probe_to_gene_2[probe_to_gene_2$symbol != '',]
#仅保留存在于表达矩阵中的探针ID
  probe_to_gene_2 <- probe_to_gene_2[probe_to_gene_2$probe_id %in% rownames(expr_matrix_2),] 
#更新表达矩阵，仅保留有基因符号的探针
  expr_matrix_2 <- expr_matrix_2[probe_to_gene_2$probe_id,]
#将表达矩阵的行名更新为基因符号
  rownames(expr_matrix_2) <- probe_to_gene_2$symbol


#合并两个平台的表达矩阵，只保留共有的基因符号
  common_genes <- intersect(rownames(expr_matrix_1), rownames(expr_matrix_2))
#合并表达矩阵
  merged_expr_matrix <- cbind(expr_matrix_2[common_genes,], expr_matrix_1[common_genes,]) 

#标准化数据
  normalized_expr_matrix <- normalizeBetweenArrays(merged_expr_matrix)  
#绘制标准化后的箱线图，检查数据分布
  boxplot(normalized_expr_matrix, las=2)  

#保存标准化后的数据
  write.table(normalized_expr_matrix, file='normalized_GSE14520_data.txt', sep='\t', quote=FALSE)  # 保存为文本文件
  
#读取临床预后数据
  clinical_data <- read.table('GSE14520_Extra_Supplement.txt',sep = '\t',header = T,row.names = 1,check.names = F)
#筛选出肿瘤样本
  clinical_data <- clinical_data[clinical_data$`Tissue Type` == 'Tumor',]  
#去除含有缺失值的行
  clinical_data <- na.omit(clinical_data)  
#获取生存状态
  clinical_data$OS <- clinical_data$`Survival status` 
  clinical_data$OS.time <- clinical_data$`Survival months`
  rownames(clinical_data) <- clinical_data$Affy_GSM
#提取生存和预后信息
  clinical <- clinical_data[,c("OS","OS.time")]
#对表达矩阵进行转置
  normalized_expr_matrix <- t(normalized_expr_matrix)
#将临床数据与表达矩阵合并
  merged_data <- cbind(clinical_data[intersect(rownames(clinical_data), rownames(normalized_expr_matrix)),], normalized_expr_matrix)  
#选取相同项
  same <- intersect(row.names(clinical),row.names(normalized_expr_matrix))
#合并获得具有临床预后信息的表达谱
  rt <- cbind(clinical[same,],normalized_expr_matrix[same,])


#加入风险评分
  rt$risk.score <- 0.173030087*rt$KIF18A+
                   0.063695246*rt$SPP1+
                   0.014878703*rt$PLK1-
                   0.037960223*rt$SLC2A2+
                   0.001804856*rt$EZH2+
                   0.121272810*rt$MMP3
#输出含风险因子且带预后的基因表达谱
  write.csv(rt,file = "含风险因子评分且带预后的GSE14520样本表达谱.csv")
  rt <- read.csv("含风险因子评分且带预后的GSE14520样本表达谱.csv")
#将生存时间由月变为日
  rt$OS.time <- rt$OS.time*30
#1年生存率,也可改为3年生存率 365*3
  ROC1<- survivalROC(Stime=rt$OS.time, 
                     status=rt$OS, 
                     marker = rt$risk.score, 
                     predict.time =365*1, #1年生存率,也可改为1年生存率 365*1
                     method="KM")
#3年生存率,也可改为3年生存率 365*3
  ROC3<- survivalROC(Stime=rt$OS.time, 
                     status=rt$OS, 
                     marker = rt$risk.score, 
                     predict.time =365*3, #3年生存率,也可改为3年生存率 365*3
                     method="KM")
#5年生存率,也可改为5年生存率 365*5
  ROC5<- survivalROC(Stime=rt$OS.time, 
                     status=rt$OS, 
                     marker = rt$risk.score, 
                     predict.time =365*5, #5年生存率,也可改为5年生存率 365*5
                     method="KM")
#查看1，3，5年对应的AUC值
  ROC1$AUC;ROC3$AUC;ROC5$AUC
#绘制1年、3年和5年的ROC曲线
  plot(ROC1$FP, ROC1$TP, xlim = c(0, 1), ylim = c(0, 1), col = "#377eb8", type = "l",
       xlab = "False Positive Rate", ylab = "True Positive Rate", 
       main = "Time-dependent ROC Curve (GSE14520)", lwd = 2, 
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)
#添加3年和5年的ROC曲线
  lines(x = ROC3$FP, y = ROC3$TP, lwd = 2, type = "l", col = "#00ba38", lty = 2)
  lines(x = ROC5$FP, y = ROC5$TP, lwd = 2, type = "l", col = "#fb8d3d", lty = 3)
#添加图例
  legend("bottomright", bty = "n",
         col = c("#377eb8", "#00ba38", "#fb8d3d"),
         legend = c(paste("AUC (1 year) =", round(ROC1$AUC, 2)), 
                    paste("AUC (3 year) =", round(ROC3$AUC, 2)), 
                    paste("AUC (5 year) =", round(ROC5$AUC, 2))),
         lty = c(1, 2, 3), lwd = 2, cex = 1, text.col = "black")
#添加对角线
  abline(0, 1, col = "gray", lty = 2)




#根据风险因子分高低组
  rt$risk.score <- ifelse(rt$risk.score>median(rt$risk.score),"High Risk Score","Low Risk Score")
#保存风险因子高低组的样本数据，用于绘制桑基图
  write.csv(rt,file = "含风险因子评分且带预后的GSE14520表达谱，用于绘制桑基图.csv")
#将生存时间转换为月  
  rt$OS.time <- rt$OS.time/30
  fit <- survfit(Surv(OS.time,OS) ~ rt$risk.score, data = rt)
#作图
  ggsurvplot(fit,
             data = rt,
             pval = T,
             conf.int = F, # 显示置信区间
             risk.table = TRUE, # 显示风险表
             risk.table.col = "strata",
             palette = "Set1", # 配色采用jco（黄蓝）
             legend.labs = c("High Risk Score","Low Risk Score"), # 图例
             size = 1,
             xlim = c(0,120), # x轴长度，一般为0-10年
             break.time.by = 20, # x轴步长为20个月
             legend.title = "Riskscore",
             surv.median.line = "hv", # 是否添加垂直和水平的中位生存
             ylab = "Survival probability (%)", # 修改y轴标签
             xlab = "Time (Months)", # 修改x轴标签
             ncensor.plot = T, # 显示删失图块
             ncensor.plot.height = 0.25,
             risk.table.y.text = FALSE)

  cox <- coxph(Surv(OS.time,OS)~risk.score,rt)
  summary(cox)




#####多因素cox回归建模####
#读取含风险因子评分，亚型分组，预后信息的表达谱矩阵
data <- read.csv("含风险因子评分及亚型分组且带预后的TCGA所有样本表达谱.csv",header = T,sep = ",")
#将生存时间变为月
data$OS.time <- data$OS.time/30
#提取样本列,预后列，亚型分组列和风险得分列
data <- data[,c("X","risk.score","Group","OS","OS.time")]
#读取预后调整数据
clinical <- read.csv("TCGA预后数据第一次调整（含stage TNM age gender 分组）.csv",header = T,sep = ",")
#将样本名中-调整为.
clinical$samples <- gsub("-",".",clinical$submitter_id.samples)
rownames(data) <- data$X
rownames(clinical) <- clinical$samples
#提取两数据中列名与行名相同的取交集
comsamples <- intersect(rownames(data),rownames(clinical))
#展示comgene第14—16位
table(substr(comsamples,14,16))
#提取两文件中共同内容  
data <- data[comsamples,]
clinical <- clinical[comsamples,]
data <- cbind(data,clinical)
#保存数据以用来进行单因素cox回归
write.csv(data,"单因素cox回归含临床信息和生存数据.csv")

data <- read.csv("单因素cox回归含临床信息和生存数据.csv",header = T,sep = ",")
#进行单因素cox回归
#风险因子
cox_risk.score <- coxph(Surv(OS.time,OS)~risk.score,data)
#年龄
cox_Age <- coxph(Surv(OS.time,OS)~Age,data)
#性别
cox_gender <- coxph(Surv(OS.time,OS)~Gender,data)
#分期
cox_Stage <- coxph(Surv(OS.time,OS)~Stage,data)
#TNM分期
cox_T <- coxph(Surv(OS.time,OS)~T,data)
cox_N <- coxph(Surv(OS.time,OS)~N,data)
cox_M <- coxph(Surv(OS.time,OS)~M,data)
#病理分期
cox_Grade <- coxph(Surv(OS.time,OS)~Grade,data)
#亚型分组
cox_Group <- coxph(Surv(OS.time,OS)~Group,data)

#查看单因素cox回归结果
summary(cox_risk.score)
summary(cox_Age)
summary(cox_gender)
summary(cox_Stage)
summary(cox_T)
summary(cox_N)
summary(cox_M)
summary(cox_Grade)
summary(cox_Group)

#对符合单因素cox回归的因子进行多因素cox回归
cox <- coxph(Surv(OS.time,OS)~data$risk.score+data$Group+
               data$Stage,data)
#查看多因素cox结果
summary(cox)
#使用逐步回归法验证
data <- na.omit(data)
cox <- coxph(Surv(OS.time,OS)~risk.score+Stage,data)
summary(cox)
fit.step <- step(cox)

##################绘制nomogram图
library(foreign)
library(survival)
library(rms)
library(Hmisc)
dd<-datadist(data)
options(datadist='dd')
coxm1 <- cph(Surv(OS.time,OS)~risk.score+Stage
               ,data=data,x=T,y=T,surv = T)
coxm1
surv <- Survival(coxm1)### 用下面函数算不同时间节点的生存率

surv1 <- function(x)surv(1*12,lp=x)
surv2 <- function(x)surv(3*12,lp=x)
surv3 <- function(x)surv(5*12,lp=x)

#创建nomogram对象
nom1 <- nomogram(
  coxm1, # Cox回归模型
  fun = list(surv1, surv2, surv3), # 需要绘制的生存曲线
  lp = F, # 是否显示线性预测变量
  funlabel = c('1-year Survival Rate', '3-year survival Rate', '5-year survival Rate'), # 生存曲线对应的标签
  maxscale = 100, # 画布最大刻度
  fun.at = c('0.9', '0.85', '0.80', '0.70', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1') # 生存曲线的刻度
)
#绘制nomogram图10x6
plot(nom1)

##################绘制校准曲线
# 定义函数
plot_calibration <- function(cal, time, col) {
  plot(cal, lwd=2, lty=1, errbar.col=col, xlim=c(0,1), ylim=c(0,1),
       xlab=paste0("Nomogram-Predicted Probabilityof ", time, " years OS"),
       ylab=paste0("Actual ", time, " years OS (proportion)"), col=col)
  lines(cal[,c("mean.predicted","KM")], type="b", lwd=2, col=col, pch=16)
}
# 定义常用参数
cmethod <- 'KM'
method <- 'boot'
m <- 100
B <- 1000
# 计算校准曲线
coxm1 <- cph(Surv(OS.time,OS)~risk.score+Stage,data=data,x=T,y=T,time.inc=12,surv=T)
cal1 <- calibrate(coxm1, cmethod=cmethod, method=method, u=12, m=m, B=B)
coxm2 <- cph(Surv(OS.time,OS)~risk.score+Stage,data=data,x=T,y=T,time.inc=36,surv=T)
cal2 <- calibrate(coxm2, cmethod=cmethod, method=method, u=36, m=m, B=B)
coxm3 <- cph(Surv(OS.time,OS)~risk.score+Stage,data=data,x=T,y=T,time.inc=60,surv=T)
cal3 <- calibrate(coxm3, cmethod=cmethod, method=method, u=60, m=m, B=B)
# 绘制校准曲线
plot(cal1, lwd=2, lty=1, errbar.col="#cb6767", xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of OS",
     ylab="Actual OS (proportion)", col="#cb6767", main="Calibration Plot for 1-, 3- and 5-Year Survival")
lines(cal1[,c("mean.predicted","KM")], type="b", lwd=2, col="#cb6767", pch=16)
plot(cal2, lwd=2, lty=1, errbar.col="#4775b4", xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of OS", add=T,
     ylab="Actual OS (proportion)", col="#4775b4")
lines(cal2[,c("mean.predicted","KM")], type="b", lwd=2, col="#4775b4", pch=16)
plot(cal3, lwd=2, lty=1, errbar.col="#308036", xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of OS", add=T,
     ylab="Actual OS (proportion)", col="#308036")
lines(cal3[,c("mean.predicted","KM")], type="b", lwd=2, col="#308036", pch=16)
# 添加图例和对角线
legend("bottomright", bty="n",
       fill=c("#cb6767","#4775b4","#308036"),legend = c("1 year", "3 years", "5 years"),
       cex=.8, border=NA, y.intersp=1.0, x.intersp=0.3)
abline(0,1,lty=3,lwd=2,col="#4d4d4d")



##################绘制DCA曲线
# 加载需要的包
library(survival)
library(rms)
library(ggDCA)
#读取数据文件
data <- read.csv("单因素cox回归含临床信息和生存数据.csv",header = T,sep = ",")
# 拟合多因素Cox模型
Stage <- cph(Surv(OS.time, OS) ~ Stage, data)
Risk.score <- cph(Surv(OS.time, OS) ~ risk.score, data)
Stage_Risk.score <- cph(Surv(OS.time, OS) ~ Stage + risk.score, data)
# 进行DCA分析
#DCA总体时间
dt <- dca(Stage, Risk.score, Stage_Risk.score)
#dt <- dca(Stage, Risk.score, Stage_Risk.score, times = 12（12月-一年）对应时间)
dt1 <- dca(Stage, Risk.score, Stage_Risk.score, times = 12)
dt3 <- dca(Stage, Risk.score, Stage_Risk.score, times = 36)
dt5 <- dca(Stage, Risk.score, Stage_Risk.score, times = 60)

# 绘制总体DCA曲线图
Fa <- ggplot(dt) + 
  labs(title = "Decision Curve Analysis", 
       subtitle = "Comparison of Three Survival Models") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5))

print(Fa)#保存8x6

# 绘制一年DCA曲线图
F1 <- ggplot(dt1) + 
  labs(title = "1-Year Decision Curve Analysis", 
       subtitle = "Comparison of Three Survival Models") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5))

print(F1)

# 绘制三年DCA曲线图
F3 <- ggplot(dt3) + 
  labs(title = "3-Year Decision Curve Analysis", 
       subtitle = "Comparison of Three Survival Models") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5))

print(F3)

# 绘制五年DCA曲线图
F5 <- ggplot(dt5) + 
  labs(title = "5-Year Decision Curve Analysis", 
       subtitle = "Comparison of Three Survival Models") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5))

print(F5)




####肿瘤突变负荷TMB####
#TCGA突变数据下载
#网址：https://portal.gdc.cancer.gov/
library(TCGAbiolinks)
#安装maftools
if (!require("BiocManager"))
  install.packages("BiocManager")
#若无法下载运行此代码
options(download.file.method = 'libcurl')
options(url.method='libcurl')
BiocManager::install("maftools")
library(maftools)
library(tidyverse)
setwd("TMB")
library(readxl)
library(readr)
#突变数据下载后，整理成数据框
mut2 <- read.maf("TCGA.LIHC.varscan.40fe9c1b-19d0-45cf-898a-f7b0cbad783e.
                   DR-10.0.somatic.maf")

a <- mut2@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:20] # 修改数字，代表TOP多少
save(mat,mat_0_1,file = "TMB-LIHC.rda")
#以下绘图代码来自解螺旋阿琛老师
oncoplot(maf = mut2,
         top = 30,   #显示前30个的突变基因信息
         fontSize = 0.6,   #设置字体大小
         showTumorSampleBarcodes = F)   #不显示病人信息
dev.off()
####计算TMB####
maf = tmb(maf = mut2,
          captureSize = 50,
          logScale = TRUE)   
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
maf$sample <- gsub("-",".",maf$sample)
rownames(maf) <- maf$sample
#手动更改后导入
write.csv(maf, file = "maf.csv")
com <- intersect(rownames(maf),rownames(group))
maf <- maf[com,]
group$x <- as.character(group$group)
group <- group %>% t() %>% as.data.frame()
group <- group[-1,]
group <- group[,com]
group <- group %>% t() %>% as.data.frame()
identical(rownames(group),rownames(maf))
tmb <- cbind(maf,group)
write.csv(tmb, file = "tmb.csv")




####NORGRAM####
#设置工作目录  
setwd("nomogram")
#加载包  
library(tidyverse)
#读取生存信息
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#将sample中_变为.  
surv$sample <- gsub("-",".",surv$sample)
#将第一列作为行名  
rownames(surv) <- surv$sample
#去除无用列  
surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,
                   stringsAsFactors = F,header = T)
#提取两数据中列名与行名相同的取交集
comgene <- intersect(colnames(expr),rownames(surv))
#展示comgene第14—16位
table(substr(comgene,14,16))
#提取两文件中共同内容  
expr <- expr[,comgene]
surv <- surv[comgene,]
#双击差异分析文件
load("LIHC_DEG.Rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)#根据自己需要
#提取变化2倍的差异基因表达谱，进行后续lasso
expr <- expr[rownames(res_deseq2),]
expr <- expr %>% t() %>% as.data.frame()
#整合表达谱与生存信息
identical(rownames(expr),rownames(surv))
#按列进行合并 <- cbind  
expr_surv <- cbind(surv,expr)
#lasso
#加载包  
library("glmnet")
library("survival")
#将第一二生存时间与生存状态统一为lasso代码中所用  
colnames(expr_surv)[1] <- 'fustat'
colnames(expr_surv)[2] <- 'futime'
#将expr_surv赋值为rt
rt <- expr_surv         
#将天数改为年数  
rt$futime=rt$futime/365   
#设计随机数种子
set.seed(3)   
#设置x为所提取第三行开始到结尾且设为矩阵形式
x=as.matrix(rt[,c(3:ncol(rt))])
#设置y为生存信息  
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
#作图  
plot(fit, xvar = "lambda", label = TRUE)
cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)
#其中两条虚线分别指示了两个特殊的λ值 <- 左侧为最佳，右侧为最简  
plot(cvfit)
dev.off()

#输出预测模型的相关系数与riskScore
#输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系数
#计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("futime", "fustat", lassoGene)
#根据风险评估中位数分为高组低组  
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
#绘制散点分布图
#读取包  
library(ggpubr)  
#作图  
p <- ggboxplot(dat, x = "fustat", y = "riskScore",
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果
#判断预测结果的准确性 <- roc根据风险评估做roc预测
#读取包  
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)
pred <- prediction(dat$riskScore, dat$fustat)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC 
#看data AUC中 y vaule中值
AUC<- 0.6873922 #改 根据结果对AUC进行赋值 <- 0.6 0.7以上最好
plot(perf,
     colorize=FALSE, 
     col="red", 
     print.auc =TRUE,
     main=paste("AUC=",AUC)) #绘制ROC曲线
#绘制对角线  
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
dev.off()

#临床信息下载整理
#网址：https://xenabrowser.net/datapages/   下载后解压缩放入当前文件夹
#导入临床信息（不要行名） TCGA-LIHC.GDC_phenotype <- 命名为行名  heading <- yes
#将样本名—改为.  
clinical$submitter_id.samples <- gsub("-",".",clinical$submitter_id.samples)
#将clinical中样本名称变为行名，并去掉样本名行
rownames(clinical) <- clinical$submitter_id.samples
clinical <- clinical[,-1]
#结合dat与clinical
#取clinical和dat行名交集  
comgene1 <- intersect(rownames(clinical),rownames(dat))
#取dat中列名为交集  
dat <- dat[comgene1,]
#取clinical中行名为交集  
clinical <- clinical[comgene1,]
#去掉无用列，只要生存状态，时间，风险高低 
dat <- dat[,-(3:13)]
#提取想要的肿瘤因素  
clinical <- clinical[,c("gender.demographic","tumor_stage.diagnoses",
                        "pathologic_T","pathologic_N","pathologic_M")]
#clinical <- clinical1
#判断顺序是否相同  
identical(rownames(dat),rownames(clinical))
#合并输出为rt2  
rt2 <- cbind(dat,clinical)
write.csv(rt2, file = "rt2.csv")
#导入修改后的rt3 <- 去掉not report ;将罗马数字改为阿拉伯数字 ;将男女改为MF
rt2 <- rt3
#加载包
library(rms)
library(foreign)
library(survival)
#设置参数，变为因子型
rt2$gender <- factor(rt2$gender,labels=c("F", "M"))#gender <- 性别
rt2$stage <- factor(rt2$stage,labels=c("Stage1", 
                                       "Stage2",
                                       "Stage3", 
                                       "Stage4")) #stage<- 分期
rt2$T <- factor(rt2$T,labels=c("T1", 
                               "T2", 
                               "T3",
                               "T4"))
rt2$M <- factor(rt2$M,labels=c("M0",
                               "M1"))
rt2$N <- factor(rt2$N,labels=c("N0",
                               "N1"))
rt2$risk <- factor(rt2$risk,labels=c("low", 
                                     "high"))
ddist <- datadist(rt2)
options(datadist='ddist')   #使用函数datadist()将数据打包

#构建Cox回归模型 <- 根据gender + stage +T + M + N + risk预测生存
f <- cph(Surv(futime, fustat) ~gender + stage +T + M + N + risk,
         x=T, y=T, surv=T, data=rt2, time.inc=1)
surv <- Survival(f)

#构建Nomogram
nom2 <- nomogram(f, fun=list(function(x) surv(1, x), 
                             function(x) surv(2, x), 
                             function(x) surv(3, x)), #对应123年
                 lp=F, funlabel=c("1-year survival", 
                                  "2-year survival", 
                                  "3-year survival"), 
                 maxscale=100, 
                 fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 
                          0.4, 0.3,0.2,0.1,0.05))
plot(nom2)
dev.off()




####免疫检查点####
#设置工作目录  
setwd("PD1")
#加载包
library(tidyverse)
#读取表达谱  
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,
                   stringsAsFactors = F,header = T)
#读取分组分析  
load("group_AY.rda")
#读取并筛选其中分组为1的命名为a 分组为2的命名为b   
a <- group %>% filter(group == "1")
b <- group%>% filter(group == "2")
#读取与ab相对应的表达谱  
expr1 <- expr[,rownames(a)]
expr2 <- expr[,rownames(b)]
#设置免疫检查点 #https://www.genecards.org/ 查询免疫检查点全名  
gene <- c("PDCD1","CD274","PDCD1LG2")
#表达谱中免疫检查点  
gene1 <- expr1[gene,]
gene2 <- expr2[gene,]
#t转换并作为数据形式保存
gene1 <- gene1 %>% t() %>% as.data.frame()
gene2 <- gene2 %>% t() %>% as.data.frame()
#保存
write.csv(gene1, file = "gene1.csv")
write.csv(gene2, file = "gene2.csv")









####X cell####
#安装devtools 右下角install
##安装Rtools：
#b站教程视频号：BV1bk4y1k7Nd 
#下载网址windows版本 https://cran.r-project.org/bin/windows/Rtools/
#验证Rtools安装是否成功
system("g++ -v")      #显示[1] 0成功
system("where make")  #显示[1] 0成功
#安装xCell包
chooseBioCmirror()
devtools::install_github('dviraran/xCell')
#加载包  
library(xCell)
library(ggpubr)
library(tidyverse)
#设置工作目录  
setwd("Xcell")
#查看细胞类型
celltypeuse<-xCell.data$spill$K
#读取表达谱  
exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,
                  stringsAsFactors = F,header = T)
#计算免疫浸润
rs<-xCellAnalysis(exp,parallel.sz=10) 
#准备分组信息
gene <- "PDCD1"#每次运行只改这个基因名
#将表达谱中基因作为数值提取并取其中位数  
med=median(as.numeric(exp[gene,]))
#提取对应基因表达谱列  
expgene <- exp[gene,]
#t转换作为数据形式  
expgene <- expgene %>% t() %>% as.data.frame()
#分高低组
expgene$group <- ifelse(expgene$PDCD1>med,"High","Low") #改基因名称
#整理免疫浸润格式  
rs <- as.data.frame(rs)
rs <- rs %>% t() %>% as.data.frame()
#取交集  
comname <- intersect(rownames(rs),rownames(expgene)) 
#提取两者共同拥有内容  
rs <- rs[comname,]
expgene <- expgene[comname,]
#判断顺序
identical(rownames(rs),rownames(expgene))
#赋予rs-group列  
rs$group <- as.factor(expgene$group)
#判断类型  
class(rs$group)
#行名转换为列并命名为sample
rs <- rs %>% rownames_to_column("sample")
#
a <- rs
#加载包
library(ggsci)
library(tidyr)
library(ggpubr)
#作图  
b <- gather(a,key=xCell,value = Expression,-c(group,sample))
ggboxplot(b, x = "xCell", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()




####WGCNA####
#安装包
install.packages("BiocManager")
BiocManager::install("preprocessCore")
BiocManager::install("impute")
install.packages("WGCNA")
#加载包  
library("tidyverse")
library("WGCNA")  
library(DESeq2)
#设置工作目录  
setwd("WGCNA")
#读取数据  
counts_01A <- read.table("LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,
                         check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("LIHC_fpkm_mRNA_01A.txt", sep = "\t",row.names = 1,
                  check.names = F,header = T)
#取两者列名交集
com <- intersect(colnames(counts_01A),colnames(exp))
#按照交集顺序排序  
exp <- exp[,com]
counts_01A <- counts_01A[,com]
#判断是否相同  
identical(colnames(counts_01A),colnames(exp))
#差异分析
gene <- "PDCD1"#每次运行只改这个基因名
#取对应基因表达中位数  
med=median(as.numeric(exp[gene,]))
#设置分组信息  
conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = 
                                     c("low","high"))) %>% 
  column_to_rownames("sample")
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
#确认结果分组比较  
resultsNames(dds)
#提取并保存结果  
res <- results(dds)
save(res,file="res_deseq2_PDCD1.Rda")
#整理输入表达谱
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)#阈值
#提取出符合阈值的表达谱 
input <- exp[rownames(DEG),]
#保证样本在行，基因在列！！！！
datExpr0 = as.data.frame(t(input))
#开始WGCNA
#检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#如果没有达标就需要筛选
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], 
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], 
                                                collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 样品聚类
# 聚类
sampleTree = hclust(dist(datExpr0), method = "average")
# 画图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#设置裁剪高度
high=67
#剪切线 <- 根据y轴位置
abline(h = high, col = "red")
#删除剪切线以下的样品 <- 注意y位置
clust = cutreeStatic(sampleTree, cutHeight = high, minSize = 10)
#展示分组  
table(clust)
#保留希望所留分组  
keepSamples = (clust==1)
#提取刚才切除所留部分  
datExpr0 = datExpr0[keepSamples, ]
dev.off()
#报错Error in hclust(dist(datExpr0), method = "average") : NaN dissimilarity value.为存在缺失值
adj = adjacency(datExpr0)
datExpr0<-datExpr0[complete.cases(datExpr0),]
#重新聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)
#记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr0)#基因数
nSamples = nrow(datExpr0)#样本数
#保存
save(datExpr0, nGenes, nSamples,file = "Step01-WGCNA_input.Rda")

#构建网络，识别模块
#power值散点图
enableWGCNAThreads()   #多线程工作!!!配置低别用这个!!!炸了不管hhhh
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9
#拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
#平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值看右侧
softPower = 2#根据右侧softpower处理
adjacency = adjacency(datExpr0, power = softPower)
#TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(TOM,file = "TOM.Rda")
#基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#动态剪切模块识别
minModuleSize = 30      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1 #剪切高度可修改
abline(h=MEDissThres, col = "red")
#相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors,
                          cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()

#整理临床信息
clinical <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",
                       row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clinical <- clinical[rownames(datExpr0),]
identical(rownames(clinical),rownames(datExpr0))
#查看临床信息
head(clinical)
#对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)
#对样本进行聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
#将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)
#样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#网络的分析
#对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")
write.table(file="Step04-modPhysiological.cor.xls",
            moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="Step04-modPhysiological.p.xls",
            moduleTraitPvalue,sep="\t",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化。
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.7,
               cex.lab=0.7,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# 不同模块与基因性状的具体分析
#矩阵一
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
#看一下目的基因和哪个模块相关性最高
a <- geneModuleMembership
a <- a %>% rownames_to_column()

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#矩阵二
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

#批量输出性状和模块散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10. 输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",
              row.names=F,col.names=F,quote=F)
}




####GEO数据库的下载####
#加载GEO下载包
  BiocManager::install("GEOquery")
library(GEOquery)
####下载GSE数据（可能存在下载过慢报错则使用网页下载）
#安装github包GEOmirror和AnnoProbe
remotes::install_github("jmzeng1314/GEOmirror")
remotes::install_github("jmzeng1314/AnnoProbe")
#三个包同时加载
library(AnnoProbe)
library(GEOmirror)
library(GEOquery) 
#下载获取GSE13507数据 
gset=AnnoProbe::geoChina('GSE14520')
gset
#提取gset中平台2测序数据（GPL571人类测序）
gpleSet=gset[2] 
eSet=gpleSet[[1]] # 提取表达矩阵 
probes_expr <- exprs(eSet)
dim(probes_expr) 
getAssayData(eSet)
probes_expr=log2(probes_expr+1) #表达矩阵行log2+1归一化处理 
# 提取表型数据信息 
phenoDat <- pData(eSet) 

#  读取GEO下载的矩阵和注释
####  加载包

library(readxl)
library(tidyverse)
library(GEOquery)

# 读取 series_matrix文件 
GSE_data <- getGEO(filename = "GSE14520-GPL571_series_matrix.txt.gz", getGPL = F)

pd_GSE_data <- pData(GSE_data) #观察临床信息中的data processing,Microarray suite,MAS 5.0，即标准化方法，如果已经经过了MAS，一般情况下不需要再标准化了，避免矫枉过正
pd_GSE_data$data_processing[1]  #查看数据处理的方式
GSE_data_targets <-  pd_GSE_data %>%  # 提取有用的信息以便后续分析
  dplyr::select(sample_id = geo_accession,
                sample_name = title,
                tissue_type = `Tissue:ch1`)  


GSE_data_expr <- exprs(GSE_data) %>% 
  as.data.frame()  # 提取表达矩阵
# 数据集的标准化结果查看
boxplot(exprs(GSE_data))  ## 箱线图查看数据是否标准化
## PCA  查看不同样本间的区别
PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  rv <- genefilter::rowVars(expr)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(expr[select, ]))#最核心的代码
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  group = group, 
                  name = colnames(expr))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
                      size = 3,
                      segment.color = "black",
                      show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}
PCA_new(exprs(GSE_data), 
        nrow(GSE_data),
        group = GSE_data_targets$tissue_type)  # 组别以颜色显示

## getGEO 读取注释信息
##　AnnotGPL = T，用的是.gz的文件
GPL_data<- getGEO(filename ="GPL571.annot.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)

## AnnotGPL = F，用的是.soft的文件
GPL_data_2 <- getGEO(filename = "GPL10687_family.soft.gz", AnnotGPL = F)
GPL_data_21 <- Table(GPL_data_2)

# 进行注释
GSE_data_expr <- GSE_data_expr %>% 
  rownames_to_column() %>% 
  reshape::rename(c(rowname = "prode_id"))

GPL <- GPL_data_11 %>% 
  dplyr::select(ID,`Gene symbol`) %>% 
  reshape::rename(c(ID = "prode_id",`Gene symbol` = "GeneSymbol")) %>% 
  add_count(GeneSymbol,name = "n_gene") %>% 
  add_count(prode_id,name = "n_prode_id") %>% 
  dplyr::filter(n_gene < 1123) %>% # 选择能有对应基因的探针
  dplyr::select(c(-3,-4))

GSE_expr <- GSE_data_expr %>% 
  inner_join(GPL, ., by = "prode_id") %>% #p2s和上一步得到的结果再取交际，p2s放在右边是以它为准
  dplyr::select(-1) %>%  #去除第一列probe_id
  as.data.frame() %>% #因为aggregate需要数据框格式
  aggregate(.~ GeneSymbol, data = ., mean) %>%  #以symbol作为因子水平，把相似的数据放在一起取均值，最大值max，中位值median
  column_to_rownames(var = "GeneSymbol")

## 可视化注释后的矩阵，因为注释过程丢失了部分探针
boxplot(GSE_expr)
PCA_new(GSE_expr, 
        nrow(GSE_expr),
        group = GSE_data_targets$tissue_type)  # 组别以颜色显示
# 以RData 格式保存信息 
GSE14520_GPL571 <- GSE_data
GSE14520_GPL571_expr <- GSE_expr
GPL571 <- GPL_data_11
GSE14520_GPL571_targets <- GSE_data_targets

picDir <- "dealed_data"
dir.create(picDir )
save(GSE14520_GPL571,  # 矩阵
     GSE14520_GPL571_targets, # metadata
     GPL571, #注释信息
     file = paste('.',picDir,"GEOqury读取的GSE14520-GPL571_data.RData",sep="/"))
save(GSE14520_GPL571_expr,
     file = paste('.',picDir,"GSE14520-GPL571_expr.RData",sep="/"))





####整理临床预后信息####
#读取整理后预后信息
LIHC_phenotype <- read.csv("预后信息.csv")
#gsub函数替换sLIHC_phenotype中sample的—换为.
LIHC_phenotype$sample <- gsub("-",".",LIHC_phenotype$submitter_id.samples)
#读取分组结果
group <- read.csv("生存大于30分组结果.csv")
#将sample换为行名
rownames(LIHC_phenotype) <- LIHC_phenotype$sample
#将患者名变为行名
rownames(group) <- group$X
#提取group和LIHC_phenotype中相同行名
com<- intersect(rownames(group),rownames(LIHC_phenotype)) 
#筛选出group中行名与LIHC_phenotype中相同的
group <- group[com,]
#筛选出LIHC_phenotype中行名与group中相同的
LIHC_phenotype <- LIHC_phenotype[com,]
#判断顺序是否相同
identical(rownames(group),rownames(LIHC_phenotype))
#LIHC_phenotype中新增group
LIHC_phenotype$group <- as.character(group$Group)
#判断数据形式  
class(LIHC_phenotype$group)
#展示分组情况  
table(LIHC_phenotype$group)
#输出为表格
write.csv(LIHC_phenotype, file = "LIHC_phenotype_筛选后大于30天含分组数据用于热图绘制.csv")




####对ICGC数据进行整理####
#载入所需R包：
library(stringr)
library(tidyverse)
#载入表达矩阵：
exp <- data.table::fread('exp_seq.tsv.gz')
#1.id串联：
#先查看一下表达矩阵数据格式：
dim(exp)
View(exp) #以宽数据形式进行存储
#提取表达矩阵中所需列：
col <- colnames(exp)
col
exp1 <- exp[,c(
  'icgc_specimen_id', #样品id
  'icgc_donor_id', #捐赠者/患者id
  'gene_id', #ENS id
  'normalized_read_count', #RPKM
  'raw_read_count' #原始read counts
)]
head(exp1)
#把样品id和患者id串联，组合成新样本名(同一个捐赠者可能存在多个癌或癌旁样品)：
exp2 <- tidyr::unite(exp1, "icgc_specimen_id _ icgc_donor_id", icgc_specimen_id, icgc_donor_id)
colnames(exp2)[1] <- c('sample_id')
head(exp2)
#可以按照需求选择RPKM或counts列：
exp_r <- exp2[,-3] #原始counts，用于差异分析

exp3 <- exp_r %>% group_by(sample_id) %>% mutate(id=1:n()) %>%
  spread(sample_id,raw_read_count) %>% select(-"id")
View(exp3) #格式转换完成
#将id作为行名：
exp4 <- as.data.frame(exp3)
rownames(exp4) <- exp4$gene_id
exp4 <- exp4[,-1]
exp_r <- exp4
View(exp_r) #表达矩阵初步整理完成
#将id作为行名：
exp4 <- as.data.frame(exp3)
rownames(exp4) <- exp4$gene_id
exp4 <- exp4[,-1]
exp_r <- exp4
View(exp_r) #表达矩阵初步整理完成
#保存所需数据并清空工作空间，准备下一步整理
save(exp_r,exp_n,file = c('ICGC_exp.Rdata'))
rm(list=ls())
#3.创建tumor/normal分组：
#载入上一步整理表达矩阵：
load('ICGC_exp.Rdata')
#读取有分组信息的specimen文件：
group <- data.table::fread('ICGC/specimen.tsv')
View(group)
#筛选列名中所需列：
colnames(group)
group1 <- group[,c(
  'icgc_specimen_id', #样品id
  'icgc_donor_id', #捐赠者/患者id
  'specimen_type' #样品类型，癌或癌旁组织
)]
head(group1)
#将样本和患者id串联(和表达矩阵清洗时相同)：
group2 <- tidyr::unite(group1, "icgc_specimen_id _ icgc_donor_id", icgc_specimen_id, icgc_donor_id)
colnames(group2)[1] <- c('sample_id')
head(group2)
table(group2$specimen_type) #可以看到normal有2个来源
#将group2样本id的顺序按照exp_r列名进行匹配：
group3 <- group2[match(colnames(exp_r),group2$sample_id),]
dim(group3) #匹配完后从原本264个样本减至136个，数量和表达矩阵相同
#检查二者是否完全匹配对应：
identical(colnames(exp_r),group3$sample_id) #若逻辑值为FALSE，需要检查上面的步骤是否有误
#简化分组信息(直接命名tumor和normal)：
group3$group <- ifelse(group3$specimen_type==c('Primary tumour - solid tissue'),'tumor','normal')
group <- group3[,-2]
table(group$group) #91个tumor，45个normal样本
head(group)
#保存所需数据并清空工作空间，准备下一步整理
save(exp_r,exp_n,group,file = c('ICGC_exp_group.Rdata'))
rm(list=ls())
#4.geneid转换：
#读取gtf注释：
library(rtracklayer)
gtf <- rtracklayer::import("gencode.v41lift37.annotation.gtf")
gtf <- as.data.frame(gtf)
View(gtf)
#从列名中提取所需列(gene_id、gene_type、gene_name)：
colnames(gtf)
gname <- subset(gtf,select = c(10:12))#提取第10,11,12列
head(gname)
#ENSid去重：
table(duplicated(gname$gene_id)) #300w+都是重复的id，需要去掉
gname <- gname[!duplicated(gname$gene_id),]
table(duplicated(gname$gene_id)) #去重完成
head(gname)
#提取mRNA和lncRNA：
lncRNA <- c(
  'lncRNA',
  '3prime_overlapping_ncRNA',
  'antisense',
  'bidirectional_promoter_lncRNA',
  'lincRNA',
  'macro_lncRNA',
  'non_coding',
  'processed_transcript',
  'sense_intronic',
  'sense_overlapping'
)
mRNA <- c("protein_coding")
#分别生成id_list：
lncRNA_list <- gname[gname$gene_type %in% lncRNA,]
mRNA_list <- gname[gname$gene_type %in% mRNA,]
head(lncRNA_list)
dim(lncRNA_list) #18146个lncRNA
head(mRNA_list)
dim(mRNA_list) #20149个mRNA
#保存所需数据并清空工作空间，准备下一步整理：
save(lncRNA_list,mRNA_list,file = c('gene_ID_list_37.Rdata'))
rm(list=ls())
#读取ID list和表达矩阵RData文件：
load('ICGC_exp_group.Rdata')
load('gene_ID_list_37.Rdata')
#合并mRNA和lncRNA(因不需要区分RNA类型所以合并，需要区分则选择所需即可)：
id <- rbind(mRNA_list,lncRNA_list)
head(id)
#ICGC表达矩阵中ENSid没有版本号，而我们的list中的ENSid是带有版本号的，因此需要去除(小数点之后)：
head(rownames(exp_r)) #表达矩阵ENSid无版本号
head(id$gene_id) #id list有版本号
id$ENSid <- unlist(str_split(id$gene_id,"[.]",simplify=T))[,1]
head(id)
#取ENSid交集:
RNA <- intersect(rownames(exp_r),id$ENSid)
length(RNA) #29511个
#表达矩阵取RNA交集:
exp2 <- exp_r[RNA,]
#表达矩阵添加symbol列：
#按表达矩阵行名的顺序来匹配id列表中ENSid的顺序，对id列表进行顺序重排:
a <- rownames(exp2)
b <- id$ENSid
id_match <- id[match(a,b),]
#检查match后是否一一匹配:
identical(a,id_match$ENSid)
#将symbol列和表达矩阵合并:
exp3 <- cbind(id_match$gene_name,exp2)
colnames(exp3)[1] <- c("gene_symbol")
View(exp3)
#5.symbol去重(去除低表达重复基因)并设为行名（行名不能重复）：
#计算每个gene平均表达量avg；
x <- exp3[,2:length(colnames(exp3))]
avg <- apply(x, 1, mean) #对每一行求平均表达量
exp4 <- as.data.frame(cbind(exp3,avg))#合并表达矩阵和均值列
#根据avg列对表格降序；
exp5 <- arrange(exp4,desc(avg))
#根据symbol列过滤重复行；
exp6 <- exp5 %>% distinct(gene_symbol, .keep_all = TRUE)
#移除均值列：
exp7 <- exp6[,1:(length(exp6)-1)]
#修改行名并去掉第一列：
rownames(exp7) <- exp7$gene_symbol
exp7 <- exp7[,-1]
dim(exp7) #symbol去重后从29511减到29496个
exp_ens <- exp3 #过滤前，ENSid为行名，含重复symbol列
exp_sym <- exp7 #过滤后表达矩阵，symbol为行名
#保存所需数据并清空工作空间：
save(group,exp_ens,exp_sym,file = c('ICGC_exp_group2.Rdata'))
write.csv(exp_sym,file = 'ICGC_exp_all.csv')
rm(list=ls())


####渐变火山图####
#加载包
library(ggplot2)
#加载数据
data <- read.csv("LIHC_deseq2_all.csv",row.names = 1)
label <- read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因 用于火山图基因标注.csv",col.names = T)
# 检查行数是否匹配
if (nrow(data) != nrow(label)) {
  # 创建一个新的标签列，初始值为NA
  data$label <- NA
  # 将匹配的标签赋值给新的标签列
  matching_rows <- match(rownames(data), label$TRUE.)
  data$label[!is.na(matching_rows)] <- label$TRUE.[matching_rows[!is.na(matching_rows)]]
} else {
  # 行数匹配，直接赋值给标签列
  data$label <- label$TRUE.
}

ggplot(data,aes(log2FoldChange, -log10(padj)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  # 调整主题和图例位置：
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.7),
        legend.justification = c(0,1)
  )+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none")+
  # 添加标签：
  geom_text(aes(label=label, color = -log10(padj)), size = 4, vjust = 2.5, hjust=2)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)") +
  geom_segment(aes(x = log2FoldChange, xend = log2FoldChange, y = -log10(padj), yend = -log10(padj) + 0.5),
               color = "black", alpha = 0.5, size = 0.5)

# 保存图片：
ggsave("vocanol_plot.pdf", height = 9, width = 10)


ggplot(data, aes(log2FoldChange, -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(padj), color = -log10(padj))) +
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none") +
  geom_text(aes(label = label, color = -log10(padj)), size = 4, vjust = 2.5, hjust = 2) +
  geom_segment(aes(x = log2FoldChange, xend = log2FoldChange, y = -log10(padj), yend = -log10(padj) + 0.5),
               color = "black", alpha = 0.5, size = 0.5) +
  xlab("Log2FC") +
  ylab("-Log10(FDR q-value)") 

  ggsave("volcano_plot.pdf", height = 9, width = 10)

  #读取全部基因coxoutput
  Coxoutput <- read.table("cox results all.txt",sep = "\t",row.names = 1,
                          check.names = F,
                          stringsAsFactors = F,header = T)
  
  
  #读取想要基因集
  gene = read.csv("GENE126新（2023.3.16）.csv",header = T)
  #改行名并取交集
  rownames(gene) <- gene$GENE
  coxgene <- intersect(rownames(gene),rownames(Coxoutput))
  #得到想要基因的cox
  Coxoutput <- Coxoutput[coxgene,]


  
  
  
  
  
  
####三维PCA(成功)####
# 加载R包，没有安装请先安装  install.packages("包名") 
  library(car)
  library(scatterplot3d)
# 读取PCA数据文件
  expression_matrix = read.table("surv大于30天的fpkm01A数据.txt",
                  header = T,    # 指定第一行是列名
                  row.names = 1  )# 指定第一列是行名
# 读取样本分组数据文件
  group = read.csv("生存大于30分组结果.csv",
                       header = T,
                       row.names = 1 )
#清洗数据去除表达为0占比百分之五十的基因
  expression_matrix = expression_matrix[apply(expression_matrix, 1, function(x) sum(x > 0) > 0.5*ncol(expression_matrix)), ]
  nrow(expression_matrix)
# 标准化数据
log_data <- log2(expression_matrix + 1)
# 对患者分组进行t转换
log_data=t(expression_matrix)
# 执行PCA
pca_result <- prcomp(log_data)
# 提取主成分得分
pca_scores <- as.data.frame(pca_result$x)
#判断是否相同  
identical(rownames(pca_scores),rownames(group))
# 添加样本分组信息
pca_scores$group <- group$Group
#设置颜色
colors <- c("#EF6262","#1D5B79")
colors <- colors[as.numeric(as.factor(group[,1]))]

scatterplot3d(x =pca_scores$PC1, y = pca_scores$PC3, z = pca_scores$PC2,
              main = "PCA", xlab='PCA1', ylab="PCA3", zlab='PC2',
              pch = 19, color = colors, cex.symbols = 1)
# 设置图例
legend("bottomright",
       legend = unique(group[,1]),
       col =  c("#EF6262","#1D5B79"),
       pch = 16,
       inset = -0.1,
       xpd = TRUE,
       horiz = TRUE)
# 绘制立体PCA图
scatter3d(pca_scores$PC1, pca_scores$PC2, pca_scores$PC3, col = pca_scores$group, pch = 16, 
          xlab = "PC1", ylab = "PC2", zlab = "PC3", main = "3D PCA Plot")







####GO/KEGG(成功)####
#下载包
  install.packages(c("clusterProfiler", "org.Hs.eg.db"))
#加载所需的包
  library(clusterProfiler)
  library(org.Hs.eg.db)
#读取数据
  DEG <- read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv", header = FALSE)
#将基因名称转换为数字形式
  genelist <- bitr(DEG$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#将转换后的基因ID加入DEG数据框
  DEG$ENTREZID <- genelist$ENTREZID
#GO分析
  ego <- enrichGO(gene = DEG$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "all",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
#KEGG分析
  kk <- enrichKEGG(gene         = DEG$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

#保存GO分析结果
  save(ego, file = "GO_cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.Rdata")
  write.csv(ego@result, file = "GO_cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv")
#保存KEGG分析结果 
  save(kk, file = "KEGG_cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.Rdata")
  write.csv(kk@result, file = "KEGG_cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv")
#可视化
  #柱状图
    barplot(ego, showCategory = 20, color = "pvalue")
  #气泡图
    dotplot(ego, showCategory = 20)
#可视化
  #柱状图
    barplot(kk, showCategory = 20, color = "pvalue")
  #气泡图
    dotplot(kk, showCategory = 20)
#分类展示
  #柱状图
    barplot(ego, drop = TRUE, showCategory = 10, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., scale = 'free')
  #气泡图
    dotplot(ego, showCategory = 10, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., scale = 'free')
    #KEGG结果注释
    #ONTOLOGY：区分是BP，MF还是CC
    #ID：具体的GO条目的ID号
    #Description：GO条目的描述
    #GeneRatio：这里是一个分数，分子是富集到这个GO条目上的gene的数目，
    #           分母是所有输入的做富集分析的gene的数目，可以是差异表达分析得到的gene
    #BgRatio：Background Ratio. 这里也是一个分数，分母是人的所有编码蛋白的
    #         基因中有GO注释的gene的数目，这里是19623个，分子是这19623个
    #         gene中注释到这个GO条目上面的gene的数目
    #pvalue：富集的p值
    #p.adjust：校正之后的p值
    #qvalue：q值
    #geneID：输入的做富集分析的gene中富集到这个GO条目上面的具体的gene名字
    #Count：输入的做富集分析的gene中富集到这个GO条目上面的gene的数目
   
    
    
     
####联合logFC进行GO/KEGG富集分析####
#加载所需的包
  library(clusterProfiler)
  library(org.Hs.eg.db)
#读取差异表达基因列表和logFC信息
#读取全部基因coxoutput
  gene <- read.csv("cox results_LIHC abs1_Anoikis0.4 p0.05 hr筛后 基因.csv")
  exp <- read.csv("LIHC_abs1.csv")
#改行名并取交集
  coxgene <- intersect(gene$Gene,exp$X)
  rownames(exp) <- exp$X
#得到想要基因的cox
  exp <- exp[coxgene,]
#提取logFC
  logFC <- exp$log2FoldChange
  exp$logFC <- exp$log2FoldChange
#将基因名称转换为数字形式
  genelist <- bitr(exp$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#将转换后的基因ID加入DEG数据框
  exp$ENTREZID <- genelist$ENTREZID
#进行GO分析
  ego <- gseGO(geneList = exp$ENTREZID,  # 假设基因ID信息在DEG数据框中的列名为"ENTREZID"
               OrgDb = org.Hs.eg.db,
               ont = "all",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE,
               geneSetID = NULL,
               universe = NULL,
               minGSSize = 1,
               maxGSSize = 500,
               by = "logFC",
               logFC = logFC)
    
#进行KEGG分析
  kegg <- gseKEGG(geneList = exp$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE,
                  geneSetID = NULL,
                  universe = NULL,
                  minGSSize = 1,
                  maxGSSize = 500,
                  by = "logFC",
                  logFC = logFC)
    
    # 保存GO和KEGG分析结果
    save(ego, kegg, file = "GO_KEGG_results.Rdata")
    
    # 保存GO和KEGG分析整体结果
    write.csv(ego@result, file = "GO_results.csv")
    write.csv(kegg@result, file = "KEGG_results.csv")
    
    
    
    
    
    
    
    
    
    
    
    
####带分组信息的热图####
# 读取基因表达量数据
    data <- read.csv("ConsensusClusterPlus聚类原始数据 surv大于30天.csv", row.names = 1)
# 读取样本的注释信息
    group <- read.csv("LIHC_phenotype_筛选后大于30天含分组数据用于热图绘制.csv", row.names = 1)
# 加载所需的库
    library(pheatmap)
# 设置颜色渐变
    color <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
# 绘制热图
    pheatmap(log2(data + 1),
             color = color,          # 图例颜色
             cluster_col = FALSE,    # 不显示样本聚类
             show_colnames = FALSE,  # 不显示样本名称
             annotation_col = group, # 添加列的注释
             fontsize_row = 10,      # 调整基因名称大小
             scale = "row"           # 对行标准化,这一步很重要
    )
    
    # 读取基因表达量数据
    data <- read.csv("ConsensusClusterPlus聚类原始数据 surv大于30天.csv", row.names = 1)
    
    # 读取样本的注释信息
    group <- read.csv("LIHC_phenotype_筛选后大于30天含分组数据用于热图绘制.csv", row.names = 1)
    
    # 设置颜色渐变
    color <- colorRampPalette(c("#EF6262", "#A2CDB0", "#FFD89C"))(80)
  #依据某一列对整体进行排序
    group <- group[order(group$group), ]  

    # 对分组数据进行排序并获取排序后的索引
    sorted_index <- rownames(group)
    
    # 根据排序后的索引重新排序基因表达量数据
    sorted_data <- data[,sorted_index ]
    
    # 绘制热图
    pheatmap(log2(sorted_data + 1),
             color = color,
             cluster_col = FALSE,
             show_colnames = FALSE,
             annotation_col = group,
             fontsize_row = 10,
             scale = "row"
    )
    
    
    
####宽数据转为长数据#####
#加载包  
  library(tidyr)
#读取宽数据
  data_kuan <- read.csv("gsva结果仅含分组.csv", check.names=FALSE)  
#预览宽数据
  data_kuan
#进行长数据转换
  data_chang <- gather(data_kuan,                 #data：待转换的宽数据
         key = "hallmark pathway",               #key：新标识变量的名称，相当于上文的“科目”；
         value = "GSVA",
         -"Group",
         na.rm = FALSE, convert = FALSE, factor_key = FALSE)
#对处理后数据进行保存
  write.csv(data_chang, file = "data_chang.csv")
  
  
  
  
#####分组箱线图####
#加载包
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(patchwork)
#创建箱线图
  p <- ggplot(data = data_chang, aes(x = `hallmark pathway`, y = GSVA)) +
    geom_boxplot(aes(fill = Group), outlier.shape = 1, outlier.colour = "gray", outlier.alpha = 0.2) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), legend.position = "bottom") +
    scale_fill_manual(values = c("C1" = "#EF6262", "C2" = "#1D5B79"), name = NULL) +
    labs(x = NULL, y = "GSVA") +
    scale_y_continuous(breaks = c(-2, -1, 0), labels = c("0.01", "0.10", "1.00"))
#提取箱线图的数据
  errorbar.df <- ggplot_build(p)$data[[1]] %>% select(x, ymin, ymax)
#添加标签和标题
  p <- p +
    labs(x = "Hallmark Pathway", y = "GSVA", title = "Comparison of GSVA across Hallmark Pathways")
#调整图例位置
  p <- p +
    theme(legend.position = "right")
#在图的底部添加误差线
  p.bottom <- p +
    geom_segment(data = errorbar.df, aes(x = x - 0.15, xend = x + 0.15, y = ymin, yend = ymin)) +
    geom_segment(data = errorbar.df, aes(x = x - 0.15, xend = x + 0.15, y = ymax, yend = ymax))
#提取unique的`hallmark pathway`信息
  group.info <- data_chang %>% pull(`hallmark pathway`) %>% unique()
#计算每个`hallmark pathway`的p-value
  pvalue.df <- tibble(x = character(), pvalue = numeric())
  for (info in group.info) {
    tmp.df <- data_chang %>% filter(`hallmark pathway` == info)
    a <- wilcox.test(`GSVA` ~ `Group`, data = tmp.df)
    pvalue.df <- pvalue.df %>% add_row(x = info, pvalue = a$p.value)
  }
#替换p-value为最小非零值
  min_p <- min(pvalue.df %>% filter(pvalue != 0) %>% pull(pvalue))
  pvalue.df <- pvalue.df %>% mutate(new_p = ifelse(pvalue == 0, min_p, pvalue))
#创建散点图
  p.top <- ggplot(data = pvalue.df, aes(x = x, y = 1)) +
    geom_point(size = 5, shape = 15, aes(color = -log10(new_p))) +
    scale_color_distiller(palette = "Greys", direction = 1, name = "-log10(Pval)") +
    theme_void() +
    theme(legend.position = "top", legend.key.size = unit(5, 'mm'), legend.text = element_text(size = 10)) +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 20))
#组合图形
  p_combined <- p.top + p.bottom + plot_layout(ncol = 1, heights = c(1, 5))
#展示图片
  p_combined
  
  
  
  
  
  
####绘制细胞系基因表达雷达图####
#加载包
  library(ggplot2)
  library(tidyverse)
#读取数据
  dt <- read.csv("Cellline.csv")
#添加Id
  dt$id <- 1:nrow(dt)
#计算标签角度
  dt$angle <- 90 - 360 * (dt$id - 0.5) / nrow(dt)
  dt$hjust <- ifelse(dt$angle < -90, 1, 0)
  dt$angle <- ifelse(dt$angle < -90, dt$angle + 180, dt$angle)
#自定义颜色
  custom_colors <- c("#EECAD5", "#F6EACB", "#D1E9F6")
#绘制雷达图
  ggplot(data = dt, aes(x = as.factor(id), y = -Gene, fill = Gene)) +
    geom_bar(stat = "identity") +
    coord_polar(start = 0) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
    geom_text(aes(y = 0.5, label = Cell.Line, hjust = hjust), color = "black", size = 3, angle = dt$angle) +
    geom_text(aes(y = -Gene / 1.5, label = Gene, hjust = hjust), color = 'black', size = 2.6, angle = dt$angle) +
    scale_fill_gradientn(colors = custom_colors, name = "KIF18A") +
    ylim(-35, 1)  # 确保数据在此范围内

  

  
  
  

  
