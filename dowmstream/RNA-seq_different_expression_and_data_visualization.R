library("edgeR")
library(org.Hs.eg.db)
library(tinyarray)
library(dplyr)
library(ggpubr)


species="hsa" 
org.db="org.Hs.eg.db"
group1=c(rep("Lx_control",3))
group2=c(rep("Lx_ip6",3))
group3=c(rep("Lx_ins",3))
group4=c(rep("Lx_ip6_ins",3))

group=c(group1,group2,group3,group4)




FPKM <- read.table("FPKM.txt",sep="\t",header=T,check.names=F)
FPKM <- FPKM[,-c(1:2)]

#pca
pdf("different_expression_gene/pca.pdf")
draw_pca(FPKM, group_list=group)
dev.off()







#loading count table and merge data with gene name
rt=read.table("count.txt",sep="\t",header=T,check.names=F)#改为自己的文件
rownames(rt)=rt[,1]
rt=rt[,-1]

rt <-aggregate(rt[,-1],by=list(rt[,1]),FUN=mean)
rownames(rt)=rt[,1]
rt=rt[,-1]

#prepair data for edgeR
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

#Filter low expression gene

data=avereps(data)
data=data[rowMeans(data)>1,] #keep genes with average count great than 1


edger_diff <- function(data,name,group,foldChange=1,padj=0.05,groupA="groupA",groupB="groupB"){
  #set model for edgeR
  design <- model.matrix(~group)
  #数据归一化
  y <- DGEList(counts=data,group=group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  
  #求差异基因，输出基因表达结果及差异基因列表
  et <- exactTest(y,pair = c(groupA,groupB))
  topTags(et)
  ordered_tags <- topTags(et, n=100000)
  allDiff=ordered_tags$table
  allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
  diff=allDiff
  newData=y$pseudo.counts
  dir.create(name)
  dir.create(paste0(name,"/different_expression_gene"))
  write.table(diff,file=paste0(name,"/different_expression_gene/edgerOut.xls"),sep="\t",quote=F)
  diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
  write.table(diffSig, file=paste0(name,"/different_expression_gene/diffSig.xls"),sep="\t",quote=F)
  diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
  write.table(diffUp,file=paste0(name,"/different_expression_gene/up.xls"),sep="\t",quote=F)
  diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
  write.table(diffDown, file=paste0(name,"/different_expression_gene/down.xls"),sep="\t",quote=F)
  normalizeExp=rbind(id=colnames(newData),newData)
  write.table(normalizeExp,file=paste0(name,"/different_expression_gene/normalizeExp.txt"),sep="\t",quote=F,col.names=F)   
  
  diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
  write.table(diffExp,file=paste0(name,"/different_expression_gene/diffmRNAExp.txt"),sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）
}
  
  #火山图
  pdf(file="different_expression_gene/vol.pdf")
  allDiff2 <- allDiff
  allDiff2[,5] <- allDiff2[,4]
  draw_volcano(allDiff2,
               pkg = 2,
               logFC_cutoff = foldChange,
               pvalue_cutoff = 0.05,
               adjust = T
  )
  dev.off()
  
  png(file="different_expression_gene/vol.png")
  allDiff2 <- allDiff
  allDiff2[,5] <- allDiff2[,4]
  draw_volcano(allDiff2,
               pkg = 2,
               logFC_cutoff = foldChange,
               pvalue_cutoff = 0.05,
               adjust = T
  )
  dev.off()
  
  #热图
  df <- newData[rownames(diffSig),]
  col <- colorRampPalette(c("blue", "white", "red"))(256)
  patient_class=RowSideColors =  c(rep("purple", length(group1)), rep("orange", length(group2)))
  #gene_class=rep(c("blue","pink"),each=16)
  png("different_expression_gene/heatmap.pdf")
  heatmap(df[rev(rownames(df)),],
          col = col,
          Colv = NA,
          ColSideColors = patient_class)
  dev.off()
  
  png("different_expression_gene/heatmap.png")
  heatmap(df[rev(rownames(df)),],
          col = col,
          Rowv = NA,
          Colv = NA,
          ColSideColors = patient_class)
  dev.off()
  
  
  #富集分析 
  #构建以Entrez ID为标准的差异基因列表
  dir.create("GO分析")
  library(clusterProfiler)
  EG2Symbol=toTable(org.Hs.egSYMBOL)
  geneLists<-diffSig
  geneLists[,5]<-rownames(diffSig)
  colnames(geneLists)=c("logFC","logCPM","PValue","FDR","symbol")
  results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  results=na.omit(results)
  id=results$gene_id
  fc=results$logFC
  names(fc)=as.character(results[,6])
  #GO分类
  ggo<-groupGO(gene=id,org.db,keyType ="ENTREZID",ont = "BP",level = 3,readable = T)
  groupgo=as.data.frame(ggo)
  write.csv(groupgo,"GO分析/groupgo.csv")
  
  
  #Go富集及制图
  #GO富集
  ALL<-enrichGO(OrgDb=org.db, gene = id,
                ont = "ALL",pvalueCutoff= 0.05,readable=T) 
  BP <-enrichGO(OrgDb=org.db, gene = id,
                ont = "BP", pvalueCutoff= 0.05,readable=T) %>% simplify()
  MF <-enrichGO(OrgDb=org.db, gene = id,
                ont = "MF", pvalueCutoff= 0.05,readable=T) %>% simplify()
  CC <-enrichGO(OrgDb=org.db, gene = id,
                ont = "CC", pvalueCutoff= 0.05,readable=T) %>% simplify()
  #条形图
  
  BP_barplot <- barplot(BP, showCategory=5,title="BPGO")
  MF_barplot <- barplot(MF, showCategory=5,title="MFGO")
  CC_barplot <- barplot(CC, showCategory=5,title="CCGO")
  
  pdf(file="GO分析/GO_barplot.pdf", bg="transparent",width = 7,height = 13)
  ggarrange(BP_barplot, MF_barplot,CC_barplot,ncol = 1,align = "v")
  dev.off()
  
  png(file="GO分析/GO_barplot.png", bg="transparent",width = 480,height = 900)
  ggarrange(BP_barplot, MF_barplot,CC_barplot,ncol = 1,align = "v")
  dev.off()
  
  
  dev.off()
  #气球图
  pdf(file="ALLGO_dotplot.pdf", bg="transparent",width = 12,height = 10)
  dotplot(ALL,showCategory=30,title="ALLGO")
  dev.off()
  pdf(file="BPGO_dotplot.pdf", bg="transparent",width = 12,height = 10)
  dotplot(BP,showCategory=30,title="BPGO")
  dev.off()
  pdf(file="MFGO_dotplot.pdf", bg="transparent",width = 12,height = 10)
  dotplot(MF,showCategory=30,title="MFGO")
  dev.off()
  pdf(file="CCGO_dotplot.pdf", bg="transparent",width = 12,height = 10)
  dotplot(CC,showCategory=30,title="CCGO")
  dev.off()
  #网状图
  pdf(file="ALLGO_cneplot.pdf", bg="transparent")
  cnetplot(ALL,categorySize="pvalue",title="ALLGO",foldChange=fc)
  dev.off()
  pdf(file="BPGO_cneplot.pdf", bg="transparent")
  cnetplot(BP,categorySize="pvalue",title="BPGO",foldChange=fc)
  dev.off()
  pdf(file="MFGO_cneplot.pdf", bg="transparent")
  cnetplot(MF,categorySize="pvalue",title="MFGO",foldChange=fc)
  dev.off()
  pdf(file="CCGO_cneplot.pdf", bg="transparent")
  cnetplot(CC,categorySize="pvalue",title="CCGO",foldChange=fc)
  dev.off()
  #网络关系图
  pdf(file="BPGO_networkplot.pdf", bg="transparent")
  plotGOgraph(BP)
  dev.off()
  pdf(file="MFGO_networkplot.pdf", bg="transparent")
  plotGOgraph(MF)
  dev.off()
  pdf(file="CCGO_networkplot.pdf", bg="transparent")
  plotGOgraph(CC)
  dev.off()
  #环状网络图
  pdf(file="ALLGO_circl_cneplot.pdf", bg="transparent",width=13,height=10)
  cnetplot(ALL,foldChange=fc, circular = TRUE, colorEdge = TRUE)
  dev.off()
  pdf(file="BPGO_circl_cneplot.pdf", bg="transparent",width=13,height=10)
  cnetplot(BP,foldChange=fc, circular = TRUE, colorEdge = TRUE)
  dev.off()
  pdf(file="CCGO_circl_cneplot.pdf", bg="transparent",width=13,height=10)
  cnetplot(CC,foldChange=fc, circular = TRUE, colorEdge = TRUE)
  dev.off()
  pdf(file="MFGO_circl_cneplot.pdf", bg="transparent",width=13,height=10)
  cnetplot(MF,foldChange=fc, circular = TRUE, colorEdge = TRUE)
  dev.off()
  #输出结果
  ALLGO<-as.data.frame(ALL)
  write.csv(ALLGO,"GO分析/ALLgo.csv")
  #go-gsea
  
  
  #GOplot绘图
  #整理表格
  library(GOplot)
  library(tidyr)
  library(dplyr)
  enrichgo2=ALLGO[,c(1,2,3,7,9,10)]
  colnames(enrichgo2)=c("category","ID","term","adj_pval","genes","count")
  enrichgo2[,5]=chartr("/",",",enrichgo2[,5])
  difflist=diffSig
  difflist[,5]=row.names(difflist)
  colnames(difflist)=c("logFC","logCPM","adj_pval","FDR","ID")
  GOplotdata=circle_dat(enrichgo2,difflist)
  #制图
  png(filename = "bubble.png",width=2000,height=1000)
  reduced_circ <- reduce_overlap(GOplotdata, overlap = 0.75)#去除重叠bubble
  GOBubble(reduced_circ, title = 'Bubble plot', colour = c('orange','darkred','gold'),
           display='multiple',labels = 3,ID=F) #labels可以改变输出的结果的多少
  dev.off()
  png(filename = "circ.png",width=1000,height=1000)
  GOCircle(GOplotdata)
  dev.off()
  #更多绘图方法参见R语言GOplot
  
  
  #kegg分析
  dir.create("KEGG")
  kk<-enrichKEGG(id,organism=species)
  enrichkegg=as.data.frame(kk)
  kkoutput=as.data.frame(kk)
  write.csv(kkoutput,"kegg.csv")
  #cnetplot(kk,showCategory = 5,categorySize=??)
  
  #制图
  #条形图
  png(file="kegg_barplot.png", bg="transparent")
  barplot(kk, showCategory=20,title="kegg")
  dev.off()
  #气球图
  png(file="kegg_dotplot.png", bg="transparent")
  dotplot(kk,showCategory=30,title="kegg")
  dev.off()
  
  #pathway view
  #整理数据
  library(pathview)
  keggdata=diff[which(diff[,3]<0.05),]
  keggdata[,5]=row.names(keggdata)
  keggdata=keggdata[,c(1,5)]
  colnames(keggdata)=c("logFC","symbol")
  keggdata=merge(keggdata,EG2Symbol,by='symbol',all.x=T)
  keggdata=na.omit(keggdata)
  row.names(keggdata)=keggdata[,3]
  a=as.matrix(keggdata[,2])
  rownames(a)=rownames(keggdata)
  
  #制图
  for (i in 1:length(enrichkegg)) {
    pathway.id = enrichkegg[i,1]
    term=enrichkegg[i,2]
    p <- pathview(gene.data =a, pathway.id = pathway.id, 
                  species = species, out.suffix = term)
    p<- pathview(gene.data =a, pathway.id = pathway.id,
                 species = species, out.suffix = term, kegg.native = F)
  }
  
  
  #GSEA
  #数据整理
  gseadata=diff
  gseadata[,5]=row.names(gseadata)
  gseadata=gseadata[,c(1,5)]
  colnames(gseadata)=c("logFC","symbol")
  gseadata=merge(gseadata,EG2Symbol,by='symbol',all.x=T)
  gseadata=na.omit(gseadata)
  gseadata <- gseadata[order(gseadata[,2],decreasing = T),]
  gseadata=gseadata[,c(2,3)]
  write.csv(gseadata,"gseadata.csv")
  gene_fc=gseadata$logFC
  names(gene_fc)=as.character(gseadata$gene_id)
  #分析及制图(制图的参数还需要改)
  gogsea <- gseGO(gene_fc,OrgDb =org.db,ont= "BP")#详细参数见gseGO
  kegggsea<-gseKEGG(gene_fc,organism=species)#详细参数见gseKEGG
  


