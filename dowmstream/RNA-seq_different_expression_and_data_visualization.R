suppressMessages(library(getopt))

options<-matrix(c(
  "help", "h", "0", "logical", "help",
  "count_file","d","1","character","data",
  "fpkm","k","1","character","data",
  "species", "s", "1", "character", "species",
  "padj", "p", "1", "double", "foldChange",
  "foldChange", "f", "1", "double", "foldChange",
  "org_db", "o", "1", "character", "org_db",
  "group","g","1","character","group",
  "groupA","a","1","character","groupA",
  "groupB","b","1","character","groupB"
  
), ncol=5, byrow=T)

opt<-getopt(options)

if(!is.null(opt$help)){
  cat(getopt(options, usage=T))
  q()
}

suppressMessages(library(edgeR))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(org.Rn.eg.db))
suppressMessages(library(tinyarray))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
suppressMessages(library(GOplot))
suppressMessages(library(tidyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pathview))
suppressMessages(library(enrichplot))
suppressMessages(library(ReactomePA))
suppressMessages(library(biomaRt))


padj=0.05
foldChange=1
count_file="count.txt"
fpkm="FPKM.txt"

if(!is.null(opt$count_file)){count_file=opt$count_file}
if(!is.null(opt$fpkm)){fpkm=opt$fpkm}
if(!is.null(opt$foldChange)){foldChange=opt$foldChange}
if(!is.null(opt$padj)){padj=opt$padj}

species <- opt$species
org_db <- opt$org_db
group_all <- c(read.table(opt$group)[,1])


print('group:')
group_all

  groupA <- opt$groupA
  groupB <- opt$groupB
  group_compair_select <- which(group_all==groupA|group_all==groupB)
  group_compair <- group_all[group_compair_select]
  
  #loading count table and merge data with gene name
  rt <- read.table(count_file,sep="\t",header=T,check.names=F)#改为自己的文件
  rownames(rt) <- rt[,1]
  rt <- rt[,-1]
  
  rt <-aggregate(rt[,-1],by=list(rt[,1]),FUN=mean)
  rownames(rt) <- rt[,1]
  rt <- rt[,-1]
  
  #prepair data for edgeR
  rt <- as.matrix(rt)
  exp <- rt
  dimnames <- list(rownames(exp),colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  
  #Filter low expression gene
  
  data <- avereps(data)
  data <- data[rowMeans(data)>1,] #keep genes with average count great than 1
  
  data <- data[,group_compair_select]
  
  name = paste0(groupA,"_",groupB)
  
  group=group_compair
  
  group1 <- group[which(group==groupA)]
  group2 <- group[which(group==groupB)]
  #set model for edgeR
  design <- model.matrix(~group)
  #data normalization and model
  y <- DGEList(counts=data,group=group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  
  #DEG
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
  write.table(diffExp,file=paste0(name,"/different_expression_gene/diffmRNAExp.txt"),sep="\t",quote=F,col.names=F) 
  
  
  FPKM <- read.table(fpkm,sep="\t",header=T,check.names=F)
  FPKM <- FPKM[,-c(1:2)]
  FPKM <- na.omit(FPKM)
  
  #pca
  pdf(paste0(name,"/different_expression_gene/pca.pdf"))
  draw_pca(FPKM, group_list=group_all)
  dev.off()
  
  #heatmap
  pdf(paste0(name,"/different_expression_gene/heatmap.pdf"))
  df <-  newData[rownames(diffSig),]
  col <- colorRampPalette(c("blue", "white", "red"))(256)
  patient_class=RowSideColors <- c(rep("purple", length(group1)), rep("orange", length(group2)))
  heatmap(df[rev(rownames(df)),],
          col = col,
          Colv = NA,
          ColSideColors = patient_class)
  dev.off()
  
  #Volcano Plot
  pdf(paste0(name,"/different_expression_gene/vol.pdf"))
  allDiff2 <- allDiff
  allDiff2[,5] <- allDiff2[,4]
  draw_volcano(allDiff2,
               pkg = 2,
               logFC_cutoff = foldChange,
               pvalue_cutoff = padj,
               adjust = T
  )
  dev.off()
  
  
  ##enrichment analysis 
  
  #Construct a list of differential genes using Entrez ID
  dir.create(paste0(name,"/GO_analysis"))
  
  if(species=="hsa"){
    EG2Symbol=toTable(org.Hs.egSYMBOL)
  }else if(species=="mmu"){
    EG2Symbol=toTable(org.Mm.egSYMBOL)
  }else if(species=="rno"){
    EG2Symbol=toTable(org.Rn.egSYMBOL)
  }else{print("only support species of hsa mmu and rno")}
  
  geneLists<-diffSig
  geneLists[,5]<-rownames(diffSig)
  colnames(geneLists)=c("logFC","logCPM","PValue","FDR","symbol")
  results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  results=na.omit(results)
  id=results$gene_id
  fc=results$logFC
  names(fc)=as.character(results[,6])
  #GO分类
  ggo<-groupGO(gene=id,org_db,keyType ="ENTREZID",ont = "BP",level = 3,readable = T)
  groupgo=as.data.frame(ggo)
  write.csv(groupgo,paste0(name,"/GO_analysis/groupgo.csv"))
  
  
  ALL<-enrichGO(OrgDb=org_db, gene = id,
                ont = "ALL",pvalueCutoff= 0.05,readable=T) 
  BP <-enrichGO(OrgDb=org_db, gene = id,
                ont = "BP", pvalueCutoff= 0.05,readable=T) %>% simplify()
  MF <-enrichGO(OrgDb=org_db, gene = id,
                ont = "MF", pvalueCutoff= 0.05,readable=T) %>% simplify()
  CC <-enrichGO(OrgDb=org_db, gene = id,
                ont = "CC", pvalueCutoff= 0.05,readable=T) %>% simplify()  
  
  
  
  #go analysis plot
  
  #barplot
  BP_barplot <- barplot(BP, showCategory=5,title="BPGO",orderBy ="p.adjust")
  MF_barplot <- barplot(MF, showCategory=5,title="MFGO",orderBy ="p.adjust")
  CC_barplot <- barplot(CC, showCategory=5,title="CCGO",orderBy ="p.adjust")
  
  pdf(file=paste0(name,"/GO_analysis/BP_barplot.pdf"), bg="transparent",width = 7,height = 7)
  BP_barplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/MF_barplot.pdf"), bg="transparent",width = 7,height = 7)
  MF_barplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/CC_barplot.pdf"), bg="transparent",width = 7,height = 7)
  CC_barplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/GO_barplot.pdf"), bg="transparent",width = 7,height = 13)
  ggarrange(BP_barplot, MF_barplot,CC_barplot,ncol = 1,align = "v")
  dev.off()
  
  #bubble plot
  BP_dotplot=dotplot(BP,showCategory=30,title="BPGO")
  MF_dotplot=dotplot(MF,showCategory=30,title="MFGO")
  CC_dotplot=dotplot(CC,showCategory=30,title="ccGO")
  
  
  
  pdf(file=paste0(name,"/GO_analysis/BP_dotplot.pdf"), bg="transparent",width = 12,height = 10)
  BP_dotplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/MF_dotplot.pdf"), bg="transparent",width = 12,height = 10)
  MF_dotplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/CC_dotplot.pdf"), bg="transparent",width = 12,height = 10)
  CC_dotplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/GO_dotplot.pdf"), bg="transparent",width = 7,height = 13)
  ggarrange(BP_dotplot, MF_dotplot,CC_dotplot,ncol = 1,align = "v")
  dev.off()
  
  
  #network
  
  BP_cnetplot=cnetplot(BP,categorySize="pvalue",title="BPGO")
  
  MF_cnetplot=cnetplot(MF,categorySize="pvalue",title="BPGO")
  
  CC_cnetplot=cnetplot(CC,categorySize="pvalue",title="BPGO")
  
  
  pdf(file=paste0(name,"/GO_analysis/BP_cnetplot.pdf"), bg="transparent")
  BP_cnetplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/MF_cnetplot.pdf"), bg="transparent")
  MF_cnetplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/CC_cnetplot.pdf"), bg="transparent")
  CC_cnetplot
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/GO_cnetplot.pdf"), bg="transparent",width = 7,height = 13)
  ggarrange(BP_cnetplot, MF_cnetplot,CC_cnetplot,ncol = 1,align = "v")
  dev.off()
  
  
  
  #networkplot
  
  
  pdf(file=paste0(name,"/GO_analysis/BOGO_networkplot.pdf"),  bg="transparent")
  plotGOgraph(BP)
  dev.off()
  pdf(file=paste0(name,"/GO_analysis/MFGO_networkplot.pdf"),  bg="transparent")
  plotGOgraph(MF)
  dev.off()
  pdf(file=paste0(name,"/GO_analysis/CCGO_networkplot.pdf"),  bg="transparent")
  plotGOgraph(CC)
  dev.off()
  
  
  
  
  #circle network
  pdf(file=paste0(name,"/GO_analysis/ALLGO_circl_cneplot.pdf"),  bg="transparent",width=13,height=10)
  cnetplot(ALL,circular = TRUE, colorEdge = TRUE)
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/BPGO_circl_cneplot.pdf"),  bg="transparent",width=13,height=10)
  cnetplot(BP,circular = TRUE, colorEdge = TRUE)
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/MFGO_circl_cneplot.pdf"),  bg="transparent",width=13,height=10)
  cnetplot(MF, circular = TRUE, colorEdge = TRUE)
  dev.off()
  
  pdf(file=paste0(name,"/GO_analysis/CCGO_circl_cneplot.pdf"),  bg="transparent",width=13,height=10)
  cnetplot(CC,circular = TRUE, colorEdge = TRUE)
  dev.off()
  
  
  ALLGO<-as.data.frame(ALL)
  write.csv(ALLGO,paste0(name,"/GO_analysis/ALLgo.csv"))
  
  
  
  #Using GOplot
  
  enrichgo2=ALLGO[,c(1,2,3,7,9,10)]
  colnames(enrichgo2)=c("category","ID","term","adj_pval","genes","count")
  enrichgo2[,5]=chartr("/",",",enrichgo2[,5])
  difflist=diffSig
  difflist[,5]=row.names(difflist)
  colnames(difflist)=c("logFC","logCPM","adj_pval","FDR","ID")
  GOplotdata=circle_dat(enrichgo2,difflist)
  
  pdf(paste0(name,"/GO_analysis/bubble.pdf"),width=2000,height=1000)
  reduced_circ <- reduce_overlap(GOplotdata, overlap = 0.75)
  GOBubble(reduced_circ, title = 'Bubble plot', colour = c('orange','darkred','gold'),
           display='multiple',labels = 3,ID=F) 
  dev.off()
  
  pdf(paste0(name,"/GO_analysis/circ.pdf"),width=1000,height=1000)
  GOCircle(GOplotdata)
  dev.off()
  
  
  
  ##kegg 暂时弃用，服务器不支持
  if(FALSE){
  R.utils::setOption("clusterProfiler.download.method",'auto')
  dir.create(paste0(name,"/KEGG"))
  kk<-enrichKEGG(id,organism=species)
  enrichkegg=as.data.frame(kk)
  kkoutput=as.data.frame(kk)
  write.csv(kkoutput,paste0(name,"/KEGG/kegg.csv"))
  
  #barplot
  pdf(file=paste0(name,"/KEGG/kegg_barplot.pdf"), bg="transparent")
  barplot(kk, showCategory=20,title="kegg")
  dev.off()
  #bubble
  pdf(file=paste0(name,"/KEGG/kegg_dotplot.pdf"), bg="transparent")
  dotplot(kk,showCategory=30,title="kegg")
  dev.off()
  
  #pathway view
  
  keggdata=diff[which(diff[,3]<0.05),]
  keggdata[,5]=row.names(keggdata)
  keggdata=keggdata[,c(1,5)]
  colnames(keggdata)=c("logFC","symbol")
  keggdata=merge(keggdata,EG2Symbol,by='symbol',all.x=T)
  keggdata=na.omit(keggdata)
  row.names(keggdata)=keggdata[,3]
  a=as.matrix(keggdata[,2])
  rownames(a)=rownames(keggdata)
  
  setwd(name)
  for (i in 1:length(enrichkegg)) {
    pathway.id = enrichkegg[i,1]
    term=enrichkegg[i,2]
    p <- pathview(gene.data =a, pathway.id = pathway.id, 
                  species = species, out.suffix = term)
    p<- pathview(gene.data =a, pathway.id = pathway.id,
                 species = species, out.suffix = term, kegg.native = F)
  }
  setwd("../")
  }
  
  #GSEA
  dir.create(paste0(name,"/GSEA"))
  aaa <- results[,c(6,2)]
  aaa$logFC<-sort(aaa$logFC,decreasing = T)
  geneList = aaa[,2]
  names(geneList) = as.character(aaa[,1])
  
  Go_gseresult <- gseGO(geneList, org_db, keyType = "ENTREZID", ont="all", 
                        nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  
  if(species=="hsa"){Go_Reactomeresult <- gsePathway(geneList, 
                                                     nPerm = 1000, minGSSize = 10, 
                                                     maxGSSize = 1000, pvalueCutoff=1)}
  
  write.table (Go_gseresult, file =paste0(name,"/GSEA/Go_gseresult.csv"), sep =",", row.names =TRUE)
  #write.table (KEGG_gseresult, file =paste0(name,"/GSEA/KEGG_gseresult.csv"), sep =",", row.names =TRUE)
  
  if(species=="hsa"){
    
    write.table (Go_Reactomeresult, file =paste0(name,"/GSEA/Go_Reactomeresult.csv"), sep =",", row.names =TRUE)
  }
  
  #go
  pdf(paste0(name,"/GSEA/GO_ridgeplot.pdf"))
  ridgeplot(Go_gseresult,10)
  dev.off()
  
  pdf(paste0(name,"/GSEA/GO_gseaplot.pdf"))
  gseaplot(Go_gseresult,1,pvalue_table = TRUE)
  dev.off()
  
  pdf(paste0(name,"/GSEA/GO_gseaplot2.pdf"))
  gseaplot2(Go_gseresult,212,pvalue_table = TRUE)
  dev.off()
  
  pdf(paste0(name,"/GSEA/GO_gseaplot_multi.pdf"))
  gseaplot2(Go_gseresult, 1:4, pvalue_table = TRUE)
  dev.off()
  
  #kegg 暂时弃用，服务器不支持
  if(FALSE){
  pdf(paste0(name,"/GSEA/kegg_ridgeplot.pdf"))
  ridgeplot(KEGG_gseresult,10)
  dev.off()
  
  pdf(paste0(name,"/GSEA/kegg_gseaplot.pdf"))
  gseaplot(KEGG_gseresult,1,pvalue_table = TRUE)
  dev.off()
  
  pdf(paste0(name,"/GSEA/kegg_gseaplot2.pdf"))
  gseaplot2(KEGG_gseresult,212,pvalue_table = TRUE)
  dev.off()
  
  pdf(paste0(name,"/GSEA/kegg_gseaplot_multi.pdf"))
  gseaplot2(KEGG_gseresult, 1:4, pvalue_table = TRUE)
  dev.off()
  }




  
  
  


