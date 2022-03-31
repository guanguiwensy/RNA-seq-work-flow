list=dir(pattern="sort.bam.count.txt$")
data1 <- read.table(list[1],header=T,sep="\t")
data1 <- data1[,c(1,7)]
for (i in c(1:(length(list)-1))) {
  data2 <- read.table(list[i+1],header=T,sep="\t")
  data2 <- data2[,c(1,7)]
  data1 <- cbind(data1,data2[,2])
  i=i+1
}

colnames(data1) <- c("gene",list)

rt=data1
gtf=read.table("human.csv",sep=",",stringsAsFactors = F,header=T)
gtf=gtf[,c(1,4)]

rt=merge(gtf,rt,by.x="Gene.stable.ID",by.y="gene",all.y=T)
rt=rt[-c(1:5),]
rt=rt[,-1]
write.table(rt,"symbol.txt",sep="\t",row.names = F)
