directory= "RNA_seq/"
setwd(directory)

sampleTable = read.table("samples.txt",h=T,sep=",")
colnames(sampleTable)=c("sampleName","Treatment","Type","batch","condition","gender")
sampleTable = sampleTable[!sampleTable$Treatment == "Ctrl",]

counts = read.table("counts.csv",h=T,row.names=1,sep=",",check.names = F)
head(counts)
counts = as.matrix(counts[,colnames(counts) %in% sampleTable$sampleName])

condition = "Ctrl"
ctr = "MDS"
cols = c("MDS" = "red","Ctrl" = "forestgreen")

qcutoff = 0.05
logfccut = 0.58

ddsHTSeq<-DESeqDataSetFromMatrix(colData= sampleTable, countData=counts,
                                       design=~ batch + condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds, contrast=c("condition",condition,ctr),cooksCutoff = FALSE)
sig <- as.data.frame(na.omit(res))
sig <- sig[which(sig$padj < qcut & abs(sig$log2FoldChange) > logcut), ]
