#!/usr/bin/env Rscript
library(DESeq2)
library(GOexpress)
library(gplots)
library(optparse)


## reading input files

readFile <- function(x) {
  conn = file(x, open = "r")
  line = readLines(conn)
  line_vec = c()
  for (i in 1:length(line)) {
    line_vec = c(line_vec, line[i])
  }
  close(conn)
  
  return(line_vec)
}


##Reading command-line arguments in R





option_list = list(
  make_option(c("--geneIDs"), type = "character", default = NULL,
              help = "List of geneIDs to be shown on the HMap"),
  make_option(c("--conditionList"), type = "character", default = NULL,
              help = "text file with sample names to be shown on the HMap")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


##Set up the working directory accordingly
#setwd("/Users/ayanmalakar/UKALLURI03162017")

fpkm_counts=read.csv("fpkm_counts_mod.2.csv", header = T, row.names = 1, check.names = F)
fpkm_counts_dummy=fpkm_counts

fpkm_counts=fpkm_counts[,c(-16,-17)]# assuming the last 2 columns have geneDescriptions which we are dropping!

fpkm_counts[]=lapply(fpkm_counts, function(x){log2(x+1)})# Log transformation of the count data

##################
file.list=c("cond_1v3.csv","cond_1v5.csv","cond_2v4.csv","cond_2v6.csv")

intsct.sig.genes.UK=character()
intsct.sig.genes.AM=character()
c=0
for (items in file.list){
  c=c+1
  cond.pair=read.csv(items,header=T,row.names = 1,stringsAsFactors =FALSE)
  temp=as.numeric(cond.pair[[1]])
  cond.pair[[1]]=temp
  temp=as.numeric(cond.pair[[2]])
  cond.pair[[2]]=temp
  cond.pair=cond.pair[complete.cases(cond.pair),]
  sig.genes.UK=cond.pair[cond.pair[[3]]==T,]####Significant genes based on Udaya's definition
  sig.genes.UK.matrix = as.matrix(fpkm_counts[rownames(sig.genes.UK), ])
  ###top20 selected /conditional pair based on average gene expression across ALL CONDITIONS
  top20.sig.genes.UK.matrix = sig.genes.UK.matrix[order(-rowMeans(sig.genes.UK.matrix)),][c(1:20),]
  #prefix = strsplit(items, split = ".", fixed = T)[[1]][1]
  #hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
  #png(paste(prefix,"UK.png",sep = "."),height = 5*300, width = 5*300, res=300, pointsize = 8)
  
  sig.genes.AM=cond.pair[((cond.pair[[1]]>=1 | cond.pair[[1]]<=-1) & cond.pair[[2]]<=0.05)==TRUE,]#our definition
  sig.genes.AM.matrix = as.matrix(fpkm_counts[rownames(sig.genes.AM), ])
  ###top20 selected /conditional pair based on average gene expression across ALL CONDITIONS
  top20.sig.genes.AM.matrix = sig.genes.AM.matrix[order(-rowMeans(sig.genes.AM.matrix)),][c(1:20),]
  #HeatMap per conditional pair AM...
  #prefix = strsplit(items, split = ".", fixed = T)[[1]][1]
  #hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
  #png(paste(prefix,"AM.png",sep = "."),height = 5*300, width = 5*300, res=300, pointsize = 8)
  #fpkm_counts[sig.genes.AM,]
  if (c==1){#i.e the gene list for first condition pair initialise the variables provided with given genes 
    intsct.sig.genes.UK=row.names(sig.genes.UK)
    union.sig.genes.UK=intsct.sig.genes.UK
    intsct.sig.genes.AM=row.names(sig.genes.AM)
    union.sig.genes.AM=intsct.sig.genes.AM
    combined.top20.sig.genes.AM = top20.sig.genes.AM.matrix
    combined.top20.sig.genes.UK = top20.sig.genes.UK.matrix
  } else{
      combined.top20.sig.genes.AM = rbind(combined.top20.sig.genes.AM, top20.sig.genes.AM.matrix)
      combined.top20.sig.genes.UK = rbind(combined.top20.sig.genes.UK, top20.sig.genes.UK.matrix)
      intsct.sig.genes.UK=intersect(intsct.sig.genes.UK,row.names(sig.genes.UK))
      union.sig.genes.UK=union(union.sig.genes.UK,row.names(sig.genes.UK))
      intsct.sig.genes.AM=intersect(intsct.sig.genes.AM, row.names(sig.genes.AM))
      union.sig.genes.AM=union(union.sig.genes.AM,row.names(sig.genes.AM))
  }
}

######Heatmap of Combined top 20s#############
# hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
# png("fpkm_counts.union.ofTop20.AM.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
# geneID = rownames(combined.top20.sig.genes.AM)
# geneNames = fpkm_counts_dummy[geneID,16]
# top80.union.AM.names.row = sub('^NA:', '', x = paste(geneNames, geneID, sep = ":"))
# heatmap.2(combined.top20.sig.genes.AM, labRow = top80.union.AM.names.row, offsetRow = -0.6, 
#           col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
# par(lend=1)
# legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
# dev.off()
# 
# hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
# png("fpkm_counts.union.ofTop20.UK.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
# geneID = rownames(combined.top20.sig.genes.UK)
# geneNames = fpkm_counts_dummy[geneID,16]
# top80.union.UK.names.row = sub('^NA:', '', x = paste(geneNames, geneID, sep = ":"))
# heatmap.2(combined.top20.sig.genes.UK, labRow = top80.union.UK.names.row, offsetRow = -0.6, 
#           col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
# par(lend=1)
# legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
# dev.off()




fpkm_counts.intsct.UK=as.matrix(fpkm_counts[intsct.sig.genes.UK,])
nrow(fpkm_counts.intsct.UK)#Total Number of genes based on itsct.UK

fpkm_counts.intsct.AM=as.matrix(fpkm_counts[intsct.sig.genes.AM,])
nrow(fpkm_counts.intsct.AM)#Total Number of genes based on intsct.AM



hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
###plotting commonly DGE based on "AM"...Top genes with high average expression across provided samples
png("fpkm_counts.intsct.AM.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
#####The list of 70 genes########
#fpkm_counts.intsct.AM <-fpkm_counts.intsct.AM [order(-rowMeans(fpkm_counts.intsct.AM)),][c(1:70),]

print(opt$conditionList)
samples_list = readFile(opt$conditionList)#sample_list must be less < length(colnames)
print(samples_list)

# df = rep(0,length(samples_list))
# dim(df) = c(1,length(samples_list))
# df = as.data.frame(df)
# colnames(df) <- samples_list
# 
# samples_list <- colnames(df)

gene_ids= readFile(opt$geneIDs) 
print(gene_ids)

print(head(fpkm_counts.intsct.AM))
fpkm_counts.intsct.AM.u <- fpkm_counts.intsct.AM[,samples_list] # subsetting the expMatrix based on samples list

##Assert if the geneIDs are SUBSET of the intsct.sig.genes.AM! ... so the user has to know the geneIDs which are "commonly" found
fpkm_counts.intsct.AM.u <- fpkm_counts.intsct.AM.u[order(-rowMeans(fpkm_counts.intsct.AM.u)),][gene_ids,]

user.intsct.AM <- rownames(fpkm_counts.intsct.AM.u) #just the gene ids

user.intsct.AM.genes <- fpkm_counts_dummy[user.intsct.AM, 16]# The ordered geneNames
user.intsct.AM.names.row <- sub('^NA:', '', x = paste(user.intsct.AM.genes, user.intsct.AM, sep = ":"))
heatmap.2(fpkm_counts.intsct.AM.u, labRow = user.intsct.AM.names.row, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
par(lend=1)
legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
dev.off()
##############################

fpkm_counts.intsct.UK.u <- fpkm_counts.intsct.UK[,samples_list]

fpkm_counts.intsct.UK.u <- fpkm_counts.intsct.UK.u[order(-rowMeans(fpkm_counts.intsct.UK.u)),][gene_ids,]
user.intsct.UK <- rownames(fpkm_counts.intsct.UK.u)

user.intsct.UK.genes <- fpkm_counts_dummy[user.intsct.UK, 16]
user.intsct.UK.names.row <- sub('^NA:', '', x = paste(user.intsct.UK.genes, user.intsct.UK, sep = ":"))
png("fpkm_counts.intsct.UK.png", height = 5*300, width = 5*300, res=300, pointsize = 8)
heatmap.2(fpkm_counts.intsct.UK.u, labRow = user.intsct.UK.names.row, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
par(lend = 1)
legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
dev.off()
################################################
#using union of all significant genes...
fpkm_counts.union.UK=as.matrix(fpkm_counts[union.sig.genes.UK,])
nrow(fpkm_counts.union.UK)

fpkm_counts.union.AM=as.matrix(fpkm_counts[union.sig.genes.AM,])     
nrow(fpkm_counts.union.AM)


#fpkm_counts.union.UK <- fpkm_counts.union.UK[order(-rowMeans(fpkm_counts.union.UK)),][c(1:70),]

fpkm_counts.union.UK.u <- fpkm_counts.union.UK[gene_ids,samples_list]
fpkm_counts.union.UK.u <- fpkm_counts.union.UK.u[order(-rowMeans(fpkm_counts.union.UK.u)),]

user.union.UK <-  rownames(fpkm_counts.union.UK.u)
user.union.UK.genes <- fpkm_counts_dummy[user.union.UK, 16]
user.union.UK.names.row <- sub('^NA:', '', x = paste(user.union.UK.genes, user.union.UK, sep = ":"))
png("fpkm_counts.union.UK.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
heatmap.2(fpkm_counts.union.UK.u, labRow = user.union.UK.names.row, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
par(lend=1)
legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
dev.off()
#######

fpkm_counts.union.AM.u <- fpkm_counts.union.AM[gene_ids,samples_list]
fpkm_counts.union.AM.u <- fpkm_counts.union.AM.u[order(-rowMeans(fpkm_counts.union.AM.u)),]

user.union.AM <- rownames(fpkm_counts.union.AM.u)
user.union.AM.genes <- fpkm_counts_dummy[user.union.AM, 16]
user.union.AM.row.names <- sub('^NA:', '', x = paste(user.union.AM.genes, user.union.AM, sep = ":"))

png("fpkm_counts.union.AM.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
heatmap.2(fpkm_counts.union.AM.u, labRow = user.union.AM.row.names, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,cexCol=1.0,ColSideColors = c(rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2)))
par(lend=1)
legend(x=0.78,y=1, legend=c("X717c.Control_A_Shoot","X717c.Control.A_Axillary.bud","IAA7.52_shoot.tip","IAA7.52_axillary.bud","IAA7.54_shoot.tip","IAA7.54_axillary.bud"), col = c("gray","blue","black","orange","darkgreen","navyblue"), lty=1:2,lwd=10, box.lty = 0,cex=0.5)
dev.off()
#########

# snames=colnames(fpkm_counts)
# c=0
# for (items in snames) {
#   c=c+1
#   fpkm_counts[[c]]=log2(fpkm_counts$items+1)
# }
