#!/usr/bin/env Rscript
library(DESeq2)
library(GOexpress)
library(gplots)
library(optparse)

##Run the script as 
##./FC_HMAP.R --geneIDs geneList.txt --conditionPair conditionPair.txt --condPair_filesMap condPair_filesMap --geneDesc geneDescriptions.tsv


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

## cutoff functions

mod_cutoff_criteria <- function(x, log_fc_up = 1, log_fc_down = -1, padj = 0.05) {
  x_rmNA = x[complete.cases(x),] #remove na values 
  cutoff = (x_rmNA[[2]] >= log_fc_up | x_rmNA[[2]] <= log_fc_down) & (x_rmNA[[3]] <= padj)
  deg_count = sum(cutoff)
  percent_deg_count = (deg_count/nrow(x)) * 100
  deg_genes = x_rmNA[[1]][cutoff] #DEG meeting the cutoff
  return (list(deg_count = deg_count, percent_deg_count = percent_deg_count, deg_genes = deg_genes))
}

mod_up_genes <- function(x, log_fc_up = 1, padj = 0.05) {
  mod_cutoff_criteria(x,log_fc_up = log_fc_up, log_fc_down = -Inf, padj = padj)
}

mod_down_genes <- function(x, log_fc_down = -1, padj = 0.05) {
  mod_cutoff_criteria(x,log_fc_up = Inf, log_fc_down = log_fc_down, padj = padj)
}



##Taking in the options
option_list = list(
  make_option(c("--geneIDs"), type = "character", default = NULL,
              help = "List of geneIDs to be shown on the HMap"),
  make_option(c("--conditionPair"), type = "character", default = NULL,
              help = "text file with condtionalPair names to be shown on the HMap"),
  make_option(c("--condPair_filesMap"),type = "character", default = NULL,
                help = "tsv file with condPair Name (1st column) and path to the condPairFile"),
  make_option(c("--geneDesc"),type = "character", default = NULL,
              help = "tsv file with Feature(gene) column and Description column")
  )

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


condition_pair = readFile(opt$conditionPair) #getiing the conditionalPair to display

condPair_file_path = read.table(opt$condPair_filesMap, header = F, sep = '\t', check.names = F, stringsAsFactors = F)
condPair = condPair_file_path[,1]
file_path = condPair_file_path[,2]

conditionalPairList = lapply(file_path, FUN = function(x){read.csv(x, header = T, stringsAsFactors = F, colClasses = c("character","numeric","numeric", "logical"))})
names(conditionalPairList) = condPair #Assigning IDENTIFIABLE Names for condPairs

##Applying the cutoff on conditionalPairList to Assess DEGs and write to a file 'intersect_genes.txt' listing all genes commonly Found in the conditional Pair
conditionalPair_result = lapply(conditionalPairList, mod_cutoff_criteria)

deg_genes = lapply(conditionalPair_result, function(x) {x["deg_genes"]})
intersect_genes = deg_genes[[1]][["deg_genes"]]
for (i in 1:length(deg_genes)) {
  print(i)
  intersect_genes = intersect(intersect_genes, deg_genes[[i]][["deg_genes"]])
}

write(intersect_genes, file = "intersect_genes.txt", ncolumns = 1)# the user should provide a subset of genes from the list provided

##intersecting list of genes  from the condPairs
intersect_genes = readFile('intersect_genes.txt')

gene_list = readFile(opt$geneIDs)

common_genes = intersect(intersect_genes, gene_list)

assertthat::assert_that(length(common_genes) > 0)##Assert if user's gene_list has common genes with intersect_genes

   





##substitute test with appropriate name 

test = lapply(conditionalPairList, FUN = function(x) {x[x$Gene.ID %in% gene_list,][,c(1,2)]})
## renaming the FC columns appropriately!
for (i in 1:length(test)){
  colnames(test[[i]])[2] <- names(test[i])
}
merge_condition_pairs = test[[1]]
for (i in 2:length(test)){
  merge_condition_pairs = merge(merge_condition_pairs , test[[i]], by = c("Gene.ID"))
}
  
merge_condition_pairs = as.matrix(merge_condition_pairs)

geneDescriptions = read.table(opt$geneDesc, header = T, row.names = 1, stringsAsFactors = F)
gene_with_names = sub('^NA:', '', x = paste(merge_condition_pairs[,1], geneDescriptions[merge_condition_pairs[,1],], sep = ":"))

## heatMaps plot

hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
png("FC_HeatMap.png",height = 5*300, width = 5*300, res=300, pointsize = 8)

legend_colors =c("gray","blue","black","orange","darkgreen","navyblue") # # assuming that is the maximum number of conditionPairs for which that many maximum colors(6) are provided
merge_condition_matrix = apply(merge_condition_pairs[,-1], 2, as.numeric)
heatmap.2(merge_condition_matrix, labRow = gene_with_names, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20),srtCol=30, offsetCol =-1.5,
          cexCol=1.0,ColSideColors = legend_colors[1:ncol(merge_condition_pairs[,-1])])
par(lend = 1)
print(colnames(merge_condition_pairs[,-1]))
print(legend_colors[1:ncol(merge_condition_pairs[,-1])])
legend(x = 0.78, y = 1, legend = colnames(merge_condition_pairs[,-1]), col = legend_colors[1:ncol(merge_condition_pairs[,-1])], lty=1:2, lwd=10, box.lty = 0, cex=0.5)
dev.off()


#rep("gray",3),rep("blue",3),rep("black",2),rep("orange",2),rep("darkgreen",3),rep("navyblue",2))
