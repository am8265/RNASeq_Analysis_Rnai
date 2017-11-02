#!/usr/bin/env Rscript
library(DESeq2)
library(GOexpress)
library(gplots)
library(optparse)

###Run the script as :

## ./Hmap_custom.R --sample2Condition <sample2condiiton.csv>  --geneIDs <geneIDs.txt> --countMatrix <countMatrix.csv>
##The script doesn't assess if the gene meets a cutoff criteria/ DE... it just provides a way to visualise. Hence the user should provide a list of DEGs previosuly found by other method
## Example Run: /Hmap_custom.R --geneIDs gene_ids.txt --sample2Condition samples2Condition.csv --countMatrix fpkm_counts.csv

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
  make_option(c("--sample2Condition"), type = "character", default = NULL,
              help = "csv file with sample names and conditions to be shown on the HMap"),
  make_option(c("--countMatrix", type = "character", default = NULL,
                help = "csv file with normalised expression values having columns as gene IDs and sample names"))
  )


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

fpkm_counts = read.csv(opt$countMatrix, header = T, row.names = 1, check.names = F) # check.names =F so that the sampels names are kept the same 
# as the sampleList
#applying the log2 transformation for visualisation
fpkm_counts[] = lapply(fpkm_counts, function(x){log2(x+1)})


gene_ids= readFile(opt$geneIDs)

sample2condition = read.csv(opt$sample2Condition, header = F, stringsAsFactors = F, check.names = F)
sample_list =as.character(sample2condition[,1]) # samplesList
condition_list = as.character(sample2condition[,2]) #conditionList

## groupings for samples ...
sampleGroup = tapply(sample2condition[,1], sample2condition[,2], length)



countMatrix = fpkm_counts[gene_ids,sample_list]
countMatrix = as.matrix(countMatrix)

color_choice = c("gray","blue", "black", "orange", "darkgreen", "navyblue")# assuming that is the maximum number of conditions for which that many maximum colors are provided 
color_choice_made = color_choice[1:length(sampleGroup)]

#for (i in 1:length(color_choice_made)){

      
      

###Plotting the HeatMap
hmcol = colorRampPalette(c("red","yellow","green"))(n=299)
png("heatmap.png",height = 5*300, width = 5*300, res=300, pointsize = 8)
heatmap.2(countMatrix, labRow = gene_ids, offsetRow = -0.6, 
          col = hmcol, trace="none", margin=c(6, 20), srtCol=30, offsetCol =-1.5, cexCol=1.0, ColSideColors = rep(color_choice_made, times = as.vector(sampleGroup)))
par(lend=1)
legend(x=0.78,y=1, legend=c(names(sampleGroup)), col = color_choice_made, lty=1:2,lwd=10, box.lty = 0,cex=0.5)
dev.off()



