#!/usr/bin/env Rscript

#####Load packages#####
if(!require(ape)){install.packages("ape")}
if(!require(dynamicTreeCut)){install.packages("dynamicTreeCut")}
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
if(!require(pheatmap)){install.packages("pheatmap")}

dynamicTreeCut_version <- packageVersion("dynamicTreeCut")
print (paste("dynamicTreeCut ", dynamicTreeCut_version, sep=""))
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
tpms_file <- args[1] #"../01_INPUT/UNvsIN_TPMs.txt"
sig_file <- args[2] #"../02_Analyses/UNvsIN_DGEA_DifferentiallyExpressedGenes.txt"
##########



#####Load parameters
experiment_name <- strsplit(basename(tpms_file), "_")[[1]][1]

txt_file.cor_gene <-"CLUSTER_corrilation_gene.txt"
txt_file.cor_sample <- "CLUSTER_corrilation_sample.txt"
tre_file.hierar_gene <- "CLUSTER_hierarchy_gene.tre"
tre_file.hierar_sample <- "CLUSTER_hierarchy_sample.tre"

txt_file.GeneModules <- "CLUSTER_modules.txt"
txt_file.HeatMap <- "CLUSTER_Heatmap.txt"

pdf_file.HeatMap <- "CLUSTER_Heatmap_ALL.pdf"
pdf_file.HeatMapBest <- "CLUSTER_Heatmap_BEST.pdf"

#Set up parameters
tpms_table <- read.table(tpms_file)
colnames(tpms_table) <- as.matrix(tpms_table[2,])
treatment <- as.factor(as.matrix(tpms_table[1,][,-1]))
row.names(tpms_table) <- as.matrix(tpms_table[,1])
tpms_table <- as.matrix(tpms_table[-1,-1][-1,])
class(tpms_table) <- "numeric"

Sig.Genes <- read.table(sig_file)
Sig.Genes <- rownames(Sig.Genes[Sig.Genes$DE==TRUE,])

Sign <- args[3] #Corrilation distance calculation "signed" or "unsigned"
MinModuleSize <- as.numeric(args[4]) #Minimum module size (eg. 30)
DS <- as.numeric(args[5]) #Degree of separation DS (0-4)
h <- args[6] #desired threshold for h, "AUTO" or fraction multiples of 0.05 (e.g. 0.65)
##########



###01: filtering DEGs###
expData_DEG <- tpms_table[rownames(tpms_table)%in%Sig.Genes,]
######



###02: Pearson's correlation###
#Pearson
distMethod <- "pearson"
TPM_cor_gene <- cor(t(expData_DEG), use="everything", method=distMethod)
TPM_cor_sample <- cor(as.matrix(expData_DEG), use="everything", method=distMethod)
if (Sign == "unsigned") {
  TPM_cor_gene <- abs(TPM_cor_gene)
  TPM_cor_sample <- abs(TPM_cor_sample)
}
print(paste("Distance method: ", distMethod, sep=""))

#Correlation distance
TPM_cordist_gene <- 1-TPM_cor_gene
TPM_cordist_sample <- 1-TPM_cor_sample

#save data
write.table(TPM_cor_gene, file = txt_file.cor_gene, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
write.table(TPM_cor_sample, file = txt_file.cor_sample, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
######



###03: Hierarchical clustering###
clusterMethod <- "complete"

#gene
TPM_hierar_gene <- hclust(as.dist(TPM_cordist_gene), method=clusterMethod)
TPM_hierar_gene$height <- (TPM_hierar_gene$height/max(TPM_hierar_gene$height))*2
write.tree(phy=ladderize(as.phylo(TPM_hierar_gene), right = TRUE), file=tre_file.hierar_gene)

#sample
TPM_hierar_sample <- hclust(as.dist(TPM_cordist_sample), method=clusterMethod)
TPM_hierar_sample$height <- (TPM_hierar_sample$height/max(TPM_hierar_sample$height))*2
write.tree(phy=ladderize(as.phylo(TPM_hierar_sample), right = TRUE), file=tre_file.hierar_sample)

print(paste("Agglomeration method for clustering: ", clusterMethod, sep=""))
######



###04: Module detection###
cuttreeMethod <- "hybrid"
for (i in 0:20) {
  t <- 0.05 * i
  c <- i + 1
  Modules <- cutreeDynamic(
    TPM_hierar_gene,
    method = cuttreeMethod,
    distM <- as.matrix(TPM_cordist_gene),
    cutHeight = t,
    pamRespectsDendro = TRUE,
    minClusterSize = MinModuleSize,
    deepSplit = DS
  )
  if (c == 1) {
    GeneModules <- data.frame(Modules)
  }
  else {
    GeneModules[,c] <- Modules
  }
  colnames(GeneModules)[c] <- paste("h=", t, sep="")
}
GeneModules[GeneModules == 0] <- NA
rownames(GeneModules) <- as.phylo(TPM_hierar_gene)$tip.label
write.table(GeneModules, file = txt_file.GeneModules, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

print(paste("dynamicTreeCut method: ", cuttreeMethod, sep=""))
######



###05: Heatmap###
TPM_heatmap <- as.data.frame(t(scale(t(expData_DEG))))
TPM_heatmap <- as.data.frame(scale(TPM_heatmap))
for (i in 1:ncol(TPM_heatmap)) {
  TPM_heatmap[which(TPM_heatmap[,i] < -3), i] <- -3
  TPM_heatmap[which(TPM_heatmap[,i] > 3), i] <- 3
}
TPM_heatmap <- TPM_heatmap[as.phylo(TPM_hierar_gene)$tip.label,]
colnames(TPM_heatmap) <- colnames(expData_DEG)
TPM_heatmap <- TPM_heatmap[,as.phylo(TPM_hierar_sample)$tip.label]
write.table(TPM_heatmap, file = txt_file.HeatMap, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
######



###06: Plot preparation###
#preparing coloring for samples
treatment_table <- data.frame(Treatment = treatment)
row.names(treatment_table) <- colnames(TPM_heatmap)

#preparing coloring for modules
module_table <- GeneModules
for (i in 1:ncol(GeneModules)) {
  x <- module_table[,i]
  x[is.na(x)] <- 0
  module_table[i] <- as.factor(as.matrix(x))
}

#preparing coloring for best threshold
if (h == "AUTO") {
  best = "h=0"
  BestScore = 0
  for(i in 1:ncol(GeneModules)) {
    x <- GeneModules[,i]
    x[is.na(x)] <- 0
    ModuleCount <- max(x, na.rm=TRUE)
    if (ModuleCount > BestScore) {
      best <- colnames(GeneModules)[i]
      BestScore <- ModuleCount
    }
  }
  h <- as.numeric(strsplit(best, split="=")[[1]][2])
} else {
  h <- as.numeric(h)
}
bestmodule <- (h/5*100) + 1
bestmodule_table <- as.data.frame(as.factor(GeneModules[,bestmodule]))
row.names(bestmodule_table) <- rownames(GeneModules)
colnames(bestmodule_table) <- colnames(GeneModules)[bestmodule]
ModuleCount <- max(as.numeric(bestmodule_table[,1]), na.rm=TRUE)

print(paste("dynamicTreeCut height threshold: ", h, sep=""))
print(paste("dynamicTreeCut module count: ", ModuleCount, sep=""))

#preparing colors
color_table <- vector(mode="list", length=ncol(GeneModules) + 1)
color_table[[1]] <- brewer.pal(n = length(levels(treatment_table[,1])), name="Accent")
names(color_table[[1]]) <- levels(treatment_table[,1])
colors <- sample(unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',]))), max(GeneModules, na.rm=TRUE) + 1)
colors[1] <- "#FFFFFF"
for (i in 1:ncol(GeneModules)) {
  color_table[[i + 1]] <- colors
  names(color_table[[i + 1]]) <- as.character(seq(from=0,to=max(GeneModules, na.rm=TRUE)))
}
names(color_table) <- c("Treatment", colnames(GeneModules))
######



###07: Generate plot###
#generate heatmap
heatmap_all <- pheatmap(
  mat = TPM_heatmap,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.5),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "z-score\n"),
  cluster_rows = TPM_hierar_gene,
  cluster_cols = TPM_hierar_sample,
  show_rownames = FALSE,
  annotation_col = treatment_table,
  annotation_row = module_table,
  annotation_colors = color_table,
  annotation_names_col = FALSE,
  annotation_legend = FALSE
)

heatmap_best <- pheatmap(
  mat = TPM_heatmap,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.5),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "z-score\n"),
  cluster_rows = TPM_hierar_gene,
  cluster_cols = TPM_hierar_sample,
  show_rownames = FALSE,
  annotation_col = treatment_table,
  annotation_row = bestmodule_table,
  annotation_colors = color_table,
  annotation_names_col = FALSE,
  annotation_names_row = TRUE,
  annotation_legend = FALSE
)

#print heatmap
pdf(pdf_file.HeatMap, width = 12, height = 8)
heatmap_all
dev.off()

pdf(pdf_file.HeatMapBest, width = 12, height = 8)
heatmap_best
dev.off()
######

