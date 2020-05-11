#!/usr/bin/env Rscript

#####Load packages#####
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
tpms_file <- args[1] #"../01_INPUT/UNvsIN_TPMs.txt"
design_file <- args[2] #"../01_INPUT/UNvsIN_Design.txt"
deg_file <- args[3] #"../02_Analyses/UNvsIN_DGEA_DifferentiallyExpressedGenes.txt"
module_file <- args[4] #"../02_Analyses/UNvsIN_CLUSTER_Modules.txt"
an_file <- args[5] #"../01_INPUT/UNvsIN_GeneAnnotation.txt"
##########






#####Load parameters
experiment_name <- strsplit(basename(tpms_file), "_")[[1]][1]

txt_file.sum <- "summary.txt"

#tpm table
tpms_table <- read.table(tpms_file)
colnames(tpms_table) <- as.matrix(tpms_table[2,])
treatment <- as.factor(as.matrix(tpms_table[1,][,-1]))
row.names(tpms_table) <- as.matrix(tpms_table[,1])
tpms_table <- as.matrix(tpms_table[-1,-1][-1,])
class(tpms_table) <- "numeric"

#unexpressed genes
unexpressed <- which(as.data.frame(apply(tpms_table, 1, max))==0)

#design table
design_table <- read.table(design_file)

#deg table
deg_table <- read.table(deg_file)

#module table
module_table <- read.table(module_file)
h <- as.numeric(args[6]) #best threshold for h, fraction multiples of 0.05 (e.g. 0.85)
module_table <- data.frame(best = module_table[,(h/5*100+1)], row.names=rownames(module_table))

#gene annotation
an_table <- read.table(an_file, fill=TRUE, sep="\t", header=TRUE)
an_table[an_table==""] <- NA
##########



#####MAIN#####
#Make table average per treatment
SUM_table <- NULL
for (i in 1:length(levels(treatment))) {
    SUM_table <- cbind(SUM_table, rowMeans(tpms_table[,which(treatment==levels(treatment)[i])]))
    colnames(SUM_table)[i] <- levels(treatment)[i]
}
SUM_table <- as.data.frame(SUM_table)

#add foldchanges
add_deg <- deg_table[,which(grepl("logFC", colnames(deg_table)))]
colnames(add_deg) <- paste("logFC", colnames(design_table), sep="_")
add_deg[unexpressed,] <- 0
add_deg <- cbind(add_deg, deg_table$adjP, deg_table$DE)
colnames(add_deg)[ncol(add_deg)-1] <- "adjP"
colnames(add_deg)[ncol(add_deg)] <- "DEG"
SUM_table <- cbind(SUM_table, add_deg)

#add modules
SUM_table$Module <- NA
for (i in 1:nrow(SUM_table)) {
    if (rownames(SUM_table)[i] %in% rownames(module_table)) {
        SUM_table$Module[i] <- module_table[which(rownames(module_table)==rownames(SUM_table)[i]),1][1]
    }
}

#add anotation
SUM_table <- cbind(SUM_table, an_table[,2:ncol(an_table)])
colnames(SUM_table)[ncol(SUM_table)] <- "GO(by Pfam)"

#print out table
write.table(x=SUM_table, file=txt_file.sum, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)
#####END#####
