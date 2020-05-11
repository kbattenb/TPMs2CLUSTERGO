#!/usr/bin/env Rscript

#####Loading Options#####
options(scipen=999)
##########

#####Load packages#####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("topGO")
require(topGO)

topGO_version <- packageVersion("topGO")
print (paste("topGO ", topGO_version, sep=""))
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
module_file <- args[1] #"../02_Analyses/UNvsIN_CLUSTER_Modules.txt"
GO_file <- args[2] #"../INPUT/Annotation.txt"
##########



#####Load parameters#####
h <- as.numeric(args[3]) #best threshold for h, fraction multiples of 0.05 (e.g. 0.85)
fdr <- as.numeric(args[4]) #p-value for go enrichment (e.g. 0.05)
m <- as.numeric(args[5]) #minimum genecount with a given go term to be kept as significant during Fisher's exact test (e.g. 3)
algorithm <- "classic"
print(paste("topGO test algorithm: ", algorithm, sep=""))
statistic <- "fisher"
print(paste("topGO test statistic: ", statistic, sep=""))

experiment_name <- strsplit(basename(module_file), "_")[[1]][1]
txt_file.go <- "GO_EnrichedPerModule.txt"

#Set up parameters
module_table <- read.table(module_file)
module_table <- data.frame(best=module_table[,h/0.05 + 1], row.names=row.names(module_table))

GO_table <- read.table(GO_file, fill=TRUE, sep="\t", header=TRUE)
GO_table <- data.frame(GO=GO_table[,4], row.names=GO_table[,1])
##########



#####MAIN#####
#set a reference for topGO
write.table(GO_table, file = "topGO_ref_temp.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
geneID2GO <- readMappings(file = "topGO_ref_temp.txt")
unlink("topGO_ref_temp.txt")

#Get all enriched GO terms
total.result.table <- NULL
for (i in 1:max(module_table[1], na.rm=TRUE)) {
    #select a set of genes of interest to analyze
    genesOfInterest <- rownames(module_table)[which(module_table[1] == i)]
    geneList <- factor(as.integer(rownames(GO_table) %in% genesOfInterest))
    names(geneList) <- rownames(GO_table)
    
    #run enrichment for 3 branches of GO terms
    topGOdataBP <- new(
        "topGOdata",
        ontology="BP", #Biological Process
        allGenes=geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO)
    result.Fisher.BP <- runTest(topGOdataBP,
        algorithm = algorithm,
        statistic = statistic)
    result.table.Fisher.BP <- GenTable(topGOdataBP,
        Fisher = result.Fisher.BP,
        orderBy = "Fisher",
        topNodes = length(usedGO(object = topGOdataBP)))
    result.table.Fisher.BP$adjP <- p.adjust(result.table.Fisher.BP[,6], method = "BH", n = length(result.table.Fisher.BP[,6]))
    result.table.Fisher.BP <- result.table.Fisher.BP[which(result.table.Fisher.BP$adjP<fdr),]
    if (nrow(result.table.Fisher.BP) > 0) {
        result.table.Fisher.BP$Ontology <- "Biological Process"
    }
    
    topGOdataMF <- new(
        "topGOdata",
        ontology="MF", #Molecular Function
        allGenes=geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO)
    result.Fisher.MF <- runTest(topGOdataMF,
                                algorithm = algorithm,
                                statistic = statistic)
    result.table.Fisher.MF <- GenTable(topGOdataMF,
                                       Fisher = result.Fisher.MF,
                                       orderBy = "Fisher",
                                       topNodes = length(usedGO(object = topGOdataMF)))
    result.table.Fisher.MF$adjP <- p.adjust(result.table.Fisher.MF[,6], method = "BH", n = length(result.table.Fisher.MF[,6]))
    result.table.Fisher.MF <- result.table.Fisher.MF[which(result.table.Fisher.MF$adjP<fdr),]
    if (nrow(result.table.Fisher.MF) > 0) {
        result.table.Fisher.MF$Ontology <- "Molecular Function"
    }
    
    topGOdataCC <- new(
        "topGOdata",
        ontology="CC", #Cellular Component
        allGenes=geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO)
    result.Fisher.CC <- runTest(topGOdataCC,
                                algorithm = algorithm,
                                statistic = statistic)
    result.table.Fisher.CC <- GenTable(topGOdataCC,
                                       Fisher = result.Fisher.CC,
                                       orderBy = "Fisher",
                                       topNodes = length(usedGO(object = topGOdataCC)))
    result.table.Fisher.CC$adjP <- p.adjust(result.table.Fisher.CC[,6], method = "BH", n = length(result.table.Fisher.CC[,6]))
    result.table.Fisher.CC <- result.table.Fisher.CC[which(result.table.Fisher.CC$adjP<fdr),]
    if (nrow(result.table.Fisher.CC) > 0) {
        result.table.Fisher.CC$Ontology <- "Cellular Component"
    }
    
    #generage a summary table
    result.table.Fisher <- rbind(result.table.Fisher.BP, result.table.Fisher.MF, result.table.Fisher.CC)
    if (nrow(result.table.Fisher) > 0) {
        result.table.Fisher <- result.table.Fisher[which(result.table.Fisher$Significant>=m),]
    }
    if (nrow(result.table.Fisher) > 0) {
        result.table.Fisher$Module <- i
        total.result.table <- rbind(total.result.table, result.table.Fisher)
    }
}
total.result.table <- total.result.table[,c(9,8,1,2,4,5,3,7)]
colnames(total.result.table)[5] <- "Observed"

#print out table
write.table(x=total.result.table, file=txt_file.go, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)
#####END#####
