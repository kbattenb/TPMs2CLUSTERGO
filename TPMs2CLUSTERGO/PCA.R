#!/usr/bin/env Rscript

#####Load packages#####
if(!require(ggplot2)){install.packages("ggplot2")}
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
tpms_file <- args[1] #"../01_INPUT/UNvsIN_TPMs.txt"
##########



#####Load function#####
compute_cv <- function(x) sd(x) / mean(x)
##########



#####Load parameters#####
experiment_name <- strsplit(basename(tpms_file), "_")[[1]][1]
pdf_file.1v2 <- "PCA_PC1vPC2.pdf"
pdf_file.2v3 <- "PCA_PC2vPC3.pdf"
pdf_file.contribution <- "PCA_PercentContribution.pdf"
txt_file.top10 <- "PCA_Top10Contributors.txt"

#Set up parameters
tpms_table <- read.table(tpms_file)
colnames(tpms_table) <- as.matrix(tpms_table[2,])
treatment <- as.factor(as.matrix(tpms_table[1,][,-1]))
row.names(tpms_table) <- as.matrix(tpms_table[,1])
tpms_table <- as.matrix(tpms_table[-1,-1][-1,])
class(tpms_table) <- "numeric"
##########



#####MAIN#####
#remove rows with only 0 values
tpms_table <- tpms_table[apply(tpms_table, 1, function(x) !all(x==0)),]

#log transform
tpms_table <- log(tpms_table + 0.01)

#select top 25% genes with the highest covariance (most variable genes)
#REFERENCE: https://jdblischak.github.io/singlecell-qtl/pca-variable.html
cutoff <- 0.25
tpms_table <- tpms_table[rank(apply(tpms_table, 1, compute_cv)) / nrow(tpms_table) > 1 - 0.25, ]

#running PCA
pca <- prcomp(t(tpms_table), scale=TRUE)
PC.contribution <- data.frame(PCs=paste("PC", seq(1, ncol(pca$x)), sep=""), Contribution=round((pca$sdev^2)/sum(pca$sdev^2)*100, 1))
PC.contribution$PCs <- factor(PC.contribution$PCs, levels = as.character(PC.contribution$PCs))

pca.data_1v2 <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
pca.data_2v3 <- data.frame(Sample=rownames(pca$x), X=pca$x[,2], Y=pca$x[,3])

lab1 <- paste("PC1(", PC.contribution$Contribution[1], "%)", sep="")
lab2 <- paste("PC2(", PC.contribution$Contribution[2], "%)", sep="")
lab3 <- paste("PC3(", PC.contribution$Contribution[3], "%)", sep="")

ContributingGenes <- data.frame(PC1=names(sort(abs(pca$rotation[,1]), decreasing=TRUE)[1:10]), PC2=names(sort(abs(pca$rotation[,2]), decreasing=TRUE)[1:10]), PC3=names(sort(abs(pca$rotation[,3]), decreasing=TRUE)[1:10]))

#plotting PCA
plot.PC1vPC2 <- ggplot(data=pca.data_1v2, aes(x=X, y=Y, label=Sample, color=treatment)) +
    geom_text() +
    xlab(lab1) +
    ylab(lab2) +
    theme_bw()

plot.PC2vPC3 <- ggplot(data=pca.data_2v3, aes(x=X, y=Y, label=Sample, color=treatment)) +
    geom_text() +
    xlab(lab2) +
    ylab(lab3) +
    theme_bw()

plot.contribution <- ggplot(PC.contribution, aes(x=PCs, y=Contribution)) +
    geom_bar(stat="identity") +
    ylab("Contribution(%)") +
    theme_bw()

#printing data
pdf(pdf_file.1v2, width = 12, height = 8)
plot.PC1vPC2
dev.off()

pdf(pdf_file.2v3, width = 12, height = 8)
plot.PC2vPC3
dev.off()

pdf(pdf_file.contribution, width = 12, height = 8)
plot.contribution
dev.off()

write.table(x=ContributingGenes, file=txt_file.top10, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)
#####END#####
