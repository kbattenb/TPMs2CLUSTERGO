#!/usr/bin/env Rscript

#####Loading Options#####
options(scipen=999)
##########



#####Load packages#####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
require(edgeR)
#BiocManager::install("DESeq2")
require(DESeq2)
#BiocManager::install("genefilter")
require(genefilter)
if(!require(statmod)){install.packages("statmod")}
if(!require(ggplot2)){install.packages("ggplot2")}

edgeR_version <- packageVersion("edgeR")
print (paste("edgeR ", edgeR_version, sep=""))
DESeq2_version <- packageVersion("DESeq2")
print (paste("DESeq2 ", DESeq2_version, sep=""))
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1] #"../01_INPUT/UNvsIN_ExpReadCounts.txt"
design_file <- args[2] #"../INPUT/ExperimentalDesign.txt"
##########



#####Load parameters#####
fdr <- as.numeric(args[3]) #Benjamini-Hochberg false discovery rate (eg. 0.05)

experiment_name <- strsplit(basename(counts_file), "_")[[1]][1]
pdf_file.BestQuantile <- "DGEA_BestQuantile_threshold.pdf"
txt_file.DEG <- "DGEA_DEGs.txt"

#Set up parameters
counts_table <- read.table(counts_file)
colnames(counts_table) <- as.matrix(counts_table[2,])
treatment <- as.factor(as.matrix(counts_table[1,][,-1]))
row.names(counts_table) <- as.matrix(counts_table[,1])
counts_table <- as.matrix(counts_table[-1,-1][-1,])
class(counts_table) <- "numeric"
counts_table <- round(counts_table)
class(counts_table) <- "integer"
experimentaldesign <- read.table(design_file)
######



#####MAIN####
###determining the the optimal CPM threshold###
#generate design matrix
for (i in 1:ncol(experimentaldesign)) {
    assign(colnames(experimentaldesign)[i], value=factor(experimentaldesign[,i]))
}
design <- model.matrix(formula(paste("~", paste(colnames(experimentaldesign), collapse = "+"), sep="")))
rownames(design) <- colnames(counts_table)

#calculate raw p-value (using edgeR)
INPUT_DGEL <- DGEList(counts=counts_table, group=treatment)
INPUT_DGEL <- calcNormFactors(INPUT_DGEL)
INPUT_DGEL <- estimateDisp(INPUT_DGEL, design)
INPUT_QLFT <- glmQLFTest(glmQLFit(INPUT_DGEL, design, robust=TRUE), coef=2:ncol(design))

#print out CPM data
cpm_table <- cpm(INPUT_DGEL)

#combining cpm table and likelihoods
INPUT_results_unfiltered <- data.frame(cpm_table, INPUT_QLFT$table)
colnames(INPUT_results_unfiltered)[1:ncol(cpm_table)] <- colnames(cpm_table)
######



###Determine the threshold for filtering genes with low CPM (implementing methods from DESeq2)###
# NOTE: The "filter" is set to be equal to the second largest CPM value for each gene, rather
# than the average CPM value for each gene. Setting it to the second largest CPM makes filter act
# like edgeR's selection where at least 2 samples have to have CPMs larger than the minimum CPM
# threshold.
filter <- apply(X=cpm_table, MARGIN=1, FUN=function(data) {data[order(rank(data), decreasing=TRUE)[2]]})

#set the range of testing
lowerQuantile <- mean(filter == 0)
upperQuantile <- 0.95
if (lowerQuantile >= 0.95) {
    upperQuantile <- 1
}
theta <- seq(lowerQuantile, upperQuantile, length=50) #specific quantile thresholds to be tested
filtPadj <- filtered_p(filter=filter, test=INPUT_results_unfiltered$PValue, theta=theta, method="BH")

#number of rejected (significant) genes at each filtering threshold
numRej <- colSums(filtPadj < fdr, na.rm = TRUE)

#smoothing the curve
lo.fit.theta <- lowess(numRej ~ theta, f=1/5)

#get the quantile threshold
k <- 1
if (max(numRej) > 10) {
    residual <- 0
    if (all(! numRej==0)) {
        residual <- numRej[numRej > 0] - lo.fit.theta$y[numRej > 0]
    }
    thresh <- max(lo.fit.theta$y) - sqrt(mean(residual^2))
    if (any(numRej > thresh)) {
        k <- which(numRej > thresh)[1]
    }
}
threshold.quantile <- as.numeric(sub("%", "", names(numRej)[k]))
print (paste("Threshold quantile ", threshold.quantile, "%", sep=""))

#print quantile data
numRej.table <- data.frame(Quantile=as.numeric(sub("%", "", names(numRej))), Sig.Genes=as.numeric(numRej), Sig.Genes.fitted=lo.fit.theta$y)
deg.count <- numRej.table[which(numRej.table$Quantile == threshold.quantile),2]
print (paste("Differentially expressed gene count ", deg.count, sep=""))

#plot quantile threshold
label_x <- "Quantile (%)"
label_y <- paste("Number of Significant Genes (fdr<", fdr, ")", sep="")
label1 <- paste("Threshold quantile=", threshold.quantile, "%", sep="")

plot.threshold <- ggplot(data=numRej.table) +
    xlab(label_x) +
    ylab(label_y) +
    ylim(ceiling(min(numRej.table$Sig.Genes.fitted))/2, NA) +
    geom_rect(
        aes(xmin=-Inf, xmax=Inf, ymin=thresh, ymax=max(Sig.Genes.fitted)),
        fill = 'red',
        alpha = 0.01,
        inherit.aes=FALSE
    ) +
    geom_hline(
        yintercept = thresh,
        color = 'red'
    ) +
    geom_vline(
        xintercept = threshold.quantile,
        color = 'blue'
    ) +
    geom_point(
        aes(x=Quantile, y=Sig.Genes),
        alpha = 0.5
    ) +
    geom_line(
        aes(x=Quantile, y=Sig.Genes.fitted),
        inherit.aes=FALSE
    ) +
    geom_text(
        aes(x=threshold.quantile, y=min(Sig.Genes.fitted), label=label1),
        hjust=0
    ) + 
    theme_bw()
pdf(pdf_file.BestQuantile, width = 12, height = 8)
plot.threshold
dev.off()
######



###Filtering according to the determined threshold###
INPUT_results_filtered <- INPUT_results_unfiltered
INPUT_results_filtered$adjP <- filtPadj[, k, drop=TRUE]
INPUT_results_filtered$DE <- FALSE
INPUT_results_filtered$DE[which(INPUT_results_filtered$adjP < fdr)] <- TRUE

#write table
write.table(INPUT_results_filtered, file=txt_file.DEG, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
######



###Generating an MA plot###
#plot.ma <- ggplot(data=INPUT_results_filtered) +
#    xlab("LogCPM") +
#    ylab("LogFC") +
#    geom_point(
#        aes(x=logCPM, y=logFC, color=DE)
#    ) +
#    scale_color_manual(values = c('grey', 'red')) +
#    theme_bw() +
#    theme(legend.position = "none")
#
#pdf(pdf_file.MA, width = 12, height = 8)
#plot.ma
#dev.off()
######
##########
