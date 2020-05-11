#!/usr/bin/env Rscript

#####Load packages#####
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(reshape)){install.packages("reshape")}
if(!require(ggpubr)){install.packages("ggpubr")}
##########



#####Load files#####
args <- commandArgs(trailingOnly=TRUE)
tpms_file <- args[1] #"../01_INPUT/UNvsIN_TPMs.txt"
module_file <- args[2] #"../02_Analyses/UNvsIN_CLUSTER_Modules.txt"
go_file <- args[3] #"../02_Analyses/UNvsIN_GO_EnrichedTerms.txt"
##########



#####Load parameters
experiment_name <- strsplit(basename(tpms_file), "_")[[1]][1]

pdf_file.module_pattern <-"GO_module_summary.pdf"

#Set up parameters
tpms_table <- read.table(tpms_file)
colnames(tpms_table) <- as.matrix(tpms_table[2,])
treatment <- as.factor(as.matrix(tpms_table[1,][,-1]))
row.names(tpms_table) <- as.matrix(tpms_table[,1])
tpms_table <- as.matrix(tpms_table[-1,-1][-1,])
class(tpms_table) <- "numeric"

module_table <- read.table(module_file)
h <- as.numeric(args[4]) #best threshold for h, fraction multiples of 0.05 (e.g. 0.85)
module_table <- data.frame(best = module_table[,(h/5*100+1)], row.names=rownames(module_table))

go_table <- read.table(go_file, sep="\t", row.names=NULL, header=TRUE)
##########



#####MAIN#####
#Make table average per treatment
Average_tpm <- NULL
for (i in 1:length(levels(treatment))) {
    Average_tpm <- cbind(Average_tpm, rowMeans(tpms_table[,which(treatment==levels(treatment)[i])]))
    colnames(Average_tpm)[i] <- levels(treatment)[i]
}
relxp <- function(x){
    x/max(x)
}
#Make table of relative expression per treatment
Relative_exp <- t(apply(Average_tpm, 1, relxp))

#Plotting
pdf(pdf_file.module_pattern, width = 12, height = 12, onefile=TRUE)
for (i in 1:max(module_table, na.rm=TRUE)) {
    total <- as.data.frame(t(Relative_exp[which(rownames(Relative_exp) %in% rownames(module_table)[which(module_table==i)]),]))
    
    exp <- total
    exp$Treatment <- as.factor(rownames(exp))
    exp <- melt(exp)
    colnames(exp)[2:3] <- c("GeneID", "RelExp")
    exp$Treatment <- as.numeric(exp$Treatment)
    exp$GeneID <- as.factor(exp$GeneID)
    
    quant <- as.data.frame(t(apply(total, 1, quantile, seq(0, 1, 0.25))))
    quant$Treatment <- rownames(quant)
    quant <- melt(quant)
    colnames(quant)[2:3] <- c("Quantile", "RelExp")
    quant$Treatment <- as.numeric(as.factor(quant$Treatment))
    
    MaxMin <- rbind(quant[which(quant$Quantile=="100%"),],quant[which(quant$Quantile=="0%"),][order(quant[which(quant$Quantile=="0%"),]$Treatment, decreasing = TRUE),])
    MaxMin$Quantile <- "MaxMin"
    
    ThirdFirst <- rbind(quant[which(quant$Quantile=="75%"),],quant[which(quant$Quantile=="25%"),][order(quant[which(quant$Quantile=="25%"),]$Treatment, decreasing = TRUE),])
    MaxMin$Quantile <- "ThirdFirst"
    
    Median <- quant[which(quant$Quantile=="50%"),]
    
    go <- go_table[which(go_table$Module==i),2:8]
    colnames(go)[c(4,7)] <- c("Observed", "p-val")
    
    #expressin pattern
    plot.pattern <- ggplot(data=exp, aes(x=Treatment, y=RelExp)) +
        geom_polygon(
            data=MaxMin,
            fill="blue",
            alpha = 0.1,
            group=1) +
        geom_polygon(
            data=ThirdFirst,
            fill="blue",
            alpha = 0.5,
            group=1) +
        geom_line(color="black", alpha=0.1, aes(group=GeneID)) +
        geom_line(data=Median, color="red") +
        xlim(1,length(levels(treatment)))+
        scale_x_continuous(
            breaks=seq(1, length(levels(treatment)), 1),
            labels=levels(treatment)) +
        labs(title=paste("Module", i, sep="_"),
             x ="Treatment",
             y = "Relative Expression Level") +
        theme_bw()
    
    if (nrow(go) > 0) {
        #go term table
        plot.goterm <- ggtexttable(go, rows = NULL, theme = ttheme("mOrange"))
        
        #Arrange the plots on the same page
        plot.total <- ggarrange(plot.pattern, plot.goterm, ncol = 1, nrow = 2, heights = c(1, 3), widths=c(1,1), align="v")
    } else {
        #Arrange the plots on the same page
        plot.total <- ggarrange(plot.pattern, NULL, ncol = 1, nrow = 2, heights = c(1, 3), widths=c(1,1), align="v")
    }
    print(plot.total)
}
dev.off()
#####END#####
