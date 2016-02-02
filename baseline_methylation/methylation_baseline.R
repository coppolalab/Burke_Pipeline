library(RnBeads)
#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Reading and writing tables
library(readr)
library(openxlsx)

#For DE analysis
#library(Biobase)
#library(marray)
#library(limma)
#library(MASS)
library(matrixStats)
#library(lumi)
#library(lumiHumanIDMapping)
#library(FDb.InfiniumMethylation.hg19)
#library(annotate)
library(R.utils)
#library(wateRmelon)

#For batch correction and PEER
library(sva)
library(peer)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(parallel)

#Functional programming
library(magrittr)
library(purrr)
library(functional)

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../../FRDA project/common_functions.R")

sample.data <- read.xlsx("../phenotypedata/2015-9037 Sample Sheet.xlsx")
sample.data$barcode <- paste(sample.data$chip.ID, sample.data$stripe, sep = "_")
colnames(sample.data)[1] <- "External.ID"
sample.data$External.ID %<>% str_replace("EPS", "EPS\\_")

pheno.read <- read.xlsx("../phenotypedata/Dana pheno filter.xlsx")
colnames(pheno.read) %<>% str_replace("@.", "")
pheno.select <- select(pheno.read, ID:GENDER, FUGL.MEYER.ADMISSION, FUGL.MEYER.DISCHARGE, ETHNICITY, high.blood.pressure, BODY.MASS.INDEX)
colnames(pheno.select) <- c("External.ID", "Age", "Sex", "FM.admit", "FM.discharge", "Ethnicity", "HBP", "BMI")
pheno.select$Delta.fm <- pheno.select$FM.discharge - pheno.select$FM.admit
pheno.select$Recovered <- (pheno.select$Delta.fm > 10)
pheno.select$External.ID %<>% sprintf(fmt = "EPS_%03d")
pheno.select$HBP %<>% factor %>% droplevels

sample.merge <- join(sample.data, pheno.select)
colnames(sample.merge)[1] <- "Sample_Name"
write_csv(sample.merge, "../phenotypedata/2015-9037 Sample Sheet.csv")

idat.dir <- "../raw_methylation/2015-9037/"
sample.annotation <- "../phenotypedata/2015-9037 Sample Sheet.csv"
report.dir <- "./reports"

rnb.initialize.reports(report.dir)
logger.start(fname = NA)
parallel.setup(7)
options(fftempdir="~/tmp/Rtmp")

data.source <- c(idat.dir, sample.annotation)
result <- rnb.run.import(data.source = data.source, data.type = "infinium.idat.dir", dir.reports = report.dir)
rnb.set <- result$rnb.set
remove.unknown <- pheno(rnb.set)$Delta.fm %>% is.na %>% which
remove.dup <- which(pheno(rnb.set)$barcode == "3998919067_R06C02")
remove.all <- c(remove.unknown, remove.dup)
rnb.known <- remove.samples(rnb.set, samples(rnb.set)[remove.all])

rnb.run.qc(rnb.known, report.dir)

rnb.filter <- rnb.execute.context.removal(rnb.known)$dataset
rnb.filter <- rnb.execute.snp.removal(rnb.filter, snp = "any")$dataset
rnb.filter <- rnb.execute.sex.removal(rnb.filter)$dataset
rnb.greedy <- rnb.execute.greedycut(rnb.filter)
filter.sites <- rnb.greedy$sites
rnb.filter <- remove.sites(rnb.filter, filter.sites)

rnb.filter <- rnb.execute.na.removal(rnb.filter)$dataset
rnb.filter <- rnb.execute.variability.removal(rnb.filter, 0.005)$dataset

rnb.norm <- rnb.execute.normalization(rnb.filter, method = "bmiq", bgcorr.method = "methylumi.noob")
saveRDS.gz(rnb.norm, "./save/rnb.norm.rda")

#Apparently there aren't any hidden factors?
sva.object <- rnb.execute.sva(rnb.norm, cmp.cols = "Recovered", columns.adj = c("Age", "Sex", "Ethnicity", "HBP", "BMI"), numSVmethod = "be")

dred.sites <- rnb.execute.dreduction(rnb.norm)
dred.promoters <- rnb.execute.dreduction(rnb.norm, target = "promoters")
dred <- list(sites = dred.sites, promoters = dred.promoters)
pca.colors <- ifelse(pheno(rnb.norm)$Recovered == TRUE, "red", "blue")

dred.plot <- data.frame(dred.promoters$mds$euclidean[,1:2])
colnames(dred.plot) <- c("PCA1", "PCA2")
dred.plot %<>% mutate(Recovered = pheno(rnb.norm)$Recovered)
#Convert this plot ggplot
CairoPDF("pca_allsites", width = 6, height = 6)
p <- ggplot(data = dred.plot, aes(x = PCA1, y = PCA2, col = Recovered)) + geom_point()
print(p)
dev.off()

HBP.fix <- pheno(rnb.norm)$HBP %>% droplevels
rnb.norm <- addPheno(rnb.norm, HBP.fix, "HBP.fix")
rnb.options(exploratory.columns = c("Age", "Sex", "Ethnicity", "HBP.fix", "BMI", "Recovered"))
assoc <- rnb.execute.batcheffects(rnb.norm, pcoordinates = dred)
assoc.qc <- rnb.execute.batch.qc(rnb.norm, pcoordinates = dred)

clustering.sites <- rnb.execute.clustering(rnb.norm, region.type = "sites")
clustering.promoters <- rnb.execute.clustering(rnb.norm, region.type = "promoters")

promoters.beta <- meth(rnb.norm, type = "promoters")
promoters.m <- lumi::beta2m(meth(rnb.norm, type = "promoters"))
sites.beta <- meth(rnb.norm)
sites.m <- lumi::beta2m(meth(rnb.norm))

gen.heatmap(promoters.beta, 1000, "promoters_beta_heatmap", clustering.promoters)
gen.heatmap(promoters.m, 1000, "promoters_mvalue_heatmap", clustering.promoters)
gen.heatmap(sites.beta, 1000, "sites_beta_heatmap", clustering.sites)
gen.heatmap(sites.m, 1000, "sites_mvalue_heatmap", clustering.sites)

gen.heatmap <- function(meths, ngenes, file.name, cluster.object)
{
    sites.ordered <-  apply(meths, 1, mad) %>% order(decreasing = TRUE)
    sites.plot <- meths[sites.ordered[1:ngenes],]
    cluster.tree <- cluster.object[[7]]@result
    attr(cluster.tree, "class") <- "hclust"
    cluster.dendro <- as.dendrogram(cluster.tree)

    CairoPDF(file.name, width = 10, height = 10)
    heatmap.2(sites.plot, Rowv = TRUE, Colv = cluster.dendro, dendrogram = "both", scale = "none", trace = "none", labRow = FALSE)
    dev.off()
}

rnb.options("covariate.adjustment.columns" = c("Age", "Sex", "Ethnicity", "HBP.fix", "BMI"))
comp.cols <- "Recovered"
reg.types <- c("genes", "promoters")
diffmeth.adj <- rnb.execute.computeDiffMeth(rnb.norm, comp.cols, region.types = reg.types)

comparison <- get.comparisons(diffmeth.adj)
tab.sites <- get.table(diffmeth.adj, comparison, "sites", return.data.frame = TRUE) 
tab.promoters <- get.table(diffmeth.adj, comparison, "promoters", return.data.frame = TRUE)

promoters.recovered <- promoters.beta[,pheno(rnb.norm)$Recovered == TRUE] 
promoters.coef <- promoters.recovered - tab.promoters$mean.mean.g1
sites.cutoff <- sort(tab.sites$combinedRank)[1000]
promoters.cutoff <- sort(tab.promoters$combinedRank)[1000]
promoters.500 <- sort(tab.promoters$combinedRank)[500]

exportDMRs2regionFile(rnb.norm, diffmeth.adj, "./sites_top1000.bed", comparison, "sites", rank.cut = sites.cutoff)
exportDMRs2regionFile(rnb.norm, diffmeth.adj, "./promoters_top1000.bed", comparison, "promoters", rank.cut = promoters.cutoff)

promoters.key <- which(tab.promoters$combinedRank <= promoters.cutoff)
promoters.top1000 <- promoters.coef[promoters.key,]

promoters.500.key <- which(tab.promoters$combinedRank <= promoters.500)
promoters.top500 <- promoters.coef[promoters.500.key,]

#Generate anova heatmaps
gen.anova.heatmap <- function(file.name, dataset)
{ 
    CairoPDF(file.name, width = 10, height = 10)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 1.0, labRow = "", keysize = 0.9)
    dev.off()
}

gen.anova.heatmap("top1000_DMR", promoters.top1000)

thresholds <- c(0.01, 0.005, 0.001)
promoters.threshold <- map(thresholds, get.sizes, tab.promoters)
names(promoters.threshold) <- thresholds
threshold.df <- melt(promoters.threshold)
colnames(threshold.df) <- c("mysum", "Direction", "Test")
threshold.df$Test %<>% factor
levels(threshold.df$Test) <- c("0.01", "0.005", "0.001")

get.sizes <- function(p.val, dataset)
{
    dataset.sig <- dataset[dataset$comb.p.val < p.val,]
    dataset.up <- dataset.sig[dataset.sig$mean.mean.quot.log2 > 0,]
    dataset.down <- dataset.sig[dataset.sig$mean.mean.quot.log2 < 0,]
    return(list(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1])))
}

gen.decideplot("threshold_selection", threshold.df)
gen.decideplot <- function(file.name, decide.plot)
{
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = max(mysum) + 60, hjust = -1.1, label = mysum), position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = min(mysum) - 60, hjust = 1.5, label = abs(mysum)), position = position_dodge(width = 1))
    #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab("Differentially Methylated Promoters")
    CairoPDF(file.name, width = 6, height = 7)
    print(p)
    dev.off()
}

rownames(promoters.beta) <- rownames(annotation(rnb.norm, type = "promoters"))
saveRDS.gz(promoters.beta, "./save/promoters_beta.rda")
promoters.annotation <- annotation(rnb.norm, type = "promoters")
saveRDS.gz(promoters.annotation, "./save/promoters_annotation.rda")

joined.table <- data.frame(tab.promoters, promoters.annotation)

top.1001 <- which(tab.promoters$combinedRank <= promoters.cutoff)
annot.1000 <- annotation(rnb.norm, type = "promoters")[top.1000,] %>% filter(!is.na(symbol)) %>% filter(!str_detect(symbol, ";"))
annot.1000$symbol %<>% str_replace("-.*$", "")
#trouble.symbols <- filter(annot.1000, str_detect(symbol, "-"))$symbol
write.xlsx(annot.1000, "annot_1000.xlsx")
#rnb.run.tnt(rnb.norm)
#annotation(rnb.norm) %>% str
