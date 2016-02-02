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
library(Biobase)
library(marray)
library(limma)
library(MASS)
library(matrixStats)
library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(R.utils)

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

#Heatmap function
gen.heatmap <- function(filename, lumi.object)
{
    intensities1.cor <- corFast(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), scale = "none", cexCol = 1.0, cexRow = 1.0)
    dev.off()
}

standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes
gen.topgenes <- function(filename, lumi.object, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    dataset <- exprs(lumi.object)
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- flashClust::hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), scale = "none", cexRow = 1.0, cexCol = 1.0)
    dev.off()
    return(top.genes.dist)
}

gen.colors <- function(diagnosis.colorscheme, targetset)
{
    diagnosis.heatmap <- data.frame("Status" = unique(as.integer(targetset$Recovered)), "Diagnosis.Color" = diagnosis.colorscheme)
    print(diagnosis.heatmap)
    colorscheme <- data.frame("Status" = targetset$Recovered) %>% join(diagnosis.heatmap)
    colorscheme <- select(colorscheme, Diagnosis.Color)
    return(colorscheme)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$External.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) 
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Run statistical cutoff tests
gen.decide <- function(test, fit.object, write.results)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] #Eliminate the row for no change in expression
    mysum[1] <- -(mysum[1])
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("negative", "positive"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot)
{
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = max(mysum) + 60, hjust = -1.1, label = mysum), position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Test, y = mysum), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Test, y = mysum, ymax = min(mysum) - 60, hjust = 1.5, label = abs(mysum)), position = position_dodge(width = 1))
    #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab("Differentially Expressed Genes")
    CairoPDF(filename, width = 6, height = 7)
    print(p)
    dev.off()
}

gen.pval.hist <- function(filename, fit.pvals)
{
    p <- ggplot(data.frame(fit.pvals), aes(x = recovered_vs_unrecovered)) + geom_histogram(binwidth = 1/80) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("P-value distribution") + theme(axis.title.y = element_blank()) + xlab("p-value")
    CairoPDF(filename, height = 7, width = 7)
    print(p)
    dev.off()
}

gen.venndiagram <- function(filename, results)
{
    sumV <- colSums(summary(results)[-2,])
    v <- paste(c("Carrier vs. Control", "Patient vs. Control", "Patient vs. Carrier"), " (", sumV, ")", sep = "")
    CairoPDF(filename, width = 6, height = 6)
    vennDiagram(results, names = v, main = "", include = c("up", "down"), counts.col=c(2,3), cex = 0.8)
    dev.off()
}

#Peer analysis
gen.peer <- function(num.factors, intensities, use.covariates, covariates)
{
    model = PEER() #Instantiate PEER object
    PEER_setNk(model, num.factors) #Specify the number of factors to find
    PEER_setPhenoMean(model, as.matrix(t(intensities))) #Provide mean values for each subject
    PEER_setAdd_mean(model, TRUE) #Enable use of mean in linear model
    if (use.covariates == TRUE)
    {
        PEER_setCovariates(model, as.matrix(covariates)) #Add additional covariates
    }
    PEER_setNmax_iterations(model, 1000) #Specify the maximum number of iterations
    PEER_update(model) #Calculate factors
    residuals.PEER = t(PEER_getResiduals(model))
    rownames(residuals.PEER) = rownames(intensities)
    colnames(residuals.PEER) = colnames(intensities)

    write.csv(data.frame(residuals.PEER), file = paste("residuals_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getX(model)), file = paste("factor_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getW(model)), file = paste("weight_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getAlpha(model)), file = paste("precision_", num.factors, sep = "", ".csv"), row.names = FALSE)

    CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    PEER_plotModel(model)
    dev.off()

    CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
    dev.off()
}

#Calculate ratios.  Really needs work!
gen.ratios <- function(lumi.object)
{
    all.unrecovered <- exprs(lumi.object[,!lumi.object$Recovered])
    all.unrecovered.means <- rowMeans(all.unrecovered)
    all.recovered <- exprs(lumi.object[,lumi.object$Recovered])

    all.coefficients <- all.recovered - all.unrecovered.means
    all.coefficients.label <- data.frame("nuID" = rownames(all.coefficients), all.coefficients)
    all.samples <- data.frame("nuID" = rownames(exprs(lumi.object)), exprs(lumi.object))
    colnames(all.samples)[2:length(all.samples)] <- paste(colnames(all.samples[2:length(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(all.coefficients.label, all.samples)
    return(ratio.exp)
}

#Generate fit object
gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(recovered_vs_unrecovered = Recovered - Not.Recovered, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef") %>% which
    dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:3, widths = "auto")
    #setColWidths(wb, 1, cols = 2, widths = 15)
    setColWidths(wb, 1, cols = 4, widths = 45)
    setColWidths(wb, 1, cols = 5:7, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("nuID" = featureNames(lumi.object), dataset)
    colnames(treat.de)[str_detect(colnames(treat.de), "Genes")] %<>% str_replace("Genes\\.", "") %>% tolower %>% capitalize
    
    fitsel.ratio.all <- merge(treat.de, ratio.exp)
    fitsel.ratio.all$RefSeq <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "REFSEQ") %>% llply(paste, collapse = ",") %>% reduce(c)
    fitsel.ratio.all$UniGene <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "UNIGENE") %>% llply(paste, collapse = ",") %>% reduce(c)
    fitsel.ratio.all$EntrezID <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "ENTREZID") %>% llply(paste, collapse = ",") %>% reduce(c)
    fitsel.return.all <- select(fitsel.ratio.all, nuID, Accession, Symbol, Definition, Coef, p.value, recovered_vs_unrecovered, t, A, RefSeq:EntrezID, matches("EPS")) %>% arrange(p.value)

    anovalist <- map_lgl(fitsel.return.all$recovered_vs_unrecovered, any) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, "_time1.xlsx", sep = ""))

    write_csv(fitsel.return.all, path = paste("./complete_genelist_time1_", suffix, ".csv", sep = ""))
    return(fitsel.return)
}

gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 3, widths = 45)
    setColWidths(wb, 1, cols = 4:6, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}
    
#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 10, height = 10)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 1.0, labRow = "", keysize = 0.9)
    dev.off()
}

enrichr.wkbk <- function(database, full.df, colname, subdir, direction)
{
    dataset <- full.df[[database]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    dir.create(file.path("./enrichr", direction), recursive = TRUE)
    filename = paste(file.path("./enrichr", direction, database), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

source("../../FRDA project/common_functions.R")

targets.read <- read.xlsx("../phenotypedata/2015-9186_Sample Key.xlsx")
targets.read$sampleID <- paste(targets.read$general.array, targets.read$genexstripe.controling.stripe, sep = "_")
targets.read$External.ID %<>% str_replace(" ", "_") 
targets.reduce <- select(targets.read, External.ID, sampleID)

rna.read <- read.xlsx("../phenotypedata/dan_burke_rna_20151207T184346.xlsx")
rna.filter <- filter(rna.read, !is.na(RIN))
duped.ids <- filter(rna.filter, Sample.Num == 2)$PIDN %>% paste(collapse = "|")
duped.rows <- filter(rna.filter, grepl(duped.ids, PIDN))
drop.rnas <- c("16622", "16626", "16653") %>% paste(collapse = "|")
rna.filter2 <- filter(rna.filter, !grepl(drop.rnas, RNA.ID))
rna.filter2$PIDN %<>% str_replace("EPS", "EPS_")
rna.reduce <- select(rna.filter2, RNA.ID, PIDN, RIN)

pheno.read <- read.xlsx("../phenotypedata/Dana pheno filter.xlsx")
colnames(pheno.read) %<>% str_replace("@.", "")
pheno.fm <- filter(pheno.read, !is.na(FUGL.MEYER.ADMISSION) & !is.na(FUGL.MEYER.DISCHARGE))
write.xlsx(pheno.fm, "pheno.fm.xlsx")
pheno.select <- select(pheno.fm, ID:GENDER, FUGL.MEYER.ADMISSION, FUGL.MEYER.DISCHARGE, ETHNICITY, high.blood.pressure, BODY.MASS.INDEX)
colnames(pheno.select) <- c("External.ID", "Age", "Sex", "FM.admit", "FM.discharge", "Ethnicity", "HBP", "BMI")
pheno.select$Delta.fm <- pheno.select$FM.discharge - pheno.select$FM.admit
pheno.select$Recovered <- (pheno.select$Delta.fm > 10)
pheno.select$External.ID %<>% sprintf(fmt = "EPS_%03d")

pheno.merge <- merge(pheno.select, targets.reduce) %>% merge(rna.reduce, by.x = "External.ID", by.y = "PIDN")

lumi.read <- lumiR("../raw_expression/2015-9186.tsv", lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, convertNuID = TRUE)
sample.key <- paste(pheno.merge$sampleID, collapse = "|")
sample.index <- grepl(sample.key, sampleNames(lumi.read))
lumi.known <- lumi.read[,sample.index]
pData(lumi.known) <- pheno.merge
sampleNames(lumi.known) <- lumi.known$External.ID
lumi.vst <- lumiT(lumi.known)
lumi.norm <- lumiN(lumi.vst, method = "rsn")
lumi.qual <- lumiQ(lumi.norm, detectionTh = 0.01)

CairoPDF("sample_clustering", height = 6, width = 6)
plot(lumi.qual, what = "sampleRelation")
dev.off()

lumi.cutoff <- detectionCall(lumi.qual, Th = 0.05)
lumi.expr <- lumi.qual[which(lumi.cutoff > 0),]
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated
saveRDS.gz(lumi.expr.annot, file = "./save/lumi.expr.annot.rda")

qcsum <- lumi.expr.annot@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum) %<>% str_replace("\\.0\\.01\\.", "")
qcsum$Sample.Name <- rownames(qcsum)
qcsum %<>% arrange(desc(distance.to.sample.mean))
#qcsum$RIN <- lumi.expr.annot$RIN

remove.samples <- qcsum[1,]$Sample.Name
outlier.index <- !(grepl(remove.samples, sampleNames(lumi.known))) 
lumi.rmout <- lumi.known[,outlier.index]

lumi.rmout.vst <- lumiT(lumi.rmout)
lumi.rmout.norm <- lumiN(lumi.rmout.vst, method = "rsn")
lumi.rmout.qual <- lumiQ(lumi.rmout.norm, detectionTh = 0.01)

CairoPDF("sample_clustering_rmout", height = 6, width = 6)
plot(lumi.rmout.qual, what = "sampleRelation")
dev.off()

lumi.rmout.cutoff <- detectionCall(lumi.rmout.qual, Th = 0.05)
lumi.rmout.expr <- lumi.rmout.qual[which(lumi.rmout.cutoff > 0),]
symbols.rmout <- getSYMBOL(rownames(lumi.rmout.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.final <- lumi.rmout.expr[!symbols.rmout,] #Drop any probe which is not annotated
lumi.final$Recovered <- !lumi.final$Recovered
saveRDS.gz(lumi.final, file = "./save/lumi.final.rda")

#Much simpler box plot code - probably not an option for complex data
CairoJPEG(filename = "baseline_intensity_norm.jpg", width = 3840, height = 2160)
par(cex.axis = 1.9, cex.lab = 1.9, cex.main = 1.9)
boxplot(lumi.final, subset = NULL, ylab = "Intensity", main = "Robust spline normalized VST signal intensity")
dev.off()

outcome.colors <- c("red", "blue")
color.vector <- gen.colors(outcome.colors, pData(lumi.final)) %>% unlist %>% as.character
gen.heatmap("baseline_heatmap_norm", lumi.final)

#Top 500 and 1000 genes
lumi.mads <- apply(exprs(lumi.final), 1, mad)
lumi.ordered <- order(lumi.mads, decreasing = TRUE)

top500.dist <- gen.topgenes("baseline_heatmap_500", lumi.final, lumi.ordered, 500)
top1000.dist <- gen.topgenes("baseline_heatmap_1000", lumi.final, lumi.ordered, 1000)

#Figure out later
#plotSampleRelation(lumi.final, subset = NULL, method = "mds", color = color.vector, pch = 20, main = "", addLegend = TRUE)

#Principal components analysis
cm1 <- cmdscale(top1000.dist, eig = TRUE)
gen.pca("baseline_mds_status", cm1, phenoData(lumi.final), outcome.colors, "Recovered")

model.continuous <- select(pData(lumi.final), Age, BMI, RIN, Recovered)
model.continuous$RIN <- model.continuous$RIN - min(model.continuous$RIN)
model.sex <- model.matrix( ~ 0 + lumi.final$Sex )[,-1]
model.ethnicity <- model.matrix( ~ 0 + lumi.final$Ethnicity )[,-1]
model.hbp <- model.matrix( ~ 0 + lumi.final$HBP )[,-1]
model.final <- cbind(model.continuous, "Male" = model.sex, "Hispanic" = model.ethnicity, "HBP" = model.hbp)

#PEER was run but deemed unnecessary
#Run PEER analysis and correlate to known covariates
#gen.peer(5, exprs(lumi.final), TRUE, model.final)

#PEER.weights <- read_csv("./weight_5.csv") %>% select(-(X1:X8))
#PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame
#PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
#colnames(PEER.weights.sums)[1] <- "Weight"

#p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
#p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
#CairoPDF("./PEER_weights", height = 4, width = 6)
#print(p)
#dev.off()

#Calculate ratios for use in tables
ratio.exp <- gen.ratios(lumi.final)
saveRDS.gz(ratio.exp, file = "./save/ratio.exp.rda")

model.df <- data.frame(model.final)
#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
model.cov <- select(model.df, -Recovered)
model.cov$RIN <- model.cov$RIN - min(model.cov$RIN)
model.design <- model.matrix( ~ 0 + lumi.final$Recovered )
colnames(model.design) <- c("Not.Recovered", "Recovered")
export.expr <- removeBatchEffect(exprs(lumi.final), covariates = model.cov, design = model.design)
export.lumi <- lumi.final

exprs(export.lumi) <- export.expr
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")
CairoJPEG(filename = "baseline_intensity_corrected.jpg", width = 3840, height = 2160)
par(cex.axis = 1.9, cex.lab = 1.9, cex.main = 1.9)
boxplot(export.lumi, subset = NULL, ylab = "Intensity", main = "Robust spline normalized VST signal intensity")
dev.off()

#Linear model fitting
model.full <- cbind(model.design, model.cov)
saveRDS.gz(model.full, file = "./save/model.full.rda")
fit.object <- gen.fit(lumi.final, model.full)
saveRDS.gz(fit.object, file = "./save/fit.object.rda")

#Generate statisical cutoff
decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01))
decide.plot <- ldply(decide, gen.decide, fit.object, FALSE) 
gen.decideplot("./threshold_selection", decide.plot)
decide.final <- gen.decide(c("none", 0.005), fit.object, TRUE) 

#Make tables
de.object <- read_tsv("./fit_none.tsv")
fit.selection <- gen.tables(de.object, lumi.final, ratio.exp, "pLess005")
saveRDS.gz(fit.selection, file = "./save/fit.selection.rda")

#P-value histogram
gen.pval.hist("./hist_pvalue", fit.object$p.value)

#Anova heatmaps
fit.selection.plot <- select(fit.selection, -contains("expr")) %>% select(contains("EPS"))
gen.anova.heatmap("./anova_heatmap", fit.selection.plot, "Recovered vs Unrecovered")

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#Submit genes to Enrichr
source('../../FRDA project/GO/enrichr.R')
enrichr.nofdr <- select(fit.selection, nuID, Symbol, recovered_vs_unrecovered)
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 

enrichr.up <- filter(enrichr.nofdr, recovered_vs_unrecovered == 1)
write.xlsx(enrichr.up, "./enrichr/up.xlsx")
enrichr.down <- filter(enrichr.nofdr, recovered_vs_unrecovered == -1)
write.xlsx(enrichr.down, "./enrichr/down.xlsx")

up.data <- lapply(enrichr.terms, get.enrichrdata, enrichr.up, FALSE)
down.data <- lapply(enrichr.terms, get.enrichrdata, enrichr.down, FALSE)
up.names <- enrichr.terms[!is.na(up.data)]
up.data <- up.data[!is.na(up.data)]
down.names <- enrichr.terms[!is.na(down.data)]
down.data <- down.data[!is.na(down.data)]

names(up.data) <- up.names
names(down.data) <- down.names

lapply(names(up.data), enrichr.wkbk, up.data, colname, subdir, "up")
lapply(names(down.data), enrichr.wkbk, down.data, colname, subdir, "down")


de.sort <- arrange(de.object, p.value)[1:1000,]
de.symbol <- unique(de.sort$Genes.SYMBOL)
write(de.symbol, "de.symbol.txt")

annot.top <- read.xlsx("../baseline_methylation/annot_1000.xlsx")
shared.genes <- de.symbol[de.symbol %in% annot.top$symbol]
write(shared.genes, "shared.genes.txt")
