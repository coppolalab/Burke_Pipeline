library(RnBeads)
#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
#library(limma)
#library(sva)
library(R.utils)
#library(Biobase)
#library(lumi)
#library(lumiHumanAll.db)
#library(annotate)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(lambda.r)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

gen.pcaplot <- function(filename, dataset, facet.bool, size.height, size.width)
{
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

gen.heatmap <- function(dataset, ME.genes)
{
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- filter(full.df, module.colors == index)
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- parLapply(cluster, enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- parLapply(cluster, names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

plot.eigencor <- function(module.color, ME.genes, col.name, pheno.vector)
{
    ME.eigengene <- ME.genes[[module.color]]
    ME.combined <- data.frame("Eigenpromoter" = ME.eigengene, pheno = pheno.vector)
    colnames(ME.combined)[2] <- col.name
    p <- ggplot(ME.combined, aes_string(x = col.name, y = "Eigenpromoter")) + geom_point(position = "jitter", color = module.color)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + stat_smooth(method = "lm", se = TRUE)

    filename <- paste(col.name, module.color, "eigengenes", sep = "_")
    CairoPDF(filename, height = 6, width = 7)
    print(p)
    dev.off()
}

rnb.final <- readRDS.gz(file = "../baseline_methylation/save/rnb.norm.rda")

promoters.beta <- meth(rnb.final, type = "promoters")# %>% lumi::beta2m()
rownames(promoters.beta) <- rownames(annotation(rnb.final, type = "promoters"))

model.continuous <- select(pheno(rnb.final), Age, BMI)
model.sex <- model.matrix( ~ 0 + pheno(rnb.final)$Sex )[,-1]
model.ethnicity <- model.matrix( ~ 0 + pheno(rnb.final)$Ethnicity )[,-1]
HBP.fix <- pheno(rnb.final)$HBP %>% droplevels
model.hbp <- model.matrix( ~ 0 + HBP.fix )[,-1]
model.final <- cbind(model.continuous, "Male" = model.sex, "Hispanic" = model.ethnicity, "HBP" = model.hbp)

#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
cleaned.beta <- removeBatchEffect(promoters.beta, covariates = model.final)
#lumi.cleaned <- lumi.final
#exprs(lumi.cleaned) <- export.expr
#boxplot(lumi.cleaned, subset = NULL, ylab = "Intensity", main = "Robust spline normalized VST signal intensity")
saveRDS.gz(cleaned.beta, file = "./save/cleaned.beta.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

meth.data.final <- t(cleaned.beta)

sft.final <- pickSoftThreshold(meth.data.final, powerVector = powers, verbose = 5, networkType = "signed")
sft.final.df <- sft.final$fitIndices
saveRDS.gz(sft.final, file = "./save/sft.final.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.final.df$multiplied <- sft.final.df$SFT.R.sq * -sign(sft.final.df$slope)
p <- ggplot(sft.final.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.final.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 16
adjacency.final <- adjacency(meth.data.final, power = softPower, type = "signed")
saveRDS.gz(adjacency.final, file = "./save/adjacency.final.rda")

TOM.final <- TOMsimilarity(adjacency.final, verbose = 5)
dissimilarity.TOM <- 1 - TOM.final
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./1-genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 50

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./1-gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(meth.data.final, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - cor(ME.genes)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")
CairoPDF(file = "./3-module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(meth.data.final, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("4-module_eigengene_clustering_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("5-eigengenes", height = 10, width = 18)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

all.degrees <- intramodularConnectivity(adjacency.final, module.colors)
fdata <- annotation(rnb.final, type = "promoters")
gene.info.join <- data.frame(Promoter = rownames(fdata), select(fdata, Chromosome:Strand))
gene.info <- mutate(gene.info.join, module.colors = module.colors, mean.count = apply(meth.data.final, 2, mean)) %>% data.frame(all.degrees)

write.csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.colors, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(cor(meth.data.final, ME.genes, use = "p"))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(meth.data.final)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Promoter:Strand), gene.module.membership, module.membership.pvalue)
write.csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Promoter <- rownames(gene.module.membership)
module.membership.pvalue$Promoter <- rownames(module.membership.pvalue)

gene.module.membership.long <- gather(gene.module.membership, module.comparison, correlation, violet:grey)
module.membership.pvalue.long <- gather(module.membership.pvalue, module.comparison, p.value, violet:grey)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Promoter, Chromosome:kscaled, module.comparison:p.value)
write.csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(meth.data.final), levels = rownames(meth.data.final))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
meth.data.plot <- data.frame(t(meth.data.final), module.colors)
cluster <- makeForkCluster(8)
split(meth.data.plot, meth.data.plot$module.colors) %>% parLapply(cl = cluster, gen.heatmap, ME.genes.plot)

modules.out <- select(gene.info, Promoter:module.colors)
gen.bedfile <- function(subset.df)
{
    color.name <- unique(subset.df$module.colors)
    write.table(select(subset.df, -Promoter, -module.colors), file = paste("./bedfiles/", color.name, ".bed", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)
}
by(modules.out, modules.out$module.colors, gen.bedfile)
write.xlsx(modules.out, "modules_out.xlsx")

#source("../../FRDA project/GO/enrichr.R")
#enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
##enrichr.submit("grey", modules.out, enrichr.terms, FALSE)
#color.names <- unique(module.colors) %>% sort
#trap1 <- l_ply(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

source("../../FRDA project/common_functions.R")
targets.final.known <- pheno(rnb.final)
colnames(targets.final.known)[1] <- "Sample.Name"
rownames(ME.genes) <- targets.final.known$Sample.Name

targets.age <- select(targets.final.known, Sample.Name, Age)
targets.age$Age %<>% as.numeric
cor.age <- gen.cor(ME.genes, targets.age)

targets.sex <- filter(targets.final.known, Sex != "UNKNOWN")
targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame 
colnames(targets.sex.m) <- c("Sex")
targets.sex.m %<>% mutate(Sample.Name = targets.final.known$Sample.Name)
cor.sex <- gen.cor(ME.genes, targets.sex.m)

targets.deltafm <- select(targets.final.known, Sample.Name, Delta.fm)
cor.deltafm <- gen.cor(ME.genes, targets.deltafm)

targets.fmadmit <- select(targets.final.known, Sample.Name, FM.admit)
cor.fmadmit <- gen.cor(ME.genes, targets.fmadmit)

targets.fmdischarge <- select(targets.final.known, Sample.Name, FM.discharge)
cor.fmdischarge <- gen.cor(ME.genes, targets.fmdischarge)

module.traits.all <- cbind(cor.deltafm, cor.fmadmit, cor.fmdischarge, cor.age, cor.sex)  %>% data.frame
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix
module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)

write_csv(module.trait.out, "module_trait_cor.csv")
text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)
gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")

#Plot interesting modules because none are currently significant
plot.eigencor("red", ME.genes, "DeltaFM", targets.final.known$Delta.fm)
plot.eigencor("orange", ME.genes, "DeltaFM", targets.final.known$Delta.fm)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

modules.join <- mutate(modules.out, Symbol = fdata$symbol)
orange.genes <- filter(modules.join, module.colors == "orange" & !is.na(Symbol) & !str_detect(Symbol, ":"))$Symbol %>% str_replace("-.*$", "")
red.genes <- filter(modules.join, module.colors == "red" & !is.na(Symbol) & !str_detect(Symbol, ":"))$Symbol %>% str_replace("-.*$", "")

gene.modules <- read.xlsx("../WGCNA/modules_out.xlsx")
pink.modules <- filter(gene.modules, module.colors == "pink" & !is.na(Symbol))$Symbol %>% str_replace("-.*$", "")
turquoise.modules <- filter(gene.modules, module.colors == "turquoise" & !is.na(Symbol))$Symbol %>% str_replace("-.*$", "")
paleturquoise.modules <- filter(gene.modules, module.colors == "paleturquoise" & !is.na(Symbol))$Symbol %>% str_replace("-.*$", "")
darkgreen.modules <- filter(gene.modules, module.colors == "darkgreen" & !is.na(Symbol))$Symbol %>% str_replace("-.*$", "")

compute.kappa <- function(methy.module, expr.module)
{
    methy.unique <- unique(methy.module)
    expr.unique <- unique(expr.module)

    shared <- which(methy.unique %in% expr.unique) %>% length
    methy.only <- which(!(methy.unique %in% expr.unique)) %>% length
    expr.only <- which(!(expr.unique %in% methy.unique)) %>%length
    
    o.num <- shared / (shared + methy.only + expr.only)
    a.num <- ((length(methy.unique) * length(expr.unique)) + ((length(methy.unique) - shared) * (length(expr.unique) - shared))) / ((shared + methy.only + expr.only) ^ 2)
    k.num <- (o.num - a.num) / (1 - a.num)

    return(shared)
}

compute.kappa(red.genes, pink.modules)
compute.kappa(red.genes, turquoise.modules)
compute.kappa(orange.genes, turquoise.modules)
compute.kappa(orange.genes, paleturquoise.modules)
