# Recommended Libraries 
library(GEOquery)
library(hgu133plus2.db)
library(dplyr)
library(stringr)
library(affy)
library(limma)
library(genefilter)
library(topGO)

# Reading in .CEL Files 
# Assumes the .CEL files and the .soft files are in a folder called data.
gse<- getGEO(filename='data/GSE50697_family.soft')
treatment <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
pd <- data.frame(treatment=as.factor(sapply(GSMList(gse), treatment)))
pd$treatment <- as.factor(pd$treatment)
levels(pd$treatment) <- c("Control", "miR")
celfiles <- paste0('data/', rownames(pd), '.CEL')

# Pre-processing 
affydata <- read.affybatch(celfiles, phenoData = new("AnnotatedDataFrame", pd))
phenoData(affydata)
eset <- rma(affydata)
pData(eset)

# Annotations 
anno_eset <- AnnotationDbi::select(hgu133plus2.db, keys = featureNames(eset), columns = c("SYMBOL", 
                                    "GENENAME"), keytype="PROBEID")
anno_eset <- subset(anno_eset, !is.na(SYMBOL))
anno_grouped <- group_by(anno_eset, PROBEID)
anno_summarized <- dplyr::summarise(anno_grouped, no_of_matches = n_distinct(SYMBOL))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
probe_stats <- anno_filtered 
ids_to_exculde <- (featureNames(eset) %in% probe_stats$PROBEID)
eset_final <- subset(eset, !ids_to_exculde)
fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- left_join(fData(eset_final), anno_eset)
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID

# Linear Models 
model <- model.matrix(~0 + eset_final$treatment)
colnames(model) <- levels(eset_final$treatment)
contrasts <- makeContrasts(miR - Control, levels=model)
fit <- lmFit(eset_final, model)
fitted.contrast <- contrasts.fit(fit, contrasts)
fitted.ebayes <- eBayes(fitted.contrast)
table_genes <- topTable(fitted.ebayes, number=Inf)
table_genes <- table_genes[complete.cases(table_genes),]

# Volcano Plot 
top_genes2 <- topTable(fitted.ebayes, number=Inf, lfc=2)
volcanoplot(fitted.ebayes, xlab = expression("log"[2]*" Fold Change"), 
            ylab = expression("-log"[10]*"(p-value)"))
points(top_genes2[['logFC']],-log10(top_genes2[['P.Value']]),col='red')

# Gene Ontology 
genes_after_fdr <- subset(table_genes, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefilter::genefinder(eset_final, as.character(genes_after_fdr), 
                                         method = "manhattan", scale="none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(eset_final)[back_genes_idx]
back_genes <- setdiff(back_genes, genes_after_fdr)
plot(density(table_genes[,"AveExpr"]))
lines(density(table_genes[genes_after_fdr, "AveExpr"]))
lines(density(table_genes[rownames(table_genes) %in% back_genes, "AveExpr"]))
gene_ids <- rownames(table_genes)
in_universe <- gene_ids %in% c(genes_after_fdr, back_genes)
in_selection <- gene_ids %in% genes_after_fdr
all_genes <- in_selection[in_universe]
all_genes <- as.factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_ids[in_universe]
top_GO_data <- new("topGOdata", ontology = "BP", allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db")
result_top_GO_elim <- 
  runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- 
  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")
res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hgu133plus2.db", geneCutOff = 1000)

helper <- function(x) {
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
}
res_top_GO$sig_genes <- sapply(genes_top_GO, helper) 
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
               useInfo = 'def')

