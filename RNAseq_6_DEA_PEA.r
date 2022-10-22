#####################################################
#################### DE analysis ####################

####################### DEseq2 ######################

library(DESeq2)

cts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, check.names=FALSE)
cts <- round(cts, digits = 0)
coldata <- read.table("sample_info.txt", header = TRUE, row.names = 1)

coldata$condition <- factor(coldata$condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds)
res <- results(dds)
res
write.table(as.data.frame(res),"DE_results.txt", sep = "\t")

vsd <- vst(dds)
pdf(file="PCA.pdf", height = 10, width = 10)
plot <- plotPCA(vsd, intgroup=c("condition"))
plot + theme_bw()
dev.off()

####################### edgeR #######################

library(edgeR)

cts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, check.names=FALSE)
cts <- round(cts, digits = 0)
coldata <- read.table("sample_info.txt", header = TRUE, row.names = 1)

y <- DGEList(counts=cts, group=coldata$condition) 

y$samples$group <- relevel(y$samples$group, ref="control")

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~y$samples$group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

tab <- topTags(qlf, n=Inf)
write.table(tab, file="DE_results.txt", sep = "\t")

################################################################
####################### Pathway analysis #######################

BiocManager::install('grimbough/biomaRt')
library(biomaRt)
packageVersion('biomaRt')
require(graph)
require(ROntoTools)

res <- tab$table

ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensg <- rownames(res)

ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
rownames(res) <- ensg.no_version

BM <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values=ensg.no_version, mart=ensembl)
BM <- na.omit(BM)	
BM <- BM[!duplicated(BM$ensembl_gene_id), ]
rownames(BM) <- BM$ensembl_gene_id
res <- res[rownames(res) %in% rownames(BM), ]
res <- unique(res)

res$entrezgene_id <- BM$entrezgene_id[match(rownames(res), rownames(BM))]
res$hgnc_symbol <- BM$hgnc_symbol[match(rownames(res), rownames(BM))]

res$entrezgene_id <- paste("hsa", res$entrezgene_id, sep=":")
res$entrezgene_id <- as.factor(res$entrezgene_id)

fc <- res$logFC[res$FDR < 0.01]
names(fc) <- res$entrezgene_id[res$FDR < 0.01]

pv <- res$PValue[res$FDR < 0.01]
names(pv) <- res$entrezgene_id[res$FDR < 0.01]

ref <- as.character(res$entrezgene_id)

kpg <- keggPathwayGraphs("hsa", verbose = FALSE)
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
edgeWeightByType = list(activation = 1, inhibition = -1,
expression = 1, repression = -1),
defaultWeight = 0)
kpn <- keggPathwayNames("hsa")
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)

peRes <- pe(x = fc, graphs = kpg, ref = ref, nboot = 200, verbose = FALSE)

pa_summary <- Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
        pAcc = TRUE, pORA = TRUE, comb.pv =  c("pAcc", "pORA"), comb.pv.func = compute.normalInv, order.by = "pComb")

write.csv(as.data.frame(pa_summary), 
          file="PA_results.csv")

pdf(file="PA.pdf", height = 3, width = 5)
plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = 0.01)
dev.off()
