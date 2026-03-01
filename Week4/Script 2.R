############################################
# 1️⃣ Install & Load Required Packages
############################################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "annotate",
                       "hgu133plus2.db", "pheatmap"))

install.packages(c("ggplot2"))

library(GEOquery)
library(limma)
library(annotate)
library(hgu133plus2.db)  # Annotation package (check platform below)
library(ggplot2)
library(pheatmap)

############################################
# 2️⃣ Download GEO Dataset
############################################

gse <- getGEO("GSE28735", GSEMatrix = TRUE)
gset <- gse[[1]]

expr <- exprs(gset)
pdata <- pData(gset)

############################################
# 3️⃣ Define Groups (Tumor vs Normal)
############################################

# Inspect metadata column names
colnames(pdata)

# Adjust column name if necessary (usually "source_name_ch1")
group <- ifelse(grepl("tumor", pdata$source_name_ch1, ignore.case=TRUE),
                "Tumor", "Normal")

group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

############################################
# 4️⃣ Differential Expression (limma)
############################################

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(Tumor - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="fdr", number=Inf)

############################################
# 5️⃣ Add Gene Annotation (DNA Gene Names)
############################################

# Map probe IDs to Gene Symbols
results$GeneSymbol <- getSYMBOL(rownames(results), "hgu133plus2.db")

# Remove probes without gene symbols
results <- results[!is.na(results$GeneSymbol), ]

############################################
# 6️⃣ Volcano Plot
############################################

results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=threshold)) +
    geom_point(alpha=0.6) +
    theme_minimal() +
    scale_color_manual(values=c("grey","red")) +
    labs(title="Volcano Plot: Tumor vs Normal",
         x="Log2 Fold Change",
         y="-Log10 Adjusted P-value")

############################################
# 7️⃣ Heatmap of Top DEGs
############################################

top_genes <- head(results[order(results$adj.P.Val), ], 50)

heatmap_data <- expr[rownames(top_genes), ]

# Replace probe IDs with gene symbols
rownames(heatmap_data) <- top_genes$GeneSymbol

pheatmap(heatmap_data,
         scale="row",
         annotation_col=data.frame(Group=group),
         show_rownames=TRUE,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")

############################################
# 8️⃣ Save DEG Results
############################################

write.csv(results, "GSE28735_DEGs_Tumor_vs_Normal.csv")
