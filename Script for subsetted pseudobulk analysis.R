#Doing UMAP analysis of the subsetted data

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@First Trimester UMAP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Reading datasets
placenta_early_3.data <- Read10X(data.dir = 'D:/Brain Placenta OCEC/Deposited placental datasets/GSM7882481_placenta_10X_early3.filtered_feature_bc_matrix/')
placenta_early_4.data <- Read10X(data.dir = 'D:/Brain Placenta OCEC/Deposited placental datasets/GSM7882482_placenta_10X_early4.filtered_feature_bc_matrix/')
placenta_early_5.data <- Read10X(data.dir = 'D:/Brain Placenta OCEC/Deposited placental datasets/GSM7882483_placenta_10X_early5.filtered_feature_bc_matrix/')
placenta_early_6.data <- Read10X(data.dir = 'D:/Brain Placenta OCEC/Deposited placental datasets/GSM7882484_placenta_10X_early6.filtered_feature_bc_matrix/')
  
#Making Seurat Objects
placenta_early_3 <- CreateSeuratObject(counts = placenta_early_3.data, project = 'First_Trim.3', min.cells = 3, min.features = 100)
placenta_early_4 <- CreateSeuratObject(counts = placenta_early_4.data, project = 'First_Trim.4', min.cells = 3, min.features = 100)
placenta_early_5 <- CreateSeuratObject(counts = placenta_early_5.data, project = 'First_Trim.5', min.cells = 3, min.features = 100)
placenta_early_6 <- CreateSeuratObject(counts = placenta_early_6.data, project = 'First_Trim.6', min.cells = 3, min.features = 100)

#Merging the six objects
first_trim <- merge(x = placenta_early_3, y = c(placenta_early_4, placenta_early_5, placenta_early_6))
first_trim[['RNA']] <- JoinLayers(first_trim[['RNA']])

#Filtering the dataset
first_trim[['percent.mt']] <- PercentageFeatureSet(first_trim, pattern = '^MT-')
VlnPlot(first_trim, features = 'percent.mt')
first_trim <- subset(first_trim, subset = percent.mt < 2)

#Preprocessing the data
first_trim <- NormalizeData(first_trim, normalization.method = 'LogNormalize', scale.factor = 10000)
first_trim <- FindVariableFeatures(first_trim, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(first_trim), 10)
LabelPoints(plot = VariableFeaturePlot(first_trim), points = top10, repel = T)
all.genes <- rownames(first_trim)
first_trim <- ScaleData(first_trim, features = all.genes)
first_trim <- RunPCA(first_trim, features = VariableFeatures(object = first_trim))
ElbowPlot(first_trim)

#Integrated Layers
first_trim <- Split_Layers(first_trim, assay = 'RNA', split.by = 'orig.ident')
first_trim <- IntegrateLayers(first_trim, method = CCAIntegration, orig.reduction = 'pca', new.reduction = 'integrated.cca')
first_trim[['RNA']] <- JoinLayers(first_trim[['RNA']])

#Running dimensional reduction
first_trim <- FindNeighbors(first_trim, reduction = 'integrated.cca', dims = 1:15)
first_trim <- FindClusters(first_trim, resolution = 0.5)
first_trim <- RunUMAP(first_trim, reduction = 'integrated.cca', dims = 1:15)
p1 <- DimPlot(first_trim, reduction = 'umap', label = T, label.size = 4) + NoLegend()
p2 <- DimPlot(first_trim, reduction = 'umap', group.by = 'orig.ident')

p3 <- DotPlot(first_trim, features = c('ITGA2', 'ITGA6', 'GATA3', 'MKI67', 'BCAM', 'MFSD2A', 'GREM2', 'ERVW-1', 'ERVFRD-1', 'SPINT2', 'GCM1', 'SDC1', 'CGB2',
                                 'CGB3', 'PAPPA', 'NOTCH1', 'NOTCH2', 'HLA-G', 'ITGB1', 'VCAM1', 'PDGFRB', 'PECAM1', 'PTPRC')) + theme(axis.text.x = 
                                                                                                                                      element_text(angle = 45, face = 'italic', vjust = 1, hjust = 1, size = 20))
#Labelling the Clusters
first_trim <- RenameIdents(object = first_trim,
                         '0' = 'ST',
                         '1' = 'pCTs',
                         '2' = 'ST',
                         '3' = 'ST',
                         '4' = 'ST',
                         '5' = 'Bi. pCTs',
                         '6' = 'Prolif. pCTs',
                         '7' = 'Stromal Cells',
                         '8' = 'EVTs',
                         '9' = 'HSCs',
                         '10' = 'Early ST',
                         '11' = 'Stromal Cells',
                         '12' = 'Fus. Comp. pCTs', 
                         '13' = 'ECs',
                         '14' = 'Hofbauer Cells',
                         '15' = 'ST')
first_trim$cluster.ident <- Idents(first_trim)

feature.cols <- colorRampPalette(c('#EBDDFB', '#C04EBD', '#8D308C'))
vln.cols <- c('ST' = '#FFA5D6', 'Early ST' = '#CED1F8', 'Fus. Comp. pCTs' = '#CDE9DC', 'Hofbauer Cells' = '#B6D8f7')

VlnPlot(first_trim, features = 'CLDN4', cols = vln.cols ) + NoLegend()

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Term UMAP@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Reading Datasets
placenta_term_3.data <- Read10X(data.dir = "J:/Juan/Pseudobulk for Protemoics Paper/Term placenta/GSM7882493_placenta_10X_late3.filtered_feature_bc_matrix")
placenta_term_4.data <- Read10X(data.dir = "J:/Juan/Pseudobulk for Protemoics Paper/Term placenta/GSM7882494_placenta_10X_late4.filtered_feature_bc_matrix")
placenta_term_5.data <- Read10X(data.dir = "J:/Juan/Pseudobulk for Protemoics Paper/Term placenta/GSM7882495_placenta_10X_late5.filtered_feature_bc_matrix")
placenta_term_6.data <- Read10X(data.dir = "J:/Juan/Pseudobulk for Protemoics Paper/Term placenta/GSM7882496_placenta_10X_late6.filtered_feature_bc_matrix")

#Making Seurat objects
placenta_term_3 <- CreateSeuratObject(counts = placenta_term_3.data, project = 'Term_3', min.cells = 3, min.features = 100)
placenta_term_4 <- CreateSeuratObject(counts = placenta_term_4.data, project = 'Term_4', min.cells = 3, min.features = 100)
placenta_term_5 <- CreateSeuratObject(counts = placenta_term_5.data, project = 'Term_5', min.cells = 3, min.features = 100)
placenta_term_6 <- CreateSeuratObject(counts = placenta_term_6.data, project = 'Term_6', min.cells = 3, min.features = 100)

#Merging the six objects
term <- merge(x = placenta_term_3, y = c(placenta_term_4, placenta_term_5, placenta_term_6))

#Filtering the dataset
term[['percent.mt']] <- PercentageFeatureSet(term, pattern = '^MT-')
VlnPlot(term, features = 'percent.mt')
term <- subset(term, subset = percent.mt < 2)

#Preprocessing the data
term <- NormalizeData(term, normalization.method = 'LogNormalize', scale.factor = 10000)
term <- FindVariableFeatures(term, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(term), 10)
LabelPoints(plot = VariableFeaturePlot(term), points = top10, repel = T)
all.genes <- rownames(term)
term <- ScaleData(term, features = all.genes)
term <- RunPCA(term, features = VariableFeatures(object = term))
ElbowPlot(term)

#Integrated Layers
term <- IntegrateLayers(term, method = CCAIntegration, orig.reduction = 'pca', new.reduction = 'integrated.cca')
term[['RNA']] <- JoinLayers(term[['RNA']])

#Running dimensional reduction
term <- FindNeighbors(term, reduction = 'integrated.cca', dims = 1:15)
term <- FindClusters(term, resolution = 0.5)
term <- RunUMAP(term, reduction = 'integrated.cca', dims = 1:15)
p1 <- DimPlot(term, reduction = 'umap', label = F) + NoLegend()
p2 <- DimPlot(term, reduction = 'umap', group.by = 'orig.ident', shuffle = T)

p3 <- DotPlot(term, features = c('ITGA2', 'ITGA6', 'GATA3', 'MKI67', 'BCAM', 'MFSD2A', 'GREM2', 'ERVW-1', 'ERVFRD-1', 'SPINT2', 'GCM1', 'SDC1', 'CGB2',
                                       'CGB3', 'PAPPA', 'NOTCH1', 'NOTCH2', 'HLA-G', 'ITGB1', 'VCAM1', 'PDGFRB', 'PECAM1', 'PTPRC')) + theme(axis.text.x = 
                                                                                                                                               element_text(angle = 45, face = 'italic', vjust = 1, hjust = 1, size = 20))
LabelClusters(p1, id = 'ident', fontface = 'bold', box = F, repel = T,  size =6)

#Labelling the Clusters
term <- RenameIdents(object = term,
                           '0' = 'ST',
                           '1' = 'ST',
                           '2' = 'ST',
                           '3' = 'ST',
                           '4' = 'ST',
                           '5' = 'SMCs',
                           '6' = 'pCTs',
                           '7' = 'ST',
                           '8' = 'ST',
                           '9' = 'ECs',
                           '10' = 'Immune Cells',
                           '11' = 'Stromal Cells')
term$cluster.ident <- Idents(term)

feature.cols <- colorRampPalette(c('#EBDDFB', '#C04EBD', '#8D308C'))
FeaturePlot_scCustom(first, features = 'TNFRSF1B', colors_use = feature.cols(3)) + theme(plot.title = element_text(size = 20))
vln.cols <- c('ST' = '#FFA5D6', 'Early ST' = '#CED1F8', 'Fus. Comp. pCTs' = '#CDE9DC', 'Hofbauer Cells' = '#B6D8f7')

VlnPlot(first, features = 'TNFRSF1B', cols = vln.cols) + NoLegend() + theme(plot.title = element_text(size = 20), axis.text.x = 
                                                                              element_text(size = 15), axis.title.x = element_text(size =20),
                                                                            axis.text.y = element_text(size = 15), 
                                                                            axis.title.y = element_text(size = 20))


#Adding labels identifiers to the term and first trimester data
placenta_merged$origin <- placenta_merged$orig.ident
placenta_merged$origin[placenta_merged$orig.ident  == 'Term_3'] <- 'Term'
placenta_merged$origin[placenta_merged$orig.ident  == 'Term_4'] <- 'Term'
placenta_merged$origin[placenta_merged$orig.ident  == 'Term_5'] <- 'Term'
placenta_merged$origin[placenta_merged$orig.ident  == 'Term_6'] <- 'Term'
placenta_merged$origin[placenta_merged$orig.ident  == 'First_Trim.3'] <- 'First_trim'
placenta_merged$origin[placenta_merged$orig.ident  == 'First_Trim.4'] <- 'First_trim'
placenta_merged$origin[placenta_merged$orig.ident  == 'First_Trim.5'] <- 'First_trim'
placenta_merged$origin[placenta_merged$orig.ident  == 'First_Trim.6'] <- 'First_trim'

#Merging both the first trimester and term dataset
placenta_merged <- merge(x= term, y= first_trim)

#Extracting counts based on sample ID and aggregating them
counts <- AggregateExpression(placenta_merged, group.by = 'origin',
                              assays = 'RNA', slot = 'counts')
counts <- counts$RNA
counts <- as.data.frame(counts)


#Creating sample info sheet
Pseudo_sample_info <- read_csv("Pseudo sample info.csv")
Pseudo_sample_info <- Pseudo_sample_info %>% remove_rownames %>% column_to_rownames(var = 'Sample')

#Quality control
all(colnames(counts) %in% rownames(Pseudo_sample_info))
all(colnames(counts) == rownames(Pseudo_sample_info))

#Constructing Deseq2 object
Pseudo_subset_dds <- DESeqDataSetFromMatrix(countData = counts, 
                                            colData = Pseudo_sample_info,
                                            design = ~Group)

#Setting Reference Level
Pseudo_subset_dds$Group <- relevel(Pseudo_subset_dds$Group, ref = 'Term')

#Filtering the dataset
keep <- rowSums(counts(Pseudo_subset_dds)) >=10
Pseudo_subset_dds <- Pseudo_subset_dds[keep,]

#running Deseq2
Pseudo_subset_dds <- DESeq(Pseudo_subset_dds)
pseudo_subset_results <- results(Pseudo_subset_dds)

#Epxloring results and seeing variance in samples
Pseudo_subset_vsd <- vst(Pseudo_subset_dds, blind = T)
nudge <- position_nudge(x = 0.5, y =2)
transcriptome_PCA <- plotPCA(Pseudo_subset_vsd, intgroup = 'Group', returnData = T)

#Making PCA plots
var <- round(100 * attr(transcriptome_PCA, 'percentVar'))
p <- ggplot(transcriptome_PCA, aes(PC1, PC2, label = rownames(transcriptome_PCA))) 
p + scale_color_manual(values = c('First_trim' = '#2AC17D', 'Term' = '#C04EBD')) +
  stat_ellipse(level = 0.95, aes(fill = Group, color = Group)) + 
  geom_point(size = 4, aes(color = factor(transcriptome_PCA$Group))) +
  geom_text_repel(label = rownames(transcriptome_PCA), aes(color = factor(transcriptome_PCA$Group)), size = 4)+ 
  labs(x = sprintf('PC1 %.0f%%', var[1]),
       y = sprintf('PC2 %.0f%%', var[2]),
       color = 'Group')  +
  theme(axis.title.y = element_text(size = 15, face = 'bold'), axis.title.x = element_text(size =15, face = 'bold'),
        legend.text = element_text(size = 10))

#getting DEGs and making GO pathways
#Up regulated genes
palette.go <- colorRampPalette(c('#c3a4cf', 'thistle2', '#24bf66'))
up_reg_genes <- rownames(pseudo_subset_results[pseudo_subset_results$log2FoldChange > 0.6 & pseudo_subset_results$padj <= 0.05,])
Go_profiler_results_up <- gost(rownames(up_genes), organism = 'hsapiens')
Go_enrich_results_up <- enrichGO(rownames(up_genes), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
Go_results_up_graph <- plot(barplot(Go_results_up, showCategory = 20))
Go_results_up_graph <- Go_results_up_graph + scale_fill_gradientn(colours = palette.go(3))
Go_results_up_graph <- Go_results_up_graph + theme(axis.text.y = element_text(size = 9))

#Down regulated genes
down_reg_genes <- rownames(pseudo_subset_results[pseudo_subset_results$log2FoldChange < -0.6 ])
Go_profiler_results_down <- gost(rownames(down_genes), organism = 'hsapiens')
Go_enrich_results_down <- enrichGO(rownames(down_genes), OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
Go_results_down_graph <- plot(barplot(Go_results_down, showCategory = 20))
Go_results_down_graph <- Go_results_down_graph + scale_fill_gradientn(colours = palette.go(3))
Go_results_down_graph <- Go_results_down_graph + theme(axis.text.y = element_text(size = 9))
