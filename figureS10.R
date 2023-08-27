library(reshape)
library(ggplot2)
library(proxy)
library(openxlsx)
library(dorothea)
library(igraph)
library(visNetwork)

inPath <- './data'
outPath <- './results'
ifelse(dir.create(outPath), return, dir.exists(outPath))
dpi <- 600

a
heatmap_top100Quant_noLegend.tiff
b
network_all_de_high_interactions.tiff
c
heatmap_go_degs.tiff
d
heatmap_top50Noble_degs.tiff

# ============================================================
# Figure S10a
# ============================================================
stats <- read.delim(gzfile(file.path('../Figure 6/data', 
	'metaseqr_all_out_KN_TGF_vs_WT_TGF.txt.gz')))

cnts <- read.delim(gzfile(file.path(inPath, 'normalized_counts_table.txt.gz')))
cnts <- merge(cnts, stats[,c('gene_id', 'gene_name')], by = 'gene_id',
	sort = FALSE, all.x = TRUE)
rownames(cnts) <- cnts$gene_id
cnts$gene_id <- NULL

# Find the DEGs
tmp <- cbind(
	stats[,1:8],
	stats[,grep('^meta_FDR|^natural_normalized_fold_change', colnames(stats))]
)
colnames(tmp)[c(ncol(tmp)-1,ncol(tmp))] <- c('FDR', 'FoldChange')
up_degs <- tmp[which(tmp$FDR < 0.05 & tmp$FoldChange >= 1.2), ]
down_degs <- tmp[which(tmp$FDR < 0.05 & tmp$FoldChange <= 1/1.2), ]
down_degs$FoldChange <- -(1/down_degs$FoldChange)
degs <- rbind(up_degs, down_degs)

degs <- degs[order(degs$FoldChange, decreasing = TRUE), ]
top100 <- c(degs[1:50, 'gene_id'], 
  degs[(nrow(degs)-49):nrow(degs), 'gene_id']
)

# Maintain normalized counts for the 100 degs
top100_df <- cnts[top100, ]
rownames(top100_df) <- top100_df$gene_name
top100_df$gene_name <- NULL
colnames(top100_df) <- c('KN_TGF_1', 'KN_TGF_2', 'KN_TGF_3',
  'WT_TGF_1', 'WT_TGF_2', 'WT_TGF_3')

# Remove genes with 0 reads
top100_df <- top100_df[which(rowSums(top100_df) != 0), ]

# Scale and cluster
top100_df <- as.data.frame(t(scale(t(top100_df))))
rowDist <- dist(top100_df, by_rows = TRUE, method = 'Euclidean')
rowClust <- hclust(rowDist)
rowOrder <-  rownames(top100_df)[rowClust$order]
colDist <- dist(top100_df, by_rows = FALSE, method = 'Euclidean')
colClust <- hclust(colDist)
colOrder <- colnames(top100_df)[colClust$order]

data <- melt(as.matrix(top100_df))
colnames(data)[1:2] <- c('Gene', 'Sample')

data$Gene <- factor(data$Gene, levels = rowOrder)
data$Sample <- factor(data$Sample, levels = colOrder)

figureS10a <- ggplot(data, aes(x = Sample, y = Gene, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colors = colorRampPalette(c('blue', 'red'))(20)) +
  guides(
      fill = guide_colourbar(
          barwidth = 0.5,
          barheight = 4
      )
  ) +
  theme_bw() +
  theme(
      plot.title=element_text(hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size= 5, angle = 45, hjust = 1),
      axis.text.y=element_text(size= 5),
      legend.position= 'right',
      legend.title= element_blank(),
      legend.text=element_text(size= 5),
      legend.key=element_blank(),
      legend.key.height=unit(4, 'cm')
  )

tiff(file.path(outPath, 'figureS10a.tiff'),
  width = 5, height = 15, units = 'cm',
  res = dpi, type='cairo')
  figureS10a
dev.off()

# ============================================================
# Figure S10c
# ============================================================
gsea <- readRDS(file.path('../Figure 6/data', 'preranked_gsea_fc.rds'))
go0062023 <- unlist(strsplit(gsea@result['GO:0062023', 'core_enrichment'], 
	split = '/'))
go0062023 <- stats[which(stats$gene_id %in% go0062023), 'gene_name']
go0062023 <- intersect(go0062023, degs$gene_name)

go_df <- cnts[which(cnts$gene_name %in% go0062023), ]
rownames(go_df) <- go_df$gene_name
go_df$gene_name <- NULL
colnames(go_df) <- c('KN_TGF_1', 'KN_TGF_2', 'KN_TGF_3',
  'WT_TGF_1', 'WT_TGF_2', 'WT_TGF_3')

# Remove genes with 0 reads
go_df <- go_df[which(rowSums(go_df) != 0), ]

# Scale and cluster
go_df <- as.data.frame(t(scale(t(go_df))))
rowDist <- dist(go_df, by_rows = TRUE, method = 'Euclidean')
rowClust <- hclust(rowDist)
rowOrder <-  rownames(go_df)[rowClust$order]
colDist <- dist(go_df, by_rows = FALSE, method = 'Euclidean')
colClust <- hclust(colDist)
colOrder <- colnames(go_df)[colClust$order]

data <- melt(as.matrix(go_df))
colnames(data)[1:2] <- c('Gene', 'Sample')

data$Gene <- factor(data$Gene, levels = rowOrder)
data$Sample <- factor(data$Sample, levels = colOrder)

# GO heatmap
figures10c <- ggplot(data, 
  aes(x = Sample, y = Gene, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colors = colorRampPalette(c('blue', 'red'))(20)) +
  guides(
      fill = guide_colourbar(
          barwidth = 0.5,
          barheight = 4
      )
  ) +
  theme_bw() +
  theme(
      plot.title=element_text(hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size= 3,
          angle= 45, hjust= 1),
      axis.text.y=element_text(size= 3),
      legend.position= 'right',
      legend.title= element_blank(),
      legend.text=element_text(size= 5),
      legend.key=element_blank(),
      legend.key.height=unit(4, 'cm'),
  )
tiff(file.path(outPath, 'figureS10c.tiff'),
  width = 5, height = 12, units = 'cm',
  res = dpi, type='cairo')
  figures10c
dev.off()

# ============================================================
# Figure S10d
# ============================================================
top50_invasive <- read.delim(file.path(inPath, 'top50_genes_suptable1_GSE118933.txt'))

# Load homology mapping (to avoid biomaRt connection problems 
# and version changes)
mart_homol <- read.delim(gzfile(file.path(inPath, 'mart_homol.txt.gz')))

top50_homol <- mart_homol[which(mart_homol$ensembl_gene_id %in% 
	top50_invasive$ENSGid), ]
top50_homol <- na.omit(top50_homol)

# Normalized expression values of the homologue genes
top50_mouse <- unique(top50_homol$mmusculus_homolog_ensembl_gene)
top50_mouse <- intersect(top50_mouse, degs$gene_id)

top50_df <- cnts[top50_mouse, ]
rownames(top50_df) <- top50_df$gene_name
top50_df$gene_name <- NULL
colnames(top50_df) <- c('KN_TGF_1', 'KN_TGF_2', 'KN_TGF_3',
  'WT_TGF_1', 'WT_TGF_2', 'WT_TGF_3')

# Remove genes with 0 reads in all samples
top50_df <- top50_df[which(rowSums(top50_df) != 0), ]

# Scale and cluster
top50_df <- as.data.frame(t(scale(t(top50_df))))
rowDist <- dist(top50_df, by_rows = TRUE, method = 'Euclidean')
rowClust <- hclust(rowDist)
rowOrder <-  rownames(top50_df)[rowClust$order]
colDist <- dist(top50_df, by_rows = FALSE, method = 'Euclidean')
colClust <- hclust(colDist)
colOrder <- colnames(top50_df)[colClust$order]

data <- melt(as.matrix(top50_df))
colnames(data)[1:2] <- c('Gene', 'Sample')

data$Gene <- factor(data$Gene, levels = rowOrder)
data$Sample <- factor(data$Sample, levels = colOrder)

# top50 Noble
figures10d <- ggplot(data, 
  aes(x = Sample, y = Gene, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colors = colorRampPalette(c('blue', 'red'))(20)) +
  guides(
      fill = guide_colourbar(
          barwidth = 0.5,
          barheight = 4
      )
  ) +
  theme_bw() +
  theme(
      plot.title=element_text(hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size= 5,
          angle= 45, hjust= 1),
      axis.text.y=element_text(size= 5),
      legend.position= 'right',
      legend.title= element_blank(),
      legend.text=element_text(size= 5),
      legend.key=element_blank(),
      legend.key.height=unit(4, 'cm')
  )
tiff(file.path(outPath, 'figureS10d.tiff'),
  width = 5, height = 5, units = 'cm',
  res = dpi, type='cairo')
  figures10d
dev.off()

# ============================================================
# Figure S10b
# ============================================================

# Function to change the relative position of
# network edges (stack overflow)
segementsDf <- function(data, shorten.end){
  data$dx = data$to.x - data$from.x
  data$dy = data$to.y - data$from.y
  data$dist = sqrt( data$dx^2 + data$dy^2 )
  data$px = data$dx/data$dist
  data$py = data$dy/data$dist

  data$to.x = data$to.x - data$px * shorten.end
  data$to.y = data$to.y - data$py * shorten.end

  return(data)
}

# All DEG targets
net <- as.data.frame(dorothea_mm)
targets <- merge(net, degs[,c('gene_name', 'FoldChange')],
  by.x = 'target', by.y = 'gene_name'
)
targets <- targets[,c('tf', 'target', 'confidence', 'mor')]

# Isolate interactions of specific confidence levels
targets <- targets[which(targets$confidence == 'A'), ]

# Which of the tfs are de?
targets <- targets[which(targets$tf %in% degs$gene_name), ]

# Remove inconsistent tf-target pairs based on mor
targets$rm <- NULL
for(i in 1:nrow(targets)){

  tmp_mor <- targets$mor[i]
  tmp_tf <- sign(degs[which(degs$gene_name == targets[i, 'tf']),
    'FoldChange'])
  tmp_target <- sign(degs[which(degs$gene_name == targets[i, 'target']),
    'FoldChange'])
  
  if(tmp_mor == 1 & tmp_tf == tmp_target) {
    targets$rm[i] <- 'keep'
  } else if(tmp_mor == -1 & tmp_tf != tmp_target) {
    targets$rm[i] <- 'keep'
  } else {
    targets$rm[i] <- 'remove'
  }
}
targets <- targets[which(targets$rm == 'keep'), ]
targets$rm <- NULL

# Create an igraph object
nodes_all <- unique(c(targets$tf, targets$target))
edges_all <- targets

g_igraph <- graph_from_data_frame(edges_all, 
  directed= FALSE, vertices= nodes_all)

# Create a layout dataframe
set.seed(21)
layout_df <- as.data.frame(layout_nicely(g_igraph))
layout_df$genes <- nodes_all

# Add direction of deregulation
layout_df <- merge(layout_df, degs[,c('gene_name', 'FoldChange')],
  by.x = 'genes', by.y = 'gene_name',
  all.x = TRUE
)
layout_df$FoldChange[is.na(layout_df$FoldChange)] <- 0
layout_df$de_direction <- as.factor(sign(layout_df$FoldChange))
layout_df <- layout_df[,c('V1', 'V2', 'genes', 'de_direction')]

# Add tf metadata
layout_df$tf <- layout_df$genes
for(i in 1:nrow(layout_df)){
  if(layout_df$tf[i] %in% targets$tf){ # is it a tf?
      layout_df$tf[i] <- '1'
  } else { # is it not a TF?
    layout_df$tf[i] <- '0'
  }
}
layout_df$tf <- factor(layout_df$tf)

# Get the edge information using the get.data.frame
g_df <- get.data.frame(g_igraph)  
# Match the from locations from the node data.frame we previously connected
g_df$from.x <- layout_df$V1[match(g_df$from, layout_df$genes)]  
g_df$from.y <- layout_df$V2[match(g_df$from, layout_df$genes)]
# Match the to locations from the node data.frame we previously connected
g_df$to.x <- layout_df$V1[match(g_df$to, layout_df$genes)]
g_df$to.y <- layout_df$V2[match(g_df$to, layout_df$genes)]

# Shorten segments end for better visualization of node labels 
g_df <- segementsDf(data = g_df, shorten.end = 0.2)

# mor to factor for different linetypes for inhibition and/or activation
g_df$mor <- as.factor(g_df$mor)

# Change names
g_df$mor <- ifelse(g_df$mor == 1, 'Activation', 'Repression')
layout_df$de_direction <- ifelse(layout_df$de_direction == 1, 'Up-regulated', 'Down-regulated')
layout_df$tf <- ifelse(layout_df$tf == 1, 'TF', 'Target')

figures10b <- ggplot() +
    geom_segment(
      data = g_df,
      aes(x = from.x, xend = to.x, y = from.y, yend = to.y, 
        linetype = mor),
      colour = rep('#D3D3D3', nrow(g_df)),
      arrow = arrow(angle = 20, length = unit(0.05, 'inches'),
        ends = 'last', type = 'open')
    ) +
    scale_linetype_manual(values=c('Repression' = 'dashed', 'Activation' = 'solid')) +
    geom_point(
      data = layout_df, 
      aes(x = V1, y = V2, color = de_direction, shape = tf), 
      size = 5
    ) +
    scale_shape_manual(values = c('Target' = 1,  'TF' = 2))+
    scale_color_manual(values = c('Down-regulated' = '#4CBB17', 'Up-regulated' = 'red')) +
    geom_text(
      data = layout_df,
      aes(x = V1, y = V2, label = genes), 
      size = 1.5
    ) + # add the node labels
    scale_x_continuous(expand = c(0, 1))+
    scale_y_continuous(expand = c(0, 1))+
    theme_bw() +
    theme(
      legend.position = 'right',
      legend.text = element_text(size = 7),
      legend.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank()
    )

tiff(file.path(outPath, 'figureS10b.tiff'),
  width = 12, height = 10, units = 'cm',
  res = dpi, type='cairo')
  figures10b
dev.off()