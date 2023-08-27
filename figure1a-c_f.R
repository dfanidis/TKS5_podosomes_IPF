library(reshape)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(openxlsx)
library(Seurat)

`%!in%` <- Negate(`%in%`)
dpi <- 600

inPath <- './data'
outPath <- './results'
ifelse(dir.exists(outPath), return, dir.create(outPath))

# ============================================================
# Figure 1A
# ============================================================
de <- read.xlsx(file.path(inPath, 'SH3PXD2A-Sh3pxd2a_DE_statistics.xlsx'))
de$Comparison <- gsub('\n', '_', de$Comparison)

# Drop non significant comparisons
de <- de[which(abs(de$log2FC) > 0.2630344 & de$FDR < 0.05), ]

# Keep lung tissue and IPF_vs_Ctrl comparison
lung <- de[which(de$Tissue == 'Lung tissue' & de$Comparison == 'IPF_vs_Ctrl'), ]

# Plot
plotDF <- rbind(
    lung[, c('Tissue', 'Comparison', 'log2FC', 'FDR')]
)
rownames(plotDF) <- NULL
plotDF$Tissue <- factor(plotDF$Tissue, levels = 'Lung tissue')

y_title <- expression(paste(italic('TKS5'), '(log2FC)'))

figure1a <- ggplot(plotDF, aes(x = Tissue, y = log2FC)) +
        geom_boxplot(color = 'black', outlier.size = 0) +
        geom_jitter(data = plotDF, aes(color = Tissue), 
            position = position_jitter(0.1), size = 1) +
        geom_hline(yintercept=0, linetype='dashed', color='black') +
        xlab(' ') +
        ylab(y_title) +
        scale_color_manual(values = '#C40000') +
        theme_bw() +
        theme(
            plot.title=element_text(size = 11, hjust = 0.5),
            axis.title.x=element_text(size= 11, face = 'bold', colour = 'black'),
            axis.title.y=element_text(size= 11, face = 'bold', colour = 'black'),
            axis.text.x=element_text(size= 8, face = 'bold', colour = 'black'),
            axis.text.y=element_text(size= 8, face = 'bold', colour = 'black'),
            axis.ticks=element_line(color = 'black'),
            legend.position= 'none',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(colour = 'black'),
            plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), 'cm')
        )

tiff(file.path(outPath, 'figure1a.tiff'), 
    width = 2.4925, height = 4.985,
    units = 'cm', res = dpi, type = 'cairo')
    figure1a
dev.off()

# ============================================================
# Figure 1B
# ============================================================
stat <- read.delim(file.path(inPath, 'GSE53845_stats_all.txt'))
stat$FDR <- p.adjust(stat$Pval, method= 'fdr')
stat$GSE <- 'GSE53845_GPL6480'

# Color points in respect to the P-value and FC thresholds of 0.05 and 1.2
stat$DEG <- rep('NS', nrow(stat))
stat$DEG[which(stat$FDR < 0.05 & stat$log2FC < -0.2630344)] <- 'Down'
stat$DEG[which(stat$FDR < 0.05 & stat$log2FC > 0.2630344)] <- 'Up'

# Which genes to name? 
stat$Name <- FALSE
stat$Name[which(stat$Symbol %in% c('COL1A1', 'SH3PXD2A'))] <- TRUE

# Plot
figure1b <- ggplot() +
    geom_point(data = stat, mapping= aes(x = log2FC, y = -log10(FDR), 
        colour= DEG), size = 0.05) +
    scale_color_manual(values = c('Down' = '#00C400', 'NS' = 'grey30',
        'Up' = '#C40000')) +
    geom_text_repel(
        data = stat, 
        mapping=aes(x= log2FC, y= -log10(FDR), 
        label=ifelse(Name==TRUE, Symbol,''),),
        fontface = 'bold.italic',
        size= 2, 
        color = 'black',                    # It controls labels' font color
        segment.colour= 'black',            # Connector line colors
        box.padding = unit(0.35, 'lines'),
        point.padding = unit(0.5, 'lines'),
        max.overlaps = Inf
    ) +
    geom_hline(yintercept=1.30103, linetype='dashed', color='black') +                      # draw Pvalue threshold line
    geom_vline(xintercept=0.2630344, linetype='dashed', color='black') +                    # draw FC threshold line
    geom_vline(xintercept=-0.2630344, linetype='dashed', color='black') +                   # draw FC threshold line
    ggtitle(unique(stat$GSE)) +
    xlab('log2FC') +
    ylab('-log10 FDR') +
    theme(
        plot.title=element_text(size = 11, hjust = 0.5),
        axis.title.x=element_text(size= 11, colour = 'black'),
        axis.title.y=element_text(size= 11, colour = 'black'),
        axis.text.x=element_text(size= 8, face = 'bold', colour = 'black'),
        axis.text.y=element_text(size= 8, face = 'bold', colour = 'black'),
        axis.ticks=element_line(color = 'black'),
        legend.position= 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), 'cm')
    )  
tiff(file.path(outPath, 'figure1b.tiff'),
  width = 4.985, height = 4.985, units = 'cm', compression = 'lzw',
  res = dpi, type='cairo')
  figure1b
dev.off()

# ============================================================
# Figure 1C
# ============================================================
data <- read.delim(file.path(inPath, 'GSE53845normalizedExpressionValues.txt'))

# Select SH3PXD2A, COL1A1
data <- data[which(data$Name %in% c('COL1A1', 'SH3PXD2A')), ]
data$Name[1] <- 'TKS5'

# Correlation
tmp <- data.frame(
    TKS5 = as.numeric(data[data$Name == 'TKS5', 2:ncol(data)]),
    COL1A1 = as.numeric(data[data$Name == 'COL1A1', 2:ncol(data)])
)

resSpear <- cor.test(
    tmp$TKS5, 
    tmp$COL1A1, 
    method= 'spearman'
) 

correlation <- data.frame(
    rho= resSpear$estimate,
    p.value_rho= resSpear$p.value
)

# Plots
plots <- data.frame(
	TKS5 = as.numeric(data[data$Name == 'TKS5', 2:ncol(data)]),
	COL1A1 = as.numeric(data[data$Name == 'COL1A1', 2:ncol(data)])
)

# Publication quality figure
x_title <- expression(italic('TKS5'))
y_title <- expression(italic('COL1A1'))
figure1c <- ggplot(plots, mapping= aes(x= TKS5, y= COL1A1)) +
	ggtitle('GSE53845_GPL6480') +
	geom_point(size = 1) +
	geom_smooth(method= 'lm') +
	theme_bw() +
	xlab(x_title) +
	ylab(y_title) +
    theme(
        plot.title=element_text(size = 11, hjust = 0.5),
        axis.title.x=element_text(size= 11, face = 'bold', colour = 'black'),
        axis.title.y=element_text(size= 11, face = 'bold', colour = 'black'),
        axis.text.x=element_text(size= 8, face = 'bold', colour = 'black'),
        axis.text.y=element_text(size= 8, face = 'bold', colour = 'black'),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        axis.ticks=element_line(color = 'black'),
        legend.position= 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black'),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), 'cm')
    ) +
    annotate(geom = 'text', x = 0.8, y = 4.5, 
    	label = 'ρ=0.72', color = 'black', size = 3) +
    annotate(geom = 'text', x = 0.8, y = 3.8, 
    	label = 'p-val=8.92E-08', color = 'black', size = 3)

tiff(file.path(outPath, 'figure1c.tiff'),
  width = 4.985, height = 4.985, units = 'cm', compression = 'lzw',
  res = dpi, type='cairo')
  figure1c
dev.off()

# ============================================================
# Figure 1F
# ============================================================
inPath <- '../Figure S1/data'
genes <- c('SH3PXD2A', 'ACTA2', 'COL1A1', 'CD274', 'PDGFRA', 
    'PDGFRB', 'PPARG', 'CD44', 'PDPN', 'AXIN2', 'HAS1')

reyf <- readRDS(file.path(inPath, 'reyfman2019.rds'))
DefaultAssay(reyf) <- 'RNA'
fibros <- subset(reyf, subset = RNA_snn_res.0.1 != '')

tks5pos <- subset(fibros, subset = SH3PXD2A > 0, slot = 'data')
tks5pos$RNA_snn_res.0.1 <- factor(tks5pos$RNA_snn_res.0.1, 
    levels = c('0', '1', '2', '3'))
Idents(tks5pos) <- tks5pos$RNA_snn_res.0.1

dot <- DotPlot(tks5pos, features = genes, 
  cols = c('#ececec', '#C40000'), col.min = 0) +
  theme(
    plot.title=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(angle = 0, vjust = 0.5, size= 10),
    axis.text.y=element_text(size= 10),
    legend.position= 'right',
    legend.title= element_text(size= 9),
    legend.text=element_text(size= 9),
    legend.key=element_blank()
  ) +
  coord_flip()

tiff(file.path(outPath, 'figure1f.tiff'),
  width = 8, height = 10, units = 'cm', compression = 'lzw',
  res = dpi, type='cairo')
  dot
dev.off()

# ============================================================
# Fig 1F - TKS5 marker genes
# ============================================================
Idents(tks5pos) <- tks5pos$RNA_snn_res.0.1
clusters <- vector("list", length = length(levels(tks5pos)))
names(clusters) <- levels(tks5pos)

for (i in names(clusters)) {
  tryCatch({
    clusters[[i]] <- FindMarkers(tks5pos, 
      ident.1 = i,
      features = genes,
      logfc.threshold = 0,
      only.pos = TRUE
    )
  }, error= function(e) {cat("ERROR :", conditionMessage(e), "\n")})

  if (is.null(clusters[[i]])) {clusters[[i]] <- NA}
  clusters[[i]]$gene <- rownames(clusters[[i]])
}
clusters[is.na(clusters)] <- NULL
clusters <- do.call("rbind", clusters)

clusters$FC <- ifelse(clusters$avg_log2FC < 0,
  paste(round(1/(2^clusters$avg_log2FC), 3), "down"),
  paste(round(2^clusters$avg_log2FC, 3), "up")
)
clusters$cluster <- gsub("\\..*", "", rownames(clusters))
rownames(clusters) <- NULL

# Significant findings
clusters[which(clusters$p_val_adj < 0.05 &
    abs(clusters$avg_log2FC) >= 0.2630344), ]
##           p_val avg_log2FC pct.1 pct.2    p_val_adj     gene        FC cluster
## 1  9.282388e-08  0.9135817  1.00 1.000 2.121490e-03 SH3PXD2A  1.884 up       0
## 7  4.976326e-15  3.6755399  1.00 0.315 1.137339e-10    ACTA2 12.778 up       1
## 11 3.901067e-16  1.3210378  0.75 0.019 8.915889e-12    PPARG  2.498 up       3