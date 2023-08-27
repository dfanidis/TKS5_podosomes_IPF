library(reshape)
library(ggplot2)
library(ggrepel)
library(Seurat)

inPath <- './data'
outPath <- './results'
ifelse(dir.create(outPath), return, dir.exists(outPath))
dpi <- 600

# ============================================================
# Figure S1A, C
# ============================================================
stats <- list(
	a = read.delim(file.path(inPath, 'GSE47460_GPL6480_stats_all.txt')),
	c = read.delim(file.path(inPath, 'GSE47460_GPL14550_stats_all.txt'))
)

stats <- lapply(stats, function(x){
	# Compute FDR
	x$FDR <- p.adjust(x$Pval, method= 'fdr')

	# Color points in respect to the P-value and FC thresholds of 0.05 and 1.2
	x$DEG <- rep('NS', nrow(x))
	x$DEG[which(x$FDR < 0.05 & x$log2FC < -0.2630344)] <- 'Down'
	x$DEG[which(x$FDR < 0.05 & x$log2FC > 0.2630344)] <- 'Up'

	# Which genes to name? 
	x$Name <- FALSE
	x$Name[which(x$Symbol %in% gn)] <- TRUE
	return(x)
})

# Plot
for(i in names(stats)) {

	title <- ifelse(i == 'a', 'GSE47460_GPL6480', 'GSE47460_GPL14550')

	figure_volcano <- ggplot() +
		geom_point(data = stats[[i]], mapping= aes(x = log2FC, y = -log10(FDR), 
			colour= DEG), size = 0.05) +
		scale_color_manual(values = c('Down' = '#00C400', 'NS' = 'grey30',
			'Up' = '#C40000')) +
		geom_text_repel(
			data= stats[[i]], 
			mapping=aes(x= log2FC, y= -log10(FDR), 
			label=ifelse(Name==TRUE, Symbol,''),),
			fontface = 'bold.italic',
			size= 2, 
			color = 'black',					# It controls labels' font color
	    	segment.colour= 'black',			# Connector line colors
	    	box.padding = unit(0.35, 'lines'),
	    	point.padding = unit(0.5, 'lines'),
	    	max.overlaps = Inf
	    ) +
		geom_hline(yintercept=1.30103, linetype='dashed', color='black') +
		geom_vline(xintercept=0.2630344, linetype='dashed', color='black') +
		geom_vline(xintercept=-0.2630344, linetype='dashed', color='black') +
		ggtitle(title) +
		xlab('log2FC') +
		ylab('-log10 FDR') +
		theme(
			legend.position = 'none',
			plot.title = element_text(size = 9, hjust = 0.5),
			axis.title = element_text(size = 9),
			axis.text.x = element_text(size = 8),  
			axis.text.y = element_text(size = 8),
			panel.grid.major = element_blank(),
		    panel.grid.minor = element_blank(),
		    panel.background = element_blank(), 
		    axis.line = element_line(colour = 'black')
		)  
	ggsave(filename=file.path(outPath, paste0('figureS1', i,'.tiff')),
		plot=figure_volcano,width = 4.985, height = 4.985, 
		dpi= dpi, units = 'cm')
}

# ============================================================
# Figure S1B, D
# ============================================================
data <- list(
	b = read.delim(file.path(inPath, 'GSE47460_A_normalizedExpressionValues.txt')),
	d = read.delim(file.path(inPath, 'GSE47460_B_normalizedExpressionValues.txt'))
)

# Select SH3PXD2A, COL1A1
data <- lapply(data, function(x){
	x <- x[which(x$Name %in% c('COL1A1', 'SH3PXD2A')), ]
	x$Name[1] <- 'TKS5'
	return(x)
})

# Correlation
cor_res <- lapply(data, function(x){

	tmp <- data.frame(
	    TKS5 = as.numeric(x[x$Name == 'TKS5', 2:ncol(x)]),
	    COL1A1 = as.numeric(x[x$Name == 'COL1A1', 2:ncol(x)])
	)

	resSpear <- cor.test(
	    tmp$TKS5, 
	    tmp$COL1A1, 
	    method= 'spearman'
	)

	out <- data.frame(
	    rho= resSpear$estimate,
	    p.value_rho= resSpear$p.value
	)

	return(out)
})

# Plots
plots_df <- lapply(data, function(x){
	data.frame(
		TKS5 = as.numeric(x[x$Name == 'TKS5', 2:ncol(x)]),
		COL1A1 = as.numeric(x[x$Name == 'COL1A1', 2:ncol(x)])
	)
})

for(i in names(plots_df)){

	title <- ifelse(i == 'b', 'GSE47460_GPL6480', 'GSE47460_GPL14550')

	title_x <- expression(italic('TKS5'))
	title_y <- expression(italic('COL1A1'))

	figure <- ggplot(plots_df[[i]], mapping= aes(x= TKS5, y= COL1A1)) +
		ggtitle(title) +
		geom_point(size = 1) +
		geom_smooth(method= 'lm') +
		theme_bw() +
		xlab(title_x) +
		ylab(title_y) +
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
	    annotate(geom = 'text', 
	    	x = ifelse(i == 'b', 11.5, 10.8), 
	    	y = ifelse(i == 'b', 16.9, 17.4), 
	    	label = ifelse(i == 'b', 
	    		paste0('ρ=', round(cor_res[['b']]$rho, 2)),
	    		paste0('ρ=', round(cor_res[['d']]$rho, 2))
	    	), 
	    	color = 'black', 
	    	size = 3
	    ) +
	    annotate(geom = 'text', 
	    	x = ifelse(i == 'b', 11.5, 10.8), 
	    	y = ifelse(i == 'b', 16, 16.5), 
	    	label = ifelse(i == 'b', 
	    		paste0('p-val=', signif(cor_res[['b']]$p.value_rho, 3)),
	    		paste0('p-val=', signif(cor_res[['d']]$p.value_rho, 3))
	    	), 
	    	color = 'black',
	    	size = 3
	    )

		ggsave(filename=file.path(outPath, paste0('figureS1', i, '.tiff')),
			plot = figure, width = 4.985, height = 4.985, 
			dpi = dpi, units = 'cm')

}

# ============================================================
# Figure S1E-F
# ============================================================
reyf <- readRDS(file.path(inPath, 'reyfman2019.rds'))
DefaultAssay(reyf) <- 'RNA'

# Remove unassigned cells
reyf <- subset(reyf, subset = integr.cell.ident != 'Unassigned')
reyf$integr.cell.ident <- droplevels(reyf$integr.cell.ident)

dot <- DotPlot(reyf, features = 'SH3PXD2A', 
  cols = c('#ececec', '#C40000'), col.min = 0) +
  theme(
    plot.title=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(angle = 45, hjust = 1,
      size= 10, color = 'black'),
    axis.text.y=element_text(size= 10, color = 'black',
      angle = 90, hjust = 0.5),
    legend.position= 'right',
    legend.title= element_text(size= 11),
    legend.text=element_text(size= 10),
    legend.key=element_blank(),
    legend.spacing.x=unit(0.2, 'cm')
  ) +
  coord_flip()

feat <- FeaturePlot(reyf, reduction = 'tsne', features = 'SH3PXD2A',
  cols = c('#ececec', '#c40000'), label = TRUE, repel = TRUE, 
  order = TRUE) + 
  theme(
    plot.title = element_blank(),
    axis.title.x=element_text(size=11),
    axis.title.y=element_text(size=11),
    axis.text.x=element_text(size=10), 
    axis.text.y=element_text(size=10),
    legend.position = 'right',
    legend.title= element_text(size= 10),
    legend.text=element_text(size= 10),
    legend.key=element_blank()
  )

tiff(file.path(outPath, 'figureS1e.tiff'),
  width = 16, height = 7.5, units = 'cm', compression = 'lzw',
  res = dpi, type='cairo')
  feat
dev.off()

tiff(file.path(outPath, 'figureS1f.tiff'),
  width = 16, height = 11, units = 'cm', compression = 'lzw',
  res = dpi, type='cairo')
  dot
dev.off()

# ============================================================
# Fig S1F - TKS5 marker gene 
# ============================================================
Idents(reyf) <- reyf$integr.cell.ident
clusters <- vector('list', length = length(levels(reyf)))
names(clusters) <- levels(reyf)

for (i in names(clusters)) {
  tryCatch({
    clusters[[i]] <- FindMarkers(reyf, 
      ident.1 = i,
      features = 'SH3PXD2A',
      logfc.threshold = 0,
      only.pos = TRUE
    )
  }, error= function(e) {cat('ERROR :', conditionMessage(e), '\n')})

  if (is.null(clusters[[i]])) {clusters[[i]] <- NA}
  clusters[[i]]$gene <- rownames(clusters[[i]])
}
clusters[is.na(clusters)] <- NULL
clusters <- do.call('rbind', clusters)

clusters$FC <- ifelse(clusters$avg_log2FC < 0,
  paste(round(1/(2^clusters$avg_log2FC), 3), 'down'),
  paste(round(2^clusters$avg_log2FC, 3), 'up')
)
clusters$cluster <- gsub('\\..*', '', rownames(clusters))
rownames(clusters) <- NULL

# Significant findings
clusters[clusters$p_val_adj < 0.05 & abs(clusters$avg_log2FC) >= 0.2630344, ]
##          p_val avg_log2FC pct.1 pct.2    p_val_adj     gene       FC
## 6 5.975007e-30  0.6143025 0.254 0.102 1.365588e-25 SH3PXD2A 1.531 up
##       cluster
## 6 Fibroblasts

# ============================================================
# Figure S1G
# ============================================================
fibros <- subset(reyf, subset = RNA_snn_res.0.1 != '') 
fibros@reductions$tsne <- readRDS(file.path(inPath, 'reyfman2019_fibros_tsne.rds'))
fibros$RNA_snn_res.0.1 <- factor(fibros$RNA_snn_res.0.1, levels = c('0', '1', '2', '3'))
Idents(fibros) <- fibros$RNA_snn_res.0.1

tsne <- DimPlot(fibros, label = TRUE,
    cols = c("#0062c4", "#00c400", "#eed202", "#c40000"),
    pt.size = 0.03
  ) +
  theme(
    plot.title=element_text(size = 9, hjust = 0.5),
    axis.title.x=element_text(size=9),
    axis.title.y=element_text(size=9),
    axis.text.x=element_text(size=8), 
    axis.text.y=element_text(size=8),
    strip.text.x=element_text(size=10),
    strip.text.y=element_text(size=10),
    legend.position = "none"
  )

tiff(file.path(outPath, "figureS1g.tiff"),
  width = 4.5972, height = 4.5972, units = "cm", compression = "lzw",
  res = dpi, type='cairo')
  tsne
dev.off()