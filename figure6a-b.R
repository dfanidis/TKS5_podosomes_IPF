library(ggplot2)
library(clusterProfiler)
library(enrichplot)

inPath <- './data'
outPath <- './results'
ifelse(dir.exists(outPath), return, dir.create(outPath))

`%!in%` <- Negate(`%in%`)
dpi <- 600

# =====================================================================
# Figure 6a
# =====================================================================

# To avoid changes in clusterProfiler versioning
gsea_res_fc <- readRDS(file.path(inPath, 'preranked_gsea_fc.rds'))

# Suppressed terms ordered by NES
supr <- gsea_res_fc
supr@result <- supr@result[which(supr@result[,'enrichmentScore'] < 0), ]
supr@result <- supr@result[order(supr@result$NES, decreasing = FALSE), ]
supr@result <- supr@result[1:7,]

# Re-order y-axis to bring terms with the higher absolute NES at the top
# of the axis
supr@result$Description <- factor(supr@result$Description, 
	levels = supr@result$Description)

dot_supr <- dotplot(
		supr,
		x = 'NES',
		color = 'p.adjust',
		showCategory = 10,
		size = NULL,
		font.size = 10,
		label_format = 100,
		title = 'Suppressed terms'
	) +
	xlab('') +
	scale_x_reverse() +
	scale_y_discrete(limits = rev(levels(supr@result$Description))) +
	theme_bw() +
	theme(
		plot.title=element_text(hjust = 0.5),
		axis.title.x=element_text(size= 10),
		axis.title.y=element_text(size= 10),
		axis.text.x=element_text(size= 9),
		axis.text.y=element_text(size= 9),
		legend.position= 'right',
		legend.title= element_text(size= 9),
		legend.text=element_text(size= 7),
		legend.key.size=unit(0.2, 'cm'),
	)
tiff(file = file.path(outPath, 'figure6a.tiff'), units = 'cm', width = 15,
	height = 8, compression = 'lzw', res = dpi, type='cairo')
	dot_supr
dev.off()

# =====================================================================
# Figure 6b
# =====================================================================
geneSetID <- gsea_res_fc@result[which(gsea_res_fc@result$Description == 
	"collagen-containing extracellular matrix"), "ID"]

collagenecm <- enrichplot::gseaplot2(gsea_res_fc, 
	geneSetID = geneSetID,
	title = "collagen-containing extracellular matrix",
	base_size = 5
)
tiff(file = file.path(outPath, "figure6b.tiff"),
	units = "cm", width = 10, height = 8, compression = "lzw",
  res = 600, type='cairo')
	collagenecm
dev.off()
