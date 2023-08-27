# Process LINCs query results

library(cmapR)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(rhdf5)

inPath <- './data'
outPath <- './results'
ifelse(dir.exists(outPath), return, dir.create(outPath))

# ===========================================================
# Figure 7a
# ===========================================================

system('tar -xvf ./data/my_analysis.tar.gz -C ./data')
inPathTmp <- paste0(inPath, '/my_analysis.sig_queryl1k_tool.63394fc6fff80600124b0b25/matrices/query')
ncs <- parse_gctx(file.path(inPathTmp, 'ncs.gctx'))
fdr <- parse_gctx(file.path(inPathTmp, 'fdr_qvalue.gctx'))

# # Alternatively, to avoid potential 'parse_gctx' bugs
# load(file.path(inPath, 'ncs_fdr.RData')) 

res <- cbind(ncs@mat[,'TAG'], fdr@mat[,'TAG'], ncs@rdesc)
colnames(res)[1:2] <- c('NCS', 'FDR')
rownames(res) <- NULL

# Filter results
res <- res[which(res$FDR < 0.05), ]
rownames(res) <- NULL

res_pos <- res[which(res$NCS > 0), ]
res_neg <- res[which(res$NCS < 0), ]

# Order
res_pos <- res_pos[order(res_pos$NCS, decreasing = TRUE), ]
res_neg <- res_neg[order(res_neg$NCS, decreasing = FALSE), ]

resList <- list(res_pos, res_neg)
resList <- lapply(resList, function(x){
	out <- subset(x, select = c(NCS, FDR, pert_iname, cell_iname,
		pert_type, pert_idose, pert_itime, moa, nsample, tas, 
		target_name))
	out[,1:2] <- round(out[,1:2], 3)
	rownames(out) <- NULL

	out <- data.frame(lapply(out, function(y){
		gsub('-666', '-', y)
	}))

	return(out)
})
names(resList) <- c('PositiveNCS', 'NegativeNCS')

# Further filtering
sel <- resList$PositiveNCS

# Keep compounds, peptides and co
sel <- sel[which(sel$pert_type %in% c('trt_cp', 'trt_lig')), ]
# Only with a known moa
sel <- sel[which(sel$moa != '-'), ]
# Tested in >1 replicates
sel <- sel[which(sel$nsample > 1), ]
# Only 24h treatments
sel <- sel[which(sel$pert_itime == '24 h'), ]
# Only with target name
sel <- sel[which(sel$target_name != '-'), ]

# Plot
sel_df <- sel[1:15,]
sel_df <- sel_df[order(sel_df$NCS, decreasing = TRUE), ]
sel_df$y <- paste(sel_df$pert_iname, sel_df$cell_iname, sep = ' - ')
sel_df$y <- factor(sel_df$y, levels = rev(sel_df$y))

sel_df$moa <- c('Aurora kinase inhibitor',
	'NTPDase inhibitor',
	'Smoothened receptor antagonist',
	'HDAC inhibitor',
	'VEGFR, KIT, PDGFR inhibitor',
	'ALK inhibitor',
	'Protein synthesis inhibitor',
	'Insulin secretagogue',
	'JAK inhibitor',
	'Phosphodiesterase inhibitor',
	'Estrogen receptor antagonist, SERM',
	'Src inhibitor',
	'Protein synthesis, BIG1 inhibitor',
	'PI3K inhibitor',
	'KIT, PDGFR, VEGFR inhibitor'
)

sel_df$moa <- factor(sel_df$moa, levels = rev(sel_df$moa))

ncols <- nrow(sel_df) + 5
cols <- colorRampPalette(brewer.pal(9, 'Reds'))(ncols)[5:ncols]
cols_y <- c('#C40000', 'grey30', 'grey30', '#C40000', rep('grey30', 6),
	'#C40000', rep('grey30', 4))

figure <- ggplot(sel_df, mapping= aes(x= moa, y= NCS, fill = NCS)) +
	geom_bar(stat = 'identity') +
	scale_fill_manual(values = cols) +
	theme_bw() +
    theme(
			plot.title=element_text(size = 9, hjust = 0.5),
			axis.title.x=element_text(size=6),
			axis.title.y=element_blank(),
			axis.text.x=element_text(size=4, angle = 45, hjust = 1), 
			axis.text.y=element_text(size=4, color = cols_y),
			legend.position = 'none',
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(colour = 'black', size = 0.1), # grey50
			axis.ticks = element_line(size = 0.1, color = 'black')
    ) + coord_flip()
tiff(file.path(outPath, 'figure7a.tiff'),
  width = 9, height = 3.5, units = 'cm', compression = 'lzw',
  res = 600, type='cairo')
  figure
dev.off()
