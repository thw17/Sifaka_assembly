# Set working directory to directory containing all .hist files
setwd()

# Imports
library("ggplot2")
library("dplyr")
library("cowplot")

# Multiplot function from Cookbook for R
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# processed using command for hg38: awk '$1 ~ /^chr/' <hg38 bed file> | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# this only looks at regions on scaffolds beginning with chr
# region lengths updated Sep 17, 2018

hg38_cds <- 35825701
hg38_exon <- 122588673
hg38_gene <- 1596285318
hg38_intron <- 1475995647
hg38_utr <- 86762972

mmul_cds <- 35561524
mmul_exon <- 87403016
mmul_gene <- 1250134352
mmul_intron <- 1167997710
mmul_utr <- 51852833

pcoq_cds <- 32220035
pcoq_exon <- 47565723
pcoq_gene <- 942704083
pcoq_intron <- 897027873
pcoq_utr <- 15357019

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

## Cds plots for NimbleGen Sifaka (4) and Rhesus (4 samples)

# Pcoq cds
cds_f249_pcoq <- mutate(read.table("F249.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F249")
cds_f406_pcoq <- mutate(read.table("F406.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F406")
cds_m288_pcoq <- mutate(read.table("M288.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M288")
cds_m418_pcoq <- mutate(read.table("M418.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M418")

cds_f249_pcoq <- mutate(cds_f249_pcoq, percent= (sum(cds_f249_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds) # pcoq_cds = 32220035
cds_f406_pcoq <- mutate(cds_f406_pcoq, percent= (sum(cds_f406_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_m288_pcoq <- mutate(cds_m288_pcoq, percent= (sum(cds_m288_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_m418_pcoq <- mutate(cds_m418_pcoq, percent= (sum(cds_m418_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)

cds_pcoq_combined <- rbind(cds_f249_pcoq, cds_f406_pcoq, cds_m288_pcoq, cds_m418_pcoq)
head(cds_pcoq_combined)

# Mmul cds
cds_WI055_mmul <- mutate(read.table("WI055.mmul.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI055")
cds_WI056_mmul <- mutate(read.table("WI056.mmul.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI056")
cds_WI057_mmul <- mutate(read.table("WI057.mmul.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI057")
cds_WI059_mmul <- mutate(read.table("WI059.mmul.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI059")

cds_WI055_mmul <- mutate(cds_WI055_mmul, percent= (sum(cds_WI055_mmul$V2) - cumsum(V2) + V2) / mmul_cds) # mmul_cds = 35561524
cds_WI056_mmul <- mutate(cds_WI056_mmul, percent= (sum(cds_WI056_mmul$V2) - cumsum(V2) + V2) / mmul_cds)
cds_WI057_mmul <- mutate(cds_WI057_mmul, percent= (sum(cds_WI057_mmul$V2) - cumsum(V2) + V2) / mmul_cds)
cds_WI059_mmul <- mutate(cds_WI059_mmul, percent= (sum(cds_WI059_mmul$V2) - cumsum(V2) + V2) / mmul_cds)

cds_mmul_combined <- rbind(cds_WI055_mmul, cds_WI056_mmul, cds_WI057_mmul, cds_WI059_mmul)
head(cds_mmul_combined)

# sifaka hg38 cds
cds_f249_hg38 <- mutate(read.table("F249.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F249")
cds_f406_hg38 <- mutate(read.table("F406.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F406")
cds_m288_hg38 <- mutate(read.table("M288.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M288")
cds_m418_hg38 <- mutate(read.table("M418.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M418")

cds_f249_hg38 <- mutate(cds_f249_hg38, percent= (sum(cds_f249_hg38$V2) - cumsum(V2) + V2) / hg38_cds) # hg38_cds = 35825701
cds_f406_hg38 <- mutate(cds_f406_hg38, percent= (sum(cds_f406_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_m288_hg38 <- mutate(cds_m288_hg38, percent= (sum(cds_m288_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_m418_hg38 <- mutate(cds_m418_hg38, percent= (sum(cds_m418_hg38$V2) - cumsum(V2) + V2) / hg38_cds)

cds_sifaka_hg38_combined <- rbind(cds_f249_hg38, cds_f406_hg38, cds_m288_hg38, cds_m418_hg38)
head(cds_sifaka_hg38_combined)

# macaque hg38 cds
cds_WI055_hg38 <- mutate(read.table("WI055.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI055")
cds_WI056_hg38 <- mutate(read.table("WI056.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI056")
cds_WI057_hg38 <- mutate(read.table("WI057.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI057")
cds_WI059_hg38 <- mutate(read.table("WI059.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="WI059")

cds_WI055_hg38 <- mutate(cds_WI055_hg38, percent= (sum(cds_WI055_hg38$V2) - cumsum(V2) + V2) / hg38_cds) # hg38_cds = 35825701
cds_WI056_hg38 <- mutate(cds_WI056_hg38, percent= (sum(cds_WI056_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_WI057_hg38 <- mutate(cds_WI057_hg38, percent= (sum(cds_WI057_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_WI059_hg38 <- mutate(cds_WI059_hg38, percent= (sum(cds_WI059_hg38$V2) - cumsum(V2) + V2) / hg38_cds)

cds_macaque_hg38_combined <- rbind(cds_WI055_hg38, cds_WI056_hg38, cds_WI057_hg38, cds_WI059_hg38)
head(cds_macaque_hg38_combined)

# Print coverage values
print("Pcoq depth of 1 or greater")
filter(cds_pcoq_combined, V1 == 1)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 1)$percent)))

print("Pcoq depth of 4 or greater")
filter(cds_pcoq_combined, V1 == 4)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 4)$percent)))

print("Pcoq depth of 8 or greater")
filter(cds_pcoq_combined, V1 == 8)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 8)$percent)))

print("Pcoq depth of 12 or greater")
filter(cds_pcoq_combined, V1 == 12)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 12)$percent)))

print("Pcoq depth of 16 or greater")
filter(cds_pcoq_combined, V1 == 16)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 16)$percent)))

print("Pcoq depth of 20 or greater")
filter(cds_pcoq_combined, V1 == 20)
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_combined, V1 == 20)$percent)))


print("MMul depth of 1 or greater")
filter(cds_mmul_combined, V1 == 1)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 1)$percent)))

print("MMul depth of 4 or greater")
filter(cds_mmul_combined, V1 == 4)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 4)$percent)))

print("MMul depth of 8 or greater")
filter(cds_mmul_combined, V1 == 8)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 8)$percent)))

print("MMul depth of 12 or greater")
filter(cds_mmul_combined, V1 == 12)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 12)$percent)))

print("MMul depth of 16 or greater")
filter(cds_mmul_combined, V1 == 16)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 16)$percent)))

print("MMul depth of 20 or greater")
filter(cds_mmul_combined, V1 == 20)
cat(sprintf("Mean = %f", mean(filter(cds_mmul_combined, V1 == 20)$percent)))


print("Sifaka (hg38) depth of 1 or greater")
filter(cds_sifaka_hg38_combined, V1 == 1)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 1)$percent)))

print("Sifaka (hg38) depth of 4 or greater")
filter(cds_sifaka_hg38_combined, V1 == 4)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 4)$percent)))

print("Sifaka (hg38) depth of 8 or greater")
filter(cds_sifaka_hg38_combined, V1 == 8)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 8)$percent)))

print("Sifaka (hg38) depth of 12 or greater")
filter(cds_sifaka_hg38_combined, V1 == 12)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 12)$percent)))

print("Sifaka (hg38) depth of 16 or greater")
filter(cds_sifaka_hg38_combined, V1 == 16)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 16)$percent)))

print("Sifaka (hg38) depth of 20 or greater")
filter(cds_sifaka_hg38_combined, V1 == 20)
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_combined, V1 == 20)$percent)))


print("Macaque (hg38) depth of 1 or greater")
filter(cds_macaque_hg38_combined, V1 == 1)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 1)$percent)))

print("Macaque (hg38) depth of 4 or greater")
filter(cds_macaque_hg38_combined, V1 == 4)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 4)$percent)))

print("Macaque (hg38) depth of 8 or greater")
filter(cds_macaque_hg38_combined, V1 == 8)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 8)$percent)))

print("Macaque (hg38) depth of 12 or greater")
filter(cds_macaque_hg38_combined, V1 == 12)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 12)$percent)))

print("Macaque (hg38) depth of 16 or greater")
filter(cds_macaque_hg38_combined, V1 == 16)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 16)$percent)))

print("Macaque (hg38) depth of 20 or greater")
filter(cds_macaque_hg38_combined, V1 == 20)
cat(sprintf("Mean = %f", mean(filter(cds_macaque_hg38_combined, V1 == 20)$percent)))

#test larger font size size
pcoq_cds_plot <- ggplot(cds_pcoq_combined, aes(x = V1, y= percent, color= Sample)) + geom_line(size = .5, alpha = 0.75) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.928) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in Pcoq1 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Sifaka\nsample") + theme(legend.title.align=0.5)

sifaka_hg38_cds_plot <- ggplot(cds_sifaka_hg38_combined, aes(x = V1, y= percent, color= Sample)) + geom_line(size = .5, alpha = 0.75) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.928) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in hg38 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Sifaka\nsample") + theme(legend.title.align=0.5)
mmul_cds_plot <- ggplot(cds_mmul_combined, aes(x = V1, y= percent, color= Sample)) + geom_line(size = .5, alpha = 0.75) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.928) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in Mmul8 for macaques") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Macaque\nsample") + theme(legend.title.align=0.5)

macaque_hg38_cds_plot <- ggplot(cds_macaque_hg38_combined, aes(x = V1, y= percent, color= Sample)) + geom_line(size = .5, alpha = 0.75) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.928) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in hg38 for macaques") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Macaque\nsample") + theme(legend.title.align=0.5)

#extract sifaka sample legend so that only one legend for both sifaka graphs displays (uses cowplot feature get_legend)
sifaka_legend <- get_legend(pcoq_cds_plot)

#extract macaque sample legend so that only one legend for both macaque graphs displays
macaque_legend <- get_legend(mmul_cds_plot)

#test changing around spacing in combined plots
sifaka_cds_plot <- ggdraw() + draw_plot(plot_grid(pcoq_cds_plot + theme(legend.position = 'none'), sifaka_hg38_cds_plot + theme(legend.position = 'none'), ncol = 2, align = 'hv'), width = 0.91) + draw_plot(sifaka_legend, x = 0.94, y = 0.275, width = 0.025, height = 0.5)
macaque_cds_plot <- ggdraw() + draw_plot(plot_grid(mmul_cds_plot + theme(legend.position = 'none'), macaque_hg38_cds_plot + theme(legend.position = 'none'), ncol = 2, align = 'hv'), width = 0.91) + draw_plot(macaque_legend, x = 0.94, y = 0.275, width = 0.025, height = 0.5)

#combine sifaka and macaque graphs
combined_cds_plot <- multiplot(sifaka_cds_plot,macaque_cds_plot)

#save image
ggsave(filename="Cds_coverage_Nimblegen_multiplot.png", plot=multiplot(sifaka_cds_plot,macaque_cds_plot), width=9, height=8, dpi = 600)

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

## Cds plots for NimbleGen Sifaka (4) and IDT Sifaka (4) comparison

# Pcoq cds
cds_f249_pcoq <- mutate(read.table("F249.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F249", Probes="NimbleGen")
cds_f406_pcoq <- mutate(read.table("F406.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F406", Probes="NimbleGen")
cds_m288_pcoq <- mutate(read.table("M288.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M288", Probes="NimbleGen")
cds_m418_pcoq <- mutate(read.table("M418.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M418", Probes="NimbleGen")
cds_137_pcoq <- mutate(read.table("137.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="137", Probes="IDT")
cds_161_pcoq <- mutate(read.table("161.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="161", Probes="IDT")
cds_184_pcoq <- mutate(read.table("184.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="184", Probes="IDT")
cds_252_pcoq <- mutate(read.table("252.pcoq.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="252", Probes="IDT")

cds_f249_pcoq <- mutate(cds_f249_pcoq, percent= (sum(cds_f249_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_f406_pcoq <- mutate(cds_f406_pcoq, percent= (sum(cds_f406_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_m288_pcoq <- mutate(cds_m288_pcoq, percent= (sum(cds_m288_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_m418_pcoq <- mutate(cds_m418_pcoq, percent= (sum(cds_m418_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_137_pcoq <- mutate(cds_137_pcoq, percent= (sum(cds_137_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_161_pcoq <- mutate(cds_161_pcoq, percent= (sum(cds_161_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_184_pcoq <- mutate(cds_184_pcoq, percent= (sum(cds_184_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)
cds_252_pcoq <- mutate(cds_252_pcoq, percent= (sum(cds_252_pcoq$V2) - cumsum(V2) + V2) / pcoq_cds)

cds_pcoq_comparison <- rbind(cds_f249_pcoq, cds_f406_pcoq, cds_m288_pcoq, cds_m418_pcoq, cds_137_pcoq, cds_161_pcoq, cds_184_pcoq, cds_252_pcoq)
head(cds_pcoq_comparison)

# sifaka hg38 cds
cds_f249_hg38 <- mutate(read.table("F249.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F249", Probes="NimbleGen")
cds_f406_hg38 <- mutate(read.table("F406.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="F406", Probes="NimbleGen")
cds_m288_hg38 <- mutate(read.table("M288.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M288", Probes="NimbleGen")
cds_m418_hg38 <- mutate(read.table("M418.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="M418", Probes="NimbleGen")
cds_137_hg38 <- mutate(read.table("137.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="137", Probes="IDT")
cds_161_hg38 <- mutate(read.table("161.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="161", Probes="IDT")
cds_184_hg38 <- mutate(read.table("184.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="184", Probes="IDT")
cds_252_hg38 <- mutate(read.table("252.hg38.downsampled.mapq20_noDup.genome_cov.cds.hist", header=FALSE, sep="\t"), Sample="252", Probes="IDT")

cds_f249_hg38 <- mutate(cds_f249_hg38, percent= (sum(cds_f249_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_f406_hg38 <- mutate(cds_f406_hg38, percent= (sum(cds_f406_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_m288_hg38 <- mutate(cds_m288_hg38, percent= (sum(cds_m288_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_m418_hg38 <- mutate(cds_m418_hg38, percent= (sum(cds_m418_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_137_hg38 <- mutate(cds_137_hg38, percent= (sum(cds_137_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_161_hg38 <- mutate(cds_161_hg38, percent= (sum(cds_161_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_184_hg38 <- mutate(cds_184_hg38, percent= (sum(cds_184_hg38$V2) - cumsum(V2) + V2) / hg38_cds)
cds_252_hg38 <- mutate(cds_252_hg38, percent= (sum(cds_252_hg38$V2) - cumsum(V2) + V2) / hg38_cds)

cds_sifaka_hg38_comparison <- rbind(cds_f249_hg38, cds_f406_hg38, cds_m288_hg38, cds_m418_hg38, cds_137_hg38, cds_161_hg38, cds_184_hg38, cds_252_hg38)
head(cds_sifaka_hg38_comparison)

print("NimbleGen (pcoq) depth of 1 or greater")
filter(cds_pcoq_comparison, V1 == 1 & Probes == "NimbleGen")
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_comparison, V1 == 1 & Probes == "NimbleGen")$percent)))

print("iDT (pcoq) depth of 1 or greater")
filter(cds_pcoq_comparison, V1 == 1 & Probes == "IDT")
cat(sprintf("Mean = %f", mean(filter(cds_pcoq_comparison, V1 == 1 & Probes == "IDT")$percent)))

print("NimbleGen (hg38) depth of 1 or greater")
filter(cds_sifaka_hg38_comparison, V1 == 1 & Probes == "NimbleGen")
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_comparison, V1 == 1 & Probes == "NimbleGen")$percent)))

print("iDT (hg38) depth of 1 or greater")
filter(cds_sifaka_hg38_comparison, V1 == 1 & Probes == "IDT")
cat(sprintf("Mean = %f", mean(filter(cds_sifaka_hg38_comparison, V1 == 1 & Probes == "IDT")$percent)))

# Pcoq
wilcox.test(filter(cds_pcoq_comparison, Probes == "NimbleGen" & V1 < 50)$V2, filter(cds_pcoq_comparison, Probes == "IDT"&  V1 < 45)$V2)
wilcox.test(filter(cds_pcoq_comparison, Probes == "NimbleGen" & 50 < V1 & V1 < 100)$V2, filter(cds_pcoq_comparison, Probes == "IDT" & 50 < V1 & V1< 100)$V2)
cat(sprintf("NimbleGen mean (<50x coverage) = %f", mean(filter(cds_pcoq_comparison, Probes == "NimbleGen" & V1 < 50)$V2)))
cat(sprintf("iDT mean (<50x coverage) = %f", mean(filter(cds_pcoq_comparison, Probes == "IDT" &  V1 < 45)$V2)))
cat(sprintf("NimbleGen mean (>50x coverage) = %f", mean(filter(cds_pcoq_comparison, Probes == "NimbleGen" & 50 < V1 & V1 < 100)$V2)))
cat(sprintf("iDT mean (>50x coverage) = %f", mean(filter(cds_pcoq_comparison, Probes == "IDT" & 50 < V1 & V1 < 100)$V2)))


#plots
pcoq_cds_plot <- ggplot(cds_pcoq_comparison, aes(x = V1, y= percent, color= Sample)) + geom_line(aes(linetype=Probes), size = .5, alpha = 0.75) + scale_linetype_manual(values=c(1,2)) + scale_color_manual(values=c( "#F8766D" , "#7CAE00" , "#00BFC4" , "#C77CFF" , "#F8766D" , "#7CAE00" , "#00BFC4" , "#C77CFF")) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.909087) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in Pcoq1 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Sifaka\nsample") + theme(legend.title.align=0.5) + guides(color = guide_legend(override.aes = list(linetype=c('solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed', 'dashed'))))
sifaka_hg38_cds_plot <- ggplot(cds_sifaka_hg38_comparison, aes(x = V1, y= percent, color= Sample)) + geom_line(aes(linetype=Probes), size = .5, alpha = 0.75) + scale_linetype_manual(values=c(1,2)) + scale_color_manual(values=c( "#F8766D" , "#7CAE00" , "#00BFC4" , "#C77CFF" , "#F8766D" , "#7CAE00" , "#00BFC4" , "#C77CFF")) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0,100,10), minor_breaks = NULL) +  scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1.0,0.1), minor_breaks = NULL) + geom_vline(xintercept=4, linetype = "dotted") + geom_vline(xintercept=8, linetype = "dotted") + geom_vline(xintercept=12, linetype = "dotted") + geom_hline(yintercept=0.909087) + theme_bw() + xlab(label = "Depth of coverage") + theme(axis.title.x = element_text(margin = margin(t = 1), size = 12)) + ylab(label = "Proportion of cds in genome\nat X coverage or greater") + theme(axis.title.y = element_text(size = 12)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Callable sites in hg38 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.text=element_text(size=10),legend.title=element_text(size=11, face="bold")) + labs(color = "Sifaka\nsample") + theme(legend.title.align=0.5) + guides(color = guide_legend(override.aes = list(linetype=c('solid', 'solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed', 'dashed'))))


#get legend
sifaka_legend <- get_legend(pcoq_cds_plot)

#combine plots
sifaka_cds_plot <- ggdraw() + draw_plot(plot_grid(pcoq_cds_plot + theme(legend.position = 'none'), sifaka_hg38_cds_plot + theme(legend.position = 'none'), nrow = 2, align = 'v'), width = 0.80) + draw_plot(sifaka_legend, x = 0.895, y = 0.27, width = 0.025, height = 0.5)

#save plots - need to specify aspect ratio so it doesn't come out really long
ggsave(filename="Cds_coverage_sifaka_kit_comparison.png", plot=sifaka_cds_plot, width=5.5, height=7.18, dpi = 600)

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

## Dot/jitter plots for counts in various genomic regions

file.names <- dir(".", pattern=".hist")
coverage_df <- data.frame(matrix(nrow=length(file.names),ncol=7))
colnames(coverage_df) <- c("sample", "genome", "region", "depth1", "depth4", "depth8", "depth12")
for(i in 1:length(file.names)){
  split <- strsplit(file.names[i], "[.]")
  temp_table <- read.table(file.names[i], header=FALSE, sep="\t")
  n1 <- sum(temp_table$V2)
  n4 <- sum(filter(temp_table, V1 >= 4)$V2)
  n8 <- sum(filter(temp_table, V1 >= 8)$V2)
  n12 <- sum(filter(temp_table, V1 >= 12)$V2)
  coverage_df[i,] <- c(split[[1]][1], split[[1]][2], split[[1]][6], n1, n4, n8, n12)
}
coverage_df
head(coverage_df, n=40)

file.names <- dir(".", pattern=".hist")
coverage_df <- data.frame(matrix(nrow=length(file.names),ncol=5))
colnames(coverage_df) <- c("sample", "genome", "region", "depth", "count")
row_num <- 1
for(i in 1:length(file.names)){
  split <- strsplit(file.names[i], "[.]")
  if(split[[1]][3]=="unsampled") next
  temp_table <- read.table(file.names[i], header=FALSE, sep="\t")
  n1 <- sum(temp_table$V2)
  n4 <- sum(filter(temp_table, V1 >= 4)$V2)
  n8 <- sum(filter(temp_table, V1 >= 8)$V2)
  n12 <- sum(filter(temp_table, V1 >= 12)$V2)
  coverage_df[row_num,] <- c(split[[1]][1], split[[1]][2], split[[1]][6], 1, n1)
  row_num = row_num + 1
  coverage_df[row_num,] <- c(split[[1]][1], split[[1]][2], split[[1]][6], 4, n4)
  row_num = row_num + 1
  coverage_df[row_num,] <- c(split[[1]][1], split[[1]][2], split[[1]][6], 8, n8)
  row_num = row_num + 1
  coverage_df[row_num,] <- c(split[[1]][1], split[[1]][2], split[[1]][6], 12, n12)
  row_num = row_num + 1
}
coverage_df$depth <- as.factor(as.character(coverage_df$depth))
coverage_df$count <- as.numeric(as.character(coverage_df$count))

sifaka_pcoq_df <- filter(coverage_df, genome == "pcoq", region != "gene", region != "exon")
sifaka_hg38_df <- filter(coverage_df, genome == "hg38", sample == "F249" | sample == "F406" | sample == "M288" | sample == "M418", region != "gene", region != "exon")
macaque_mmul_df <- filter(coverage_df, genome == "mmul", region != "gene", region != "exon")
macaque_hg38_df <- filter(coverage_df, genome == "hg38", sample == "WI055" | sample == "WI056" | sample == "WI057" | sample == "WI059", region != "gene", region != "exon")

#To make the plots clustered by region rather than depth
sifaka_pcoq_df$depth <- factor(sifaka_pcoq_df$depth, levels = c(1,4,8,12))
sifaka_hg38_df$depth <- factor(sifaka_hg38_df$depth, levels = c(1,4,8,12))
macaque_mmul_df$depth <- factor(macaque_mmul_df$depth, levels = c(1,4,8,12))
macaque_hg38_df$depth <- factor(macaque_hg38_df$depth, levels = c(1,4,8,12))

#box plots
sifaka_pcoq.box <- ggplot(sifaka_pcoq_df, aes(x=region, y=count, color=depth, shape=depth)) + geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_dodge(0.8)) + scale_x_discrete(limits=c("cds", "intergenic", "intron", "utr")) + scale_y_continuous(breaks=c(0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)) + ylab("Number of sites (Mb)") + theme(axis.title.y = element_text(size = 12)) + background_grid(major="y", minor="y") + geom_vline(xintercept=c(1.5, 2.5, 3.5), color="gray") + ggtitle("Callable sites in Pcoq1 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.title=element_text(size=11, face="bold")) + labs(shape = "Minimum\ndepth") + labs(color = "Minimum\ndepth") + theme(legend.title.align=0.5) + theme(axis.text.y=element_text(size=10)) + theme(axis.title.x=element_blank())
sifaka_hg38.box <- ggplot(sifaka_hg38_df, aes(x=region, y=count, color=depth, shape=depth)) + geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_dodge(0.8)) + scale_x_discrete(limits=c("cds", "intergenic", "intron", "utr")) + scale_y_continuous(breaks=c(0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000, 130000000, 140000000, 150000000, 160000000, 170000000, 180000000, 190000000, 200000000), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)) + ylab("Number of sites (Mb)") + theme(axis.title.y = element_text(size = 12)) + background_grid(major="y", minor="y") + geom_vline(xintercept=c(1.5, 2.5, 3.5), color="gray")+ ggtitle("Callable sites in hg38 for sifakas") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.title=element_text(size=11, face="bold")) + labs(shape = "Minimum\ndepth") + labs(color = "Minimum\ndepth") + theme(legend.title.align=0.5) + theme(axis.text.y=element_text(size=10)) + theme(axis.title.x=element_blank())
macaque_mmul.box <- ggplot(macaque_mmul_df, aes(x=region, y=count, color=depth, shape=depth)) + geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_dodge(0.8)) + scale_x_discrete(limits=c("cds", "intergenic", "intron", "utr")) + scale_y_continuous(breaks=c(0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000, 130000000, 140000000, 150000000, 160000000, 170000000, 180000000, 190000000, 200000000), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)) + ylab("Number of sites (Mb)") + theme(axis.title.y = element_text(size = 12)) + background_grid(major="y", minor="y") + geom_vline(xintercept=c(1.5, 2.5, 3.5), color="gray")+ ggtitle("Callable sites in Mmul8 for macaques") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.title=element_text(size=11, face="bold")) + labs(shape = "Minimum\ndepth") + labs(color = "Minimum\ndepth") + theme(legend.title.align=0.5) + theme(axis.text.y=element_text(size=10)) + theme(axis.title.x=element_blank())
macaque_hg38.box <- ggplot(macaque_hg38_df, aes(x=region, y=count, color=depth, shape=depth)) + geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_dodge(0.8)) + scale_x_discrete(limits=c("cds", "intergenic", "intron", "utr")) + scale_y_continuous(breaks=c(0, 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000, 130000000, 140000000, 150000000, 160000000, 170000000, 180000000, 190000000, 200000000), labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)) + ylab("Number of sites (Mb)") + theme(axis.title.y = element_text(size = 12)) + background_grid(major="y", minor="y") + geom_vline(xintercept=c(1.5, 2.5, 3.5), color="gray")+ ggtitle("Callable sites in hg38 for macaques") + theme(plot.title = element_text(hjust = .5, size = 13, face="bold")) + theme(legend.title=element_text(size=11, face="bold")) + labs(shape = "Minimum\ndepth") + labs(color = "Minimum\ndepth") + theme(legend.title.align=0.5) + theme(axis.text.y=element_text(size=10)) + theme(axis.title.x=element_blank())

#extract legend for all
depth_legend <- get_legend(sifaka_pcoq.box)

#test changing around spacing in combined plots
sifakas_box_plot <- ggdraw() + draw_plot(plot_grid(sifaka_pcoq.box + theme(legend.position = 'none'), sifaka_hg38.box + theme(legend.position = 'none'), ncol = 2, align = 'hv'), width = 0.91) + draw_plot(depth_legend, x = 0.94, y = 0.275, width = 0.025, height = 0.5)
macaques_box_plot <- ggdraw() + draw_plot(plot_grid(macaque_mmul.box + theme(legend.position = 'none'), macaque_hg38.box + theme(legend.position = 'none'), ncol = 2, align = 'hv'), width = 0.91) + draw_plot(depth_legend, x = 0.94, y = 0.275, width = 0.025, height = 0.5)

#combine sifaka and macaque graphs
combined_cds_plot <- multiplot(sifakas_box_plot,macaques_box_plot)

#save image
ggsave(filename="combined_cov_boxplots.png", multiplot(sifakas_box_plot,macaques_box_plot), height = 8.7, width = 7.8, dpi = 600)
