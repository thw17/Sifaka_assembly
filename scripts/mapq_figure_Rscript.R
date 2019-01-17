# Set working directory
setwd()

# Set table name
input_table <- ""

# Imports
library("cowplot")
library("dplyr")

# Read table
full_df <- read.table(input_table, sep="\t", header=TRUE)

# Quick dplyr processing to isolate downsampled and add info for error bars
downsampled_df <- filter(full_df, Downsampled == "Yes")

#adds in mean minus stdev to the part of the dataframe as ymin, ymax set to 60
downsampled_df <- mutate(downsampled_df, ymin = Mean_MAPQ - MAPQ_SD, ymax = 60)

# Plot data
#p <- ggplot(downsampled_df, aes(x=Sample, y=Mean_MAPQ, fill=Reference)) + geom_bar(position="dodge", stat="identity") + geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(0.9), width=.2) + theme_bw() + theme(panel.grid.major.x =element_blank(), axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold"), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + labs(fill="Reference Genome") + scale_fill_discrete(labels=c("Human", "Rhesus", "Sifaka")) + ylab("MAPQ")
p <- ggplot(downsampled_df, aes(x=Sample, y=Mean_MAPQ, fill=Reference)) + geom_bar(position="dodge", stat="identity", width = .7) +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(0.7), width=.2) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x =element_blank(), axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold"), axis.text.x = element_text(size=12, angle=45, hjust=1), axis.text.y = element_text(size=12)) +
    labs(fill="Reference Genome") + scale_fill_discrete(labels=c("Human (hg38)", "Rhesus Macaque (Mmul8)", "Coquerel's Sifaka (Pcoq1)")) + ylab("\nMapping Quality\n") +
    xlab("\n\nSample") + theme(axis.title.x = element_text(vjust=4)) + theme(axis.title.y = element_text(vjust=4)) + (scale_y_continuous(expand = c(0, 0))) +
    geom_vline(xintercept=c(4.5, 8.5), linetype="dotted") + theme(legend.text=element_text(size=12),legend.title=element_text(size=14, face="bold"))
p
ggsave(filename="mapq_bargraph.pdf", p, height = 8.7, width = 10, dpi = 600)
