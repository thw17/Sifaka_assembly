# Set working directory
setwd()

# Libraries not needed for analysis
# library("ggplot2")
# library("dplyr")
# library("cowplot")


# Get NimbleGen Data

f249 <- mutate(read.table("F249.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
    header=FALSE, sep="\t"), Sample="F249", Probe = "NimbleGen")

f406 <- mutate(read.table("F406.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="F406", Probe = "NimbleGen")

m288 <- mutate(read.table("M288.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="M288", Probe = "NimbleGen")

m418 <- mutate(read.table("M418.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="M418", Probe = "NimbleGen")

nimblegen_pcoq_df <- rbind(f249, f406, m288, m418)

mean(f249$V7)
mean(f406$V7)
mean(m288$V7)
mean(m418$V7)

head(table(f249$V7))
head(table(f406$V7))
head(table(m288$V7))
head(table(m418$V7))

head(table(nimblegen_pcoq_df$V7))

# Get IDT Data

s137 <- mutate(read.table("137.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="137", Probe = "iDT")

s161 <- mutate(read.table("161.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="161", Probe = "iDT")

s184 <- mutate(read.table("184.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="184", Probe = "iDT")

s252 <- mutate(read.table("252.pcoq.downsampled.mapq20_noDup.genome_cov.cdscoverage.bed",
                          header=FALSE, sep="\t"), Sample="252", Probe = "iDT")

idt_pcoq_df <- rbind(s137, s161, s184, s252)

mean(s137$V7)
mean(s161$V7)
mean(s184$V7)
mean(s252$V7)

head(table(s137$V7))
head(table(s161$V7))
head(table(s184$V7))
head(table(s252$V7))

head(table(idt_pcoq_df$V7))

# Quick Stats
combined_df <- rbind(nimblegen_pcoq_df, idt_pcoq_df)
mean(nimblegen_pcoq_df$V7)
mean(idt_pcoq_df$V7)
wilcox.test(V7 ~ Probe, data = combined_df)

combined_df_no_zero <- filter(combined_df, V7 > 0)
mean(filter(nimblegen_pcoq_df, V7 > 0)$V7)
mean(filter(idt_pcoq_df, V7 > 0)$V7)

wilcox.test(V7 ~ Probe, data = combined_df_no_zero)

# # Plot data
# ggplot(filter(idt_pcoq_df, Sample == "252"), aes(V7, color = Sample)) +
#   geom_freqpoly(bins = 50)
