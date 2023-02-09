egg_df <- read.csv("ecoli.emapper.annotations", header= T, sep = '\t')
egg_new <- snakemake.input[1]
egg_new <- read.csv("bsubtilis.emapper.annotations", header= T, sep = '\t')
#read in subset.csv 
subset_df <- read.csv("subset_vchol", header=TRUE)

colnames(egg_df) <- c('query_name','seed_eggNOG_ortholog','evalue','score',
                      'OGs','max_annot','COGs','Description','Preferred_name','GOs','EC',
                      'KEGG_KO','KEGG_pathway','KEGG_module','KEGG_reaction','KEGG_rclass','BRITE',
                      'KEGG_TC','CAZy','BIGG_reaction','PFAMs')

library(ggplot2)

#log the data
log_df <- -log10(egg_df$evalue)
log_df <- -log10(meth$evalue)

#plot1
plot1 <- ggplot(data= egg_df, aes(x=log_df)) + geom_histogram(fill='dodgerblue4', bins=15, col='white')
plot1 + xlab("E -value (-log10)") + ylab("Frequency") + theme(axis.title = element_text(size = 14)) + 
  theme(title =element_text(size=18)) + ggtitle('Distribution of E-value of EggNOG matches in B.subtilis') + theme_light()
ggsave("plot1bsub.png")

ggplot(data= meth, aes(x=log_df)) + geom_histogram(fill='dodgerblue4', bins=10, col='white')

#plot2
plot2 <- ggplot(egg_df, aes(x=(max_annot_lvl), y=(-log10(evalue)), col=max_annot_lvl)) + geom_boxplot() + coord_flip() + theme_classic() +
  ylab("E-value (-log10)") + xlab("Taxonomy level of gene") + ggtitle("Comparing taxonomy level with E-value B.subtilis")
plot2 + theme(legend.position="none") + theme(axis.text = element_text(size = 6)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 8)) + theme(title =element_text(size=8))

ggplot(data= egg_df, aes(x=max_annot_lvl, fill= max_annot_lvl)) + geom_bar() + coord_flip() + theme(legend.position = "none")
ggsave("plotbsub2.png")

#plot3
methanoegg <- egg_df[grep("methano", egg_df$max_annot_lvl),]
plot3 <- ggplot(methanoegg, aes(x=(max_annot_lvl), y=(-log10(evalue)), fill=max_annot_lvl, col=max_annot_lvl)) + geom_boxplot() + coord_flip() + theme_classic() +
  ylab("E-value (-log10)") + xlab("Taxonomy level of gene") + ggtitle("Figure  . Comparing 'Methano-' taxonomy levels with E-value")
plot3 +theme(legend.position="none") + scale_fill_manual(values = c( "lightskyblue2",  "lightskyblue2",  "lightskyblue2"))+
  scale_colour_manual(values = c("deepskyblue4", "deepskyblue4", "deepskyblue4")) + theme(text = element_text(size= 12)) +
  theme(axis.text = element_text(size = 14))
ggsave("plotbsub3.png")

#plot4
plot4 <- ggplot(methanoegg, aes(x=(-log10(evalue)), y=(max_annot_lvl), col=max_annot_lvl,
                                fill=max_annot_lvl)) + geom_violin() + coord_flip() + theme_classic() +
  ylab("Taxonomy level of gene") + xlab("E-value (-log10)") + ggtitle("Figure  . Illustrating 'Methano-' taxonomy levels compared with E-value")
plot4 + scale_colour_manual(values = c("#8c96c6", "#8856a7", "#810f7c"))  + geom_jitter(size=1) + scale_fill_manual(values = c("lavender", "mediumpurple1", "mediumorchid1")) +
  theme(axis.text = element_text(size = 14)) + theme(legend.position = "none") +
  theme(axis.title = element_text(size = 14)) + theme(title =element_text(size=12))
ggsave("plotbsub4.png")

#extracting genes from methanomicrobia taxonomy level
methanomicrobia <- methanoegg[grep("Methanomicrobia", methanoegg$best_tax_level),]

#plot5
plot5 <- ggplot(methanomicrobia, aes(x=(-log10(evalue)), y=(max_annot_lvl), col=max_annot_lvl)) + geom_jitter(size=1.2) + coord_flip() + theme_classic() +
  ylab("E-value (-log10)") + xlab("Taxonomy level of gene") + ggtitle("Figure  . Extracting genes of interest based on E-value") + xlim(c(0,300)) + theme(legend.position = "none") + scale_colour_manual(values = "#810f7c") +
  theme(text = element_text(size= 16)) + theme(title=element_text(size=16))
plot5
ggsave("plotbsub5.png")

#plot6
plot6 <- ggplot(data= egg_df, aes(COG_category)) + geom_bar(fill='coral2', col='white')
plot6 + xlab("Functional COG") + ylab("Frequency") +
  ggtitle("Figure  . EggNOG Functional COGs") + theme_light() + coord_flip()
ggsave("plotbsub6.png")


## confidence rank plot (working )
x <- subset_df$confidence_value 
as.factor(x)

conplot <- ggplot(data= subset_df, aes(x = x, fill= confidence_value)) +
  geom_bar(col= 'grey', fill= c('#282AD0', '#5B9EF8','#BCF5FD','#FFFDC6','#F4B07B','#C6393A')) +  
  theme_light() + coord_flip() + ggtitle("(e) Confidence Value Output Vibrio cholerae") + ylab("Frequency") + xlab("Confidence Value") + scale_x_continuous(breaks= c(1:6), n.breaks=6)
conplot 
#install.packages('colorBlindness')
library(colorBlindness)
cvdPlot(conplot)


conplot <- ggplot(data= subset_df, aes(x = confidence_value, fill= x)) +
  geom_bar(fill= c('darkolivegreen2','gold','darkorange1','firebrick3', 'hotpink2','darkorchid3')) + 
  theme_light() + coord_flip() + ggtitle("Confidence Value Output ") + ylab("Frequency") + xlab("Confidence Value")
conplot

displayAvailablePalette(color="white")


plot_con <- ggplot(data= subset_df, aes(x = confidence_value, fill= x)) + stat_count(geom = "bar", aes(fill= confidence_value))
plot_con + coord_flip() + xlab("Confidence Value") + ylab("Frequency") + scale_fill_manual(c("red", "orange", "yellow", "green", "blue"))
con <- ggplot(subset_df, aes(x = confidence_value)) + geom_histogram(bins =12) + coord_flip()
con + scale_fill_manual(c("red", "orange", "yellow", "green", "blue"))

conplot <- ggplot(data= subset_df, aes(x = confidence_value, fill= x)) +
  geom_bar(fill= c('darkolivegreen2','khaki2','darkorange1','indianred2')) + theme_light()
conplot + coord_flip() + ggtitle("Confidence Value Output ") + ylab("Frequency") + xlab("Confidence Value")


ggsave("plotbsubcon.png")
ggplot(subse, aes(x = class)) + geom_bar(aes(fill = drv))
ggsave("plot6.png")
