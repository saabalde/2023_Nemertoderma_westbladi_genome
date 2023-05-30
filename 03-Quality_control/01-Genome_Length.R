############################
## Compare genome lengths ##
############################

## Load libraries
library(ggplot2)
library(reshape2)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/00-Figures/")


## Load the data
Genomes <- read.table(file = "Genome_stats.csv", header = TRUE, sep = ",", dec = ".", na.strings = "")
Genomes$Genome <- factor(Genomes$Genome, levels = c("Illumina", "HiFi", "Pnaikaiensis", "Sroscoffensis"))
Genomes$Category <- factor(Genomes$Category, levels = c("Raw Assembly", "Decontaminated"))

################################################################################

#########################
## Plot genome lengths ##
#########################

Lineplot.Lengths <- ggplot(Genomes, aes(x = Index, y = Length, col = Genome, linetype = Category)) + 
                           geom_line(size = 1.5) + 
                           scale_y_continuous(limits = c(0, 1150), breaks = seq(0, 1150, 50)) +
                           scale_x_continuous(labels = c("0", "5,000", "10,000", "15,000", "20,000", "25,000", "30,000"), 
                                              breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) + 
                           theme_light(base_size = 8) + 
                           scale_color_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                           scale_linetype_manual(values=c("solid", "twodash")) +
                           ggtitle("Cumulative Genome Length") + 
                           xlab("Contig index") +  ylab("Genome Length (Mbps)") + 
                           theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                           theme(legend.position = "top",legend.title = element_blank()) + 
                           theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                           theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                           theme(panel.grid.minor = element_blank()) + 
                           theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                           theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                           theme(axis.line = element_line(linewidth = 1, colour = "black")) + 
                           theme(axis.ticks.length = unit(.85, "cm")) + 
                           theme(axis.ticks.y = element_line(colour = "#222222")) + 
                           theme(axis.ticks.x = element_line(colour = "#222222")) + 
                           theme(axis.ticks.length = unit(0.4, "cm")) + 
                           theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                           theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                           guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                           guides(fill = guide_legend(nrow = 1, byrow = TRUE))
Lineplot.Lengths

################################################################################

####################
## Save all plots ##
####################

ggsave(Lineplot.Lengths, file="01-Cumulative_genome_length.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")

