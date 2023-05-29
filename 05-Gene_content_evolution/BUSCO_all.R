#############################################
## Compare the BUSCO Content among genomes ##
#############################################

## Load libraries
library(ggplot2)
library(reshape2)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/00-Figures/")


## Load the data
BUSCO_genes <- read.table(file = "BUSCO_all.txt", header = TRUE, sep = ",", dec = ".")
BUSCO_genes$Clade <- factor(BUSCO_genes$Clade, levels = c("Acoelomorpha", "Deuterostomia", "Protostomia", "Cnidaria"))
BUSCO_genes$BUSCO <- factor(BUSCO_genes$BUSCO, levels = c("Complete", "Fragmented", "Missing"))

################################################################################

########################
## Plot missing stats ##
########################

Barplot.BUSCO <- ggplot(BUSCO_genes, aes(x = Clade, y = Percentage, fill = BUSCO)) + 
                        geom_bar(stat = "identity", width = 0.8) + 
                        scale_y_continuous(labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
                                           breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
                        theme_light(base_size = 8) + 
                        scale_fill_manual(values = c("#56B4E9", "#F0E442", "#F04442")) + 
                        ggtitle("Percentage of BUSCO Genes Clade") + 
                        xlab("") +  ylab("% BUSCO Genes") + 
                        theme(plot.title = element_text(family = "sans", colour = "black", 
                                                        size = rel(2.2), face = "bold")) + 
                        theme(legend.position = "top", legend.justification='left', legend.title = element_blank()) + 
                        theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                        theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                        theme(panel.grid.minor = element_blank()) + 
                        theme(panel.grid.major = element_blank()) +
                        theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                        theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                        theme(axis.line = element_line(linewidth = 1, colour = "black")) + 
                        theme(axis.ticks.length = unit(.85, "cm")) + 
                        theme(axis.ticks.y = element_line(colour = "#222222")) + 
                        theme(axis.ticks.x = element_line(colour = "#222222")) + 
                        theme(axis.ticks.length = unit(0.4, "cm")) + 
                        theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                        guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                        guides(fill = guide_legend(nrow = 1, byrow = TRUE))

Barplot.BUSCO

################################################################################

####################
## Save all plots ##
####################

ggsave(Barplot.BUSCO, file="05-BUSCO_genes.tiff", device = "tiff", dpi = 300, 
       bg = NULL, width = 25, height = 20, units = "cm")

