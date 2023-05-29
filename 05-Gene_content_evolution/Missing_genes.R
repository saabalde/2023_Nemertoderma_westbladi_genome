########################################
## Compare Gene Content among genomes ##
########################################

## Load libraries
library(ggplot2)
library(reshape2)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/16-GeneFamily_evolution/06-Shared_orthogroups/03-PFAM/")


## Load the data
Missing_genes <- read.table(file = "Orthogroups.SharedGenes.csv", header = TRUE, sep = ",", dec = ".")
Missing_genes$Taxa <- factor(Missing_genes$Taxa, levels = c("Acoelomorpha", "Deuterostomia", "Protostomia"))
Missing_genes$Group <- factor(Missing_genes$Group, levels = c("Metazoa", "Bilateria"))

################################################################################

########################
## Plot missing stats ##
########################

Barplot.Missing <- ggplot(Missing_genes, aes(x = Taxa, y = P_Missing, fill = Group)) + 
                          geom_bar(stat = "identity", width = 0.8, position = position_dodge2()) + 
                          scale_y_continuous(labels = c("0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60"), 
                                             breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)) + 
                          theme_light(base_size = 8) + 
                          scale_fill_manual(values = c("grey24", "grey54")) + 
                          ggtitle("Percentage of Genes Missing per Taxa and Higher Clade") + 
                          xlab("") +  ylab("% Missing Genes") + 
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

Barplot.Missing

################################################################################

####################
## Save all plots ##
####################

ggsave(Barplot.Missing, file="Missing_genes.tiff", device = "tiff", dpi = 300, 
       bg = NULL, width = 3590, units = "px")

