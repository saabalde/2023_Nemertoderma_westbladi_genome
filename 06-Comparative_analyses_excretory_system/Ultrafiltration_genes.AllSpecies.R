##################################################
## Generate stats for the ultrafiltration genes ##
##################################################

## Load libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(FSA)


## Set the working directory, if necessary
#setwd("")


## Load the data
Ultrafiltration_genes <- read.table(file = "Ultrafiltration_genes.AllSpecies.txt", header = TRUE, sep = ",")
Ultrafiltration_genes$Main_clade <- factor(Ultrafiltration_genes$Main_clade, 
                                           levels = c("Cnidaria", "Deuterostomia", "Protostomia", "Acoelomorpha"))
Ultrafiltration_genes$Blast_modified <- factor(Ultrafiltration_genes$Blast_modified,
                                               levels = c("Nephrin", "Kirrel", "ZO1", 
                                                          "Eya", "Six", "POU3", "Sall", "Lhx", "Osr"))

# Reduce the dataset by removing genes shorter than half the average length
genes <- unique(Ultrafiltration_genes$Blast_modified)

Ultrafiltration_genes_filtered <- data.frame()
for(i in 1:length(genes)){
    average <- mean(Ultrafiltration_genes[grep(genes[i], Ultrafiltration_genes$Blast_modified), "Prot_length"])
    
    Extract_rows <- Ultrafiltration_genes[grep(genes[i], Ultrafiltration_genes$Blast_modified), ]
    Remove_short_sequences <- Extract_rows[Extract_rows$Prot_length > (average / 2), ]
    Ultrafiltration_genes_filtered <- rbind(Ultrafiltration_genes_filtered, Remove_short_sequences)
}

################################################################################

## Set the names of the labels in the facet plot
Category <- c("Structural proteins", "Transcription factors")
names(Category) <- c("Structural_protein", "Transcription_factor")

## Protein length
BarPlot.GeneLength <- ggplot(Ultrafiltration_genes_filtered, aes(x = Blast_modified, y = Prot_length, fill = Main_clade)) + 
                             geom_boxplot() + 
                             scale_y_continuous(limits = c(0, 2000), n.breaks = 10) + 
                             guides(fill=guide_legend(title="Clades")) + 
                             theme_light(base_size = 8) + 
                             scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                             ggtitle("Protein Length per Gene") + 
                             xlab("") +  ylab("Protein length (aa)") + 
                             theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                             theme(legend.position = "top",legend.title = element_blank()) + 
                             theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                             theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                             theme(panel.grid.minor = element_blank()) + 
                             theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                             theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                             theme(axis.line = element_line(linewidth = 1, colour = "black")) + 
                             theme(axis.ticks.length = unit(0, "cm")) + 
                             theme(axis.ticks.y = element_line(colour = "#222222")) + 
                             theme(axis.ticks.x = element_line(colour = "#222222")) + 
                             theme(axis.ticks.length = unit(0.4, "cm")) + 
                             theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                             theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                             guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                             guides(fill = guide_legend(nrow = 1, byrow = TRUE))

BarPlot.GeneLength <- BarPlot.GeneLength + facet_grid(~ Structural_or_Transcription, scales = "free_x", space = "free_x", 
                                                      labeller = labeller(Structural_or_Transcription = Category)) +
                                           theme(strip.text.x = element_text(size = 15, color = "black", face = "bold"),
                                                 strip.background = element_rect(color = "black", fill = "white", size = 2)) 

BarPlot.GeneLength
                     

## Exons per gene
BarPlot.ExonsPerGene <- ggplot(Ultrafiltration_genes_filtered, aes(x = Blast_modified, y = N_exons, fill = Main_clade)) + 
                               geom_boxplot() + 
                               scale_y_continuous(limits = c(0, 100), n.breaks = 10) + 
                               guides(fill=guide_legend(title="Clades")) + 
                               theme_light(base_size = 8) + 
                               scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                               ggtitle("Number of Exons per Gene") + 
                               xlab("") +  ylab("Number of exons") + 
                               theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                               theme(legend.position = "top",legend.title = element_blank()) + 
                               theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                               theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                               theme(panel.grid.minor = element_blank()) + 
                               theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                               theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                               theme(axis.line = element_line(linewidth = 1, colour = "black")) + 
                               theme(axis.ticks.length = unit(0, "cm")) + 
                               theme(axis.ticks.y = element_line(colour = "#222222")) + 
                               theme(axis.ticks.x = element_line(colour = "#222222")) + 
                               theme(axis.ticks.length = unit(0.4, "cm")) + 
                               theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                               theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                               guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                               guides(fill = guide_legend(nrow = 1, byrow = TRUE))

BarPlot.ExonsPerGene <- BarPlot.ExonsPerGene + facet_grid(~ Structural_or_Transcription, scales = "free_x", space = "free_x", 
                                                          labeller = labeller(Structural_or_Transcription = Category)) +
                                               theme(strip.text.x = element_text(size = 15, color = "black", face = "bold"),
                                                     strip.background = element_rect(color = "black", fill = "white", size = 2)) 

BarPlot.ExonsPerGene


## Exon length
BarPlot.ExonLength <- ggplot(Ultrafiltration_genes_filtered, aes(x = Blast_modified, y = Av_exon_length, fill = Main_clade)) + 
                             geom_boxplot() + 
                             scale_y_continuous(limits = c(0, 3500), n.breaks = 10) + 
                             guides(fill=guide_legend(title="Clades")) + 
                             theme_light(base_size = 8) + 
                             scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                             ggtitle("Average Exon Length per Gene") + 
                             xlab("") +  ylab("Exon length (bp)") + 
                             theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                             theme(legend.position = "top",legend.title = element_blank()) + 
                             theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                             theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                             theme(panel.grid.minor = element_blank()) + 
                             theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                             theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                             theme(axis.line = element_line(linewidth = 1, colour = "black")) + 
                             theme(axis.ticks.length = unit(0, "cm")) + 
                             theme(axis.ticks.y = element_line(colour = "#222222")) + 
                             theme(axis.ticks.x = element_line(colour = "#222222")) + 
                             theme(axis.ticks.length = unit(0.4, "cm")) + 
                             theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                             theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                             guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                             guides(fill = guide_legend(nrow = 1, byrow = TRUE))

BarPlot.ExonLength <- BarPlot.ExonLength + facet_grid(~ Structural_or_Transcription, scales = "free_x", space = "free_x", 
                                                      labeller = labeller(Structural_or_Transcription = Category)) +
                                           theme(strip.text.x = element_text(size = 15, color = "black", face = "bold"),
                                                 strip.background = element_rect(color = "black", fill = "white", size = 2)) 

BarPlot.ExonLength

################################################################################

####################
## Save all plots ##
####################

ggsave(BarPlot.GeneLength, file="00-Gene_length.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(BarPlot.ExonsPerGene, file="01-Exons_per_gene.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(BarPlot.ExonLength, file="02-Average_exon_length.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
