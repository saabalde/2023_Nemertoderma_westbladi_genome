##################################################
## Generate stats for the ultrafiltration genes ##
##################################################

## Load libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(FSA)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/Ultrafiltration_genes/00-Filtered_hits/")


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

## For each gene, check the normality of the data
for(i in 1:length(genes)){
    subset <- Ultrafiltration_genes[grep(genes[i], Ultrafiltration_genes$Blast_modified), ]
    
    # Print the gene name
    print(genes[i])
    
    # Calculate the normality for protein length, number of exons, and average 
    # exon length. In that order
    print(shapiro.test(subset$Prot_length))
    print(shapiro.test(subset$N_exons))
    print(shapiro.test(subset$Av_exon_length))
    
    # Print a blank line to separate genes
    print("   ")
    print("   ")
}

## All but three of the metrics are not normal, so instead of running a for loop
## I will simply run the homocedasticity test over them. I will use the Barlett
## test, since it is robust if the data are normal and it allows comparisons 
## among more than two populations

bartlett.test(x = Ultrafiltration_genes_filtered$Prot_length[grep("POU3", Ultrafiltration_genes_filtered$Blast_modified, )],
              g = Ultrafiltration_genes_filtered$Main_clade[grep("POU3", Ultrafiltration_genes_filtered$Blast_modified, )])

bartlett.test(x = Ultrafiltration_genes_filtered$Av_exon_length[grep("Sall", Ultrafiltration_genes_filtered$Blast_modified, )],
              g = Ultrafiltration_genes_filtered$Main_clade[grep("Sall", Ultrafiltration_genes_filtered$Blast_modified, )])

bartlett.test(x = Ultrafiltration_genes_filtered$Av_exon_length[grep("Osr", Ultrafiltration_genes_filtered$Blast_modified, )],
              g = Ultrafiltration_genes_filtered$Main_clade[grep("Osr", Ultrafiltration_genes_filtered$Blast_modified, )])

## The three are also homocedastic (equal variances), so I can run an ANOVA on them.
## For all others, I will run a Kruskal-Wallis

for(i in 1:length(genes)){
    subset <- Ultrafiltration_genes[grep(genes[i], Ultrafiltration_genes$Blast_modified), ]
    
    # Print the gene name
    print(genes[i])
    
    # Calculate the normality for protein length, number of exons, and average 
    # exon length. In that order. Do this with an if statement to separate
    # the ANOVA from Kruskal-Wallis
    if(genes[i] == "POU3"){
        anova <- aov(Prot_length ~ Main_clade, data = subset)
        print(summary(anova))
        print(kruskal.test(N_exons ~ Main_clade, data = subset))
        print(kruskal.test(Av_exon_length ~ Main_clade, data = subset))
        
    } else if(genes[i] == "Sall" || genes[i] == "Osr"){
        print(kruskal.test(Prot_length ~ Main_clade, data = subset))
        print(kruskal.test(N_exons ~ Main_clade, data = subset))
        anova <- aov(Av_exon_length ~ Main_clade, data = subset)
        print(summary(anova))
    } else{
        print(kruskal.test(Prot_length ~ Main_clade, data = subset))
        print(kruskal.test(N_exons ~ Main_clade, data = subset))
        print(kruskal.test(Av_exon_length ~ Main_clade, data = subset))
    }

    # Print a blank line to separate genes
    print("   ")
    print("   ")
}

## There are several comparisons that have returned significant differences. Only
## two of them are normal (Sall and Osr for Average exon length), for which I will
## use Bonferroni For all others, I will use the Dunn test adjusted with the Holm
## method. Because only two samples need to be analysed with Bonferroni, I will
## run Dunn over all of them and then run Bonferroni individually

## Note: Osr is failing, so I skip it now. It only has two categories anyway:
##       Protostomia and Deuterostomia

for(i in c(1:4,6:length(genes))){
    subset <- Ultrafiltration_genes[grep(genes[i], Ultrafiltration_genes$Blast_modified), ]
    
    # Print the gene name
    print(genes[i])
    
    # Run multiple comparisons to check which pairs are significantly different
    print(dunnTest(Prot_length ~ Main_clade, data = subset))
    print(dunnTest(N_exons ~ Main_clade, data = subset))
    print(dunnTest(Av_exon_length ~ Main_clade, data = subset))
    
    # Print a blank line to separate genes
    print("   ")
    print("   ")
}
subset <- Ultrafiltration_genes[grep("Sall", Ultrafiltration_genes$Blast_modified), ]
pairwise.t.test(subset$Av_exon_length, subset$Main_clade, p.adjust.method = "bonferroni")

################################################################################

## Import the table summarising the signifcance results
Significance_summary <- read.table(file = "Ultrafiltration_genes.AllSpecies.Significances_Summary.txt", header = TRUE, sep = ",")
Significance_summary$Main_clade <- factor(Significance_summary$Main_clade, 
                                          levels = c("Cnidaria", "Deuterostomia", "Protostomia", "Acoelomorpha"))
Significance_summary$Blast_modified <- factor(Significance_summary$Blast_modified,
                                               levels = c("Nephrin", "Kirrel", "ZO1", 
                                                          "Eya", "Six", "POU3", "Sall", "Lhx", "Osr"))
Significance_summary$Metric <- factor(Significance_summary$Metric, 
                                      levels = c("Prot_length", "N_exons", "Av_exon_length"))

# Reduce the dataset by removing genes shorter than half the average length
genes <- unique(Significance_summary$Blast_modified)

Significance_summary_filtered <- data.frame()
for(i in 1:length(genes)){
    average <- mean(Significance_summary[grep(genes[i], Significance_summary$Blast_modified), "Prot_length"])
    
    Extract_rows <- Significance_summary[grep(genes[i], Significance_summary$Blast_modified), ]
    Remove_short_sequences <- Extract_rows[Extract_rows$Prot_length > (average / 2), ]
    Significance_summary_filtered <- rbind(Significance_summary_filtered, Remove_short_sequences)
}

################################################################################

## Set the names of the labels in the facet plot
Category <- c("Protein length (aa)", "Number of exons", "Av. exon length (bp)")
names(Category) <- c("Prot_length", "N_exons", "Av_exon_length")

## Plot the significance results
BarPlot.Significance <- ggplot(Significance_summary_filtered, aes(x = Blast_modified, y = Results, fill = Main_clade)) + 
                             geom_boxplot() + 
                             guides(fill=guide_legend(title="Clades")) + 
                             theme_light(base_size = 8) + 
                             scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                             ggtitle("Summary of the Statistically Significant Differences in Gene Metrics") + 
                             xlab("") +  ylab("") + 
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

BarPlot.Significance <- BarPlot.Significance + ggh4x::facet_grid2(~ Metric, scales = "free", space = "free_x",
                                                                  labeller = labeller(Metric = Category),
                                                                  independent = "y") +
                                               theme(strip.text.x = element_text(size = 15, color = "black", face = "bold"),
                                                     strip.background = element_rect(color = "black", fill = "white", size = 2))

################################################################################

####################
## Save all plots ##
####################

ggsave(BarPlot.Significance, file="03-Statistically_Different_Metrics.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 40, height = 20, units = "cm")
