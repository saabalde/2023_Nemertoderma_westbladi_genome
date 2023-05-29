########################################
## Compare Gene Content among genomes ##
########################################

## Load libraries
library(ggplot2)
library(reshape2)
library(UpSetR)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/16-GeneFamily_evolution/06-Shared_orthogroups/03-PFAM/")


## Load the data to a list
Gene_content <- read.table(file = "Orthogroups.UpSetPlot.csv", header = TRUE, sep = ",", 
                                dec = ".")
Gene_content.list <- list(`Cnidaria (n = 3)` = Gene_content$Orthogroup[Gene_content$Cnidaria..n...3. == 1], 
                          `Acoelomorpha (n = 3)` = Gene_content$Orthogroup[Gene_content$Acoelomorpha..n...3. == 1],
                          `Protostomia (n = 8)` = Gene_content$Orthogroup[Gene_content$Protostomia..n...8. == 1],
                          `Deuterostomia (n = 4)` = Gene_content$Orthogroup[Gene_content$Deuterostomia..n...4. == 1])

################################################################################

################
## UpSet plot ##
################

## Make the plot
UpSet.plot <- upset(fromList(Gene_content.list), nsets = 4, order.by = "degree", decreasing = FALSE, 
                    mb.ratio = c(0.7, 0.3),
                    number.angles = 0, scale.intersections = "identity",
                    point.size = 5, line.size = 1, shade.alpha = 0.4,
                    text.scale = c(1.5, 1.25, 1.25, 1, 2, 1.5),
                    mainbar.y.label = "Orthogroups Per Intersection", sets.x.label = "Orthogroups Per Group", 
                    matrix.color = "grey24", main.bar.color = "grey24", sets.bar.color = "grey24")
UpSet.plot

################################################################################

####################
## Save all plots ##
####################

tiff("UpSet_plot.tiff", units="px", width=3590, height = 2286, res =300)
UpSet.plot <- upset(fromList(Gene_content.list), nsets = 4, order.by = "degree", decreasing = FALSE, 
                    mb.ratio = c(0.7, 0.3),
                    number.angles = 0, scale.intersections = "identity",
                    point.size = 5, line.size = 1, shade.alpha = 0.4,
                    text.scale = c(1.5, 1.25, 1.25, 1, 2, 1.5),
                    mainbar.y.label = "Orthogroups Per Intersection", sets.x.label = "Orthogroups Per Group", 
                    matrix.color = "grey24", main.bar.color = "grey24", sets.bar.color = "grey24")
UpSet.plot
dev.off()
