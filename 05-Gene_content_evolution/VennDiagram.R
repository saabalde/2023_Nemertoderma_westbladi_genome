##################
## Venn Diagram ##
##################

## Set working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/16-GeneFamily_evolution/06-Shared_orthogroups/03-PFAM/")


## Load the libraries
library(VennDiagram)
library(ggplot2)


## Load the data
Acoelomorpha <- read.table(file = "Orthogroups.Acoelomorpha.csv", header = TRUE, sep = ",")

# Extract the genes present in each group
Acoela <- Acoelomorpha[Acoelomorpha$Acoela..n...2. == 1, 1]
Nemertodermatida <- Acoelomorpha[Acoelomorpha$Nemertodermatida..n...1. == 1, 1]


## Make a Venn Diagram
venn.diagram(x = list(Acoela, Nemertodermatida), category.names = c("Acoela (n = 2)", "Nemertodermatida (n = 1)"),
             output = TRUE, 
             # Appearance
             lwd = 2 , # width of the circle's circumference 
             rotation.degree = 90,
             cex = 1, # Size of the numbers
             cat.cex = 1, # Size of the labels
             fontfamily = "sans",
             cat.default.pos = "outer", # Position of the labels
             cat.pos = c(180, 0), # Position of the labels around the circle in degrees
             cat.dist = c(0.015, 0.015), # Distance of the labels from the circumference

             # Colors
             col = c("#6CACE4", "#D50032"), # Circumferences
             fill = c(alpha("#6CACE4", 0.5), alpha("#D50032", 0.5)), # Inside circles (color, opacity)
             
             # Output image
             imagetype = "tiff",
             height = 1080,
             width = 1080,
             resolution = 300,
             filename= "Venn_diagram.tiff")

