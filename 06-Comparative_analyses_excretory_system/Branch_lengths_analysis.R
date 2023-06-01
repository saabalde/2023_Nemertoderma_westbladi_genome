################################################
## Extract the branch length of each sequence ##
## and compare them among phyla and genes     ##
################################################


## Define my working environment, if necessary
#setwd("")


## Load modules
library(ggtree)
library(phytools)
library(TreeTools)
library(ggplot2)
library(FSA)


## Load the trees
Nephrin <- read.tree(file = "Nephrin.fas.treefile")
Kirrel <- read.tree(file = "Kirrel.fas.treefile")
ZO1 <- read.tree(file = "ZO1.fas.treefile")
Eya <- read.tree(file = "Eya.fas.treefile")
Lhx <- read.tree(file = "Lhx.fas.treefile")
Osr <- read.tree(file = "Osr.fas.treefile")
POU3 <- read.tree(file = "POU3.fas.treefile")
Sall <- read.tree(file = "Sall.fas.treefile")
Six <- read.tree(file = "Six.fas.treefile")


## Useful variables
phyla <- c("Annelida", "Brachiopoda", "Bryozoa", "Chordata", "Cnidaria", "Echinodermata", "Mollusca", "Xenacoelomorpha")
genes <- c("Nephrin", "Kirrel", "ZO1", "Eya", "Lhx", "Osr", "POU3", "Sall", "Six")
clades <- c("Acoelomorpha", "Cnidaria", "Deuterostomia", "Protostomia")


## Collapse all nodes
# Nephrin
Nephrin_collapsed <- CollapseNode(Nephrin, c((length(Nephrin$tip.label) + 2):(length(Nephrin$edge.length) + 1)))
plot(Nephrin_collapsed)

# Kirrel
Kirrel_collapsed <- CollapseNode(Kirrel, c((length(Kirrel$tip.label) + 2):(length(Kirrel$edge.length) + 1)))
plot(Kirrel_collapsed)

# ZO1
ZO1_collapsed <- CollapseNode(ZO1, c((length(ZO1$tip.label) + 2):(length(ZO1$edge.length) + 1)))
plot(ZO1_collapsed)

# Eya
Eya_collapsed <- CollapseNode(Eya, c((length(Eya$tip.label) + 2):(length(Eya$edge.length) + 1)))
plot(Eya_collapsed)

# Lhx
Lhx_collapsed <- CollapseNode(Lhx, c((length(Lhx$tip.label) + 2):(length(Lhx$edge.length) + 1)))
plot(Lhx_collapsed)

# Osr
Osr_collapsed <- CollapseNode(Osr, c((length(Osr$tip.label) + 2):(length(Osr$edge.length) + 1)))
plot(Osr_collapsed)

# POU3
POU3_collapsed <- CollapseNode(POU3, c((length(POU3$tip.label) + 2):(length(POU3$edge.length) + 1)))
plot(POU3_collapsed)

# Sall
Sall_collapsed <- CollapseNode(Sall, c((length(Sall$tip.label) + 2):(length(Sall$edge.length) + 1)))
plot(Sall_collapsed)

# Six
Six_collapsed <- CollapseNode(Six, c((length(Six$tip.label) + 2):(length(Six$edge.length) + 1)))
plot(Six_collapsed)


## Extract the branch length of each sequence
# Make one dataframe per gene
Nephrin_BranchLengths <- data.frame(Gene = rep("Nephrin", times = length(Nephrin_collapsed$tip.label)),
                                    Sequence = Nephrin_collapsed$tip.label,
                                    Branch_length = Nephrin_collapsed$edge.length)
Kirrel_BranchLengths <- data.frame(Gene = rep("Kirrel", times = length(Kirrel_collapsed$tip.label)),
                                    Sequence = Kirrel_collapsed$tip.label,
                                    Branch_length = Kirrel_collapsed$edge.length)
ZO1_BranchLengths <- data.frame(Gene = rep("ZO1", times = length(ZO1_collapsed$tip.label)),
                                    Sequence = ZO1_collapsed$tip.label,
                                    Branch_length = ZO1_collapsed$edge.length)
Eya_BranchLengths <- data.frame(Gene = rep("Eya", times = length(Eya_collapsed$tip.label)),
                                    Sequence = Eya_collapsed$tip.label,
                                    Branch_length = Eya_collapsed$edge.length)
Lhx_BranchLengths <- data.frame(Gene = rep("Lhx", times = length(Lhx_collapsed$tip.label)),
                                    Sequence = Lhx_collapsed$tip.label,
                                    Branch_length = Lhx_collapsed$edge.length)
Osr_BranchLengths <- data.frame(Gene = rep("Osr", times = length(Osr_collapsed$tip.label)),
                                    Sequence = Osr_collapsed$tip.label,
                                    Branch_length = Osr_collapsed$edge.length)
POU3_BranchLengths <- data.frame(Gene = rep("POU3", times = length(POU3_collapsed$tip.label)),
                                    Sequence = POU3_collapsed$tip.label,
                                    Branch_length = POU3_collapsed$edge.length)
Sall_BranchLengths <- data.frame(Gene = rep("Sall", times = length(Sall_collapsed$tip.label)),
                                    Sequence = Sall_collapsed$tip.label,
                                    Branch_length = Sall_collapsed$edge.length)
Six_BranchLengths <- data.frame(Gene = rep("Six", times = length(Six_collapsed$tip.label)),
                                    Sequence = Six_collapsed$tip.label,
                                    Branch_length = Six_collapsed$edge.length)

# Merge them all
Ultrafiltration_BranchLengths <- rbind(Nephrin_BranchLengths, Kirrel_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, ZO1_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, Eya_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, Lhx_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, Osr_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, POU3_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, Sall_BranchLengths)
Ultrafiltration_BranchLengths <- rbind(Ultrafiltration_BranchLengths, Six_BranchLengths)

# Remove the GenBank sequences
Ultrafiltration_BranchLengths <- Ultrafiltration_BranchLengths[-grep("QRF7", Ultrafiltration_BranchLengths$Sequence), ]

# Add the clade and phylum information
dummy_dataframe <- data.frame(Main_clade = rep(NA, times = nrow(Ultrafiltration_BranchLengths)),
                              Phylum = rep(NA, times = nrow(Ultrafiltration_BranchLengths)))
Ultrafiltration_BranchLengths <- cbind(dummy_dataframe, Ultrafiltration_BranchLengths)

for(phylum in 1:length(phyla)){
    Ultrafiltration_BranchLengths$Phylum[grep(phyla[phylum], Ultrafiltration_BranchLengths$Sequence)] <- phyla[phylum]
    
    if(phyla[phylum] == "Annelida" || phyla[phylum] == "Brachiopoda" || phyla[phylum] == "Bryozoa" || phyla[phylum] == "Mollusca"){
        Ultrafiltration_BranchLengths$Main_clade[grep(phyla[phylum], Ultrafiltration_BranchLengths$Sequence)] <- "Protostomia"
    } else if(phyla[phylum] == "Chordata" || phyla[phylum] == "Echinodermata"){
        Ultrafiltration_BranchLengths$Main_clade[grep(phyla[phylum], Ultrafiltration_BranchLengths$Sequence)] <- "Deuterostomia"
    } else if(phyla[phylum] == "Cnidaria"){
        Ultrafiltration_BranchLengths$Main_clade[grep(phyla[phylum], Ultrafiltration_BranchLengths$Sequence)] <- "Cnidaria"
    } else if(phyla[phylum] == "Xenacoelomorpha"){
        Ultrafiltration_BranchLengths$Main_clade[grep(phyla[phylum], Ultrafiltration_BranchLengths$Sequence)] <- "Acoelomorpha"
    }
}

# Assign factors to the dataframe
Ultrafiltration_BranchLengths$Main_clade <- factor(Ultrafiltration_BranchLengths$Main_clade, 
                                                   levels = c("Cnidaria", "Deuterostomia", "Protostomia", "Acoelomorpha"))
Ultrafiltration_BranchLengths$Gene <- factor(Ultrafiltration_BranchLengths$Gene, 
                                             levels = c("Nephrin", "Kirrel", "ZO1", 
                                                        "Eya", "Lhx", "Osr", "POU3", "Sall", "Six"))


## Write this data frame to a file
write.table(Ultrafiltration_BranchLengths, file = "Branch_lengths_analysis.txt", 
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
#Ultrafiltration_BranchLengths <- read.table(file = "Branch_lengths_analysis.txt", header = TRUE, sep = ",")


## Calculate the average branch length per clade and gene
Ultrafiltration_BranchLengths_Average <- data.frame(Main_clade = rep(NA, times = (length(clades) * length(genes))),
                                                    Gene = rep(NA, times = (length(clades) * length(genes))),
                                                    Branch_length = rep(NA, times = (length(clades) * length(genes))),
                                                    Standard_error = rep(NA, times = (length(clades) * length(genes))))

count <- 0
for(clade in 1:length(clades)){
    subset <- Ultrafiltration_BranchLengths[Ultrafiltration_BranchLengths$Main_clade == clades[clade], ]
    
    for(gene in 1:length(genes)){
        count <- count + 1
        subset_v2 <- subset[subset$Gene == genes[gene],]
        print(subset_v2[1,1])
        Ultrafiltration_BranchLengths_Average$Main_clade[count] <- as.character(subset_v2[1, 1])
        Ultrafiltration_BranchLengths_Average$Gene[count] <- as.character(subset_v2[1, 3])
        Ultrafiltration_BranchLengths_Average$Branch_length[count] <- mean(subset_v2$Branch_length)
        Ultrafiltration_BranchLengths_Average$Standard_error[count] <- sd(subset_v2$Branch_length)/sqrt(length(subset_v2$Branch_length))
    }
}
Ultrafiltration_BranchLengths_Average <- na.omit(Ultrafiltration_BranchLengths_Average)

Ultrafiltration_BranchLengths_Average$Main_clade <- factor(Ultrafiltration_BranchLengths_Average$Main_clade, 
                                                           levels = c("Cnidaria", "Deuterostomia", "Protostomia", "Acoelomorpha"))
Ultrafiltration_BranchLengths_Average$Gene <- factor(Ultrafiltration_BranchLengths_Average$Gene, 
                                                     levels = c("Nephrin", "Kirrel", "ZO1", 
                                                                "Eya", "Lhx", "Osr", "POU3", "Sall", "Six"))

## Plot these results
Barplot.BranchLengths <- ggplot(Ultrafiltration_BranchLengths_Average, aes(x = Gene, y = Branch_length, fill = Main_clade)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_y_continuous(limits = c(0, 6), n.breaks = 12) + 
  guides(fill=guide_legend(title="Clades")) + 
  theme_light(base_size = 8) + 
  scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
  ggtitle("Average Branch Length per Gene and Clade") + 
  xlab("") +  ylab("Branch length") + 
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

Barplot.BranchLengths <- Barplot.BranchLengths + geom_errorbar(aes(ymin=Branch_length-Standard_error, ymax=Branch_length+Standard_error), width=.2,
                                                               position=position_dodge(.9)) 

Barplot.BranchLengths

################################################################################

## Calculate if the differences are statistically significant

## For each gene, check the normality of the data
for(i in 1:length(genes)){
  subset <- Ultrafiltration_BranchLengths[grep(genes[i], Ultrafiltration_BranchLengths$Gene), ]
  
  # Print the gene name
  print(genes[i])
  
  # Calculate the normality for protein length, number of exons, and average 
  # exon length. In that order
  print(shapiro.test(subset$Branch_length))

  # Print a blank line to separate genes
  print("   ")
  print("   ")
}

## All but Lhx of the metrics are not normal. I will use the Barlett
## test, since it is robust if the data are normal and it allows comparisons 
## among more than two populations

bartlett.test(x = Ultrafiltration_BranchLengths$Branch_length[grep("Lhx", Ultrafiltration_BranchLengths$Gene, )],
              g = Ultrafiltration_BranchLengths$Main_clade[grep("Lhx", Ultrafiltration_BranchLengths$Gene, )])

## The branch lengths of this gene are also homocedastic (equal variances), so I 
## can run an ANOVA on it. For all others, I will run a Kruskal-Wallis

for(i in 1:length(genes)){
  subset <- Ultrafiltration_BranchLengths[grep(genes[i], Ultrafiltration_BranchLengths$Gene), ]
  
  # Print the gene name
  print(genes[i])
  
  # Calculate the normality for protein length, number of exons, and average 
  # exon length. In that order. Do this with an if statement to separate
  # the ANOVA from Kruskal-Wallis
  if(genes[i] == "Lhx"){
    anova <- aov(Branch_length ~ Main_clade, data = subset)
    print(summary(anova))
  } else{
    print(kruskal.test(Branch_length ~ Main_clade, data = subset))
  }
  
  # Print a blank line to separate genes
  print("   ")
  print("   ")
}

## Branch lengths are statistically difference in all genes. Only Lhx is normal, 
## for which I will use Bonferroni For all others, I will use the Dunn test 
## adjusted with the Holm method. Because only one sample needs to be analysed 
## with Bonferroni, I will run Dunn over all of them and then run Bonferroni 
## individually

## Note: Osr is failing, so I skip it now. It only has two categories anyway:
##       Protostomia and Deuterostomia

for(i in c(1:4, 7:length(genes))){
  subset <- Ultrafiltration_BranchLengths[grep(genes[i], Ultrafiltration_BranchLengths$Gene), ]
  
  # Print the gene name
  print(genes[i])
  
  # Run multiple comparisons to check which pairs are significantly different
  print(dunnTest(Branch_length ~ Main_clade, data = subset))

  # Print a blank line to separate genes
  print("   ")
  print("   ")
}
subset <- Ultrafiltration_BranchLengths[grep("Lhx", Ultrafiltration_BranchLengths$Gene), ]
pairwise.t.test(subset$Branch_length, subset$Main_clade, p.adjust.method = "bonferroni")

################################################################################
################################################################################

####################
## Save all plots ##
####################

ggsave(Boxplot.BranchLengths, file="Branch_lengths_analysis.Boxplot.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(Barplot.BranchLengths, file="Branch_lengths_analysis.Barplot.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
