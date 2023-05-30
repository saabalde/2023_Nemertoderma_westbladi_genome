############################################################
## Generate annotation stats from the N. westbladi genome ##
############################################################

## Load libraries
library(ggplot2)
library(reshape2)
library(tidyr)


## Set the working directory
setwd("/home/saabalde/Escritorio/Nwestbladi_genome/00-Figures/")


## Load the data
# Cumulative length of the genome
Illumina_genomeLength <- read.table(file = "Illumina_CleanAssembly.length.txt", header = TRUE, sep = ",")
HiFi_genomeLength <- read.table(file = "pt_087_CleanAssembly.nomtDNA.length.txt", header = TRUE, sep = ",")
Pnaikaiensis_genomeLength <- read.table(file = "Pnaikaiensis.CleanAssembly.length.txt", header = TRUE, sep = ",")
Sroscoffensis_genomeLength <- read.table(file = "Sroscoffensis.CleanAssembly.length.txt", header = TRUE, sep = ",")

# Genomes annotation
Illumina_annotations <- read.table(file = "Illumina_CleanAssembly.annotations.gff3", header = FALSE, sep = "\t")
names(Illumina_annotations) <- c("Contig", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Atributes")
HiFi_annotations <- read.table(file = "pt_087_CleanAssembly.annotations.gff3", header = FALSE, sep = "\t")
names(HiFi_annotations) <- c("Contig", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Atributes")
Pnaikaiensis_annotations <- read.table(file = "Pnaikaiensis.annotations.gff3", header = FALSE, sep = "\t")
names(Pnaikaiensis_annotations) <- c("Contig", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Atributes")
Sroscoffensis_annotations <- read.table(file = "Sroscoffensis.annotations.gff3", header = FALSE, sep = "\t")
names(Sroscoffensis_annotations) <- c("Contig", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Atributes")

# Create subsets for gene, exon and intron positions
Illumina_fullgene <- Illumina_annotations[grep("gene", Illumina_annotations$Type),]
Illumina_exons <- Illumina_annotations[grep("exon", Illumina_annotations$Type),]
Illumina_intron <- Illumina_annotations[grep("intron", Illumina_annotations$Type),]

HiFi_fullgene <- HiFi_annotations[grep("gene", HiFi_annotations$Type),]
HiFi_exons <- HiFi_annotations[grep("exon", HiFi_annotations$Type),]
HiFi_intron <- HiFi_annotations[grep("intron", HiFi_annotations$Type),]

Pnaikaiensis_fullgene <- Pnaikaiensis_annotations[grep("gene", Pnaikaiensis_annotations$Type),]
Pnaikaiensis_exons <- Pnaikaiensis_annotations[grep("exon", Pnaikaiensis_annotations$Type),]
Pnaikaiensis_intron <- Pnaikaiensis_annotations[grep("intron", Pnaikaiensis_annotations$Type),]

Sroscoffensis_fullgene <- Sroscoffensis_annotations[grep("gene", Sroscoffensis_annotations$Type),]
Sroscoffensis_exons <- Sroscoffensis_annotations[grep("exon", Sroscoffensis_annotations$Type),]
Sroscoffensis_intron <- Sroscoffensis_annotations[grep("intron", Sroscoffensis_annotations$Type),]

# Calculate the intron length of S. roscoffensis
Sroscoffensis_exons_modified <- separate(Sroscoffensis_exons, Atributes, into = c("dump", "Parent"), sep = ";")
Sroscoffensis_genes_names <- unique(Sroscoffensis_exons_modified$Parent)

count <- 0
for(i in 1:length(Sroscoffensis_genes_names)){
    gene <- Sroscoffensis_genes_names[i]
    
    if(length(grep(gene, Sroscoffensis_exons_modified$Parent)) > 1){
        dataframe <- Sroscoffensis_exons[grep(gene, Sroscoffensis_exons$Atributes), c(4:5)]
        
        for(j in 2:nrow(dataframe)){
            count <- count + 1
            Sroscoffensis_intron[count, c(1, 3:5, 9)] <- c(Sroscoffensis_exons[grep(gene, Sroscoffensis_exons$Atributes)[1], 1], "intron", dataframe$End[(j -1)] + 1, dataframe$Start[j] - 1, gene)
        }
    }
}
write.table(Sroscoffensis_intron, file = "Sroscoffensis_intron.tsv", sep = "\t", dec = ".")
#Sroscoffensis_intron <- read.table(file = "Sroscoffensis_intron.tsv", header = TRUE, sep = "\t", dec = ".")

################################################################################

######################
## Stats per contig ##
######################

## Obtain the number of genes, exons and introns per contig in the genomes.
# Illumina genome
Illumina_Stats_Per_Contig <- data.frame(matrix(NA, nrow = length(Illumina_genomeLength$Contig), ncol = 4))
names(Illumina_Stats_Per_Contig) <- c("Contig", "N_Genes", "N_Exons", "N_Introns")

Illumina_Stats_Per_Contig$Contig <- Illumina_genomeLength$Contig
for(i in 1:nrow(Illumina_Stats_Per_Contig)){
  contig <- paste("^", Illumina_Stats_Per_Contig$Contig[i], "$", sep = "")
  Illumina_Stats_Per_Contig$N_Genes[i]   <- length(grep("gene", Illumina_annotations[grep(contig, Illumina_annotations$Contig), 3]))
  Illumina_Stats_Per_Contig$N_Exons[i]   <- length(grep("exon", Illumina_annotations[grep(contig, Illumina_annotations$Contig), 3]))
  Illumina_Stats_Per_Contig$N_Introns[i] <- length(grep("intron", Illumina_annotations[grep(contig, Illumina_annotations$Contig), 3]))
}
write.table(Illumina_Stats_Per_Contig, file = "Illumina_Stats_Per_Contig.tsv", sep = "\t", dec = ".")
#Illumina_Stats_Per_Contig <- read.table(file = "Illumina_Stats_Per_Contig.tsv", header = TRUE, sep = "\t", dec = ".")

# HiFi genome
HiFi_Stats_Per_Contig <- data.frame(matrix(NA, nrow = length(HiFi_genomeLength$Contig), ncol = 4))
names(HiFi_Stats_Per_Contig) <- c("Contig", "N_Genes", "N_Exons", "N_Introns")

HiFi_Stats_Per_Contig$Contig <- HiFi_genomeLength$Contig
for(i in 1:nrow(HiFi_Stats_Per_Contig)){
  contig <- paste("^", HiFi_Stats_Per_Contig$Contig[i], "$", sep = "")
  HiFi_Stats_Per_Contig$N_Genes[i]   <- length(grep("gene", HiFi_annotations[grep(contig, HiFi_annotations$Contig), 3]))
  HiFi_Stats_Per_Contig$N_Exons[i]   <- length(grep("exon", HiFi_annotations[grep(contig, HiFi_annotations$Contig), 3]))
  HiFi_Stats_Per_Contig$N_Introns[i] <- length(grep("intron", HiFi_annotations[grep(contig, HiFi_annotations$Contig), 3]))
}
write.table(HiFi_Stats_Per_Contig, file = "pt_087_Stats_Per_Contig.tsv", sep = "\t", dec = ".")
#HiFi_Stats_Per_Contig <- read.table(file = "pt_087_Stats_Per_Contig.tsv", header = TRUE, sep = "\t", dec = ".")

# Pnaikaiensis genome
Pnaikaiensis_Stats_Per_Contig <- data.frame(matrix(NA, nrow = length(Pnaikaiensis_genomeLength$Contig), ncol = 4))
names(Pnaikaiensis_Stats_Per_Contig) <- c("Contig", "N_Genes", "N_Exons", "N_Introns")

Pnaikaiensis_Stats_Per_Contig$Contig <- Pnaikaiensis_genomeLength$Contig
for(i in 1:nrow(Pnaikaiensis_Stats_Per_Contig)){
  contig <- paste("^", Pnaikaiensis_Stats_Per_Contig$Contig[i], "$", sep = "")
  Pnaikaiensis_Stats_Per_Contig$N_Genes[i]   <- length(grep("gene", Pnaikaiensis_annotations[grep(contig, Pnaikaiensis_annotations$Contig), 3]))
  Pnaikaiensis_Stats_Per_Contig$N_Exons[i]   <- length(grep("exon", Pnaikaiensis_annotations[grep(contig, Pnaikaiensis_annotations$Contig), 3]))
  Pnaikaiensis_Stats_Per_Contig$N_Introns[i] <- length(grep("intron", Pnaikaiensis_annotations[grep(contig, Pnaikaiensis_annotations$Contig), 3]))
  print(paste(i, " out of ", nrow(Pnaikaiensis_Stats_Per_Contig), sep = ""))
}
write.table(Pnaikaiensis_Stats_Per_Contig, file = "Pnaikaiensis_Stats_Per_Contig.tsv", sep = "\t", dec = ".")
#Pnaikaiensis_Stats_Per_Contig <- read.table(file = "Pnaikaiensis_Stats_Per_Contig.tsv", header = TRUE, sep = "\t", dec = ".")

# Sroscoffensis genome
Sroscoffensis_Stats_Per_Contig <- data.frame(matrix(NA, nrow = length(Sroscoffensis_genomeLength$Contig), ncol = 4))
names(Sroscoffensis_Stats_Per_Contig) <- c("Contig", "N_Genes", "N_Exons", "N_Introns")

Sroscoffensis_Stats_Per_Contig$Contig <- Sroscoffensis_genomeLength$Contig
for(i in 1:nrow(Sroscoffensis_Stats_Per_Contig)){
  contig <- paste("^", Sroscoffensis_Stats_Per_Contig$Contig[i], "$", sep = "")
  Sroscoffensis_Stats_Per_Contig$N_Genes[i]   <- length(grep("gene", Sroscoffensis_annotations[grep(contig, Sroscoffensis_annotations$Contig), 3]))
  Sroscoffensis_Stats_Per_Contig$N_Exons[i]   <- length(grep("exon", Sroscoffensis_annotations[grep(contig, Sroscoffensis_annotations$Contig), 3]))
  Sroscoffensis_Stats_Per_Contig$N_Introns[i] <- length(grep("intron", Sroscoffensis_annotations[grep(contig, Sroscoffensis_annotations$Contig), 3]))
  print(paste(i, " out of ", nrow(Sroscoffensis_Stats_Per_Contig), sep = ""))
}
write.table(Sroscoffensis_Stats_Per_Contig, file = "Sroscoffensis_Stats_Per_Contig.tsv", sep = "\t", dec = ".")
#Sroscoffensis_Stats_Per_Contig <- read.table(file = "Sroscoffensis_Stats_Per_Contig.tsv", header = TRUE, sep = "\t", dec = ".")


# Compare some basic stats to have an overview of the genomes
mean(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 0])                # 1.137114
mean(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 0])                        # 2.504068
mean(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 0])        # 3.983324
mean(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 0])      # 16.50812

mean(Illumina_Stats_Per_Contig$N_Exons[Illumina_Stats_Per_Contig$N_Exons > 0])                # 1.740685
mean(HiFi_Stats_Per_Contig$N_Exons[HiFi_Stats_Per_Contig$N_Exons > 0])                        # 7.62197
mean(Pnaikaiensis_Stats_Per_Contig$N_Exons[Pnaikaiensis_Stats_Per_Contig$N_Exons > 0])        # 25.81048
mean(Sroscoffensis_Stats_Per_Contig$N_Exons[Sroscoffensis_Stats_Per_Contig$N_Exons > 0])      # 70.06155

mean(Illumina_Stats_Per_Contig$N_Introns[Illumina_Stats_Per_Contig$N_Introns > 0])            # 1.892656
mean(HiFi_Stats_Per_Contig$N_Introns[HiFi_Stats_Per_Contig$N_Introns > 0])                    # 6.296576
mean(Pnaikaiensis_Stats_Per_Contig$N_Introns[Pnaikaiensis_Stats_Per_Contig$N_Introns > 0])    # 20.91962
mean(Sroscoffensis_Stats_Per_Contig$N_Introns[Sroscoffensis_Stats_Per_Contig$N_Introns > 0])  # NA

max(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 0])                 # 33
max(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 0])                         # 89
max(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 0])         # 37
max(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 0])       # 280

sum(Illumina_Stats_Per_Contig$N_Exons) / sum(Illumina_Stats_Per_Contig$N_Genes)               # 1.530792
sum(HiFi_Stats_Per_Contig$N_Exons) / sum(HiFi_Stats_Per_Contig$N_Genes)                       # 3.043834
sum(Pnaikaiensis_Stats_Per_Contig$N_Exons) / sum(Pnaikaiensis_Stats_Per_Contig$N_Genes)       # 6.479634
sum(Sroscoffensis_Stats_Per_Contig$N_Exons) / sum(Sroscoffensis_Stats_Per_Contig$N_Genes)     # 4.244065

length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 0])              # 20049
length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 1])              # 2223
length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 2])              # 375
length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 3])              # 74
length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 4])              # 22
length(Illumina_Stats_Per_Contig$N_Genes[Illumina_Stats_Per_Contig$N_Genes > 5])              # 11

length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 0])                      # 11798
length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 1])                      # 6799
length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 2])                      # 3884
length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 3])                      # 2277
length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 4])                      # 1397
length(HiFi_Stats_Per_Contig$N_Genes[HiFi_Stats_Per_Contig$N_Genes > 5])                      # 893

length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 0])      # 5097
length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 1])      # 3524
length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 2])      # 2546
length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 3])      # 1894
length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 4])      # 1473
length(Pnaikaiensis_Stats_Per_Contig$N_Genes[Pnaikaiensis_Stats_Per_Contig$N_Genes > 5])      # 1165

length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 0])    # 2031
length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 1])    # 1608
length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 2])    # 1367
length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 3])    # 1219
length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 4])    # 1107
length(Sroscoffensis_Stats_Per_Contig$N_Genes[Sroscoffensis_Stats_Per_Contig$N_Genes > 5])    # 1002

# Make two simple plots comparing the % of contigs with X number of genes
Illumina.Piechart.data <- data.frame(NGenes = c("1", "2", "3 or more"),
                                     Count = c(length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 1, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 2, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes > 2, 1])))

HiFi.Piechart.data <- data.frame(NGenes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                     Count = c(length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 1, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 2, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 3, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 4, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 5, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 6, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 7, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 8, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 9, 1]),
                                               length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes > 9, 1])))

Pnaikaiensis.Piechart.data <- data.frame(NGenes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                 Count = c(length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 1, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 2, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 3, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 4, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 5, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 6, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 7, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 8, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 9, 1]),
                                           length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes > 9, 1])))

Sroscoffensis.Piechart.data <- data.frame(NGenes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                         Count = c(length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 1, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 2, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 3, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 4, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 5, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 6, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 7, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 8, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 9, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes > 9, 1])))

pie(Illumina.Piechart.data$Count, labels = Illumina.Piechart.data$NGenes, clockwise = TRUE)
pie(HiFi.Piechart.data$Count, labels = HiFi.Piechart.data$NGenes, clockwise = TRUE)
pie(Pnaikaiensis.Piechart.data$Count, labels = Pnaikaiensis.Piechart.data$NGenes, clockwise = TRUE)
pie(Sroscoffensis.Piechart.data$Count, labels = Sroscoffensis.Piechart.data$NGenes, clockwise = TRUE)


Illumina.BarPlot.data <- data.frame(Count = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                    NGenes = c(length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 1, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 2, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 3, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 4, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 5, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 6, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 7, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 8, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes == 9, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Genes > 9, 1])),
                                    NExons = c(length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 1, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 2, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 3, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 4, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 5, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 6, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 7, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 8, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons == 9, 1]),
                                               length(Illumina_Stats_Per_Contig[Illumina_Stats_Per_Contig$N_Exons > 9, 1])))


HiFi.BarPlot.data <- data.frame(Count = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                NGenes = c(length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 1, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 2, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 3, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 4, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 5, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 6, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 7, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 8, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes == 9, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Genes > 9, 1])),
                                NExons = c(length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 1, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 2, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 3, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 4, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 5, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 6, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 7, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 8, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons == 9, 1]),
                                           length(HiFi_Stats_Per_Contig[HiFi_Stats_Per_Contig$N_Exons > 9, 1])))

Pnaikaiensis.BarPlot.data <- data.frame(Count = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                        NGenes = c(length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 1, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 2, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 3, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 4, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 5, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 6, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 7, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 8, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes == 9, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Genes > 9, 1])),
                                        NExons = c(length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 1, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 2, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 3, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 4, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 5, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 6, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 7, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 8, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons == 9, 1]),
                                                   length(Pnaikaiensis_Stats_Per_Contig[Pnaikaiensis_Stats_Per_Contig$N_Exons > 9, 1])))

Sroscoffensis.BarPlot.data <- data.frame(Count = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"),
                                        NGenes = c(length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 1, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 2, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 3, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 4, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 5, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 6, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 7, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 8, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes == 9, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Genes > 9, 1])),
                                        NExons = c(length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 1, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 2, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 3, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 4, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 5, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 6, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 7, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 8, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons == 9, 1]),
                                                   length(Sroscoffensis_Stats_Per_Contig[Sroscoffensis_Stats_Per_Contig$N_Exons > 9, 1])))

NGenes.Percentage <- c()
NExons.Percentage <- c()
for(i in 1:nrow(Illumina.BarPlot.data)){
  percentage.genes <- Illumina.BarPlot.data[i, 2] / sum(Illumina.BarPlot.data[, 2]) * 100
  percentage.exons <- Illumina.BarPlot.data[i, 3] / sum(Illumina.BarPlot.data[, 3]) * 100
  NGenes.Percentage <- c(NGenes.Percentage, percentage.genes)
  NExons.Percentage <- c(NExons.Percentage, percentage.exons)
}
Illumina.BarPlot.data$NGenes.Percentage <- NGenes.Percentage
Illumina.BarPlot.data$NExons.Percentage <- NExons.Percentage

NGenes.Percentage <- c()
NExons.Percentage <- c()
for(i in 1:nrow(HiFi.BarPlot.data)){
  percentage.genes <- HiFi.BarPlot.data[i, 2] / sum(HiFi.BarPlot.data[, 2]) * 100
  percentage.exons <- HiFi.BarPlot.data[i, 3] / sum(HiFi.BarPlot.data[, 3]) * 100
  NGenes.Percentage <- c(NGenes.Percentage, percentage.genes)
  NExons.Percentage <- c(NExons.Percentage, percentage.exons)
}
HiFi.BarPlot.data$NGenes.Percentage <- NGenes.Percentage
HiFi.BarPlot.data$NExons.Percentage <- NExons.Percentage

NGenes.Percentage <- c()
NExons.Percentage <- c()
for(i in 1:nrow(Pnaikaiensis.BarPlot.data)){
  percentage.genes <- Pnaikaiensis.BarPlot.data[i, 2] / sum(Pnaikaiensis.BarPlot.data[, 2]) * 100
  percentage.exons <- Pnaikaiensis.BarPlot.data[i, 3] / sum(Pnaikaiensis.BarPlot.data[, 3]) * 100
  NGenes.Percentage <- c(NGenes.Percentage, percentage.genes)
  NExons.Percentage <- c(NExons.Percentage, percentage.exons)
}
Pnaikaiensis.BarPlot.data$NGenes.Percentage <- NGenes.Percentage
Pnaikaiensis.BarPlot.data$NExons.Percentage <- NExons.Percentage

NGenes.Percentage <- c()
NExons.Percentage <- c()
for(i in 1:nrow(Sroscoffensis.BarPlot.data)){
  percentage.genes <- Sroscoffensis.BarPlot.data[i, 2] / sum(Sroscoffensis.BarPlot.data[, 2]) * 100
  percentage.exons <- Sroscoffensis.BarPlot.data[i, 3] / sum(Sroscoffensis.BarPlot.data[, 3]) * 100
  NGenes.Percentage <- c(NGenes.Percentage, percentage.genes)
  NExons.Percentage <- c(NExons.Percentage, percentage.exons)
}
Sroscoffensis.BarPlot.data$NGenes.Percentage <- NGenes.Percentage
Sroscoffensis.BarPlot.data$NExons.Percentage <- NExons.Percentage


NGenes_Per_Contig.BarPlot <- data.frame(NGenes = Illumina.BarPlot.data[, 1],
                                        Illumina = Illumina.BarPlot.data[, c(4)], 
                                        HiFi = HiFi.BarPlot.data[, 4], 
                                        Pnaikaiensis = Pnaikaiensis.BarPlot.data[, 4],
                                        Sroscoffensis = Sroscoffensis.BarPlot.data[, 4])
NGenes_Per_Contig.BarPlot.reshape <- melt(NGenes_Per_Contig.BarPlot[, c('NGenes', 'Illumina', 'HiFi', 'Pnaikaiensis', 'Sroscoffensis')], id.vars = 1)
NGenes_Per_Contig.BarPlot.reshape$NGenes <- factor(NGenes_Per_Contig.BarPlot.reshape$NGenes, 
                                                   levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"))

BarPlot.NGenes <- ggplot(NGenes_Per_Contig.BarPlot.reshape, aes(x = NGenes, y = value, fill = variable)) + 
                         geom_bar(stat = "identity", position = "dodge") + 
                         scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 5)) + 
                         guides(fill=guide_legend(title="Genome")) + 
                         theme_light(base_size = 8) + 
                         scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                         ggtitle("Number of Genes per Contig") + 
                         xlab("# Genes") +  ylab("% Contigs") + 
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
BarPlot.NGenes


# Make another plot with the number of exons per contig
NExons_Per_Contig.BarPlot <- data.frame(NExons = Illumina.BarPlot.data[, 1],
                                         Illumina = Illumina.BarPlot.data[, 5], 
                                         HiFi = HiFi.BarPlot.data[, 5], 
                                         Pnaikaiensis = Pnaikaiensis.BarPlot.data[, 5],
                                         Sroscoffensis = Sroscoffensis.BarPlot.data[, 5])
NExons_Per_Contig.BarPlot.reshape <- melt(NExons_Per_Contig.BarPlot[, c('NExons', 'Illumina', 'HiFi', 'Pnaikaiensis', 'Sroscoffensis')], id.vars = 1)
NExons_Per_Contig.BarPlot.reshape$NExons <- factor(NExons_Per_Contig.BarPlot.reshape$NExons, 
                                                   levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"))

BarPlot.NExons <- ggplot(NExons_Per_Contig.BarPlot.reshape, aes(x = NExons, y = value, fill = variable)) + 
                         geom_bar(stat = "identity", position = "dodge") + 
                         scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 5)) + 
                         guides(fill=guide_legend(title="Genome")) + 
                         theme_light(base_size = 8) + 
                         scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                         ggtitle("Number of Exons per Contig") + 
                         xlab("# Exons") +  ylab("% Contigs") + 
                         theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                         theme(legend.position = "top",legend.title = element_blank()) + 
                         theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                         theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                         theme(panel.grid.minor = element_blank()) + 
                         theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                         theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                         theme(axis.line = element_line(size = 1, colour = "black")) + 
                         theme(axis.ticks.length = unit(0, "cm")) + 
                         theme(axis.ticks.y = element_line(colour = "#222222")) + 
                         theme(axis.ticks.x = element_line(colour = "#222222")) + 
                         theme(axis.ticks.length = unit(0.4, "cm")) + 
                         theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                         theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                         guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                         guides(fill = guide_legend(nrow = 1, byrow = TRUE))
BarPlot.NExons

# Make one plot with the number of exons per gene
Illumina_Exons_Per_Gene <- data.frame(matrix(NA, nrow = length(Illumina_fullgene$Type), ncol = 2))
names(Illumina_Exons_Per_Gene) <- c("Gene_name", "N_Exons")

for(i in 1:nrow(Illumina_Exons_Per_Gene)){
  Illumina_Exons_Per_Gene[i, 1] <- substring(Illumina_fullgene$Atributes[i], 4, nchar(Illumina_fullgene$Atributes[i])-1)

  gene <- paste("=", Illumina_Exons_Per_Gene[i, 1], ".t", sep = "")
  Illumina_Exons_Per_Gene[i, 2] <- length(grep(gene, Illumina_exons$Atributes))
}
write.table(Illumina_Exons_Per_Gene, file = "Illumina_Exons_Per_Gene.tsv", sep = "\t", dec = ".")
#Illumina_Exons_Per_Gene <- read.table(file = "Illumina_Exons_Per_Gene.tsv", header = TRUE, sep = "\t", dec = ".")


HiFi_Exons_Per_Gene <- data.frame(matrix(NA, nrow = length(HiFi_fullgene$Type), ncol = 2))
names(HiFi_Exons_Per_Gene) <- c("Gene_name", "N_Exons")

for(i in 1:nrow(HiFi_Exons_Per_Gene)){
  HiFi_Exons_Per_Gene[i, 1] <- substring(HiFi_fullgene$Atributes[i], 4, nchar(HiFi_fullgene$Atributes[i])-1)
  
  gene <- paste("=", HiFi_Exons_Per_Gene[i, 1], ".t", sep = "")
  HiFi_Exons_Per_Gene[i, 2] <- length(grep(gene, HiFi_exons$Atributes))
}
write.table(HiFi_Exons_Per_Gene, file = "pt_087_Exons_Per_Gene.tsv", sep = "\t", dec = ".")
#HiFi_Exons_Per_Gene <- read.table(file = "pt_087_Exons_Per_Gene.tsv", header = TRUE, sep = "\t", dec = ".")


Pnaikaiensis_Exons_Per_Gene <- data.frame(matrix(NA, nrow = length(Pnaikaiensis_fullgene$Type), ncol = 2))
names(Pnaikaiensis_Exons_Per_Gene) <- c("Gene_name", "N_Exons")

for(i in 1:nrow(Pnaikaiensis_Exons_Per_Gene)){
  Pnaikaiensis_Exons_Per_Gene[i, 1] <- substring(Pnaikaiensis_fullgene$Atributes[i], 4, nchar(Pnaikaiensis_fullgene$Atributes[i]))
  
  gene <- paste("=", Pnaikaiensis_Exons_Per_Gene[i, 1], ".t", sep = "")
  Pnaikaiensis_Exons_Per_Gene[i, 2] <- length(grep(gene, Pnaikaiensis_exons$Atributes))
}
write.table(Pnaikaiensis_Exons_Per_Gene, file = "Pnaikaiensis_Exons_Per_Gene.tsv", sep = "\t", dec = ".")
#Pnaikaiensis_Exons_Per_Gene <- read.table(file = "Pnaikaiensis_Exons_Per_Gene.tsv", header = TRUE, sep = "\t", dec = ".")


Sroscoffensis_Exons_Per_Gene <- data.frame(matrix(NA, nrow = length(Sroscoffensis_fullgene$Type), ncol = 2))
names(Sroscoffensis_Exons_Per_Gene) <- c("Gene_name", "N_Exons")

for(i in 1:nrow(Sroscoffensis_Exons_Per_Gene)){
  Sroscoffensis_Exons_Per_Gene[i, 1] <- substring(Sroscoffensis_fullgene$Atributes[i], 4, (nchar(Sroscoffensis_fullgene$Atributes[i]) - 5))
  
  gene <- paste("=", Sroscoffensis_Exons_Per_Gene[i, 1], ".exon", sep = "")
  Sroscoffensis_Exons_Per_Gene[i, 2] <- length(grep(gene, Sroscoffensis_exons$Atributes))
}
write.table(Sroscoffensis_Exons_Per_Gene, file = "Sroscoffensis_Exons_Per_Gene.tsv", sep = "\t", dec = ".")
#Sroscoffensis_Exons_Per_Gene <- read.table(file = "Sroscoffensis_Exons_Per_Gene.tsv", header = TRUE, sep = "\t", dec = ".")


max(nrow(Illumina_Exons_Per_Gene), nrow(HiFi_Exons_Per_Gene), nrow(Pnaikaiensis_Exons_Per_Gene), nrow(Sroscoffensis_Exons_Per_Gene)) # Get the number of rows: 33528
NExons_Per_Gene.Barplot <- data.frame(Index = c(1:33528),
                                      Illumina = c(Illumina_Exons_Per_Gene$N_Exons, rep(NA, 33528 - length(Illumina_Exons_Per_Gene$N_Exons))),
                                      HiFi = c(HiFi_Exons_Per_Gene$N_Exons, rep(NA, 33528 - length(HiFi_Exons_Per_Gene$N_Exons))),
                                      Pnaikaiensis = c(Pnaikaiensis_Exons_Per_Gene$N_Exons, rep(NA, 33528 - length(Pnaikaiensis_Exons_Per_Gene$N_Exons))),
                                      Sroscoffensis = c(Sroscoffensis_Exons_Per_Gene$N_Exons, rep(NA, 33528 - length(Sroscoffensis_Exons_Per_Gene$N_Exons))))
NExons_Per_Gene.Barplot.reshape <- melt(NExons_Per_Gene.Barplot[, c('Index', 'Illumina', 'HiFi', 'Pnaikaiensis', 'Sroscoffensis')], id.vars = 1, na.rm = TRUE)
NExons_Per_Contig.BarPlot.reshape$NExons <- factor(NExons_Per_Contig.BarPlot.reshape$NExons, 
                                                   levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more"))

BarPlot.NExons_Per_Gene <- ggplot(NExons_Per_Gene.Barplot.reshape, aes(x = variable, y = value, fill = variable)) +
                                  geom_violin() + scale_y_log10(breaks=c(1, 10, 50, 100, 250, 500)) + 
                                  theme_light() + 
                                  guides(fill=guide_legend(title="Genome")) + 
                                  theme_light(base_size = 8) + 
                                  scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                                  ggtitle("Number of Exons per Gene") + 
                                  xlab("") +  ylab("# Exons (log10)") + 
                                  theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                                  theme(legend.position = "top",legend.title = element_blank()) + 
                                  theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                                  theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                                  theme(panel.grid.minor = element_blank()) + 
                                  theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                                  theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                                  theme(axis.line = element_line(size = 1, colour = "black")) + 
                                  theme(axis.ticks.length = unit(0, "cm")) + 
                                  theme(axis.ticks.y = element_line(colour = "#222222")) + 
                                  theme(axis.ticks.x = element_line(colour = "#222222")) + 
                                  theme(axis.ticks.length = unit(0.4, "cm")) + 
                                  theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                                  theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                                  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                                  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
BarPlot.NExons_Per_Gene


## Calculate the intron length for all the genomes
Intron_length <- data.frame(Index = c(1:(max(nrow(HiFi_intron), nrow(Illumina_intron), nrow(Pnaikaiensis_intron), nrow(Sroscoffensis_intron)) - 1)),
                            HiFi = rep(NA, times = (max(nrow(HiFi_intron), nrow(Illumina_intron), nrow(Pnaikaiensis_intron), nrow(Sroscoffensis_intron)) - 1)),
                            Illumina = rep(NA, times = (max(nrow(HiFi_intron), nrow(Illumina_intron), nrow(Pnaikaiensis_intron), nrow(Sroscoffensis_intron)) - 1)),
                            Pnaikaiensis = rep(NA, times = (max(nrow(HiFi_intron), nrow(Illumina_intron), nrow(Pnaikaiensis_intron), nrow(Sroscoffensis_intron)) - 1)),
                            Sroscoffensis = rep(NA, times = (max(nrow(HiFi_intron), nrow(Illumina_intron), nrow(Pnaikaiensis_intron), nrow(Sroscoffensis_intron)) - 1)))

for(i in 1:nrow(Intron_length)){
    if(i <= nrow(HiFi_intron)){
        Intron_length$HiFi[i] <- abs(HiFi_intron$End[i] - HiFi_intron$Start[i] + 1)
    }
    if(i <= nrow(Illumina_intron)){
        Intron_length$Illumina[i] <- abs(Illumina_intron$End[i] - Illumina_intron$Start[i] + 1)
    }
    if(i <= nrow(Pnaikaiensis_intron)){
        Intron_length$Pnaikaiensis[i] <- abs(Pnaikaiensis_intron$End[i] - Pnaikaiensis_intron$Start[i] + 1)
    }
    if(i <= nrow(Sroscoffensis_intron)){
        Intron_length$Sroscoffensis[i] <- abs(Sroscoffensis_intron$End[i] - Sroscoffensis_intron$Start[i] + 1)
    }
}
write.table(Intron_length, file = "Intron_length.tsv", sep = "\t", dec = ".")
#Intron_length <- read.table(file = "Intron_length.tsv", header = TRUE, sep = "\t", dec = ".")


# Calculate the average, min, and max
min(Intron_length$HiFi, na.rm = TRUE)
min(Intron_length$Illumina, na.rm = TRUE)
min(Intron_length$Pnaikaiensis, na.rm = TRUE)
min(Intron_length$Sroscoffensis, na.rm = TRUE)

mean(Intron_length$HiFi, na.rm = TRUE)
mean(Intron_length$Illumina, na.rm = TRUE)
mean(Intron_length$Pnaikaiensis, na.rm = TRUE)
mean(Intron_length$Sroscoffensis, na.rm = TRUE)

max(Intron_length$HiFi, na.rm = TRUE)
max(Intron_length$Illumina, na.rm = TRUE)
max(Intron_length$Pnaikaiensis, na.rm = TRUE)
max(Intron_length$Sroscoffensis, na.rm = TRUE)

# Plot this
Intron_length.reshape <- melt(Intron_length[, c('Index', 'Illumina', 'HiFi', 'Pnaikaiensis', 'Sroscoffensis')], id.vars = 1, na.rm = TRUE)

BarPlot.Intron_length <- ggplot(Intron_length.reshape, aes(x = variable, y = (value / 1000), fill = variable)) +
                                geom_violin() + scale_y_log10(breaks = c(0, 0.01, 0.1, 1, 10, 100, 200),
                                                              labels = c("0", "0.01", "0.1", "1", "10", "100", "200")) + 
                                theme_light() + 
                                guides(fill=guide_legend(title="Genome")) + 
                                theme_light(base_size = 8) + 
                                scale_fill_manual(values = c("#EB9F17", "#EB3717", "#2877D1", "#00B366")) + 
                                ggtitle("Intron length") + 
                                xlab("") +  ylab("Intron length in kbps (log10)") + 
                                theme(plot.title = element_text(family = "sans", colour = "black", size = rel(2.2), face = "bold")) + 
                                theme(legend.position = "top",legend.title = element_blank()) + 
                                theme(legend.text = element_text(family = "sans", size = rel(1.5))) + 
                                theme(panel.background = element_rect(color = "#FFFFFF", fill = "white")) + 
                                theme(panel.grid.minor = element_blank()) + 
                                theme(axis.text.y = element_text(family = "sans", colour = "black", size = rel(2))) + 
                                theme(axis.text.x = element_text(family = "sans", colour = "black", size = rel(2))) + 
                                theme(axis.line = element_line(size = 1, colour = "black")) + 
                                theme(axis.ticks.length = unit(0, "cm")) + 
                                theme(axis.ticks.y = element_line(colour = "#222222")) + 
                                theme(axis.ticks.x = element_line(colour = "#222222")) + 
                                theme(axis.ticks.length = unit(0.4, "cm")) + 
                                theme(axis.title.x = element_text(family = "sans", size = rel(2))) + 
                                theme(axis.title.y = element_text(family = "sans", size = rel(2))) + 
                                guides(fill = guide_legend(override.aes = list(colour = NULL))) +
                                guides(fill = guide_legend(nrow = 1, byrow = TRUE))
BarPlot.Intron_length

################################################################################

####################
## Save all plots ##
####################

ggsave(BarPlot.NGenes, file="02-NGenes_Per_Contig.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(BarPlot.NExons, file="02-NExons_Per_Contig.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(BarPlot.NExons_Per_Gene, file="02-NExons_Per_Gene.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
ggsave(BarPlot.Intron_length, file="02-Intron_length.tiff", device = "tiff", dpi = 300, bg = NULL,
       width = 25, height = 20, units = "cm")
