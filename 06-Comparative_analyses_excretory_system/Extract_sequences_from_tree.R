setwd("/home/saabalde/Escritorio/Nwestbladi_genome/Ultrafiltration_genes/")

library(ape)

my_tree <- read.tree("AllSpecies_and_Genes.diamond.fas.treefile")

getMRCA(my_tree, tip = c("Cnidaria_Scyphozoa_Rhopilema_esculentum_Proteins_mRNA.RE14024", 
                      "Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf11609_2_9"))
clade_01 <- extract.clade(my_tree, node = 868)

getMRCA(my_tree, tip = c("Arthropoda_Crustacea_Idotea_baltica_Genome_Proteins_MCL4144777_1", 
                         "Mollusca_Bivalvia_Pecten_maximus_Proteins_XP_033727763_1"))
clade_02 <- extract.clade(my_tree, node = 1492)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf8278_4_9", 
                         "Xenacoelomorpha_Acoela_Praesagittifera_naikaiensis_Proteins_g26370_t1"))
clade_03 <- extract.clade(my_tree, node = 1572)

getMRCA(my_tree, tip = c("Arthropoda_Crustacea_Idotea_baltica_Genome_Proteins_MCL4130465_1", 
                         "Platyhelminthes_Rhabditophora_Macrostomum_lignano_Proteins_PAA87012_1"))
clade_04 <- extract.clade(my_tree, node = 1623)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf11952_18_6", 
                         "Cnidaria_Scyphozoa_Rhopilema_esculentum_Proteins_mRNA.RE05182"))
clade_05 <- extract.clade(my_tree, node = 1072)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf4976_1_4", 
                         "Echinodermata_Asteroidea_Plazaster_borealis_Proteins_KPB_00001379-RA"))
clade_06 <- extract.clade(my_tree, node = 1143)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf11327_5_2", 
                         "Xenacoelomorpha_Nemertodermatida_Nemertoderma_westbladi_Proteins_jg6601_t1"))
clade_07 <- extract.clade(my_tree, node = 1204)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf9524_0_17", 
                         "Cnidaria_Scyphozoa_Rhopilema_esculentum_Proteins_mRNA.RE16409"))
clade_08 <- extract.clade(my_tree, node = 1248)

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf1584_4_9", 
                         "QRF78303.1"))
clade_09 <- extract.clade(my_tree, node = 1300)

all_branches <- c(clade_01$tip.label, clade_02$tip.label, clade_03$tip.label, clade_04$tip.label,
                  clade_05$tip.label, clade_06$tip.label, clade_07$tip.label, clade_08$tip.label,
                  clade_09$tip.label, "Arthropoda_Crustacea_Idotea_baltica_Genome_Proteins_MCL4127324_1",
                  "Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf8763_1_6",
                  "Mollusca_Bivalvia_Pecten_maximus_Proteins_XP_033727046_1", 
                  "Mollusca_Bivalvia_Pecten_maximus_Proteins_XP_033727047_1", 
                  "Mollusca_Bivalvia_Pecten_maximus_Proteins_XP_033727048_1", 
                  "Mollusca_Bivalvia_Pecten_maximus_Proteins_XP_033727049_1")

all_branches_but_QRF <- all_branches[-grep("QRF78", all_branches)]

#################

eya <- clade_07$tip.label
six1_2 <- clade_06$tip.label
pou3 <- clade_05$tip.label
sall <- clade_04$tip.label
lhx1_5 <- clade_08$tip.label
osr <- clade_02$tip.label

getMRCA(my_tree, tip = c("Annelida_Polychaeta_Paraescarpia_echinospica_Proteins_PE_Scaf11609_2_9",
                         "Chordata_Craniata_Homo_sapiens_Proteins_NP_001273278_1"))
clade_nephrin <- extract.clade(my_tree, 877)

nephrin <- clade_nephrin$tip.label
kirre <- setdiff(clade_01$tip.label, clade_nephrin$tip.label)
zo1 <- clade_03$tip.label

