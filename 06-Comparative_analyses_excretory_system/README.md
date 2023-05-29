## Evolution of the excretory system
Xenacoelomorphs are morphologically simple worms, lacking several structures typical of bilaterian animals such as a circulatory system, a through-gut, or an excretory system. Regardless if this is because they are an [early bilaterian](https://www.nature.com/articles/nature16520) or because of [secondary losses](https://www.sciencedirect.com/science/article/pii/S0960982219304075), they are an important model for understanding the evolution of animals' body plans. Interestingly, [Andrikou et al. (2019)](https://doi.org/10.1371/journal.pbio.3000408) found that, unlike what was originally thought, xenacoelomorphs do have some sort of active excretion through the digestive system. thus far, they annotated several genes known to participate in the excretory function of nephrozoan animals (i.e., [the clade including protostomes and deuterostomes, animals with an excretory system](https://doi.org/10.1046/j.1463-6409.2002.00090.x)). A couple of years later, [Gąsiorowski et al., (2021)](https://doi.org/10.1016/j.cub.2021.05.057) described the molecular machinery underlying the development of the excretory organs in both protostomes and deuterostomes. 
Building on these findings, we decided to investigate the presence of these genes in Xenacoelomorpha, as well as to compare their general structure in the genomes of these animals.

As in other analyses, the set of genomes used are listed in 'SRA_genomes.xlsx', but in these cases those whose GFF annotation file is available were considered. the reason is simple: we can study presence/absence in all of them, but not the features of their exons.

Gąsiorowski et al. (2021) focused on nine genes: three structural proteins (_Nephrin_, _Kirrel_, and _ZO1_) and six transcription factors (_eya_, _Six1/2_, _POU3_, _Sall_, _Lhx1/5_, and _Osr_). They gene _Hunchback_ was also studied, but it was not found to participate in nephridiogenesis and thus was not considered here. 
We will need:
1. The [proteins](https://www.ncbi.nlm.nih.gov/protein/?term=A+single+origin+of+animal+excretory+organs) annotated by Gąsiorowski et al., as reference. I named this file 'Ultrafiltration_genes.fasta'.
2. The proteomes, which can be downloaded from the NCBI.
3. The GFF annotation files, which are also available online.

The first thing we need to do is to prepare the files for [diamond](https://github.com/bbuchfink/diamond).

    # First, concatenate the proteomes into a single file. I will just copy them from the gene content directory. Remember that we modified the headers
    # to include the name of the species
    cat ../05-Gene_content_evolution/02-CD_hit/* > AllSpecies.faa
    
    # Create a diamond database.
    diamond makedb -p 16 --in Ultrafiltration_genes.fasta -d Ultrafiltration_genes.dmnd
    
    # Run diamond
    diamond blastp --query AllSpecies.fasta --db Ultrafiltration_genes.dmnd --outfmt 6 \
                   --sensitive --max-target-seqs 3 --evalue 1e-5 --threads 4 > AllSpecies.diamond.out

This file will contain all hits, but I want to separate them per gene:

    # For each gene, create a list with the accesion numbers. For instance:
    grep 'Osr' Ultrafiltration_genes.fasta | sed 's/>//g' | awk '{print $1}' > Ultrafiltration_genes.Osr.list
    
    # Now, use these lists to separate the genes. Note that the ooutput will include both the query and the hit sequences
    for gene in $( ls Ultrafiltration_genes.*.list | sed 's/\.list//g' )
        do
        while read id
            do
            grep "$id" AllSpecies.diamond.out >> ${gene}.diamond.out
            grep "$id" AllSpecies.diamond.out | awk '{print $1} | sort | uniq >> ${gene}.diamond.list
            grep "$id" AllSpecies.diamond.out | awk '{print $2} | sort | uniq >> ${gene}.diamond.list
        done < ${gene}.list
    done

Now, create a fasta file for each gene using [seqtk](https://github.com/lh3/seqtk):

    # First, concatenate all sequences into a single file
    cat AllSpecies.fasta Ultrafiltration_genes.fasta > AllSpecies_and_Genes.fasta
    
    # Now, extract the sequences using seqtk
    for gene in $( ls Ultrafiltration_genes.*.diamond.list | sed 's/\.list//g' )
        do
        seqtk subseq AllSpecies_and_Genes.fasta ${gene}.list > ${gene}.fasta
    done

Align the sequences in each fasta file with [MAFFT](https://mafft.cbrc.jp/alignment/server/). Note that I will align all sequences from all genes together. I decided to separate them per gene first so they are sorted in the final file, which makes it easier to rename them to include the proposed gene name. Then, infer a phylogenetic tree from this alignment using [IQ-TREE](http://www.iqtree.org/doc/Tutorial):

    cat Ultrafiltration_genes.*.diamond.fasta > Ultrafiltration_genes.AllSpecies.fasta
    
    # And align them all
    mafft-linsi --reorder --thread 12 Ultrafiltration_genes.AllSpecies.fasta > Ultrafiltration_genes.AllSpecies.fas
    
    # Infer the tree
    

The reason I inferred this tree is simple: all sequences not clustered with any of the ultrafiltration genes will likely be false positives that need to be removed. I wish I had a script to automate this, but unfortunately this step involves some manual labour. Open the tree, draw it as you prefer, and locate the GenBank sequences (names QRF...). You need to focus on these nodes, go back in the tree, blasting the sequences from preceding nodes, until you stop getting the right hits (i.e. ultrafiltration genes). Once you have defined the nodes of interest, select two species whose most recent common ancestor (MRCA) is the ancestor of the full node. I provide the R script 'Extract_sequences_from_tree.R', which uses the R library [ape](https://cran.r-project.org/web/packages/ape/index.html) to find the MRCA of the provided sequences and print tall the tip labels included in that node.
Finally, I put all these sequences in a file, which I used to: (1) create a fasta file with seqtk, as above; and (2) blast each of these sequences, one by one, against the full NCBI database. There are still *some* false positives in this list, and by blasting against the full NCBI database we will be able to separate them.

At this point, you will have a complete list of all the ultrafiltration genes present in your genomes. If you want, you can make a new tree from this filtered fasta file to make sure all genes are monophyletic:

    # Align
    mafft-linsi --reorder --thread 12 AllSpecies_and_Genes.Ultrafiltration_v2.fasta > AllSpecies_and_Genes.Ultrafiltration_v2.fas
    
    # Infer the tree
    iqtree -s AllSpecies_and_Genes.Ultrafiltration_v2.fas -m LG -bb 1000 -alrt 1000 -bnni -nt AUTO -safe

This is the result. You will see some sequences are misplaced, but this is likely an artifact during tree inference. the tree looks generally fine:

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/06-Comparative_analyses_excretory_system/AllSpecies_and_Genes.Ultrafiltration_v2.fas.treefile.png)

Now, we can use the names of the genes to generate all the metrics we want to analyse. This is fairly easy. First, calculate the length of all sequences

    awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../AllSpecies_and_Genes.Ultrafiltration.fasta | tr '\n' ';' | tr '>' '\n' | tr ';' '\t' > AllSpecies_and_Genes.Ultrafiltration.seq_length

Second, calculate the number of exons per gene:

    # Create a list with the names of the sequences
    grep '>' AllSpecies_and_Genes.Ultrafiltration_v2.fasta | sed 's/>//g' > AllSpecies_and_Genes.Ultrafiltration.list
    
    # Concatenate all GFF files into one
    cat *gff* > all_species.gff
    
    # Find the name of the name of the sequence on the GFF file and count how many exons it has
    while read id
        do
        echo ${id}
        species=$( echo ${id} | sed 's/_Proteins_.*//g' )
        gene=$( echo ${id} | sed 's/.*_Proteins_//g' )
        n_exons=$( grep "${gene}" all_species.gff | grep -c "CDS" )
        printf "%s\t%s\t%s\n" "${species}" "${gene}" "${n_exons}" >> AllSpecies_and_Genes.Ultrafiltration.N_Exons
    done < AllSpecies_and_Genes.Ultrafiltration.list

Third, use the coordinates to extract the length of each exon:

    while read id
        do
        echo ${id}
        species=$( echo ${id} | sed 's/_Proteins_.*//g' )
        gene=$( echo ${id} | sed 's/.*_Proteins_//g' )
        grep "${gene}" all_species.gff | grep "CDS" > tmp; awk -v var="$id" '{print var,$1,$4,$5,$9}' tmp >> AllSpecies_and_Genes.Ultrafiltration.All_Exons
        rm tmp
    done < AllSpecies_and_Genes.Ultrafiltration.list

Note that the 'AllSpecies_and_Genes.Ultrafiltration.All_Exons' list is much, much longer than the others. This is because we have calculated the length of all exons, use it to calculate the average exon length for each gene (I did this in a spreadsheet, but I'm sure it should be relatively easy to do on the command line). With all three tables, I built the attached spreadsheet 'Ultrafiltration_genes.AllSpecies.txt', with 13 columns: Main_clade, Acoelomorph_or_not, Phylum, Species, Sequence, Blast, Blast_modified, Structural_or_Transcription, Prot_length, N_exons, Av_exon_length, Av_intron_length, Gene_length.
I do not include it here becuse there is too much variance in the data, with thousands of base pairs difference between an intron and the other, but I also calculated the average intron length (the number of base pairs between exons) and the gene length (the number of base pairs between the firsr nucleotide of the first exon and the last nucleotide of the last exon). you can ignore them.

To analyze this table, I created two R scripts:
First, I used 'Ultrafiltration_genes.AllSpecies.R' to visualize the results, and have a glimpse of the variation in the data. The script is commented so I won't delve into careful explanations here, but it basically plots a boxplot for each of the metric of interest, separating the data for clade (Acoelomorpha, Cnidaria, Deuterostomia, and Protostomia) and type of gene (Structural gene or transcription factor). These are the final plots:

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/06-Comparative_analyses_excretory_system/Ultrafiltration_genes.AllSpecies.png)

Second, given the huge variation observed I decided to check the differences in this gene and just plot those significantly different. Very briefly, I used the the Shapiro-Wilk’s method and the Barlett test were used to check if they follow a normal distribution and the homogeneity of their variances, respectively. For each gene, the differences among clades were tested with either an ANOVA or a Kruskal-Wallis test, depending on the result of the normality and homoscedasticity tests. Finally, the Bonferroni correction (ANOVA) and the Dunn test (Kruskal-Wallis) were selected to run pairwise comparisons in all cases identified as statistically different.

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/06-Comparative_analyses_excretory_system/Ultrafiltration_genes.AllSpecies.Significances_Summary.png)

As you can see, many of the statistically significant differences (summarize in brackets) are related to Acoelomorpha. This is interesting, because it indicates the relativley conserve structure that these genes have in animals with an excretory system. However, we need to run a final test. Some authors have proposed 





Ultrafiltration_genes.AllSpecies.Significances_Summary.R
