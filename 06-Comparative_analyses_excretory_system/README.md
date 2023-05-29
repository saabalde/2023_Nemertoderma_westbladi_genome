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

This is the result. You will see some sequences are misplaced, but this is likely an artifact during tree inference:

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/06-Comparative_analyses_excretory_system/AllSpecies_and_Genes.Ultrafiltration_v2.fas.treefile.png)




