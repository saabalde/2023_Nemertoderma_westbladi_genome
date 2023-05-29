## Analysis of gene content
Analysis of gene content is an interesting approach to understanding animal evolution. In particular, it would be interesting to know whether the morphological simplicity of Xenacoelomorpha is somehow correlated with some form of genomic simplicity. Previous studies (see [this](https://www.biorxiv.org/content/10.1101/2022.06.24.497508v2.full) and [this](https://doi.org/10.1126/science.aau6173)) have shown that both _Xenoturbella_ and Acola have a similar number of genes than other bilaterians, so this might not be true. Likewise, we see a similar pattern in [our annotation](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/tree/main/04-Annotation) of the _N. westbladi_ genome. However, it would be interesting to see if, beyond the number of genes, we find the same gene families present in other animals. To test this, we have compared the gene content of _N. westbladi_ to the genomes listed in the attached file 'SRA_genomes.xlsx'. You can download their proteomes from GenBank.

### 1. Filter redundancies in the proteomes
It is common to find duplicates of the same protein or recently duplicated paralogs (with very little variation from the original gene) in the annotation of a genome. Hence, one of the first things we can do is to remove this redundancy. We will be looking at presence/absence of gene families, so this is not particularly important, but reducing the number of proteins can signifcantly reduce computing time. I have used [CD-HIT](https://sites.google.com/view/cd-hit) to cluster sequences that are more than **XXXX**% similar.

    CODE HERE

### 2. Functional annotation
I have tried to cluster the proteins in orthogroups in several ways, but working just with functionally annotated proteins seems to return the best results. There are no big differences in the number of genes shared among phyla, but by comparing functionally annotated proteins the number of unique proteins (i.e. those present only in a given phylum or clade) is drastically reduced. To annotate the proteins, I have used [pfam_scan](https://anaconda.org/bioconda/pfam_scan) to compare all proteins to the [PFAM 35 database](http://pfam.xfam.org/).
Please, note that PFAM is discontinued. Now the go-to database is [InterPro](https://www.ebi.ac.uk/interpro/).

To annotate the sequences:

    for species in $( ls *cdhit.faa | sed 's/\.Proteins\.faa\.cdhit\.faa//g' )
        do
        pfam_scan.pl -fasta ${species}.Proteins.faa.cdhit.faa -dir DATABASE -outfile ${species}.pfam.out
    done

And to extract the annotated sequences:

    # This loop will: (1) read the PFAM output, (2) extract the names of the sequences with hit into a list, and 
    # (3) use the list to extract the sequences using seqtk
    for i in $( ls *pfam.out | sed 's/\.pfam\.out//g' )
        do
        echo ${i}
    
        awk '{print $1}' ${i}.pfam.out | grep -v '#' | sort | uniq > ${i}.pfam.out.list
        awk '{print $1}' ${i}.pfam.out | grep -v '#' | sort | uniq | wc -l
    
        seqtk subseq ../02-CD-HIT/${i}.Proteins.faa.cdhit.faa ${i}.pfam.out.list > ${i}.pfam.faa
        grep -c '>' ${i}.pfam.faa
    
        echo ""
    done

### 3. Infer gene families
From the set of annotated proteins, I will use [OrthoFinder](https://github.com/davidemms/OrthoFinder) to cluster them into orthogroups, which are thought to be conformed by proteins with a common origin together (ortholog sequences).

    # Create a directory for OrthoFinder
    cd 01-Create_orthogroups
    cp *pfam.aa 01-Create_orthogroups/
    
    # Run OrthoFinder
    CODE HERE

### 4. 




