## Decontamination
Unlike any of the other xenacoelomorph genomes, this genome was not sequenced from a cultured specimen. This means it is likely many contaminants were sequenced alongside _N. westbladi_, such as epidermic symbionts, the microbiome, or DNA suspended in the seawater. To filter out these contigs, we have used [BlobTools2](https://blobtoolkit.genomehubs.org/blobtools2/), which uses BLAST searches, coverage information, and the results from BUSCO to produce several metrics.

Herein I will assume all the relevant files are in the working directory:

### 1. Download the taxonomic information from the NCBI
We first need to download the taxonomic information from GenBank. We need this information to format the Uniprot database and later on to tell BlobTools the organisms we are working with:

    mkdir 01-NCBI_taxonomy
    cd 01-NCBI_taxonomy/
    curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
    cd ../

### 2. Download and format the Uniprot database
Now, we download the latest version of the Uniprot database and create a diamond database with the taxonomic information. I followed [this tutorial](https://github.com/blobtoolkit/blobtoolkit/issues/41) to make sure everything works:

    # First, set the location of the taxonomic information
    TAXDUMP=./01-NCBI_taxonomy/
    
    # Download Uniprot
    wget -q -O reference_proteomes.tar.gz \
         ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
         -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
         awk '/tar.gz/ {print $9}')
    
    # Extract the sequences
    tar xf reference_proteomes.tar.gz
    
    # Add the taxonomic information
    touch reference_proteomes.fasta.gz
    find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz
    printf "accession\taccession.version\ttaxid\tgi\n" > reference_proteomes.taxid_map
    zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map
    
    # Create the diamond database
    diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map \
                   --taxonnodes $TAXDUMP/nodes.dmp --taxonnames $TAXDUMP/names.dmp -d reference_proteomes.dmnd

### 3. Run the analyses to create the input files
BlobTools is a very flexible program, accepting many different files and formats and summarising everything is very visual, easy to understand figures. The first file I generated was the coverage information, mapping the filtered reads to the assembly using [minimap2](https://github.com/lh3/minimap2):

    mkdir 02-minimap2
    cd 02-minimap2/
    
    # This is a two step process: first, map the reads to the genome; then sort and convert the SAM file to BAM
    minimap2 -ax map-pb -t 16 ../pt_087_001_flye20211205meta.fasta ../pt_087_001_trimmed_reads.fastq.gz | samtools sort -@16 -O BAM -o assembly.reads.bam -
    
    cd ../

Second, I identified the origin of each contig by blasting them to the Uniprot database using [diamond](https://github.com/bbuchfink/diamond):

    mkdir 03-DIAMOND
    cd 03-DIAMOND/
    
    diamond blastx --query ../pt_087_001_flye20211205meta.fasta --db ../reference_proteomes.dmnd \
                   --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                   --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > diamond.out
    
    cd ../

Third, I used [BUSCO v.5](https://busco.ezlab.org/) and the Metazoa_odb10 database to assess the completeness of the assembly. Note that this is not necessary to decontaminate the genome, but it will include this information in the summary and will work and I think it is interesting:

    mkdir 04-BUSCO
    cd 04-BUSCO/
    
    run_BUSCO.py -i ../pt_087_001_flye20211205meta.fasta -l $BUSCO_LINEAGE_SETS/metazoa_odb10 -o pt_087_001_flye20211205meta --long -m genome -c 8
    
    cd ../
    
Now that all files are ready, run BlobTools:

    # I installed Blobtools in a conda environment that I need to load
    conda activate blobtools2
    
    # Create a Blobdir, called 'Nwestbladi', where to merge all the information
    blobtools create --fasta pt_087_001_flye20211205meta.fasta --taxid 172109 --taxdump ./01-NCBI_taxonomy/ Nwestbladi
    
    # Add the Diamond, Minimap2 and BUSCO information
    blobtools add --cov 02-minimap2/assembly.reads.bam Nwestbladi
    blobtools add --hits 03-DIAMOND/diamond.out --taxrule bestsumorder --taxdump ./01-NCBI_taxonomy/ Nwestbladi
    blobtools add --busco 04-BUSCO/pt_087_001_flye20211205meta/run_metazoa_odb10/full_table.tsv Nwestbladi
    
    # Connect to the webserver to visualize the results
    blobtools host `pwd`
    
    # When you are done, close everything and deactivate the environment
    conda deactivate

### 4. How clean is the genome?
Exloring the webserver you will realize BlobTools has generated several figures and tables, full of interesting information. We all recognise this figure, which summarises the taxonomic information and coverage of all contigs in the assembly: the 'blobplot'

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/02-Genome_decontamination/Nwestbladi_blobplot.png)

As you can see, Uniprot has identified a lot of different things in the genome. This was expected, so now we need to separate all the contigs that are not of interest. From the '70517384-88fb-4a33-acee-82461b0bef4d.csv' table, I will create a list with the names of all the contigs that are not 'Metazoa' (almost 20% of the assembled genome). It is likely that some of these will also be contaminants, but it is hard to tell just from this analysis. Xenacoelomorpha is not a very well represented phylum in the databases and this results might just be an artifact. I will check this in a saparate analysis.

From the '70517384-88fb-4a33-acee-82461b0bef4d.csv' table I created two files: Metazoa_contigs.txt and no_Metazoa_contigs.txt. Just use [seqtk](https://github.com/lh3/seqtk) to extract the contigs of interest:

    seqtk subseq pt_087_001_flye20211205meta.fasta Metazoa_contigs.txt > pt_087_001_flye20211205meta.Metazoa.fasta
    seqtk subseq pt_087_001_flye20211205meta.fasta no_Metazoa_contigs.txt > pt_087_001_flye20211205meta.no_Metazoa.fasta

### 5. Identify the sources of contamination
Using the no_Metazoa_contigs.txt file we can also find the relevant hits in the diamond output and identify, to the finest taxonomic level possible, the source of the contmainations:

    # Extract the relevant lines from the diamond output
    while read seq
        do
        grep -w "${seq}" ../03-DIAMOND/diamond.out >> diamond.out
    done < Relevant_contigs.txt
    
    # Extract the protein ID
    awk '{print $5}' diamond.out | sort | uniq | tr '|' '\t' | awk '{print $3}' > diamond.GeneID
    
    # Now, this is a bit tricky and, admitedly, dirty. The link between the gene ID and the taxonomic information is stored in 
    # many thousands of files. Each file corresponds to one gene, so we need to locate them.
    # First, put them all together in a directory and uncompress them
    cp */*/*idmapping.gz uncompress_files/
    gzip -d uncompress_files/*
    
    # Second, find the file about oyur genes of interes
    while read id
        do
        echo $id
        grep -wl "$id" uncompress_files/*
    done < diamond.GeneID
    
    # This can take a long time (>2 days if you do not parallelize). Now, since you have the file names, you can copy them to a 
    # new directory called 00-TaxID
    cat 00-TaxID/* | grep 'NCBI_TaxID' | awk '{print $3}' | sort | uniq > diamond.TaxID
    
    # And fourth, create a new file including the full taxonomic information of each TaxID
    while read id
        do
        classification=$( grep -w "^$id" fullnamelineage.dmp | sed -E 's/\t//g' | sed 's/\ //g' | tr '|' '\t' | awk '{print $3}' | sed 's/cellularorganisms\;//g' )
        printf "%s\t%s\n" "$id" "$classification"
    done < diamond.TaxID > diamond.TaxID_classification
    
    # Don't forget to delete the uncompress_files directory
    rm -r uncompress_files

Now, you can open 'diamond.TaxID_classification' in an excel file and parse it as you wish.

### 6. Filter mitochondrial contigs
Besides contaminants from other species, the assembly will also contain mitochondrial genomes, either native from _N. westbladi_ or from the contaminants themshelves. In either case, it is a good idea to remove them. To do so, create a custom database with the set of mitochondrial genes that you want. I limited it to Xenacoelomorpha, because all animal genomes have the same genes and it will work anyway, but you can also add other species if you think it will be more efficient.

    # First, create a database from the fasta file
    makeblastdb -in Xenacoelomorpha_mtDNA.fasta -out Xenacoelomorpha_mtDNA.fasta -dbtype nucl
    
    # Run blast
    blastn -query pt_087_001_flye20211205meta.Metazoa.fasta -db Xenacoelomorpha_mtDNA.fasta \
           -num_threads 16 -outfmt 6 -evalue 1e-5 -out Nwestbladi_mtDNA.txt
           
    # Extract the names of the mitochondrial contigs
    awk '{print $1}' Nwestbladi_mtDNA.txt | sort | uniq > Nwestbladi_mtDNA.Contigs.txt
    
    # Create a list with the names of the non-mitochondrial contigs
    grep '>' pt_087_001_flye20211205meta.Metazoa.fasta | sed 's/>//g' > Contig_names.txt
    
    while read contig
        do
        count=$( grep -c "${contig}" no_Metazoa_contigs.txt )
        if [ $count -lt 1 ]; then
            ${contig} >> Contig_names.Clean.txt
        fi
    done < Contig_names.txt
    
    # Extract the non-mitochondrial contigs
    seqtk subseq pt_087_001_flye20211205meta.Metazoa.fasta Contig_names.Clean.txt > pt_087_001_flye20211205meta.Metazoa.Clean.fasta

### 7. Clean the HiFi reads
Now that our genome is clean, we can use it as reference to separate contaminant from native reads using minimap2, [samtools](https://github.com/samtools/samtools), and [bam2fastq](https://github.com/jts/bam2fastq):

    # Map the reads to the clean genome
    minimap2 -ax map-pb -t 16 pt_087_001_flye20211205meta.Metazoa.Clean.fasta pt_087_001_trimmed_reads.fastq.gz | samtools sort -@16 -O BAM -o Assembly.reads.bam -
    
    # Extract the mapped reads
    samtools view -h -u -F 4 Assembly.reads.bam > Assembly.mapped_reads.bam
    
    # Extract the unmapped reads
    samtools view -h -u -f 4 Assembly.reads.bam > Assembly.unmapped_reads.bam
    
    # Sort the two BAM files
    samtools sort -n Assembly.mapped_reads.bam -o Assembly.mapped_reads.bam.sort
    samtools sort -n Assembly.unmapped_reads.bam -o Assembly.unmapped_reads.bam.sort
    
    # Convert the BAM files to fastq
    bam2fastq -o pt_087_CleanAssembly.mapped_reads#.fastq Assembly.mapped_reads.bam.sort
    bam2fastq -o pt_087_CleanAssembly.unmapped_reads#.fastq Assembly.unmapped_reads.bam.sort
    
    # Compress the fastq files
    gzip -9 pt_087_CleanAssembly.mapped_reads_M.fastq
    gzip -9 pt_087_CleanAssembly.unmapped_reads_M.fastq


## Make sure the metazoan contigs come from _Nemertoderma westbladi_
As outlined above, BlobTools identified the most likely source of all contigs in the assembly, but that analysis did not exclude metazoan contigs. Hence, we cannot know if our genome actually comes from _N. westbladi_ or if other animals were sequenced instead. To test this, I have compared our genome to published genomes from other phyla. The list of reference genomes can be found in the attached spreadsheet 'SRA_genomes.xlsx'. I have basically followed two approaches: first, I used kmer to calculate the distances among all genomes; second, I used IQ-TREE to infer a phylogenetic tree. In both cases, the goal is to make sure we recover our genome in a clade with the two Acoela genomes.

### 1. Genomic distances: Skmer
[Skmer](https://github.com/shahab-sarmashghi/Skmer) is a tool designed to extract all kmer from a set of genomes, compare them and calculate the genomic distances. Although this tool was originally designed to identify species from skimming data, a sort of barcoding but with genome-wide information, the same principle can be applied here. 


### 2. Phylogenetic inference
After annotating the _N. westbladi_ genome, we used this proteome to infer a phylogenetic tree, including the same species as above. We have not seen how to annotate the genome, yet, but you can check this step [here: genome annotation](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/tree/main/04-Annotation).
This was a relatively quick test, so I will not delve in depth into the reasoning of each analysis. If you are interested in learning a bit more about this, you can follow [this much more detailed explanation](https://github.com/saabalde/2023_Lventricosus_Species_Complex/tree/main/Phylogenomics). It's the same idea.

[Orthofinder](https://github.com/davidemms/OrthoFinder)

[MAFFT](https://mafft.cbrc.jp/alignment/software/)

[BMGE](https://doi.org/10.1186/1471-2148-10-210)

[FASconCat](https://github.com/PatrickKueck/FASconCAT-G)

[IQ-TREE tutorial](http://www.iqtree.org/doc/Assessing-Phylogenetic-Assumptions)

[ASTRAL](https://github.com/smirarab/ASTRAL)

[IQ-TREE](http://www.iqtree.org/doc/Tutorial)



---
