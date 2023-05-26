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

### How clean is the genome?
Exloring the webserver you will realize BlobTools has generated several figures and tables, full of interesting information. We all recognise this figure, which summarises the taxonomic information and coverage of all contigs in the assembly: the 'blobplot'






I used the '70517384-88fb-4a33-acee-82461b0bef4d.csv' to explore the taxonomy of all contigs and identify which ones I need to remove




## This is the final tree inferred with IQ-TREE

All nodes have maximum statistic support unless otherwise specified.
![image](../Phylogenomics/IQtree_437genes.png)

Modified version of Fig. 3.

---

## Identify the contaminant contigs


---
