## Genome annotation
The annotation of the genome was performed with [BRAKER](https://github.com/Gaius-Augustus/BRAKER), following standard practices. About the same time I was running this protocol, the genome of the chiton [Hanleya hanleyi](https://f1000research.com/articles/11-555) came out, which led me to find [this script](https://github.com/kmkocot/GenomeAnnotation/blob/main/BRAKER2_pipeline.sh) by Kevin M. Kocot, from the University of Alabama. It helped me understand some of the steps and the parameters a litte better, so I share it here in case you find it useful as well.

Moving on, genome annotation consisted on four main steps: (1) repeat masking, generating (2) RNA and (3) protein evidences, and (4) the annotation itself.

### 1. Repeat masking
I first used [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/tree/master) to identify the repeats in the genome, and then [RepeatMasker](https://github.com/rmhubley/RepeatMasker) to locate them and soft-mask them (i.e. convert the repeats from upper to lowercase, so the nucleotides remain there but the programs know they should avoid these regions if possible). Installing these programs can be tricky because they have a lot of requisites, so I hope you have access to a pre-installed version somewhere.

The first step is to create a database of repeats specific to your genome:

    BuildDatabase -engine ncbi -name "Nwestbladi_repeats" pt_087_001_flye20211205meta.Metazoa.Clean.fasta
    
RepeatModeler will use this database to locate all repeats in the genome. Note that I use the -LTRStruct option to activate the search of LTR retrotransposons:

    RepeatModeler -pa 10 -engine ncbi -LTRStruct -database Nwestbladi_repeats 2>&1 | tee repeatmodeler.log

RepeatModeler will take some time to run, but it will create two interesting files: (1) Nwestbladi_repeats-families.fa, a fasta file containing the repeat sequences, and (2) Nwestbladi_repeats-families.stk, a similar file but in Stockholm format. RepeatMasker will use the fatsa file to locate these sequences in the genome and mask them:

    RepeatMasker -parallel 20 -engine rmblast -lib Nwestbladi_repeats-families.fa pt_087_001_flye20211205meta.Metazoa.Clean.fasta -s -xsmall -alignments
    
That is basically it. The _N. westbladi_ is relatively repetitive, as almost 60% of the genomes was masked. You can check this and other stats in the file 'pt_087_CleanAssembly.nomtDNA.fasta.tbl'.
With this masker genome, we can now map the transcriptomes and the proteins to generate the evidences that BREAKER will use. The masked genome has ben stored in 'pt_087_CleanAssembly.nomtDNA.fasta.masked'.

### 2. RNA mapping
Ideally, we should have sequenced a transcriptome from the same organism that was used for genome sequencing. Unfortunately, we had to use the whole body to obtain enough DNA. Instead, we used the two published transcriptome from this species ([SRX1343819](https://www.ncbi.nlm.nih.gov/sra/SRX1343819[accn]) and [SRX5296516](https://www.ncbi.nlm.nih.gov/sra/SRX5296516[accn])). After downloading these transcriptomes, the first step is to quality trim the sequences, which I did in a two step approach (I only show one of the two files):

    # First, a light trim and adapter removal (if there is any) with Trimmomatic
    java -jar $PATH_TO_TRINITYv2.6.6/trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads 16 SRR2682004_1.fastq.gz SRR2682004_2.fastq.gz \
          SRR2682004_1P.fastq.gz SRR2682004_1U.fastq.gz SRR2682004_2P.fastq.gz SRR2682004_2U.fastq.gz \
          ILLUMINACLIP:TruSeq_adapters.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:35
    
    # Remove the unpaired reads and uncompress the paired ones
    rm SRR2682004_*U.fastq.gz
    gzip -d SRR2682004_*P.fastq.gz
    
    # Reformat the fastq files so they don't return errors
    for i in SRR2682004_*P.fastq; do sed -E '/\+/ s/\+SRR2682004.*/\+/g' ${i} > ${i}.tmp; mv ${i}.tmp ${i}; done
    for i in SRR2682004_*P.fastq; do sed -i 's/_forward//g' ${i}; sed -i 's/_reverse//g' ${i}; done
    
    # Finally, filter low quality reads with Prinseq. Briefly, remove reads with a mean quality below 25, with more than 25% N's, with low
    # complexity (minimum entropy 50), remove nucleotides from the 3' and 5' ends whose quality is under 30, and then remove all sequences
    # shorter than 75 nucleotides
    prinseq-lite -fastq SRR2682004_1P.fastq -fastq2 SRR2682004_2P.fastq \
                 -out_good SRR2682004_good -out_bad SRR2682004_bad \
                 -min_qual_mean 25 -ns_max_p 25 \
                 -lc_method entropy -lc_threshold 50 \
                 -trim_qual_left 30 -trim_qual_right 30 \
                 -min_len 75
                
    # Remove the unnecessary files and compress everything
    rm *bad* *singleton*
    gzip -9 *fastq

With the reads clean, we can map them to the genome. After a careful screening of the literature and tutorials, it seems that the best mapping algorithms are [STAR](https://github.com/alexdobin/STAR) and [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml). STAR is generally considered as faster and more accurate, but if RAM is a limiting factor then TopHat is the wat to go. In my case, I could use STAR.

Besides the transcriptomes and the genome, STAR can use alternative sources of evidence to improve the quality of the mapping. In this case, I used the complete genes identified by BUSCO in the previous step, following [this tutorial](https://darencard.net/blog/2020-07-23-augustus-optimization/#:~:text=BUSCO%20can%2C%20therefore%2C%20be%20co,it%20also%20utilizes%20multiple%20cores.) to conver the results to a GTF file. I share the script 'Prepare_GTF.sh', which is fully commented and is easy to follow but, please, note that the paths to the files might change. The output file is called ''Nwestbladi.133genes.gtf'.

Now, we can prepare the reference genome. Here it is important to note two things: I stored the reference genome in the "Genome" directory, and the parameters "--genomeSAindexNbases" and "--genomeChrBinNbits" are genome specific, so you should follow the instructions in the manual to calculate the best values for your genome.

    STAR --runThreadN 10 --runMode genomeGenerate --genomeDir Genome --genomeFastaFiles Genome/pt_087_001_flye20211205meta.Metazoa.Clean.fasta.masked \
         --sjdbGTFfile Nwestbladi.133genes.gtf --sjdbOverhang 100 --genomeSAindexNbases 13 --genomeChrBinNbits 15

Once the genome is indexed, we can map the reads to the genome:

    STAR --runThreadN 10 --genomeDir Genome --readFilesIn SRR2682004_good_1.fastq.gz SRR2682004_good_2.fastq.gz --readFilesCommand zcat \
         --outFileNamePrefix SRR2682004 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 \
         --chimSegmentMin 40 --twopassMode Basic --limitBAMsortRAM 14500000000

### 3. Protein mapping
Selecting the right set of proteins is crucial for this step. the closer they are to your genome the better, since they will be easier to map and you will have more evidence for BRAKER. Unfortunately, the closest published genome to _N. westbladi_ is either _P. naikaiensis_ or _S. roscoffensis_, which diverged from _Nemertoderma_ more than 500 millions years ago. Thus, I will use three sources of proteins to maximize the changes of mapping the major amount of genes possible: (1) the Metazoa dataset from [OrthoDB](https://www.orthodb.org/?level=33208&species=33208), (2) the annotation of [_P. naikaiensis_](https://doi.org/10.1093/gigascience/giz023), and (3) a custom dataset of orthogroups inferred from [published transcriptomes](https://www.ncbi.nlm.nih.gov/sra/?term=Xenacoelomorpha).

After merging the three datasets into one fasta, I used [ProtHint](https://github.com/gatech-genemark/ProtHint) to map them to the genome:

    # The concatenated fasta file is called 'All_references.faa'
    
    # Remove redundancies
    cd-hit-est -i All_references.faa -o All_references.cdhit.0.90.faa -T 2 -M 10000 -c 0.95
    
    # Map them to the genome
    prothint.py --threads 10 --cleanup --evalue 1e-25 pt_087_001_flye20211205meta.Metazoa.Clean.fasta.masked All_references.cdhit.0.90.faa

However, the output of this analysis kept failing on BRAKER, so I decided to feed the 'All_references.faa' file directly. Let BRAKER run the mapping.

### 4. Annotation
Now that all files are prepared, it is as easy as running the BRAKER pipeline. The only thing that needs some consideration is that the installation might need some tuning, because the different versions require different prerequisites. For instance, from Augustus 3.3.3 to 3.4.0 changes the version of perl, which changes the packages that the program uses to parallelise the analysis. In the version 3.4.0, the one I used, the CRF training algorithm (in theory more accurate than HMM) does not work. In my case it is not a problem, as I have run several tests and haven't noticed any improvement from HMM to CRF, but it is worth to keep this in mind.

To run the program, simply type:

    braker.pl --cores 10 --etpmode --verbosity=4 --softmasking --filterOutShort  \
              --gff3 --min_contig=20000 --species=Nemertoderma_westbladi \
              --genome pt_087_001_flye20211205meta.Metazoa.Clean.fasta.masked \
              --bam=SRR2682004Aligned.sortedByCoord.out.bam,SRR8491950Aligned.sortedByCoord.out.bam \
              --prot_seq=All_references.cdhit.0.90.faa

And that's it.

---
