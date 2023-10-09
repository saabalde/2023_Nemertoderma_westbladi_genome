## Read filtering and assembly
After library prep, the _Nemertoderma_ genome was sequenced on HiFi. The cleaning and assembly processes were pretty straightforward. First, the ‘Trim gDNA Amplification Adapters’ pipeline from [SMRT Link v11](https://www.pacb.com/wp-content/uploads/SMRT_Link_User_Guide_v11.0.pdf) was used to remove the adapter sequences. Three genome assembly strategies were attempted and compared: the [IPA HiFi Genome Assembler](https://github.com/PacificBiosciences/pbipa), [Hifiasm v.0.7](https://github.com/chhylp123/hifiasm)⁠, and [Flye v.2.8.3](https://github.com/fenderglass/Flye). Default parameters were used in all three approaches, but with the “--meta” option activated in Flye. We selected the Flye assembly in downstream analyses because of its better completeness report.

More importantly, after the first BUSCO assembly, we detected a lot of duplication. To remove redundant contigs, we used [kmerDedup](https://github.com/xiekunwhy/kmerDedup). First, we extracted all the K-mers from the assembly and the cleaned reads using [Jellyfish](https://github.com/gmarcais/Jellyfish):

    # Count the kmers from the genome
    jellyfish count -C -m 21 -s 10G -t 10 pt_087_001_flye20211205meta.fasta -o pt_087_001_flye20211205meta.fasta.count
    # Count the kmers from the reads
    jellyfish count -C -m 21 -s 10G -t 10 pt_087_001_trimmed_reads.fastq -o pt_087_001_trimmed_reads.count

    # Merge them
    jellyfish merge -o kmer_genome_reads.jf *.count 

Before doing anything, we checked te merging was successful by comparing some stats

    # Calculate some stats from the kmer counts
    jellyfish stats -o pt_087_001_flye20211205meta.fasta.count.stats pt_087_001_flye20211205meta.fasta.count
    jellyfish histo -t 8 -o pt_087_001_flye20211205meta.fasta.count.histo pt_087_001_flye20211205meta.fasta.count
    
    jellyfish stats -o pt_087_001_trimmed_reads.count.stats pt_087_001_trimmed_reads.count
    jellyfish histo -t 8 -o pt_087_001_trimmed_reads.count.histo pt_087_001_trimmed_reads.count
    
    jellyfish stats -o kmer_genome_reads.jf.stats kmer_genome_reads.jf
    jellyfish histo -t 8 -o kmer_genome_reads.jf.histo kmer_genome_reads.jf

Create a fasta file with all the kmers:

    # Dump the kmers
    jellyfish dump -c -t -o kmer_genome_reads.jf.dump kmer_genome_reads.jf

    # Witj kmerDedup, remove kmers present less than three times and create a fasta file with the others
    perl kmerFilter.pl -d kmer_genome_reads.jf.dump -o kmer_genome_reads.jf.fa -l 3 

Now, map the K-mers back to the assembly using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml#:~:text=Bowtie%202%20is%20an%20ultrafast,long%20(e.g.%20mammalian)%20genomes.) and calculate coverage statistics with [BamDeal](https://github.com/BGI-shenzhen/BamDeal). Reformat the genome first, just in case:

    # Reformat the genome with kmerDedup
    perl kmerDedupfa2fa.pl -f pt_087_001_flye20211205meta.fasta -o pt_087_001_flye20211205meta.fasta.reformat.fasta -c F -n F
    
    # Map the kmers to the genome (pipe it to samtools to convert the SAM to BAM)
    bowtie2 --very-sensitive -k 1000 --score-min L,-0.6,-0.2 --end-to-end --reorder -L 21 -p 10 -f kmer_genome_reads.jf.fa -x pt_087_001_flye20211205meta.fasta.indices | samtools view -@ 4 -F 4 -b -S -o pt_087_001_flye20211205meta.fasta.bam
    
    ## Run BamDeal
    BamDeal_Linux statistics Coverage -i pt_087_001_flye20211205meta.fasta.bam -r pt_087_001_flye20211205meta.fasta.reformat.fasta -q 0 -o pt_087_001_flye20211205meta.fasta.cov

Finally, run kmerdedup to remove all redundancies. To find the best parameters, we ran three jobs in parallel changing the max duplication percentage (-mpr):

    ## Run kmerDedup
    # Remove redundancies from the genome
    perl kmerDedup.pl -k pt_087_001_flye20211205meta.fasta.Dedup -o 02.kmerDedup \
                      -mpr 0.2 -mcv 30 -kmer 21 \
                      -f pt_087_001_flye20211205meta.fasta.reformat.fasta -bam pt_087_001_flye20211205meta.fasta.bam -cov pt_087_001_flye20211205meta.fasta.cov.stat \
                      -s /home/bin/samtools
    
    perl kmerDedup.pl -k pt_087_001_flye20211205meta.fasta.Dedup -o 03.kmerDedup \
                      -mpr 0.3 -mcv 30 -kmer 21 \
                      -f pt_087_001_flye20211205meta.fasta.reformat.fasta -bam pt_087_001_flye20211205meta.fasta.ba -cov pt_087_001_flye20211205meta.fasta.cov.stat \
                      -s /home/bin/samtools
    
    perl kmerDedup.pl -k pt_087_001_flye20211205meta.fasta.Dedup -o 05.kmerDedup \
                      -mpr 0.5 -mcv 30 -kmer 21 \
                      -f pt_087_001_flye20211205meta.fasta.reformat.fasta -bam pt_087_001_flye20211205meta.fasta.ba -cov pt_087_001_flye20211205meta.fasta.cov.stat \
                      -s /home/bin/samtools

Once finished, we used BUSCO to find the best assembly: the one with the minimum duplicates while keeping the best completeness. We decided to continue with a max duplication percentage of 30% (-mpr 0.3):

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/01-Read_filtering_and_assembly/03-BUSCO_scores-Metazoa.png)

---
