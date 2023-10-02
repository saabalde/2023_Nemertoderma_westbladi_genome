## Quality control
Once decontaminated, it is time to evaluate the quality of the genome. I used several tools to get different metrics.

First, I used [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) to get a quick idea about the contiguity of the genome. It basically returns some basic metrics such as number of contigs, N50, longest contig size, and so on. I find it very useful to very quickly evaluate genome assemblies and it is very easy to run:

    assembly-stats pt_087_001_flye20211205meta.Metazoa.Clean.fasta > pt_087_001_flye20211205meta.Metazoa.Clean.log
    # Note that you can simply print the result on the screen if you remove the "> pt_087_001_flye20211205meta.Metazoa.Clean.log"

Second, I used [QUAST](https://github.com/ablab/quast) to obtain a more detailed analysis. I tried to include the _P. naikaiensis_ genome as a reference, but I did not see any remarkable difference (they are too divergent) so I won't include it here. The no-reference script looks like this:

    time quast.py pt_087_001_flye20211205meta.Metazoa.Clean.fasta --pacbio pt_087_001_trimmed_reads.fastq.gz \
                  --gene-finding --eukaryote \
                  --large --circos --k-mer-stats \
                  --threads 8 -o QUAST_output

These two analyses circle around genome contiguity, size, coverage, etc. Unfortunately, it looks like the genome is quite fragmented. There are several explanations for this, such as low coverage and low DNA quality, but it is hard to guess just from these results. However, contig size is not everything. Having a good, complete genome is also important and can be very useful for many analyses. In this case, we used [BUSCO](https://busco.ezlab.org/) with the Metazoa and Eukaryota databases to test genome completeness:

    run_BUSCO.py -i pt_087_001_flye20211205meta.Metazoa.Clean.fasta -l $BUSCO_LINEAGE_SETS/metazoa_odb10 -o run_Nwestbladi_clean_busco_odb10.metazoa --long -m genome -c 4
    run_BUSCO.py -i pt_087_001_flye20211205meta.Metazoa.Clean.fasta -l $BUSCO_LINEAGE_SETS/eukaryota_odb10 -o run_Nwestbladi_clean_busco_odb10.eukaryota --long -m genome -c 4

It looks like we got a quite complete genome! As you can see in the figure, the number of missing genes is relatively low in both datasets. It could be better, of course, but it is similar to any other published acoelomorph genome, so this might be inherent to these animals.

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/03-Quality_control/BUSCO_results.png)

These are pretty basic analyses that you can find in any genome report. Moving into a little more complexity, we also analyzed the kmer spectra of the genome. First, we inferred the kmer spectra with [Jellyfish](https://github.com/gmarcais/Jellyfish) and then used [KAT](https://github.com/TGAC/KAT) to analysed it:

    # Calculate the kmer spectra from the genome and the reads
    jellyfish count -C -m 21 -s 2000000000 -t 8 pt_087_001_flye20211205meta.Metazoa.Clean.fasta -o kmers_genome.jf
    jellyfish count -C -m 21 -s 2000000000 -t 8 pt_087_001_trimmed_reads.fastq.gz -o kmers_reads.jf

    # Convert the kmer hashes to histograms
    jellyfish histo -t 8 kmers_reads.jf > kmers_reads.histo
    jellyfish histo -t 8 kmers_genome.jf > kmers_genome.histo
    
    # Compare the two spectra with KAT
    kat comp -t 8 -o QC-Reads_vs_Assembly kmers_reads.jf kmers_genome.jf
    
    # Plot the individual spectra
    kat plot spectra-hist -o QC-Reads_spectra kmers_reads.histo
    kat plot spectra-hist -o QC-Genome_spectra.png kmers_genome.histo

As you can see in the figure below, the results are not fantastic. The curve looks basically like an exponential function, with too many unique (or, at least, repeated too few times) K-mers. This could be due to low coverage in the assembly, which we also suspected from other analyses.

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/03-Quality_control/KAT-Reads_vs_Assembly.png)

Finally, we compared our genome to the already published _P. naikaiensis_ and _S. roscoffensis_. Thus far, we have included an old _M. westbladi_ genome sequenced with Illumina to compare the performance of PacBio with respect to short-read approaches. I created the Rscripts 'Genome_Length.R' and 'Calculate_stats.R' for that matter. The former compares the assembly and contig lengths (before and after decontamination), whereas the former uses the annotation files to compare several metrics such as the number of genes per contig, number of exons per gene, or intron length. The two scripts are fully commented, but you can see in the figures below how, despite the worse contiguity of our genome, all the metrics analyzed are comparable.

Note: unfortunately, the annotation files are too big for GitHub so you will have to download them from somewhere else. I include all the results in the next figure and hope you still find the script useful.

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/03-Quality_control/Comparisons-Genome_stats.png)

Overall, we are satisfied with our results. It is true that we could have a much more contiguos genomes, such fragmentation precludes analyses of synteny, but we have a pretty complete genome and, more importantly, the annotated proteins are complete. It will definitely help us understand better the evolution of acoelomorph genomes.

---
