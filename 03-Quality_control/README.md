## Quality control
Once decontaminated, it is time to evaluate the quality of the genome. I used several tools to get different metrics.

First, I used [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) to have a quick idea about the contiguity of the genome. It basically returns some basic metrics such as number of contigs, N50, longest contig size, and so on. I find it very useful to very quickly evaluate genome assemblies and it is very easy to run:

    assembly-stats pt_087_001_flye20211205meta.Metazoa.Clean.fasta > pt_087_001_flye20211205meta.Metazoa.Clean.log
    # Note that you can simply print the result on the screen if you remove the "> pt_087_001_flye20211205meta.Metazoa.Clean.log"

Second, I used [QUAST](https://github.com/ablab/quast) to obtain a more detailed analysis. I tried to include the _P. naikaiensis_ genome as a reference, but I did not see any remarkable difference (they are too divergent) so I won't include it here. The no-reference script looks like:

    time quast.py pt_087_001_flye20211205meta.Metazoa.Clean.fasta --pacbio pt_087_001_trimmed_reads.fastq.gz \
                  --gene-finding --eukaryote \
                  --large --circos --k-mer-stats \
                  --threads 8 -o QUAST_output

These two analyses circle around genome contiguity, size, coverage, etc. Unfortunately, it looks the genome is quite fragmented. There are several explanations for this, such as low coverage and low DNA quality, but it is hard to guess just from these results. However, contig size is not everything. Having a good, complete genome is also important and can be very useful for many analyses. In this case, we used [BUSCO](https://busco.ezlab.org/) with the Metazoa and Eukaryota databases to test genome completeness:

    run_BUSCO.py -i pt_087_001_flye20211205meta.Metazoa.Clean.fasta -l $BUSCO_LINEAGE_SETS/metazoa_odb10 -o run_Nwestbladi_clean_busco_odb10.metazoa --long -m genome -c 4
    run_BUSCO.py -i pt_087_001_flye20211205meta.Metazoa.Clean.fasta -l $BUSCO_LINEAGE_SETS/eukaryota_odb10 -o run_Nwestbladi_clean_busco_odb10.eukaryota --long -m genome -c 4

It looks like we got a quite complete genome! As you can see in the figure, the numer of missing genes is relatively low in both datasets. It could be better, of course, but it is similar to any other published acoelomorph genome, so this might be inherent to these animals.

![image](https://github.com/saabalde/2023_Nemertoderma_westbladi_genome/blob/main/03-Quality_control/BUSCO_results.png)

XXX
    
