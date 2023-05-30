## Script based on: 
## https://darencard.net/blog/2020-07-23-augustus-optimization/#:~:text=BUSCO%20can%2C%20therefore%2C%20be%20co,it%20also%20utilizes%20multiple%20cores.


## Create a text file containing the information for all the genes. This script will:
# 1) Create a txt file containing the columns "group, gene_name, orthodb_url, etc.
# 2) Open the file containing the ids and read them, one by pne (second column)
# 3) Call OrthoDB to download the information of each gene, and save it into the text file
cat \
<(echo -e "#group\tgene_name\torthodb_url\tevolutionary_rate\tmedian_exon_count\tstdev_exon_count\tmedian_protein_length\tstdev_protein_length") \
<(cat ./metazoa_odb10/scores_cutoff | cut -f 1 | \
while read id; \
do \
wget -qO - https://dev.orthodb.org/group?id=${id} | \
jq -r '. | [.data.id, .data.name, .url, .data.evolutionary_rate, .data.gene_architecture.exon_median_counts, .data.gene_architecture.exon_stdev_counts, .data.gene_architecture.protein_median_length, .data.gene_architecture.protein_stdev_length] | @tsv'; \
sleep 1s; \
done) > metazoa_odb10.info.txt


## Sort the list based on the average and standard deviation of the number of exons and the gene length
cat <(cat metazoa_odb10.info.txt | head -1) \
<(cat metazoa_odb10.info.txt | grep -v "^#" | sort -t $'\t' -k5,5nr -k7,7nr -k4,4n) > metazoa_odb10.info.ranked.txt


## We run BUSCO, returning this result: C:77.3%[S:55.6%,D:21.7%],F:4.1%,M:18.6%,n:954


## Now, from the BUSCO output we need to extract the Complete and Single-Copy genes
cat metazoa_odb10.info.ranked.txt | grep -v "^#" | cut -f 1 | \
while read id; \
do \
status=`cat ../../05-QC_BUSCO_CleanAssembly/05-BUSCO_odb10/run_Nwestbladi_clean_busco_odb10/full_table_Nwestbladi_clean_busco_odb10.tsv | awk -v id="${id}" '{ if ($1 == id) print $2 }' | sort | uniq`; \
file="../../05-QC_BUSCO_CleanAssembly/05-BUSCO_odb10/run_Nwestbladi_clean_busco_odb10/augustus_output/extracted_proteins/"${id}".faa.1"; \
count=`cat ${file} | grep -c ">"`; \
echo -e "${status}\t${count}\t${file}" | awk '{ if ($1 == "Complete" && $2 == 1) print $3 }'; \
done > Nwestbladi_complete_singlecopy.aafiles.txt


## 133 genes were extracted during the previous step. I was expecting 530 complete and single-copy genes, but some of them contained
## more than one sequence. This dataset is smaller, but we are certain that we are dealing with complete, single-copy genes.


## Now we rename the sequences in the fasta files. Augustus includes a "g1" symbol by default, and we want to remove it.
cat Nwestbladi_complete_singlecopy.aafiles.txt | \
while read file; \
do \
id=`basename ${file} .faa.1`; \
cat ${file} | seqkit fx2tab | awk -v id="${id}" -v OFS="\t" '{ print id, $2 }' | seqkit tab2fx; \
done > Nwestbladi_complete_singlecopy.seqs.rename.fasta


## Now that we have a fasta file with our sequences, we can run cd-hit to remove sequences too similar (>80%)
cd-hit -o Nwestbladi_complete_singlecopy.seqs.rename.cdhit -c 0.8 -i Nwestbladi_complete_singlecopy.seqs.rename.fasta -p 1 -d 0 -b 3 -T 0 -M 10000

# let's rename the default FASTA output to provide an extension
mv Nwestbladi_complete_singlecopy.seqs.rename.cdhit Nwestbladi_complete_singlecopy.seqs.rename.cdhit.fasta


## No sequences were clustered, which means we still have 133 genes to work with.


## Prepare a list with the names of the 133 genes.
grep '>' Nwestbladi_complete_singlecopy.seqs.rename.cdhit.fasta | sed 's/>//g' | sort | uniq > Nwestbladi_complete_singlecopy.seqs.rename.cdhit.names

## Prepare the gff file with the annotation of these genes.
counter=1; cat metazoa_odb10.info.ranked.txt | grep -v "^#" | cut -f 1 | \
while read id; \
do \
grep -w "${id}" Nwestbladi_complete_singlecopy.seqs.rename.cdhit.fasta; \
done | \
tr -d ">" | head -323 | \
while read busco; \
do \
cat ../../05-QC_BUSCO_CleanAssembly/05-BUSCO_odb10/run_Nwestbladi_clean_busco_odb10/augustus_output/predicted_genes/${busco}.out.1 | sed "s/g1/g${counter}/g"; \
let counter+=1; \
done | \
grep -v "^#" | sort -k1,1 -k4,4n > Nwestbladi.133genes.gff


## Convert the GFF file to GTF.
/home/saabalde/Escritorio/software/gffread/gffread -g pt_087_CleanAssembly.nomtDNA.fasta -T -o Nwestbladi.133genes.gtf Nwestbladi.133genes.gff

## gffread failed to find some of the sequences and the output contains only 132 transcripts.


