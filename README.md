# Code & data compilation information for paper in review

Sarah K. Hu, Zhenfeng Liu, Harriet Alexander, Victoria Campbell, Paige Connell, Karla B. Heidelberg, Sonya Dyhrman, & David A. Caron. Shifting metabolic priorities among key protistan taxonomic groups within and below the euphotic zone. [Submitted]

Brief description of methods for sequence QC, assembly, annotation, and transcript abundance estimates, followed by review of the R and python scripts used for data compilation, normalization, and code for data visualization. All figures for this paper were made in R.

Contact: sarah.hu[at]usc.edu

## Available data
Raw fastq files available at SRA BioProject: PRJNA391503
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.846380.svg)](https://doi.org/10.5281/zenodo.846380)

* Combined_reference.fa - Custom database (nucleotide)
* Combined_reference_pep.fa - Custom database (amino acid)
* alldepths.fa.gz - contigs
* results.clstr - ortholog cluster output
* contigID_key_depth.txt - list of contigs for each depth
* RawCounts_byContigs.csv - Output count files from Salmon
* RawCounts_bycontig_annotated.csv - Output count file with annotation information
* TaxaName_bycontig.txt - taxonomy annotation by contig
* miTags_alldepths.txt - miTag taxonomy table
* Ortho_groups_by_taxa_binary.csv - (UpSetR input) ortholog groups by taxa
* Ortho_by_depth_Binary.csv - (UpSetR input) ortholog groups by depth
* Annotation_withKO_ids.txt - Annotations for each read
* AnnotationInfo_bycontig.txt - Annotation information for each contig
* Raw_Counts_ByTaxa_CommonKO.csv - Counts of shared/common by taxa (CCA input)
* test_file_05022017.csv - test dataset

## Scripts in this repository
* Parse_and_Compile_Data_SPOTmetaT.ipynb - compiles raw data from annotation and transcript abundance estimation
* MetaT_data_compile_10172017.ipynb .r - R script that continues compilation of raw data, performs normalization, and generates all plots
* miTag_R.ipynb .r - manual taxonomic curation and generation of summary tables and figures from miTag results

## Abbreviated methods (see Hu et al. 2017)
### Sequence QC
Trim raw fastq sequences using Trimmomatic. TruSeq3-PE-2.fa was customized to include string of “AAA’s” and “TTT’s” to removal residual fragments. This is a known artifact from some polyA tail selection library preparation practices (personal communication KAPA Biosystems).
```
java -jar /usr/local/bioinf/Trimmomatic-0.32/trimmomatic-0.32.jar PE [list of fastq files, R1, R2, SE] ILLUMINACLIP:/usr/local/bioinf/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
####Trinity (2.1.1) and in-house script to count & remove ERCC
```
align_and_estimtate_abundance.pl --transcripts ERCCA92 --est_method RSEM --align_method bowtie -- trinity_mode --prep_reference --output_dir
#INSERT in-house script to remove excess ERCC seqs
```
###Separate "bleed-through" rRNA and mRNA using SortMeRNA
```
merge-paired-reads.sh [R1 PE fast#q file] [R1 PE fastq file] [output fastq file for merged reads]
# Example code for paired end reads:
sortmerna --ref /usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s --reads [merged fastq file from above] —sam --fastx --aligned [output file of rRNA reads aligned to databases] --other [output for reads that did not hit database] --log -v --paired_in --best 1

# Example code for single end reads:
sortmerna --ref /usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\/usr/local/bioinf/sortmerna-2.0-linux-64/rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s --reads [SE fastq file] --sam --fastx --aligned [output file of rRNA reads aligned to databases] --other [output for reads that did not hit database] --log -v --paired_in --best 1
```
### miTag analysis

Unmerge rRNA PE reads. Assign taxonomy using QIIME-1's parallel assign taxonomy script that uses uclust.

```
#Unmerge
bash /usr/local/bioinf/sortmerna-2.0-linux-64/unmerge-paired-reads.sh [output file of rRNA reads aligned to databases] [rRNA R1 fastq] [rRNA R2 fastq] 

#Use R1 and SE reads to assign taxonomy. Merged paired end reads is ideal, but sequencing was not designed for merged PE reads in this experiment. I am using split libraries in QIIME so I can get the sequences in a format to assign taxonomy and run a non-stringent quality check.

split_libraries_fastq.py -i [R1 or SE fastq] -m map.txt --barcode_type 'not-barcoded' --sample_ids [name of sample] -q 15 -n 5 -o out_dir
mv out_dir/seqs.fna [output fasta file]

parallel_assign_taxonomy_uclust.py -i [output fasta file] -o [out directory with assignment results] -r /beleriand/SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta -t /beleriand/SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_all_levels.txt --similarity 0.97 --min_consensus_fraction 0.51 -O 4
```

To plot miTag results, see miTag_R.ipynb. R script that takes output from 'parallel_assign_taxonomy.py' to generating  taxonomic composition plots. Since taxonomy assignment included all domains and taxonomic levels, Rscript includes manual curation of taxonomic identities to simplify plotting.


```
### Assembly and ortholog clustering
Performed de novo sequence assembly (mRNA reads) using MEGAHIT. Ortholog clustering was performed at 75% similarity using uclust.
```
#Assembly
megahit -t 16 -1 [list of R1 PE reads] -2 [list of R2 PE reads] -r [list of SE reads]

#Ortholog clustering
uclust --input alldepths.fa --uc results.clstr --id 0.75

#Output contigs and ortholog clustering available (10.5281/zenodo.846380)
```
### Transcript abundance estimates using Salmon
Index assembled reads and perform transcript quantification.
```
#index
salmon index -t [contig file] -i alldepths_index --type quasi -k 31
#transcript abundance estimation
salmon quant -i alldepths_index -l A -1 [R1 mRNA fastq reads] -2 [R2 mRNA fastq reads] -o output.quant
```
## Annotation
### Taxonomy
First ran blastx against custom database. Any reads that did not hit a reference database were re-run with blastn.
```
#Example code:
diamond blastx -d [db] -q [non-rRNA reads] --sensitive -e 0.01 -o [output.blastx]

#Any non-assigned reads were re-run using blastn
blastn -query [non-rRNA reads left unannotated] -evalue 1e-5 -db Combined_reference.fa -task dc-megablast -outfmt 6 -num_threads 16 -out [output.blastn]
```

### Protein annotation
Predicted amino acid sequences using GeneMark, which created .fnn and .faa files from assembled contigs.
```
#Example code:
GeneMarkS-T/gmst.pl --fnn -faa final.contigs.fa
```
Concatenate output .faa sequences and run through GhostKOALA for KEGG annotation (http://www.kegg.jp/ghostkoala/).

## Data compilation
Zenodo (10.5281/zenodo.846380) houses raw data that can be downloaded and then run to re-create all figures.

### Parse_and_Compile_Data_SPOTmetaT.ipynb
Code required to run ahead of R script. Takes raw count information and annotation data and compiles. Output files can then be input into the R script: MetaT_data_compile_plot.ipynb.

### MetaT_data_compile_10172017.ipynb
Script used to generate all figures. Many figures were processed further using InkScape. See "test_file_05022017.csv" to run through small subset of data.

Required R packages
```
install.packages("edgeR")
install.packages("plyr")
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("ggtern")
install.packages("UpSetR")
install.packages("vegan")
```
Imports raw count data and performs EdgeR normalization. Resulting transcript counts per million (CPM) tables are joined with annotation data. 

## Figures
All figures generated in R, then downstream customization performed in InkScape (.svg files).

### Barplot (Fig 2)
CPM was averages over replicates at each sample. Manual taxonomic assignments were added to aid data visualization (see supplementary material).

![Barplot](https://github.com/shu251/figs/blob/master/barplot_SPOTmetaT.png)

### UpSetR plot (Fig 3)
UpsetR package generates easier to digest venn diagrams or "intersecting sets" (http://www.biorxiv.org/content/early/2017/03/25/120600). Results from ortholog group clustering (uclust 75%) were used here to find the distribution of contigs from all depths.

![UpSetR](https://github.com/shu251/figs/blob/master/upsetR_SPOTmetaT.png)

### Ternary plots (Figs 4 and 6)
Using "ggplot2" and "ggtern" R packages (http://www.ggtern.com/), I created the 3 axes "triangle plots". This was performed at both the whole community level and for each individual taxonomic group.

![Ternaryplot](https://github.com/shu251/figs/blob/master/ternary_SPOTmetaT.png)

Statistical analyses are also included at the end of this script. ANOVAs and post-hoc Tukey HSD (to obtain p-values) were performed to determine significance. 

### CCA plot (Fig 5)
Output file from the python jupyter notebook included: "Raw_Counts_ByTaxa_CommonKO.csv". This file includes KEGG identifiers (K0s) that were COMMON across all taxonomic groups. Common K0 IDs are normalized (using edgeR again) and then CCA was computed using the "vegan" R package.

![CCA](https://github.com/shu251/figs/blob/master/CCA_SPOTmetaT.png)

## Contributors
- Jay Liu
- Harriet Alexander

## Citations

Bolger AM, Lohse M, Usadel B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30: 2114-2120.

Buchfink B, Xie C, Huson DH. (2014). Fast and sensitive protein alignment using DIAMOND. Nat Methods 12: 59-60.

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics 10: 421-429.

Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK et al. (2010). QIIME allows analysis of high-throughput community sequencing data. Nat Methods 7: 335-336.

Conway JR, Lex A, Gehlenborg N. (2017). UpSetR: An R Package for the visualization of intersecting sets and their properties. Bioinformatics: 1-3. DOI: 10.1093/bioinformatics/btx1364.

Edgar RC. (2010). Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26: 2460-2461.

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I et al. (2011). Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol 29: 644-652.

Hamilton N. (2016). ggtern: An Extension to 'ggplot2', for the Creation of Ternary Diagrams. R package version 2.1.5. https://CRAN.R-project.org/package=ggtern.

Kanehisa M, Sato Y, Morishima K. (2016). BlastKOALA and GhostKOALA: KEGG Tools for Functional Characterization of Genome and Metagenome Sequences. J Mol Biol 428: 726-731.

Kopylova E, Noe L, Touzet H. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics 28: 3211-3217.

Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 31: 1674-1676.

Li W, Jaroszewski L, Godzik A. (2001). Clustering of highly homologous sequences to reduce the size of large protein databases. Bioinformatics 17: 282-283.

Lund SP, Nettleton D, McCarthy DJ, Smyth GK. (2012). Detecting Differential Expression in RNA-sequence Data Using Quasi-likelihood with Shrunken Dispersion Estimates. Stat Appl Genet Molec Biol 11: 1-44.

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods 14: 1-10.

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P et al. (2012). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res 41: D590-D596.

R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, http://www.R-project.org/ edn: Vienna, Austria.

Tang S, Lomsadze A, Borodovsky M. (2015). Identification of protein coding regions in RNA transcripts. Nucleic Acids Res 43: e78-e78.

Wickham H (2009). ggplot2: Elegant Graphics for Data Analysis. Springer Publishing Company, Incorporated.


## Last Updated 10-23-2017

