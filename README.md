# dsCAGE: direct cDNA CAGE for paired-end patterned flowcell sequencing. 

A complete analytical framework tailored for the analysis of direct cDNA CAGE data.

Authors: [TO ADD]

Table of content

1.	Introduction
2.	Installation
3.	Input
4.	[Usage]
5.	Main steps
  1. Raw base distribution and quality score distribution calculation	5
  2. Sequence trimming	6
  3. rRNA sequences removal	7
  4. Paired-end reads matching	8
  5. Mapping with STAR	9
  6. BAM to CTSS with G correction	11
  7. Hierarchical Intersection	14
  8. Analysis of Paired BAM files	15
  9. Generating peak file from single-nucleotide CTSS	16

## 1.	Introduction

  This document presents the analysis pipeline that processes data generated from the direct cDNA CAGE library protocol. It provides a comprehensive description of the pipeline and paths of each output file generated throughout the analysis.
  This pipeline supports paired-end sequence data only.
  This pipeline takes demultiplexed paired-end FASTQ files as input, and includes steps to perform: read quality control, read trimming for low quality bases, filtering of short reads, rRNA removal, mapping, conversion to BED file containing CAGE transcription start sites (CTSS) with G correction, annotations of these CTSSs and clustering for peak calling.

  This analysis pipeline handles jobs using Slurm [ref] allowing for simultaneous processing of multiple samples and comes in a Singularity [ref] and Docker containers [ref], available via <https://github.com/hSueki/dscage-pe2>.


## 2.	Installation
  1.	Install any of Singularity [ref] or Docker [ref].
  2.	Download the container image file containing CAGE pipeline.

_With Docker :_

		Pull the docker image from ghcr.io
```
$ docker pull ghcr.io/hsueki/dscage-pe2:latest
```

>[!WARNING] 
>This docker image cannot be used with singularity.
>When you are using singularity, please use sif file as below.


        
_With Singularity :_

			Download the singularity image file from Github and save on your local directory.
  		Download file: dscage-pe2.sif


  3.	Prepare reference

    [In case of Human and Mouse data]
      Download the archived reference from Github and save on your local directory.
      Extract the reference using tar command.
```
　$ tar -zxvf hg38.tar.gz
```
    [In case of other species]
      Prepare the fasta file of genome and rDNA.

## 3.	Input
     The input files required are paired FASTQ files. The accepted naming convention for the FASTQ files is <sample>.R1.fq.gz and <sample>.R2.fq.gz.
     The the input files must be copied within a single input directory: </path/input_fastq_dir/>

## 4.	Usage

    How to Set up and start Docker or Singularity

  With Docker 
  1.	Run a Docker container using the loaded image (hsueki/dscage-pe2). Mount the reference  and your local data directory to the container.
```
  $ docker run -it |
    --mount type=bind,source=/path/to/reference,target=/usr/local/reference |
    --mount type=bind,source=/path/to/your/data,target=/root/data |  
    hsueki/dscage-pe2
```

The path of target directories should not be changed.

The Slurm (a job scheduling and management system) is started automatically within the container. You can verify this by running sinfo.


2.	Change directory into the target directory 
```
(docker)$ cd /root/data
```

3.	Start the pipeline
```
(docker)$ CAGE_PE_pipeline.sh -s 8 -c 4 -t 8
```

4.	Once the pipeline is finished the Docker can be closed.


With Singularity

1.	Shell in the Singularity container image file *.sif containing the entire environment:

$ singularity shell --writable |
--bind /path/to/reference:/usr/local/reference |
dscage-pe2.sif


2.	Set up and start a Slurm cluster

(singularity)$ /etc/start-slurm-services.sh

3.	Change the current directory to the directory containing the input fastq files

(singularity)$ cd </path/input_fastq_dir/>

4. Start the pipeline
   (singularity)$ CAGE_PE_pipeline.sh -s 8 -c 4 -t 8

5.	Once the pipeline is finished, the Docker can be closed, but stop the local-slurm service before exit container.
   (singularity)$ /etc/stop-slurm-services.sh
   (singularity)$ exit



CAGE_PE_pipeline.sh -s <number of samples> -c <concurrent samples to process> -t <number of threads>

Options:
•	-s <number of samples> (mandatory): Specify the number of samples (numeric value).
•	-c <concurrent jobs to be processed> (mandatory): Specify the number of concurrent jobs to be processed (numeric value).
•	-t <number of threads> (optional): Specify the number of threads (numeric value). If not provided, it defaults to 8.


After submitting jobs using this pipeline, you can check their status and cancel them if necessary, using the following commands: squeue, scancel.
## For detailed information and usage, please refer to the Slurm manual available at (https://slurm.schedmd.com/documentation.html)


[In case of other species]
!! When you'd like to use genome other than human and mouse, at first you have to prepare some files. !!

1. genome.fa
2. rDNA.fa
3. bed files for hierarchical intersect

1.	genome.fa
Save the fasta file of genome as /path/to/reference/genome.fa
And you have to make STAR index.

Start docker or singularity container with mount reference directory.
$ docker run -it |
--mount type=bind,source=/path/to/reference,target=/usr/local/reference |
hsueki/dscage-pe2

(docker)$ cd /usr/local/reference
(docker)$ mkdir STAR 
(docker)$ STAR --runMode genomeGenerate \
--genomeDir STAR \
--genomeFastaFiles genome.fa \
--sjdbGTFfile genome.gtf \
--limitGenomeGenerateRAM 3400000000

(docker)$ mv genome.fa STAR

# It takes time to make STAR index...
# Please refer to STAR docs how to prepare STAR index.
# After the STAR index is generaged, please exit the container. 
# Then run the container as written  as above, with mount your data and reference directories..


2.	rDNA.fa
Save the rDNA.fa as /path/to/reference/ribosomalRNA/rDNA.fa


3.	bed files for hierarchical intersect
We prepare annotation bed files using UCSC table browser.
(https://genome.ucsc.edu/cgi-bin/hgTables)

Select these datasets and get output knownGene as BED files.

upstream100.bed
5UTR_exon.bed
coding_exon.bed
3UTR_exon.bed
downstream100.bed
intron.bed

Save these bed files to /path/to/reference/hierarchical_intersect/

If you cannot prepare these annotation bed files, you can run CAGE pipeline and get results without annotation of CAGE tags. 


<reference directory>
The directory and file name under the reference directory should be the same for all species.

 ----- hg38/
   |    |-- reference/
   |          |-- STAR/STAR_INDEX_FILES, genome.fa
   |          |-- ribosomalRNA/rDNA.fa
   |          |-- hierarchical_intersect/BED_FILES
   |
   |-- mm10/
   |    |-- reference/
   |          |-- STAR/STAR_INDEX_FILES, genome.fa
   |          |-- ribosomalRNA/rDNA.fa
   |          |-- hierarchical_intersect/BED_FILES
   |
   |-- other_species/
   |    |-- reference/
   |          |-- STAR/STAR_INDEX_FILES, genome.fa
   |          |-- ribosomalRNA/rDNA.fa
   |          |-- hierarchical_intersect/BED_FILES
   |


5.	Main steps

The pipeline is structured to perform the following steps in sequence:

1. Raw base distribution and quality score distribution calculation
2. Sequence trimming
3. rRNA sequences removal
4. Paired-end reads matching
5. STAR Mapping
6. BAM to CTSS with G correction
7. Hierarchical Intersection and Clustering

1. Raw base distribution and quality score distribution calculation

This step generates QC plots to examine the raw base distribution and raw quality score distribution for read 1 (R1) and read 2 (R2).

Inputs:
<sample>.R1.fq.gz 
<sample>.R2.fq.gz

Tools:
•	baseSeq: A tool for calculating raw base distribution [ref].
•	Graph.sh: A script for generating graphical plots [ref].
•	qvSeq: A tool for calculating raw quality score distribution [ref].

Outputs:
Figures presenting the ｒaw base and raw quality score distribution for read 1 (R1) and read 2 (R2) are generated.
figs/QC/ 
 <sample>.R1.raw_base_distribution.png
 <sample>.R2.raw_base_distribution.png      <sample>.R1.raw_QV_distribution.png
 <sample>.R2.raw_QV_distribution.png　
      

2. Sequence trimming

This process consists in trimming of low-quality base calls from the end of the reads (Phred score > 30), adapter removal (based on a comprehensive list of standard adapters), filtering of reads contains any “Ns” and short sequence filtering (< 15 bp).


Inputs:
<sample>.R1.fq.gz 
<sample>.R2.fq.gz

Tool:
•	trim_galore version 0.6.10 [ref]: A tool for quality-based trimming and adapter removal. TrimGalore is wrapper that applies Cutadapt to trim perform adapter and quality trimming to FASTQ files.

Parameters:
trim_galore main parameters used:
•	--paired: Input data are paired-end reads.
•	--phred33: Phred quality score offset.
•	-q 30: Quality threshold for trimming.
•	--length 15: Minimum length of reads to be kept after trimming.
•	--trim-n: Trims low-quality ends of reads based on quality scores.
•	--max_n 0: Discards reads containing more than a specified number of 'N' bases.

Outputs:
|Directory|File|Description|
|-----|----|----|
|figs/QC/|sample.R1.baseN_distribution.png sample.R2.baseN_distribution.png|distribution of base N (uncertain or ambiguous bases) across the sequences in the input FASTQ files
|tmp/|sample.R1_val_1.fq, sample.R2_val_2.fq|Trimmed and filtered FASTQ files|

 

3. rRNA sequences removal

This portion of the script is responsible for removing sequences that match ribosomal RNA (rRNA) sequences.

Inputs:
•	tmp/<sample>.R1_val_1.fq
•	tmp/<sample>.R2_val_2.fq
•	Sequences of ribosomal RNAs for human or mouse:
•	Human ribosomal DNA complete repeating unit: (Accession U13369; Version U13369.1 [ref])
•	TPA_exp: Mus musculus ribosomal DNA, complete repeating unit (Accession BK000964; Version BK000964.1 [ref])

Tools:
•	RNAdust 1.06 [ref] 
•	SampleTopSeq

Outputs:

|Directory|File|Description|
|----|----|----|
|summary/|sample.R1.top100_rRNA.txt, sample.R2.top100_rRNA.txt|Top100 sequences matched to rRNA.
|summary/|sample.R1.top100_extracted.txt, sample.R2.top100_extracted.txt|Top100 sequences extracted
|tmp/|sample.R1.matchrRNA.fq, sample.R2.matchrRNA.fq|Pair of FASTQ files containing the filtered out rRNA data
|tmp/|sample.R1.norRNA.fq, sample.R2.norRNA.fq|pair of FASTQ files contains the data without rRNA
|log/|sample.R1.matchrRNA.fq.gz, sample.R2.matchrRNA.fq.gz|Pair of FASTQ files containing the filtered out rRNA data

  

4. Paired-end reads matching

This step matches and selects paired sequences by comparing the original FASTQ file pairs to the filtered FASTQ file pairs.

Inputs:
•	<sample>.R1.fq.gz: Input R1 FASTQ file (original).
•	<sample>.R2.fq.gz: Input R2 FASTQ file (original).
•	tmp/ <sample>.R1.norRNA.fq: R1 FASTQ file after rRNA removal.
•	tmp/ <sample>.R2.norRNA.fq: R2 FASTQ file after rRNA removal.


Tool:
•	MatchPairedEndSeq [ref]

Output:
Directory	File	Description
/tmp	<sample>.R1.paired.fq
<sample>.R2.paired.fq	Matched paired FASTQ files


5. Mapping with STAR

Sequences are mapped with STAR (2.7.4a). Mapped data is then sorted with read name, ofiltered to rto retain properly mapped reads (MAPQ > 10).. OReads from the filtered BAM file that are the first and second mates of a pair (R1 and R2 reads) are selected, sorted by coordinates and indexed.

Inputs:
•	Genome: Specifies the genome being used (either “hg38” or “mm10”)
o	hg38: /usr/local/reference/STAR/hg38analysisset
o	mm10: /usr/local/reference/STAR/mm10
•	tmp/<sample>.R1.paired.fq and tmp/<sample>.R2.paired.fq: Paired-end FASTQ files after rRNA removal and read matching. 

Tools:
•	STAR (2.7.4a) [ref]
•	Samtools 1.19 (Using htslib 1.19) [ref]
•	SampleTopSeq

Parameters:
STAR
•	--outFilterMultimapNmax 20: Maximum number of multiple alignments allowed per read.
•	--outSAMtype BAM Unsorted: Format is output is unsorted BAM.
•	--outReadsUnmapped Fastx: unmapped reads outputted in Fastx format.
•	--alignIntronMin 20 and --alignIntronMax 1000000: Minimum and maximum intron lengths for alignment.

Samtools
•	-bSq 10: Reads with a MAPQ score greater than 10 are kept
•	-f 0x2, -F 0x104, -f 0x40, -f 0x80: Filters reads based on flags to select properly mapped reads and exclude non-primary alignments and unmapped reads.

Outputs:
Directory	File	Description
/map	<sample>.bam	Aligned BAM file
/map	<sample>.mapq10.bam	Aligned BAM file filtered with Samtools
/map	<sample>.R1.mapq10.bam
<sample>.R1.mapq10.bam.bai
<sample>.R2.mapq10.bam	Sorted and indexed BAM files for R1 and R2 reads
/map	<sample>.Unmapped.out.mate1
<sample>.Unmapped.out.mate2	Unmapped reads in Fastx format
summary/	<sample>.R1.top100_unmapped.txt
<sample>.R2.top100_unmapped.txt
<sample>.R1.top100_mapped.txt
<sample>.R2.top100_mapped.txt
	Top 100 sequences from mapped and unmapped files
/log	<sample>.Log.final.out	Statistics from STAR: Number of input reads, % Uniquely mapped reads, of reads mapped to multiple loci, % of reads mapped to too many loci, % of reads unmapped: too many mismatches, % of reads unmapped: too short, % of reads unmapped: other
/log	<sample>.Log.out	Description of the execution of the STAR, including command line parameters, initial and final user parameters, effective command line, genome generation parameters and genome information.
/log
	<sample>.Log.progress.out	time-stamped record of progress during a process
/log	<sample>.SJ.out.tab	List of the splice junctions detected during the alignment process. The columns are as follow:   Chromosome: Chromosomal location of the splice junction; Start: Starting position of the splice junction; End: Ending position of the splice junction; Strand: Orientation of the splice junction; Intron Motif: Represents the splice junction's motif, with "0" indicating unknown motif; Support: Number of reads supporting the splice junction; Annotated: Indicates whether the splice junction is annotated or novel.; Unique: Indicates whether the splice junction is unique or shared; Overhang: Length of the overhang.
Of note, since in STAR, paired-end reads are treated as a single read,  the paired-end reads are counted as a single read in the uniquely mapped reads number. The uniquely mapped reads number also includes reads in which only one of the paired-end reads is uniquely mapped.  The file containing the top 100 sequences from mapped and unmapped files was generated.


 				 		
6. BAM to CTSS with G correction
 This step performs a conversion of the BAM file into a CTSS BED file while applying a G correction to account for incorrect G in first base, based on the supplementary note 3e of Nature genet 38:626-35.


Input:
•	BAM file: map/<sample>.R1.mapq10.bam
•	Chromosome name length file as used by STAR: chrNameLength.txt (related to hg38 or mm10)
•	Reference genome FASTA file as used by STAR (hg38analysisset.fa or mm10.fa)

Tools:
•	starBam2GcorrectedCtss_nbg.sh Version 0.1 [ref]

Parameters:
•	-q 10: Only reads with a mapping quality of 10 or higher will be considered for CTSS generation.

Outputs:
Directory	File	Description
map/	<sample>.R1.ctss.bed	1 G-corrected CTSS BED file. The output is formatted as BED (for CTSS), where names represent internal scores and the score represent the corrected counts.
The columns are as follow: 
Chromosome
Start position
End position
name: contains multiple sub-fields separated by commas. Each sub-field provides specific information:
•	X: Observed read counts corresponding to the CTSS.
•	A0: Counts of reads that are observed in this CTSS, with an extra G mismatching to the genome.
•	Nuc: The nucleotide observed at the CTSS position.
•	State: Indicates the state of the CTSS (Start, Generic, End, or Other).
•	A: Counts of reads that are observed in this CTSS, with an extra G.
•	N: Corrected read counts corresponding to the CTSS.
•	U: Counts of reads that are observed in this CTSS, without an extra G.
•	F: The counts of reads that are observed in this CTSS but expected to belong to the next (1bp downstream) CTSS.
Corrected Counts: Represents the corrected read counts corresponding to the CTSS after applying corrections.
Strand
summary/	mapping.txt	mapping summary for all samples which includes the following counts of reads: raw, trimmed_adapters_Nsequences,paired, mapq10, multimap, unmapped, and rRNA.
tmp/	<sample>.mapsummary.txt	Making a Mapping summary file which includes the following counts of reads: raw, trimmed_adapters_Nsequences,paired, mapq10, multimap, unmapped, and rRNA.
summary/	STARsummary_perc.txt	mapping summary for all samples which includes the Statistics from STAR, including input: total number of input reads processed by the aligner; mapq10: % of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; multi: % of reads that mapped to multiple genomic locations; too_many: % of reads that mapped to too many genomic locations, possibly indicating a problem with mapping specificity; unmap_mismatch: % of reads that were unmapped due to mismatches with the reference genome; unmap_short: % of reads that were too short to be mapped.; unmap_other: % of unmapped reads due to reasons other.

tmp/	<sample>.STARsummary_perc.txt	Statistics from STAR, including input: total number of input reads processed by the aligner; mapq10: % of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; multi: % of reads that mapped to multiple genomic locations; too_many: % of reads that mapped to too many genomic locations, possibly indicating a problem with mapping specificity; unmap_mismatch: % of reads that were unmapped due to mismatches with the reference genome; unmap_short: % of reads that were too short to be mapped.; unmap_other: % of unmapped reads due to reasons other.
figs/	STARmap_percent.pdf	Bar plot with stacked bars, where each bar represents the percentage distribution of different STAR-specific mapping statistics for different samples.
tmp/	<sample>.mapsummary.txt	Statistics from STAR, including: the total number of raw sequencing reads; the number of sequences that were trimmed due to the presence of adapters; the number of reads that remained paired after preprocessing and filtering steps; the number of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; the number of reads that mapped to multiple genomic locations; the number of reads that were not successfully mapped to the reference genome; the number of reads that were identified as ribosomal RNA (rRNA) sequences. 

		

7. Hierarchical Intersection
The CTSS files are filtered and annotated using hierarchical_intersect. This tool intersects with each annotation bed file in a hierarchical manner the several regions stored as BED files in /usr/local/reference/hierarchical_intersect : upstream100, 5UTR_exon, coding_exon, 3UTR_exon, intron and downstream100. This step is specific to the reference genome specified: hg38 or mm10.

Input:
•	CTSS file: tmp/<sample>.R1.ctss.bed

Tools:
•	hierarchical_intersect.sh [ref]
•	bedtools intersect (aka intersectBed) v2.31.1

Parameters:
•	-u : (unique) Reporting the mere presence of any overlapping features
•	-s: Enforcing same strandedness

Outputs:
Directory	File	Description
tmp/	<sample>.R1.ctss.bed	CTSS file containing TSSs with Corrected Counts > 0.
tmp/intersect	<sample>.R1.ctss.3UTR_exon.bed
<sample>.R1.ctss.5UTR_exon.bed
<sample>.R1.ctss.coding_exon.bed
<sample>.R1.ctss.downstream100.bed
<sample>.R1.ctss.intergenic.bed
<sample>.R1.ctss.intron.bed
<sample>.R1.ctss.upstream100.bed	Each BED file contains the TSSs that overlap with regions stored as BED files in /usr/local/reference/hierarchical_intersect.
All remaining TSSs not found in any of these annotation BED files, will be stored in <sample>.R1.ctss.intergenic.bed
summary/	promoter_counts.txt	Sum of the read counts for each region: upstream100, 5UTR_exon, coding_exon, 3UTR_exon, intron, downstream100, intergenic and total for each sample of the dataset.
figs/	promoter_counts.pdf	Bar plot with stacked bars, where each bar represents the count of promoters for different types of regions (e.g.,Intergenic, Downstream 100bp, Intron, etc.) across all samples.
figs/	promoter_percent.pdf	Bar plot with stacked bars, where each bar represents the percentage distribution of promoters for different types of regions.


8. Analysis of Paired BAM files

BED12 format provides a more comprehensive representation of genomic features, particularly useful for analyzing spliced alignments. This step converts 'properly paired' BAM alignments to BED12 format and create CAGEscan clusters from the BED12 file.

Caveat:
•	This step does not perform G correction.

Tools:
•	PairedBamToBed12:https://github.com/Population-Transcriptomics/pairedBamToBed12
Input:
o	BAM file: map/<sample>.mapq10.bam
•	
Parameters:
o	-u : (unique) Reporting the mere presence of any overlapping features
o	-s: Enforcing same strandedness

Outputs:
Directory	File	Description
map/	<sample>.bed	BED12 file where the start position is the 5’' end of Read 1 and the end position is the 3' end of read 2. Columns 7-8 (thick coordinates) indicate where the contribution of read 1.
Column 10-12 (BED12 blocks) indicate positions where the reads match.
map/	<sample>.CAGEscan.bed	Column 4 (cluster name): Lx_chr_strand_start_end with the chr,strand,start and end of the TSS single-linkage derived clusters
Column5 (score): number of input paired-end tags of the CAGEscan cluster


9. Generating peak file from single-nucleotide CTSS
This step performs generates a peak file using all single-nucleotide CTSS files of a given project performs peak calling using the Paraclu algorithm, which is commonly employed for clustering and peak detection in genomic data analysis pipelines. Counts for each peak is detailed per sample.


Input:
•	All CTSS files in the same project: map/<sample>.R1.ctss.bed, which are concatenated and formatted into a four-column peak BED file.

Tools:
•	Paraclu and paraclu-cut.sh (https://gitlab.com/mcfrith/paraclu)

Parameters:
•	minValue = 10: omits clusters where the sum of the data values in the cluster is less than minValue
•	d = 2 (default) minimum density increase 
•	l = 200 (default) maximum cluster length

Outputs:
Directory	File	Description
map/	peak.bed	8-column file containing information about the peaks detected by the Paraclu algorithm. The columns are: Chromosome, start position, end position, peak name, the sum of the data values in the cluster, strand, the number of positions with data in the cluster, the cluster's "minimum density", and the cluster's "maximum density".
map/	<sample>.count.bed	7-column file containing the peak details (column 1-6) and the count number for this sample (column 7).
map/	<sample>.txt	Peak name and count number for this sample. The first line contains the total counts
summary/	totals.txt	Number of raw counts of CTSSs, annotated ones and number of peaks per sample.
figs/	ctss_counts.pdf	Bar plot where each bar represents the counts of clustered and raw CTSS for different samples.
figs/	mapping_counts.pdf	Bar plot with stacked bars, where each bar represents the counts of different mapping statistics such as the counts of rRNA, unmatched pairs, removed post-trimming reads, unmapped reads, multimapped reads, and reads with MAPQ10 scores, in millions across different samples.
figs/	mapping_percent.pdf	Bar plot with stacked bars, where each bar represents the percentage distribution of different mapping statistics for different samples.






6.	Brief descriptions of all output file: 
log/
sample.Log.final.out		Log file of STAR
sample.Log.out			Log file of STAR
sample.Log.progress.out          Log file of STAR
sample.R1.matchrRNA.fq.gz	removed rRNA sequence by rRNAdust (R1)
sample.R2.matchrRNA.fq.gz	removed rRNA sequence by rRNAdust (R2)
sample.SJ.out.tab		        Log file of STAR

map/
peak.bed			        peak file from single-nucleotide CTSS using paraclu (from all ctss files)
sample.bed			BED12 read file converted from sample.mapq10.bam (R1+R2 reads)
sample.mapq10.bam		After STAR mapping, sorted and filtered BAM file
sample.R1.count.bed		counts for each peak
sample.R1.count.txt		Output total counts
sample.R1.ctss.bed		R1 G-corrected ctss file (starBam2GcorrectedCtss_nbg.sh)
sample.R1.mapq10.bam	Select R1 reads from sample.mapq10.bam
sample.R1.mapq10.bam.bai	Indexed R1 reads
sample.R2.mapq10.bam	Select R2 reads from sample.mapq10.bam
sample.Unmapped.out.mate1	Unmapped reads, R1
sample.Unmapped.out.mate2	Unmapped reads, R2

summary/mapping.txt			Summary of mapping reads taken from STAR's Log filepromoter_counts.txt		Summary of hierarchical intersectsample.R1.top100_extracted.txt	containing top100 sequences passing rRNA removalsample.R1.top100_mapped.txt		containing top100 sequences from mapped(mapq10)sample.R1.top100_rRNA.txt		containing top100 sequences rRNAsample.R1.top100_unmapped.txt	containing top100 sequences from unmappedSTARsummary_perc.txt			STAR% map file ( from STAR log files )total.txt				summary file for ctss-intersect-peak graph


