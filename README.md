# dsCAGE: direct cDNA CAGE for paired-end patterned flowcell sequencing. 

A complete analytical framework tailored for the analysis of direct cDNA CAGE data.

**Authors:** Diane Delobel<sup>1</sup>, Hiromi Nishiyori-Sueki<sup>1</sup>, Ilaria Nisoli<sup>2</sup>, Hideya Kawaji<sup>3</sup>,  Pauline Robbe<sup>1</sup> ,Piero Carninci<sup>1,2</sup>, Hazuki Takahashi<sup>1</sup>
<br />

1: RIKEN Center for Integrative Medical Sciences (IMS), Yokohama, 230-0045, Japan <br/>
2: Human Technopole Research Center for Genomics, Milan, 20157, Italy<br/>
3: Tokyo Metropolitan Institute of Medical Science, Research Center for Genome & Medical Sciences, Tokyo, 156-8506, Japan <br/>



### Table of contents

1.	Introduction
2.	Installation
3.	Input
4.	Usage
5.	Main steps
    1. Raw base distribution and quality score distribution calculation
    2. Sequence trimming
    3. rRNA sequences removal
    4. Paired-end reads matching
    5. Mapping with STAR
    6. BAM to CTSS with G correction
    7. Hierarchical Intersection
    8. Analysis of Paired BAM files
    9. Generating a peak file from single-nucleotide CTSS
6.	Brief descriptions of all output files

## 1.	Introduction

  This document presents the analysis pipeline that processes data generated from the direct cDNA CAGE library protocol. It provides a comprehensive description of the pipeline and paths of each output file generated throughout the analysis.
  This pipeline supports paired-end sequence data only.
  This pipeline takes demultiplexed paired-end FASTQ files as input, and includes steps to perform: read quality control, read trimming for low quality bases, filtering of short reads, rRNA removal, mapping, conversion to BED file containing **CAGE transcription start sites (CTSS)** with G correction, annotation of these CTSSs and clustering for peak calling.

  This analysis pipeline handles jobs using [Slurm](https://slurm.schedmd.com/documentation.html) allowing for simultaneous processing of multiple samples and comes in a [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)  and [Docker container](https://docs.docker.com/).


## 2.	Installation
1.	Install either [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) or [Docker](https://docs.docker.com/).
2.	Download the container image file containing the CAGE pipeline.

	**With Docker**: 
	Pull the Docker image from ghcr.io 
	```
	docker pull ghcr.io/hsueki/dscage_pe2:latest
	```

	
	**With Singularity**: 
	Pull the Singularity image file from Sylabs. 
	```
	singularity pull --arch amd64 library://hsueki/dscage/dscage_pe2:latest
	```


>[!WARNING] 
>This docker image cannot be used with Singularity.<br/>
>When you are using Singularity, please use Singularity image file.


<br/>

3. Prepare fasta and bed files
   
	**In the case of Human and Mouse data**<br />
	- Download the archived reference and save on your local directory.<br />
		  Extract the reference using the tar command.　<br />
		  This archive contains annotation bed files and rDNA.fa for Human/Mouse.
	```shell
		tar -zxvf hg38.tar.gz
	```
	- Prepare the FASTA file of [human](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)/[mouse genome](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/) and save it to extracted reference directory as `/path/to/reference/STAR/genome.fa`. <br/>
 	- Prepare the GTF file of the [human](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz)/ [mouse](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz) and save it to extracted reference directory as `/path/to/reference/STAR/genome.gtf`.
 	

	**In the case of other species**<br/>
	- Prepare the FASTA files and BED files:

	   + genome.fa:	Save the FASTA file of genome as `/path/to/reference/STAR/genome.fa`
   	   + genome.gtf: Save the GTF file as `/path/to/reference/STAR/genome.gtf`.
	   + rDNA.fa:	Save the FASTA file of ribosomal RNA as `/path/to/reference/ribosomalRNA/rDNA.fa`
	   + annotation BED files for hierarchical intersection <br/>
		Prepare annotation BED files using [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables).<br />
		Select the following datasets and output the knownGene as BED files.
   
				- upstream100.bed
				- 5UTR_exon.bed
				- coding_exon.bed
				- 3UTR_exon.bed
				- downstream100.bed
				- intron.bed

	        Save these BED files to `/path/to/reference/hierarchical_intersect/`<br/>
>[!NOTE]
>When you cannot prepare these annotation bed files, you can still run the pipeline and obtain simple results without annotating the CAGE tags. 
<br/>

4. Create the STAR index

	Save the fasta file of genome as `/path/to/reference/STAR/genome.fa`. <br/>
   You will also need to create a STAR index.

	Start Docker or Singularity container with mounted reference directory.

	```
	docker run -it --mount type=bind,source=/path/to/reference,target=/usr/local/reference ghcr.io/hsueki/dscage_pe2
	```
	OR
	```
	singularity shell --writable --bind /path/to/reference:/usr/local/reference dscage_pe2_latest.sif
	```

Within the container, run the following commands to generate the STAR index.
```
	cd /usr/local/reference

	STAR --runThreadN 8 \
             --runMode genomeGenerate \
	     --genomeDir STAR \
	     --genomeFastaFiles STAR/genome.fa \
	     --sjdbGTFfile STAR/genome.gtf 
```

>[!NOTE]
> It takes time to create the STAR index.<br/>
> Please refer to STAR documentation on how to prepare the STAR index.<br/>


>[!CAUTION]
> The STAR index should be generated by the same version of STAR used for mapping. <br/>
> It is recommended to use STAR provided in the container.

After the STAR index is generated, please exit the container. 

 Ensure the directory and file names under the reference directory are the same for all species, like this...<br/><br/>
 ![Reference directory structure](https://github.com/user-attachments/assets/76811d84-029f-4f19-ba7a-3a752595308b)



<br/>

## 3.	Input
- The input files required are paired FASTQ files. The accepted naming convention for the FASTQ files is **_sample_.R1.fq.gz** and **_sample_.R2.fq.gz**.
- The input files must be copied within a single input directory: `/path/input_fastq_dir/` <br/>

## 4.	Usage

How to Set up and start Docker or Singularity

**With Docker** 
1. Run a Docker container using the loaded image (hsueki/dscage-pe2). Mount the reference and your local data directory to the container.
```
  docker run -it \
    --mount type=bind,source=/path/to/reference,target=/usr/local/reference \
    --mount type=bind,source=/path/input_fastq_dir,target=/root/data \  
    ghcr.io/hsueki/dscage_pe2
```
>[!IMPORTANT]
>The path of target directories should not be changed.

>[!NOTE]
>The Slurm (a job scheduling and management system) is started automatically within the container. You can verify this by running `sinfo`.
>Please refer to [Slurm manual](https://slurm.schedmd.com/documentation.html).

<br/>

_Within the docker container_ 

2. Change directory into the target directory 
```
 cd /root/data
```

3. Start the pipeline. (About options, please see the last part of this section.) 
```
 CAGE_PE_pipeline.sh -s 8 -c 4 -t 8
```

4. Once the pipeline is finished the Docker can be closed.
<br/>

**With Singularity**

1. Shell in the Singularity container image file `dscage-pe2.sif` containing the entire environment:
```
 singularity shell --writable --bind /path/to/reference:/usr/local/reference \
             dscage_pe2_latest.sif
```

_Within the singularity container_

2. Set up and start a Slurm cluster
```
 /etc/start-slurm-services.sh
```

3. Change the current directory to the directory containing the input fastq files
```
 cd /path/input_fastq_dir/
```

4. Start the pipeline. (About options, please see the last part of this section.) 
```
 CAGE_PE_pipeline.sh -s 8 -c 4 -t 8
```

5. Once the pipeline is finished, the Docker can be closed, but stop the local-slurm service before exit container.
```
/etc/stop-slurm-services.sh

exit
```

<br/>
<br/>

**CAGE_PE_pipeline.sh -s number_of_samples -c concurrent_samples_to_process -t number_of_threads**

Options:
- `-s` number_of_samples (mandatory): Specify the number of samples (numeric value).
- `-c` concurrent_samples_to_process (mandatory): Specify the number of concurrent jobs to be processed (numeric value).
- `-t` number_of_threads (optional): Specify the number of threads (numeric value). If not provided, it defaults to 8.


After submitting jobs using this pipeline, you can check their status and cancel them if necessary, using the following commands: `squeue`, `scancel`.
For detailed information and usage, please refer to the [Slurm manual](https://slurm.schedmd.com/documentation.html).

<br/>

## 5.	Main steps

The pipeline is structured to perform the following steps in sequence:

	   i. Raw base distribution and quality score distribution calculation
	  ii. Sequence trimming
	 iii. rRNA sequences removal
	  iv. Paired-end reads matching
	   v. Mapping with STAR
	  vi. BAM to CTSS with G correction
	 vii. Hierarchical Intersection
	viii. Analysis of Paired BAM files
	  iX. Generating peak file from single-nucleotide CTSS

### i. Raw base distribution and quality score distribution calculation

This step generates QC plots to examine the raw base distribution and raw quality score distribution for read 1 (R1) and read 2 (R2).

Inputs:
- _sample_.R1.fq.gz
- _sample_.R2.fq.gz

Tools:
- `baseSeq`: A tool for calculating raw base distribution[^1].
- `Graph.sh`: A script for generating graphical plots[^1].
- `qvSeq`: A tool for calculating raw quality score distribution[^1].

Outputs:
Figures presenting the raw base and raw quality score distribution for read 1 (R1) and read 2 (R2) are generated.
|Directory|File|
|----|----|
|figs/QC/|_sample_.R1.raw_base_distribution.png<br/>_sample_.R2.raw_base_distribution.png|
|figs/QC/|_sample_.R1.raw_QV_distribution.png<br/>_sample_.R2.raw_QV_distribution.png|　
      
[^1]: These scripts are part of [Moirai](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-144)

### ii. Sequence trimming

This process consists in trimming of low-quality base calls from the end of the reads (Phred score > 30), adapter removal (based on a comprehensive list of standard adapters), filtering of reads contains any “Ns” and short sequence filtering (< 15 bp).


Inputs:
- _sample_.R1.fq.gz
- _sample_.R2.fq.gz

Tool:
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) version 0.6.10 : A tool for quality-based trimming and adapter removal. TrimGalore is wrapper that applies Cutadapt to trim perform adapter and quality trimming to FASTQ files.

Parameters:
trim_galore main parameters used:
- `--paired`: Input data are paired-end reads.
- `--phred33`: Phred quality score offset.
- `-q 30`: Quality threshold for trimming.
- `--length 15`: Minimum length of reads to be kept after trimming.
- `--trim-n`: Trims low-quality ends of reads based on quality scores.
- `--max_n 0`: Discards reads containing more than a specified number of 'N' bases.

Outputs:
|Directory|File|Description|
|-----|----|----|
|figs/QC/|_sample_.R1.baseN_distribution.png<br />_sample_.R2.baseN_distribution.png|distribution of base N (uncertain or ambiguous bases) across the sequences in the input FASTQ files
|tmp/|_sample_.R1_val_1.fq<br />_sample_.R2_val_2.fq|Trimmed and filtered FASTQ files|

 

### iii. rRNA sequences removal

This portion of the script is responsible for removing sequences that match ribosomal RNA (rRNA) sequences.

Inputs:
- tmp/_sample_.R1_val_1.fq
- tmp/_sample_.R2_val_2.fq
- Sequences of ribosomal RNAs (../reference/ribosomalRNA/rDNA.fa):
    - Human ribosomal DNA complete repeating unit: ([Accession U13369; Version U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/555853))
    - TPA_exp: Mus musculus ribosomal DNA, complete repeating unit ([Accession BK000964; Version BK000964.1](https://www.ncbi.nlm.nih.gov/nuccore/BK000964.1))

Tools:
- [rRNAdust](https://fantom.gsc.riken.jp/5/sstar/Protocols:rRNAdust)[^1] 1.06 
- `SampleTopSeq`[^1]

Outputs:

|Directory|File                                                          |Description                                              |
|----     |----                                                          |----                                                     |
|summary/ |_sample_.R1.top100_rRNA.txt<br />_sample_.R2.top100_rRNA.txt          |Top100 sequences matched to rRNA.                        |
|summary/ |_sample_.R1.top100_extracted.txt<br />_sample_.R2.top100_extracted.txt|Top100 sequences extracted                               |
|tmp/     |_sample_.R1.matchrRNA.fq<br />_sample_.R2.matchrRNA.fq                |Pair of FASTQ files containing the filtered out rRNA data|
|tmp/     |_sample_.R1.norRNA.fq<br />_sample_.R2.norRNA.fq                      |pair of FASTQ files contains the data without rRNA       |
|log/     |_sample_.R1.matchrRNA.fq.gz<br />_sample_.R2.matchrRNA.fq.gz          |Pair of FASTQ files containing the filtered out rRNA data|

  

### iv. Paired-end reads matching

This step matches and selects paired sequences by comparing the original FASTQ file pairs to the filtered FASTQ file pairs.

Inputs:
- _sample_.R1.fq.gz: Input R1 FASTQ file (original).
- _sample_.R2.fq.gz: Input R2 FASTQ file (original).
- tmp/ _sample_.R1.norRNA.fq: R1 FASTQ file after rRNA removal.
- tmp/ _sample_.R2.norRNA.fq: R2 FASTQ file after rRNA removal.


Tool:
- MatchPairedEndSeq[^1]

Output:
|Directory|File|Description|
|----|----|----|
|tmp/|_sample_.R1.paired.fq<br />_sample_.R2.paired.fq|Matched paired FASTQ files|


### v. Mapping with STAR

Sequences are mapped with STAR (2.7.4a). Mapped data is then sorted with read name, filtered to retain properly mapped reads (MAPQ > 10). Reads from the filtered BAM file that are the first and second mates of a pair (R1 and R2 reads) are selected, sorted by coordinates and indexed.

Inputs:
- Genome: Specifies the genome being used (genome.fa)
- tmp/_sample_.R1.paired.fq and tmp/_sample_.R2.paired.fq: Paired-end FASTQ files after rRNA removal and read matching. 

Tools:
- [STAR (2.7.4a)](https://github.com/alexdobin/STAR/releases/tag/2.7.4a)
- [Samtools 1.19](https://github.com/samtools/samtools/releases/tag/1.19)
- `SampleTopSeq`[^1]

Parameters:
STAR
- `--outFilterMultimapNmax 20`: Maximum number of multiple alignments allowed per read.
- `--outSAMtype BAM Unsorted`: Format output is unsorted BAM.
- `--outReadsUnmapped Fastx`: unmapped reads output in Fastx format.
- `--alignIntronMin 20 and --alignIntronMax 1000000`: Minimum and maximum intron lengths for alignment.

Samtools
- `-bSq 10`: Reads with a MAPQ score greater than 10 are kept
- `-f 0x2`, `-F 0x104`, `-f 0x40`, `-f 0x80`: Filters reads based on flags to select properly mapped reads and exclude non-primary alignments and unmapped reads.

Outputs:
|Directory|File                                                                                                                   |Description|
|----       |----                                                                                                                   |----                                            |
|map/       |_sample_.bam                                                                                                             |Aligned BAM file
|map/       |_sample_.mapq10.bam                                                                                                      |Aligned BAM file filtered with Samtools
|map/       |_sample_.R1.mapq10.bam<br />_sample_.R1.mapq10.bam.bai<br />_sample_.R2.mapq10.bam                                           |Sorted and indexed BAM files for R1 and R2 reads
|map/       |_sample_.Unmapped.out.mate1<br />_sample_.Unmapped.out.mate2                                                                   |Unmapped reads in Fastx format
|summary/   |_sample_.R1.top100_unmapped.txt<br />_sample_.R2.top100_unmapped.txt<br />_sample_.R1.top100_mapped.txt<br />_sample_.R2.top100_mapped.txt    |Top 100 sequences from mapped and unmapped files
|log/       |_sample_.Log.final.out|Statistics from STAR: Number of input reads<br/> Number and % of : <br/>- Uniquely mapped reads<br/> - reads mapped to multiple loci<br/>- reads mapped to too many loci<br/> - reads unmapped: too many mismatches<br/> - reads unmapped: too short<br> - reads unmapped: other|
|log/      |_sample_.Log.out|Description of the execution of the STAR, including command line parameters, initial and final user parameters, effective command line, genome generation parameters and genome information.
|log/      |_sample_.Log.progress.out|time-stamped record of progress during a process
|log/      |_sample_.SJ.out.tab|List of the splice junctions detected during the alignment process. <br />The columns are as follow:   <br />Chromosome: Chromosomal location of the splice junction; <br />Start: Starting position of the splice junction; <br />End: Ending position of the splice junction; <br />Strand: Orientation of the splice junction; <br />Intron Motif: Represents the splice junction's motif, with "0" indicating unknown motif; <br />Support: Number of reads supporting the splice junction; <br />Annotated: Indicates whether the splice junction is annotated or novel.; <br />Unique: Indicates whether the splice junction is unique or shared; <br />Overhang: Length of the overhang.

Of note, since in STAR, paired-end reads are treated as a single read,  the paired-end reads are counted as a single read in the uniquely mapped reads number. The uniquely mapped reads number also includes reads in which only one of the paired-end reads is uniquely mapped.  The file containing the top 100 sequences from mapped and unmapped files was generated.


 				 		
### vi. BAM to CTSS with G correction
 This step performs a conversion of the BAM file into a CTSS BED file while applying a G correction to account for incorrect G in first base, based on the supplementary note 3e of Nature genet 38:626-35.


Input:
- BAM file: map/_sample_.R1.mapq10.bam
- Chromosome name length file as used by STAR: chrNameLength.txt 
- Reference genome FASTA file as used by STAR: genome.fa

Tools:
- [starBam2GcorrectedCtss_nbg.sh](https://github.com/hkawaji/2016-starBam2GcorrectedCtss) Version 0.1

Parameters:
- `-q 10`: Only reads with a mapping quality of 10 or higher will be considered for CTSS generation.

Outputs:
|Directory|File|Description|
|----|----|----|
|map/|_sample_.R1.ctss.bed| 1G-corrected CTSS BED file. The output is formatted as BED (for CTSS), where names represent internal scores and the score represent the corrected counts. The columns are as follow: <br />col1: Chromosome <br /> col2: Start position <br />col3: End position <br />col4: name: contains multiple sub-fields separated by commas. Each sub-field provides specific information:<br />- **X**: Observed read counts corresponding to the CTSS.<br />- **A0**: Counts of reads that are observed in this CTSS, with an extra G mismatching to the genome.<br />- **Nuc**: The nucleotide observed at the CTSS position.<br />- **State**: Indicates the state of the CTSS (Start, Generic, End, or Other).<br />- **A**: Counts of reads that are observed in this CTSS, with an extra G.<br />- **N**: Corrected read counts corresponding to the CTSS.<br />- **U**: Counts of reads that are observed in this CTSS, without an extra G.<br />- **F**: The counts of reads that are observed in this CTSS but expected to belong to the next (1bp downstream) CTSS.<br />col5: Corrected Counts: Represents the corrected read counts corresponding to the CTSS after applying corrections.<br />col6: Strand|
|summary/|mapping.txt|mapping summary for all samples which includes the following counts of reads: raw, trimmed_adapters_Nsequences,paired, mapq10, multimap, unmapped, and rRNA.|
|tmp/|_sample_.mapsummary.txt|Making a Mapping summary file which includes the following counts of reads: raw, trimmed_adapters_Nsequences,paired, mapq10, multimap, unmapped, and rRNA.|
|summary/|STARsummary_perc.txt|mapping summary for all samples which includes the Statistics from STAR, including input: total number of input reads processed by the aligner; mapq10: % of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; multi: % of reads that mapped to multiple genomic locations; too_many: % of reads that mapped to too many genomic locations, possibly indicating a problem with mapping specificity; unmap_mismatch: % of reads that were unmapped due to mismatches with the reference genome; unmap_short: % of reads that were too short to be mapped.; unmap_other: % of unmapped reads due to reasons other.|
|tmp/|sample.STARsummary_perc.txt|Statistics from STAR, including input: total number of input reads processed by the aligner; mapq10: % of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; multi: % of reads that mapped to multiple genomic locations; too_many: % of reads that mapped to too many genomic locations, possibly indicating a problem with mapping specificity; unmap_mismatch: % of reads that were unmapped due to mismatches with the reference genome; unmap_short: % of reads that were too short to be mapped.; unmap_other: % of unmapped reads due to reasons other.|
|figs/|STARmap_percent.pdf|Bar plot with stacked bars, where each bar represents the percentage distribution of different STAR-specific mapping statistics for different samples.|
|tmp/|sample.mapsummary.txt|Statistics from STAR, including: the total number of raw sequencing reads; the number of sequences that were trimmed due to the presence of adapters; the number of reads that remained paired after preprocessing and filtering steps; the number of reads that were mapped with a mapping quality (MAPQ) score of 10 or higher; the number of reads that mapped to multiple genomic locations; the number of reads that were not successfully mapped to the reference genome; the number of reads that were identified as ribosomal RNA (rRNA) sequences.| 

		

### vii. Hierarchical Intersection
The CTSS files are filtered and annotated using hierarchical_intersect. This tool intersects with each annotation bed file in a hierarchical manner the several regions stored as BED files in /usr/local/reference/hierarchical_intersect : upstream100, 5UTR_exon, coding_exon, 3UTR_exon, intron and downstream100. When the annotation bed files are not found, this step will be skipped.

Input:
- CTSS file: tmp/_sample_.R1.ctss.bed

Tools:
- `hierarchical_intersect.sh`[^1]
- `bedtools intersect` (aka intersectBed) v2.31.1

Parameters:
- `-u` : (unique) Reporting the mere presence of any overlapping features
- `-s`: Enforcing same strandedness

Outputs:
|Directory|File|Description|
|----|----|----|
|tmp/|_sample_.R1.ctss.bed|CTSS file containing TSSs with Corrected Counts > 0.
|tmp/intersect|_sample_.R1.ctss.3UTR_exon.bed<br />_sample_.R1.ctss.5UTR_exon.bed<br />_sample_.R1.ctss.coding_exon.bed<br />_sample_.R1.ctss.downstream100.bed<br />_sample_.R1.ctss.intergenic.bed<br /> _sample_.R1.ctss.intron.bed<br /> _sample_.R1.ctss.upstream100.bed|Each BED file contains the TSSs that overlap with regions stored as BED files in /usr/local/reference/hierarchical_intersect. All remaining TSSs not found in any of these annotation BED files, will be stored in _sample_.R1.ctss.intergenic.bed
|summary/|promoter_counts.txt|Sum of the read counts for each region: upstream100, 5UTR_exon, coding_exon, 3UTR_exon, intron, downstream100, intergenic and total for each sample of the dataset.
|figs/|promoter_counts.pdf|Bar plot with stacked bars, where each bar represents the count of promoters for different types of regions (e.g.,Intergenic, Downstream 100bp, Intron, etc.) across all samples.
|figs/|promoter_percent.pdf|Bar plot with stacked bars, where each bar represents the percentage distribution of promoters for different types of regions.


### viii. Analysis of Paired BAM files

BED12 format provides a more comprehensive representation of genomic features, particularly useful for analyzing spliced alignments. This step converts 'properly paired' BAM alignments to BED12 format and create CAGEscan clusters from the BED12 file.

Caveat:
- This step does not perform G correction.

Tools:
- [PairedBamToBed12](https://github.com/Population-Transcriptomics/pairedBamToBed12)

Input:
- BAM file: map/_sample_.mapq10.bam

Parameters:
- `-u` : (unique) Reporting the mere presence of any overlapping features
- `-s`: Enforcing same strandedness

Outputs:
|Directory|File|Description|
|----|----|----|
|map/|_sample_.bed|BED12 file where the start position is the 5’' end of Read 1 and the end position is the 3' end of read 2. Columns 7-8 (thick coordinates) indicate where the contribution of read 1.Column 10-12 (BED12 blocks) indicate positions where the reads match.|
|map/|_sample_.CAGEscan.bed|Column 4 (cluster name): Lx_chr_strand_start_end with the chr,strand,start and end of the TSS single-linkage derived clusters, Column5 (score): number of input paired-end tags of the CAGEscan cluster
<br/>

### ix. Generating peak file from single-nucleotide CTSS
This step generates a peak file using all single-nucleotide CTSS files of a given project, calling peaks using the Paraclu algorithm, which is commonly employed for clustering and peak detection in genomic data analysis pipelines. Counts for each peak is detailed per sample.


Input:
- All CTSS files in the same project: map/_sample_.R1.ctss.bed, are concatenated and formatted into a four-column peak BED file.

Tools:
- [Paraclu and paraclu-cut.sh](https://gitlab.com/mcfrith/paraclu)

Parameters:<br/>
- `minValue = 10`: omits clusters where the sum of the data values in the cluster is less than minValue
- `d = 2` (default) minimum density increase 
- `l = 200` (default) maximum cluster length

Outputs:
|Directory|File|Description|
|----|----|----|
|map/|peak.bed|8-column file containing information about the peaks detected by the Paraclu algorithm. The columns are: Chromosome, start position, end position, peak name, the sum of the data values in the cluster, strand, the number of positions with data in the cluster, the cluster's "minimum density", and the cluster's "maximum density".|
|map/|_sample_.count.bed|7-column file containing the peak details (column 1-6) and the count number for this sample (column 7).|
|map/|_sample_.txt|Peak name and count number for this sample. The first line contains the total counts.|
|summary/|totals.txt|Number of raw counts of CTSSs, annotated ones and number of peaks per sample.|
|figs/|ctss_counts.pdf|Bar plot where each bar represents the counts of clustered and raw CTSS for different samples.|
|figs/|mapping_counts.pdf|Bar plot with stacked bars, where each bar represents the counts of different mapping statistics such as the counts of rRNA, unmatched pairs, removed post-trimming reads, unmapped reads, multimapped reads, and reads with MAPQ10 scores, in millions across different samples.|
|figs/|mapping_percent.pdf|Bar plot with stacked bars, where each bar represents the percentage distribution of different mapping statistics for different samples.|

<br/>
<br/>




## 6.	Brief descriptions of all output file:
   
|log/| |
|----|----|
|_sample_.Log.final.out|Log file of STAR|
|_sample_.Log.out|Log file of STAR|
|_sample_.Log.progress.out|Log file of STAR|
|_sample_.R1.matchrRNA.fq.gz|removed rRNA sequence by rRNAdust (R1)|
|_sample_.R2.matchrRNA.fq.gz|removed rRNA sequence by rRNAdust (R2)|
|_sample_.SJ.out.tab|Log file of STAR|

|map/|  |
|----|----|
|peak.bed|peak file from single-nucleotide CTSS using paraclu (from all ctss files)|
|_sample_.bed|BED12 read file converted from sample.mapq10.bam (R1+R2 reads)|
|_sample_.mapq10.bam|After STAR mapping, sorted and filtered BAM file|
|_sample_.R1.count.bed|counts for each peak|
|_sample_.R1.count.txt|Output total counts|
|_sample_.R1.ctss.bed|R1 G-corrected ctss file (starBam2GcorrectedCtss_nbg.sh)|
|_sample_.R1.mapq10.bam|Select R1 reads from sample.mapq10.bam|
|_sample_.R1.mapq10.bam.bai|Indexed R1 reads|
|_sample_.R2.mapq10.bam|Select R2 reads from sample.mapq10.bam|
|_sample_.Unmapped.out.mate1|Unmapped reads, R1|
|_sample_.Unmapped.out.mate2|Unmapped reads, R2|

|summary/|  |
|----|----|
|mapping.txt|Summary of mapping reads taken from STAR's Log file|
|promoter_counts.txt|Summary of hierarchical intersect|
|_sample_.R1.top100_extracted.txt|containing top100 sequences passing rRNA removal|
|_sample_.R1.top100_mapped.txt|containing top100 sequences from mapped(mapq10)|
|_sample_.R1.top100_rRNA.txt|containing top100 sequences rRNA|
|_sample_.R1.top100_unmapped.txt|containing top100 sequences from unmapped|
|STARsummary_perc.txt|STAR% map file ( from STAR log files )|
|total.txt|summary file for ctss-intersect-peak graph|


