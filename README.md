# iCLIP-AtGRP7-pipeline
Step by step workflow of AtGRP7 iCLIP sequencing analysis and binding site determination

## Software and versions
- FastQC (v 0.11.5)
- cutadapt (v 1.16)
- flexbar (v 3.4.0)
- STAR (v 2.6.0a)
- PureCLIP (v 1.0.4)
- bedtools (v2.27.1)
- bioawk (conda)
- Python (v 3.5) 
- R (v 3.5)
  
## Data
The data was generated in the following publication:


Meyer, K., Köster, T., Nolte, C. et al. Adaptation of iCLIP to plants determines the binding landscape of the clock-regulated RNA-binding protein AtGRP7. Genome Biol 18, 204 (2017). https://doi.org/10.1186/s13059-017-1332-x

Demultiplexed reads are downloadable via the SRA by the accession number SRP108277.

## pipeline workflow
The follwing paragraphs describe the pipeline workflow to determine reproducible binding sites from AtGRP7 iCLIP reads. The quality control, adapter trimming and demultiplexing steps can be skipped if the reads have been downloaded in demultiplexed files. The sequencing files were renamed to the follwing naming scheme for better readability:

```
SRR5628249.fastq.gz -> AtGRP7-GFP-LL36_rep1.fastq.gz
SRR5628250.fastq.gz -> AtGRP7-GFP_LL36_rep2.fastq.gz
SRR5628251.fastq.gz -> AtGRP7-GFP_LL36_rep3.fastq.gz
SRR5628252.fastq.gz -> AtGRP7-GFP_LL36_rep4.fastq.gz
SRR5628253.fastq.gz -> AtGRP7-GFP_LL36_rep5.fastq.gz
```

### quality control
Create a subfolder for the results from FastQC and check the sequencing quality of the reads using FastQC:
```
mkdir -p fastqc
fastqc -oc fastqc/ *.fastq.gz
```


### trimming of the 3' adapter
As displayed in the quality control results from FastQC, the provided sequences are containing 3' sequencing adapter which will be trimmed using cutadapt: 
```
mkdir -p adapterless
ADAPTER="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
cutadapt -j 3 -a $ADAPTER -o adapterless/AtGRP7-GFP_LL36_trimmed.fastq.gz  > adapterless/AtGRP7-GFP_LL36.log
```

### demultiplexing
The sequenced library is provided multiplexed, i.e. the sequencing file contains reads with a barcode adapter sequence at the 5' end to distinguish the individual samples, the reads need to be assigned to the corresponding experiment and saved into separate files:
```
mkdir -p demultiplexed
flexbar -r adapterless/AtGRP7-GFP_LL36_trimmed.fastq.gz -b Barcodes_AtGRP7-GFP_LL36_rc.fa \
  -n 2 \
  -z GZ \
  -bk \
  -m 24 \
  -o \
  -t demultiplexed/ > demultiplexed/flexbar_AtGRP7-GFP_LL36.log
```

The barcode file has to be provided in FASTA format, where the experimental barcode is defined and random barcode positions are marked with N: 
```
>AtGRP7-GFP-LL36_rep1
NNTCCGNNN
>AtGRP7-GFP-LL36_rep2
NNGCCANNN
>AtGRP7-GFP-LL36_rep3
NNAACCNNN
>AtGRP7-GFP-LL36_rep4
NNCGCCNNN
>AtGRP7-GFP-LL36_rep5
NNGCCANNN
```

### quality trimming

After demultiplexing, the reads undergo a quality and length trim with flexbar:

```
for FQ in demultiplexed/*.gz;
do
  OUTPREFIX="${FQ#demultiplexed/}"
  OUTPREFIX="${OUTPREFIX%.fastq.gz}"
  flexbar -r $FQ  --zip-output GZ -t $OUTDIR/$OUTPREFIX  -q WIN -qf sanger  -qt 24 --min-read-length 15 -n 3
done
```

To conserve the random barcode part for later PCR duplicate removal, it is written to the FASTQ *read_id* field using bioawk:

```
OUTDIR=demultiplexed/correct_IDs
mkdir -p $OUTDIR
for FQ in demultiplexed/*.gz;
do
  OUTPREFIX="${FQ#demultiplexed/}"
  bioawk -c fastx '{print "@"$name"#"substr($comment,8) "\n"$seq"\n+\n"$qual}' $FQ | gzip >$OUTDIR/$OUTPREFIX.fastq.gz
done
```

### genomic read mapping
Before mapping reads to the Arabidopsis thalian genome, download the [TAIR10 genome FASTA](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz) and create an index using STAR:

```
REF=TAIR10_chr_all.fas.gz
mkdir -p STARindex
STAR --runMode genomeGenerate \
  --runThreadN 6 \
  --genomeDir STARindex \
  --genomeFastaFiles $REF \
  --genomeSAindexNbases 12 \
  --outFileNamePrefix STARindex_
```

After index generation uniquely map iCLIP reads using STAR and creating a separate folder for each replicate:

```
prefix=demultiplexed/correct_IDs/
outdir=mapped_unique

mkdir -p $outdir

for fastq in demultiplexed/correct_IDs/*.gz; do
 outname="${fastq#$prefix}"
 outname="${outname%.fastq.gz}"
 mkdir -p $outdir/$outname
 STAR --genomeDir STARindex \
  --readFilesIn $fastq \
  --readFilesCommand zcat \
  --runThreadN 12 \
  --alignIntronMin 11 \
  --alignIntronMax 28000 \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 1 \
  --outFileNamePrefix $outdir/$outname/ \
  --outSAMprimaryFlag AllBestScore \
  --outSAMtype BAM SortedByCoordinate \
  --outStd BAM_SortedByCoordinate \
  --alignEndsType Extend5pOfRead1 > $outdir/$outname".bam"
done
```

### PCR duplicate removal
Mapped reads are deduplicated using the custom Python3 script that groups Alignments at the Start position and removes reads with identical random barcode tags:

```
prefix=mapped_unique
outdir=deduplicated

mkdir -p $outdir

for bam in 06_mapped_A11_protein_coding/*.bam; do
 outname="${bam#$prefix}"
 outname="${outname%.bam}"
 mkdir -p $outdir/$outname
 python3.6 python/remove_PCR_duplicate_from_bam.py $bam $outdir/$outname.bam
 echo $SECONDS
done
```

### peak calling
### binding site definition
### reproducibility filtering
