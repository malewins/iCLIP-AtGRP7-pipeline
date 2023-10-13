# iCLIP-AtGRP7-pipeline
Step by step workflow of AtGRP7 iCLIP sequencing analysis and binding site determination

## Software and versions
- FastQC (v 0.11.5)
- cutadapt (v 1.16)
- flexbar (v 3.4.0)
- STAR (v 2.6.0a)
- PureCLIP (v 1.0.4)
- bedtools (v2.27.1)
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
flexbar -r adapterless/AtGRP7-GFP_LL36_trimmed.fastq.gz -b Barcodes_AtGRP7-GFP_LL36_rc.fa -n 2 -z GZ -bk -m 24 -o -t demultiplexed/ > demultiplexed/flexbar_AtGRP7-GFP_LL36.log
```

The barcode file has to be provided in FASTA format, where the experimental barcode is defined and random barcode positions are marked with N: 
```
>rep1
NNNAAAANNN
...
```

### quality trimming

### genome mapping
### demultiplexing
### peak calling
### binding site definition
### reproducibility filtering
