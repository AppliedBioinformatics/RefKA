# RefKA: A fast and efficient long-read genome assembly approach for large and complex genomes

## Introduction

Recent advances in long-read sequencing have the potential to produce more complete genome assemblies using sequence reads which can span repetitive regions. However, overlap based assembly methods routinely used for this data require significant computing time and resources. We have developed RefKA, a reference-based approach for long read genome assembly. This approach relies on breaking up a closely related reference genome into bins, aligning k-mers unique to each bin with PacBio reads, and then assembling each bin in parallel followed by a final bin-stitching step. 

## Requirements

* Jellyfish
* Soap Aligner
* Soap3: GPU aligner
* Canu

## Inputs

* Long genome sequencing reads. 
* Short Illumina Reads
* Reference Genome


## Process:


**Determine Unique Kmers** 

Use Jellyfish to find unique canonical kmers in the reference genome. 

**Align Unique Kmers to the Reference Genome**

Use Soap to align kmers to the reference genome.

**Split Kmers into Bins**

Use the python script to split the kmers into 500Kbp bins with 100kbp overlaps. 

Input: 
* Chromosome fasta file
* Chromosome Kmer file (Soap output format)

```
python split_soap_file.py -f chromosome.fasta -S chromosome_kmer_file.soap -o output_directory
```

**Align Kmers to Sequencing reads**

Use Soap3 to align kmers to long read sequences. 

**Sort Reads into Bins**

Use Kmer alignments to the reference genome and long sequencing reads to sort reads into bins. 

Input: 
* Split chromosome list: txt file containing complete file paths to all split chromosome files in numerical order. 
* Soap3 alignment list: txt file containing complete file paths to all kmer alignments to reads
* Read Index: this is a fastq SeqIO index of all long sequencing reads, and needs to be in the same location as a fasta file containings all of the long sequencing reads. 

```
python makeFastaPerBucket_edited.py split_chromosome_list.txt soap3_alignment_list.txt read_index.idx

```

**Assemble Bins**

Use Canu to assemble bins separately. 

**Stitch Chromosomes**

Stitch the bins together according to kmer content and overlap.

Input: 
* Split chromosome list: txt file containing complete file paths to all split chromosomes in numerical order
* Bin contig list: txt file containing complete file paths to all assembled bins in numerical order

```
python LiloAndStitch.py split_chromosome_list.txt bin_contig_list.txt output_log.log stitched_output.fasta
```

**Polish with Illumina reads**

Polish stitched genome with ntEdit and Illumina reads. 

## Citations

## Other Resources

## Issues 

## Acknowledgements
