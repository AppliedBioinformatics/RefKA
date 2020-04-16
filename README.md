# RefKA: A fast and efficient long-read genome assembly approach for large and complex genomes

## Introduction

Recent advances in long-read sequencing have the potential to produce more complete genome assemblies using sequence reads which can span repetitive regions. However, overlap based assembly methods routinely used for this data require significant computing time and resources. We have developed RefKA, a reference-based approach for long read genome assembly. This approach relies on breaking up a closely related reference genome into bins, aligning k-mers unique to each bin with PacBio reads, and then assembling each bin in parallel followed by a final bin-stitching step. 

## Requirements

### Jellyfish

###  Soap Alignner

### Soap3: GPU aligner


## Inputs

Long genome sequencing reads. 
Short Illumina Reads
Reference Genome


## Process

### Determine Unique Kmers 

### Align Unique Kmers to the Reference Genome

### Split Kmers into Bins

### Align Kmers to Sequencing reads

###Sort Reads into Bins

### Assemble Bins

### Stitch Chromosomes

### Polish with Illumina reads

## Citations

##Other Resources

## Issues 

## Acknowledgements
