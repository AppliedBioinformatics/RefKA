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

**Kmer counting**

jellyfish count -t $threads -m $kmer_size -C -s 200M -o $outFile $myFasta

**Create unikmer fasta file**


jellyfish dump -U 1 $outFile > $outFile.fasta

Assign an ID to each unikmer

**Option 1: pretty slow**

    n=0; while read line; do if [[ "$line" == ">1" ]]; then n=$((++n)) && echo $line | sed -e 's/>.*/>'$(($n))'/' >> uniq_kmer_ids.fasta; else echo $line >> uniq_kmer_ids.fasta; fi; done <$outFile.fasta

**Option 2: preferred**


#/usr/bin/env python2
```python
__author__= "Andy Yuan"
__contact__="<yuxuan.yuan@outlook.com>"

from Bio import SeqIO
import argparse
import sys
import os

def assignID(fasta, path):
    fasta = os.path.abspath(fasta)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print >> sys.stderr, '\nSomething is wrong with your output directory! Please check!'
        sys.exit(1)
    count=0
    with open (fasta, 'r') as fa:
        for seqs in SeqIO.parse (fa, 'fasta'):
            seq = str(seqs.seq.upper())
            seq_id =seqs.id
            if seq_id == "1":
                count+=1
                fh=open('%s/uniq_kmer_ids.fasta'%path, 'a')
                fh.write('>%s\n'%count)
                fh.write('%s\n'%seq)

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='\nAssign an id to unikmer: start from 1.')
    parser.add_argument('-v', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='unikmer fasta file', type = str)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.output]:
        try:
            assignID(args.fasta, args.output)
        except:
            print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
    else:
        print
        parser.print_help()
        print
        sys.exit(0)
```

**Mapping unikmers back to the reference and get position for each unikmer**
soap2 index the reference fasta

    2bwt-builder $myFasta

**soap2 Mapping**


    soap -D $myFasta.index -a uniq_kmer_ids.fasta -o uniq_kmer_ids.fasta.soap -M 0 -v 0 -r 0 -p $threads -u uniq_kmer_ids.fasta.unmapped

Make an optimised soap file

    sort -k9n uniq_kmer_ids.fasta.soap | awk 'BEGIN{OFS="\t"}''{print $1, $2, $7, $8, $9}' > uniq_kmer_id.fasta.optimised.soap

**Step 2: Splitting the optimised soap file into pieces**

The script used is:

```python
#/usr/bin/env python2

__author__= "Andy Yuan"
__contact__="<yuxuan.yuan@outlook.com>"

from Bio import SeqIO
import argparse
import sys
import os

def splitSOAP(fasta, soap, windowsize, stepsize, path):
    fasta = os.path.abspath(fasta)
    soap = os.path.abspath(soap)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print >> sys.stderr, '\nSomething is wrong with your output directory! Please check!'
        sys.exit(1)
    windowsize=int(windowsize)
    stepsize=int(stepsize)
    with open (fasta, 'r') as fa:
        for seqs in SeqIO.parse (fa, 'fasta'):
            seq = str(seqs.seq.upper())
            seq_id = seqs.id
            piece_counter = 0
            position = 1
            while position < len(seq):
                start, end = position, position+windowsize-1
                position += windowsize - stepsize
                with open (soap, 'r') as Soap:
                    for line in Soap:
                        line=line.strip()
                        pos=int(line.split()[-1])
                        chr=line.split()[3]
                        if chr==seq_id and pos>=start and pos <=end:
                            if not os.path.exists ('%s/%s'%(path,seq_id)):
                                os.makedirs('%s/%s'%(path,seq_id))
                            fh = open('%s/%s/uniq_kmer_id_%s.soap.%s'%(path,seq_id,seq_id,piece_counter), 'a')
                            fh.write('%s\n'%line)
                        else:
                            pass
                piece_counter += 1

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='\nSplit optimised soap file into pieces.')
    parser.add_argument('-v', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='a fasta format file that contains all genome sequences', type = str)
    parser.add_argument('-S', dest='soap', help='optimised soap file', type = str)
    parser.add_argument('-w', dest='windowsize', help='size of the target pieces (bp).Default: 500000', default=500000, type = int)
    parser.add_argument('-s', dest='stepsize', help='stepsize to move the window (bp). Default: 50000', default=50000, type = int)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.soap, args.output]:
        if args.windowsize is None:
            args.windowsize = 500000
        if args.stepsize is None:
            args.stepsize = 50000
        splitSOAP(args.fasta, args.soap, args.windowsize, args.stepsize, args.output)
    else:
        print
        parser.print_help()
        print
        sys.exit(0)
```

This script can create chromosome based directory and chop the input reference into pieces. However, the running time may be long. To save the time, chromosome specified fasta files and soap files are preferred.


**Step 3: Splitting long reads**

The script used to split long reads dataset into chunks is: makeChunks.py

```python
#/usr/bin/env python2

__author__= "Andy Yuan"
__contact__="<yuxuan.yuan@outlook.com>"

# Created date: 01/09/2017
# Last modification date: 03/09/2017

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from math import trunc

def makeChunks(fasta_file,nreads,ns,path):
    """
    This script is used to split a given fasta file into tiny pieces.
    """
    fasta_file = os.path.abspath(fasta_file)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print '\nSomething is wrong with your output directory! Please check!'
        sys.exit()
    nreads=float(nreads)
    ns=float(ns)
    if nreads%ns == 0:
        step=int(nreads/ns)
    else:
        step=int(trunc(nreads/ns)+1)
    i=1
    j=1

    for seqs in SeqIO.parse(fasta_file, 'fasta'):
        ID=seqs.id
        if ID[0]=="@":
            ID=ID[1:]
        seq=str(seqs.seq.upper())
        #start=(i-1)*step+1
        #end=i*step
        #if end >nreads:
        #    end = int(nreads)
        fd=open('%s/chunk_%s.fa'%(path,i),'a')
        if j%step!=0:
            fd.write('>%s\n'%ID)
            fd.write('%s\n'%seq)
        else:
            fd.write('>%s\n'%ID)
            fd.write('%s\n'%seq)
            i+=1
        j+=1

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='\nSplit the raw nucleotide data to chunks.')
    parser.add_argument('-v', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='a fasta format file that contains all genome sequences', type = str)
    parser.add_argument('-t', dest='nreads', help='total number of reads in the fasta file', type = int)
    parser.add_argument('-p', dest='npieces', help='number of chunks to split to ', type=int)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.nreads, args.npieces, args.output]:
        try:
            makeChunks(args.fasta, args.nreads, args.npieces, args.output)
        except:
            print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
    else:
        print
        parser.print_help()
        print
        sys.exit()
```

Accoding to SOAP3-dp, the chunks cannot contain more than 65,000 reads and the total length of all sequences is at most 4 billion. So you may adjust your parameter settings in using ‘makeChunks.py’


```

#========== Modify the settings below ========
myFasta="/path/to/your/fasta"
outDir="/path/to/your/outputDir"
nChunks="" # num of total chunks
#========= make chunks ===========
nReads=`grep -c ">" "$myFasta"`
python makeChunks.py -f "$myFasta" -t $nReads -p $nChunks -o "$outDIR"
```

**Step 4: SOAP3-dp indexing read chunks**

When using SOAP3-dp, please ensure “.ini” files in the same directory as your software. An example is



    /path/to/SOAP3-dp/bin/soap3-dp-builder /path/to/chunk_1.fa

**Step 5: SOAP3-dp converting normal index to GPU index**

The command line used is:


    /path/to/SOAP3-dp/bin/BGS-Build /path/to/chunk_1.fa.index

**Step 6: SOAP3-dp mapping**

The command lines are:
```
#========== Modify the settings below ========
nChunks="" # num of total chunks
nPieces="" # num of pieces from the reference
index="/path/to/soap3_index_Dir"
outDir="/path/to/outputDir"
SOAP3-dp="/path/to/SOAP3-dp_Dir"
unikmers="/path/to/unikmer_fasta" # see below the format
#========= start mapping ===========
for i in `seq 1 $nChunks`; do mkdir -p "$outDir/chunk.$i"; cd "$outDir/chunk.$i"; "$SOAP3-dp/soap3-dp" single "$index/chunk_${i}.fa.index" "$unikmers" -o uniq_kmer.soap.chunk.$i -L 41 -h 1 -s 0; cat uniq_kmer.soap.chunk.$i.gout.* | grep -v "@" | awk '{if($3!="*") print $1, $3}' >uniq_kmer.soap.chunk.$i.ids; rm -rf uniq_kmer.soap.chunk.$i.gout.*; rm -rf uniq_kmer.soap.chunk.$i.dpout.*; rm -rf uniq_kmer.soap.chunk.$i.done; done
for i in `seq 1 $nChunks`; do cd "$outDir/chunk.$i"; for j in `seq 0 $nPieces`; do grep soap.${j}_kmerID uniq_kmer.soap.chunk.$i.ids | cut -d" " -f2 | sort | uniq >soap.$j.chunk${i}.reads; done; done
```
The format of the unikmer_fasta file is:

```
>soap.0_kmerID_2396792851_pos_761
AAACTATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACA
>soap.0_kmerID_3338392134_pos_762
AACTATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAA
>soap.0_kmerID_4867289698_pos_763
ACTATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAA
>soap.0_kmerID_2852094694_pos_764
CTATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAAA
>soap.0_kmerID_2664500755_pos_765
TATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAAAC
>soap.0_kmerID_4135776802_pos_766
ATTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAAACG
>soap.0_kmerID_1626108915_pos_767
TTGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAAACGT
>soap.0_kmerID_2534315346_pos_768
TGTTTTTCATCCTGTAGTCCCATTTAGAATTACAAAACGTC
```
How to create such a file

If you have a soap file like (lets’s say uniq_kmer_id_chr1A.soap.0):

```
3072848642      CTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAAC       +       chr1A   1
4903828995      TAAACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAACC       -       chr1A   2
4938805162      AAACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAACCC       +       chr1A   3
3182373650      AACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAACCCT       +       chr1A   4
2291210296      ACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAACCCTA       +       chr1A   5
3413536015      CCCTAAACCCTAAACCCTAAACCCTAAACCCCTAACCCTAA       +       chr1A   6
2603724574      CCTAAACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAA       +       chr1A   7
2933780398      CTAAACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAAC       +       chr1A   8
3046773362      TAAACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAACC       -       chr1A   9
4335988878      AAACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAACCC       +       chr1A   10
2503614940      AACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAACCCT       +       chr1A   11
254313026       ACCCTAAACCCTAAACCCTAAACCCCTAACCCTAAACCCTA       +       chr1A   12
1662286652      CCCTAAACCCTAAACCCTAAACCCCTAACCCTAAACCCTAA       +       chr1A   13
```
You may change the file into a fasta using:

```bash
i=0
cat uniq_kmer_id_chr1A.soap.0 | awk -v id=$i '{print ">soap."id"_kmerID_"$1"_pos_"$5"\n"$2}' > uniq_kmer_id_chr1A.soap.$id.fasta
```
So in the whole chromosome, you may use:

```bash
#========== Modify the settings below ========
soapFileDir="/path/to/soap_file_Dir"
chr="" #name of chromosome, such as chr1A
outDir="/path/to/outputDir"
#====== start converting =======
cd $soapFileDir
for i in `ls -l | grep -v fasta | awk '{print $9}' | grep soap | rev | cut -d"." -f1 | rev | sort -n`; do echo $i; cat uniq_kmer_id_${chr}.soap.$i | awk -v id=$i '{print ">soap."id"_kmerID_"$1"_pos_"$5"\n"$2}' >> "$outDir/uniq_kmer_id_${chr}.soap.fasta"; done
```
The whole chromsome based fasta file is the one used in SOAP3-dp mapping
**Step 7: Getting long reads id for each piece**

The command lines are:

```bash
#========== Modify the settings below ========
nChunks="" # num of total chunks
nPieces="" # num of pieces from the reference
outDir="/path/to/outputDir"
#====== getting ids =======
for i in `seq 0 $nPieces`; do for j in `seq 1 $nChunks`; do cd "$outDir/chunk.$j"; cat soap.$i.chunk$j.reads | awk -v id=$i 'BEGIN{OFS="\t"}{print $1, id}'>> "$outDir/all.txt"; done; done
```
**Step 8: Pulling out long reads to each piece**

A script used to pull out all long reads to each piece is:

```python

#!/usr/bin/env python2

__author__= "Andy Yuan"
__contact__="<yuxuan.yuan@outlook.com>"

import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

def extract(fa, npieces, info, path):
    """
    fa: the whole fasta dataset
    npieces: total number of pieces from the reference
    info: the read id and corresponsing soap id, such as all.txt file.
    """
    fa = os.path.abspath(fa)
    info = os.path.abspath(info)
    path = os.path.abspath(path)
    npieces=int(npieces)
    faDict={}
    infoDict=defaultdict(list)
    with open (fa, 'r') as fasta:
        for seqs in SeqIO.parse(fasta, 'fasta'):
            seqID=seqs.id
            seq=str(seqs.seq)
            faDict[seqID]=seq
    with open (info, 'r') as info:
        for line in info:
            line=line.strip()
            try:
                ID=line.split()[0]
                soap=line.split()[1]
                infoDict[int(soap)].append(ID)
            except IndexError:
                pass
    for i in range(npieces):
        with open ('soap.%s.fasta'%i, 'a') as myResult:
            for j in infoDict[i]:
                myResult.write('>%s\n'%j)
                myResult.write('%s\n'%faDict[j])

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='\nExtract sequences using the id from the whole dataset')
    parser.add_argument('-v', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='the whole dataset', type = str)
    parser.add_argument('-n', dest='npieces', help='total number of pieces from specific chromosome', type = int)
    parser.add_argument('-i', dest='info', help='the read id and corresponsing soap id from one chromosome, such as all.txt file', type = str)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.npieces, args.info, args.output]:
        try:
             extract(args.fasta, args.npieces, args.infpo, args.output)
        except:
            print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
    else:
        print
        parser.print_help()
        print
        sys.exit(0)
```

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
