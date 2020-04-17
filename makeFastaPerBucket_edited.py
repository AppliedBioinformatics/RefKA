# We need a mapping from split chromosomes file to soap3 file first
import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Make a fasta file of PacBio reads per bucket.')
parser.add_argument('split_chrom_file', help='Looks like this: "13228 GATTTCATGAAAGAATAGGAACGTGTGATTGATTCCAGACG   NC_035433.1 131293". These are the buckets too!!! Each of these files is one bucket of k-mers we have to assemble.') 
parser.add_argument('soap3_output', help='Looks like this: "194 SRR5637818.4131 9253 + 0". This stores which PacBio read has which k-mer aligning.')
parser.add_argument('idx', help='This is the .idx file created by Biopython')

args = parser.parse_args()

# first, we load the k-mers per bucket
buckets_dict = {} # key: filename, value: list of k-mer IDs contained therein
kmer_to_bucket = {} # key: kmer, value: bucketfile
bucket_fh_dict = {} # key: bucket name, value: file handle to write read to
bucket_seen_dict = {} # key: bucket name, value: set of read names that we have already written to that bucket

print('Loading split_chrom_file with collections of k-mers')
with open(args.split_chrom_file) as kmer_fh:
    for bucketline in kmer_fh:
        # each line is a file containing a collection of k-mers
        this_bucket = bucketline.rstrip()
        buckets_dict[this_bucket] = []
        with open(this_bucket) as bucket_fh:
            for line in bucket_fh:
                ll = line.split()
                # ['322174363', 'TTTTTAATATGGCGTCACGCTGTTATGAATGAATTCATAAC', 'NC_035433.1', '31428']
                kmerid, kmerseq, hit_chrom, hit_pos = ll
                buckets_dict[this_bucket].append(kmerid)
                kmer_to_bucket[kmerid] = this_bucket
                bucket_fh_dict[this_bucket] = open(os.path.basename(this_bucket) + '_READS.fasta', 'w')
                bucket_seen_dict[this_bucket] = set()

print('Finished loading split_chrom_file')

print('Now reading in reads index')
index = SeqIO.index_db(args.idx)
print('Now iterating over the soap3 output')

with open(args.soap3_output) as soap_fh:
    for soapline in soap_fh:
        # each line is a file with soap3 results
        this_soap = soapline.rstrip()
        for line in open(this_soap):
            ll = line.split()
            # ['319903679', 'SRR5637818.49517', '19682', '+', '0']
            kmerid, pacbioread, hitpos, strand, offset = ll
            if kmerid not in kmer_to_bucket:
                # TODO: DELETE ME
                #print('WARNING: Could not find bucket for %s'%kmerid)
                continue
            # now get the bucket this kmer is in
            kmer_bucket = kmer_to_bucket[kmerid]
            # now get the file handle for this bucket
            this_bucket_fh = bucket_fh_dict[kmer_bucket]
            # and get the set of all seen reads for this bucket
            this_bucket_read_set = bucket_seen_dict[kmer_bucket]
            if pacbioread in this_bucket_read_set:
                # skip if we've already written this read for this bucket
                continue
            this_bucket_read_set.add(pacbioread)
            bucket_seen_dict[kmer_bucket] = this_bucket_read_set

            # now get the read we're writing
            read = index[pacbioread]
            this_bucket_fh.write(read.format('fasta'))

