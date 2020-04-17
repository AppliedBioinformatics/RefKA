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
