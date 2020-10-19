#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import re

from pyfaidx import Fasta, Sequence
#from scipy.stats import fisher_exact

def locate(args):
    kmer, fd, fo = args.kmer, args.db, args.out
    seqs = Fasta(fd)

    fho = open(fo, 'w')
    #fho.write('kmer\tsid\tsrd\tstart\n')
    i = 1
    for seqid in seqs.keys():
        seq = seqs[seqid][0:].seq
        kseq = kmer
        kseq2 = Sequence(name='kmer',seq=kseq).reverse.complement.seq
        for m in re.finditer(kseq, seq):
            fho.write("%s\t%d\t%s\t+\t%d\t%d\n" % (kseq, i, seqid, m.start()+1, m.start()+len(kmer)))
            i += 1
        for m in re.finditer(kseq2, seq):
            fho.write("%s\t%d\t%s\t-\t%d\t%d\n" % (kseq, i, seqid, m.start()+1, m.start()+len(kmer)))
            i += 1
    fho.close()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "kmer utilities"
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("locate",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "find given kmers in sequence database and report locations")
    sp1.add_argument('kmer', help = 'kmer sequence')
    sp1.add_argument('db', help = 'database (fasta file) to search')
    sp1.add_argument('out', help = 'output file (*.tsv)')
    #sp1.add_argument('--envdir', default='~/git/nf/configs/environments', help = 'config folder')
    sp1.set_defaults(func = locate)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

