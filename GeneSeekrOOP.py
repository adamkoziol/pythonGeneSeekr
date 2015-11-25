#!/usr/bin/env python
from ARMI_Lt import GeneSeekr
from Bio.Blast.Applications import NcbiblastnCommandline
__author__ = 'mike knowles'


class GeneSeeker(GeneSeekr):

    def _blast(self, (fasta, db)):
        self.yeah()
        blastn = NcbiblastnCommandline('/usr/local/bin/blastn',
                                       query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen'",
                                       perc_identity=self.cutoff,
                                       num_descriptions=10000,
                                       num_alignments=10000)
        stdout, stderr = blastn()
        if stdout != '':
            return [[fasta, aln[0].split('_')[0], int(aln[0].split('_')[1])]
                    for aln in [hsp.split('\t')
                                for hsp in stdout.rstrip().split("\n")]
                    if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]


def helper(genes, targets, out, cuttoff, threads):
    from glob import glob
    from json import dump
    import time
    import os
    assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
    assert os.path.isfile(genes), u'rMLST fasta not valid {0!r:s}'.format(genes)
    assert isinstance(threads, int)
    ispath = (lambda x: glob(x + "/*.fa*") if os.path.isdir(x) else [x])
    genes = ispath(genes)
    targets = ispath(targets)
    result = GeneSeekr(genes, targets, threads)
    result.mpblast(cuttoff)
    dump(result.plus, open("%s/MLST-gene_results_%s.json" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w'),
         sort_keys=True, indent=4, separators=(',', ': '))


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Multilocus Seqeunce Typing Assay with BLAST:\n'
                                        'Use to find markers for any bacterial genome')
    parser.add_argument('--version', action='version', version='%(prog)s v0.5')
    parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
    parser.add_argument('-m', '--marker', required=True, help='Specify MLST markers folder')
    parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
    parser.add_argument('-c', '--cutoff', type=int, default=100, help='Threshold for maximum unique bacteria'
                                                                    ' for a single MLST allele')
    parser.add_argument('-t', '--threads', type=int, default=12, help='Specify number of threads')
    args = vars(parser.parse_args())
    helper(args['marker'], args['input'], args['output'], args['cutoff'], args['threads'])