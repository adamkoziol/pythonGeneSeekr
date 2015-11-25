#!/usr/bin/env python
import os
import sys
import time
import json
import copy_reg
import types
from collections import defaultdict
from Bio.Application import _Option, AbstractCommandline, _Switch
from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing import Pool
from ARMICARD import decipher
__author__ = 'mike knowles'

__doc__ = 'The purpose of this set of modules is to improve upon earlier development of ARMISeekr.py and eventually' \
          'to include generalized functionality for with OOP for GeneSeekr'


class MakeBlastDB(AbstractCommandline):
    """Base makeblastdb wrapper"""
    def __init__(self, cmd='makeblastdb', **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file prefix for db.",
                    filename=True,
                    equate=False),
            _Option(["-in", "db"],
                    "The sequence create db with.",
                    filename=True,
                    equate=False),  # Should this be required?
            _Option(["-dbtype", "dbtype"],
                    "Molecule type of target db (string, 'nucl' or 'prot').",
                    equate=False)]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):  # deal with mangled names
        cls_name = cls.__name__.lstrip('_')
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)




copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class GeneSeekr(object):

    def yeah(self, reset=None):
        """
        :type reset: int
        :rtype: yeah
        """
        if reset is not None:
            self.count = 1
        if self.count == 1:
            sys.stdout.write('[{}] 1 ( \xE2\x80\xA2_\xE2\x80\xA2)'.format(time.strftime("%H:%M:%S")))
        elif self. count % 3 == 0:
            sys.stdout.write('\r[{}] {} (\xE2\x8C\x90\xE2\x96\xA0_\xE2\x96\xA0) #Yeeeaaahhhh'
                             .format(time.strftime("%H:%M:%S"), self.count))
        elif self.count % 2 == 0:
            sys.stdout.write('\r[{}] {} ( \xE2\x80\xA2_\xE2\x80\xA2)>\xE2\x8C\x90\xE2\x96\xA0-\xE2\x96\xA0'
                             .format(time.strftime("%H:%M:%S"), self.count))
        else:
            sys.stdout.write('\r[{}] {} ( \xE2\x80\xA2_\xE2\x80\xA2)'.format(time.strftime("%H:%M:%S"), self.count))
        self.count += 1

    def makeblastdb(self, (fasta, db)):
        if not os.path.isfile('{}.nhr'.format(db)):  # add nhr for searching
            assert os.path.isfile(fasta) # check that the fasta has been specified properly
            MakeBlastDB('/usr/local/bin/makeblastdb', db=fasta, out=db, dbtype='nucl')()  # Use MakeBlastDB above
            self.yeah()
        return 0

    def __init__(self, subject, query, threads=12):
        """:type subject: list of genes
           :type query: list of target genomes"""
        assert isinstance(subject, list), 'Subject is not a list "{0!r:s}"'.format(subject)
        assert isinstance(query, list), 'Query is not a list"{0!r:s}"'.format(query)
        self.count, self.subject, self.query, self.threads = 0, subject, query, threads
        self.cutoff, self.genelist = 70, []
        self.db = map((lambda x: os.path.splitext(x)[0]), subject)  # remove the file extension for easier globing
        self.plus = dict((target, defaultdict(list)) for target in self.query)  # Initialize :return dict
        print '[{}] GeneSeekr input is path with {} files'.format(time.strftime("%H:%M:%S"), len(query))
        print "[{}] Creating necessary databases for BLAST".format(time.strftime("%H:%M:%S"))
        Pool(self.threads).map(self.makeblastdb, zip(self.subject, self.db))
        print "\r[{0}] BLAST database(s) created".format(time.strftime("%H:%M:%S"))

    def _blast(self, (fasta, db)):
        self.yeah()
        blastn = NcbiblastnCommandline('/usr/local/bin/blastn',
                                       query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen'",
                                       perc_identity=self.cutoff)
        stdout, stderr = blastn()
        if stdout != '':
            return [[fasta, aln[0][4:], abs(float(aln[1]) / float(aln[2]))]
                    for aln in [hsp.split('\t')
                                for hsp in stdout.rstrip().split("\n")]
                    if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]

    def mpblast(self, cutoff=70):
        assert isinstance(cutoff, int), u'Cutoff is not an integer {0!r:s}'.format(cutoff)
        self.cutoff = cutoff
        print "[{}] Now performing and parsing BLAST database searches".format(time.strftime("%H:%M:%S"))
        # self.yeah(0)
        start = time.time()
        p = Pool(12)
        for genes in self.db:
            mapblast = p.map(self._blast, [(genome, genes) for genome in self.query])
            for fastaline in mapblast:
                if fastaline is not None:
                    for fasta, k , v in fastaline:
                        if k not in self.genelist:
                            self.genelist.append(k)
                        self.plus[fasta][k].append(v)

        end = time.time() - start
        print "\n[{0:s}] Elapsed time for GeneSeekr is {1:0.2f} mins with {2:0.2f}s per genome".format(
            time.strftime("%H:%M:%S"), end / 60, end / float(len(self.query)))



    def csvwriter(self, out, name):
        assert isinstance(out, str), u'Output location is not a string {0!r:s}'.format(out)
        assert isinstance(name, str), u'Output name is not a string {0!r:s}'.format(name)
        assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
        self.genelist.sort()
        row, rowcount, csvheader = '', 0, 'Strain'
        for genomerow in self.plus:
            row += "\n" + genomerow.split('/')[-1].split('.')[0]
            rowcount += 1

            for generow in self.genelist:
                genename = generow
                if rowcount <= 1:
                    csvheader += ', ' + genename
                if genename in self.plus[genomerow]:
                    allelenum = ""
                    self.plus[genomerow][genename].sort()
                    for allele in self.plus[genomerow][genename]:
                        allelenum += str(allele) + ' '
                    row += ',' + allelenum
                else:
                    row += ',N'
        with open("%s/%s_results_%s.csv" % (out, name, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
            csvfile.write(csvheader)
            csvfile.write(row)


def helper(genes, targets, out, cuttoff, aro, threads):
    from glob import glob
    assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
    assert os.path.isfile(genes), u'ARMI-genes.fa not valid {0!r:s}'.format(genes)
    assert os.path.isfile(aro), u'Antibiotic JSON not valid {0!r:s}'.format(aro)
    assert isinstance(threads, int)
    ispath =(lambda x: glob(x + "/*.fa*") if os.path.isdir(x) else [x])
    genes = ispath(genes)
    targets = ispath(targets)
    result = GeneSeekr(genes, targets, threads)
    result.mpblast(cuttoff)
    json.dump(result.plus, open("%s/ARMI-gene_results_%s.json" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w'),
              sort_keys=True, indent=4, separators=(',', ': '))
    decipher(result.plus, json.load(open(aro)), out)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                                 'Use to find markers for any bacterial genome')
    parser.add_argument('--version', action='version', version='%(prog)s v0.5')
    parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
    parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
    parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
    parser.add_argument('-a', '--anti', type=str, required=True, help='JSON file location')
    parser.add_argument('-c', '--cutoff', type=int, default=70, help='Threshold for maximum unique bacteria'
                                                                    ' for a single antibiotic')
    parser.add_argument('-t', '--threads', type=int, default=12, help='Specify number of threads')
    args = vars(parser.parse_args())
    helper(args['marker'], args['input'], args['output'], args['cutoff'], args['anti'], args['threads'])
