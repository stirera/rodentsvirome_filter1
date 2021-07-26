#!/usr/bin/python3
# -*-coding:Latin-1 -*

import argparse
import datetime
import os

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

import Contig
import CsvIO
import SubjectProt


def opts_and_args():
    description = 'Postprocessing of blastx csv files'
    usage = """usage : filter_blast [-h] -f fastafile -c input
            csvfile -d shift -l length"""
    parser = argparse.ArgumentParser(description=description,
                                     epilog=usage)
    parser.add_argument('-f',
                        '--fasta',
                        type=str,
                        help='fasta file used as BLAST input',
                        dest='fastafile',
                        metavar='fastafile',
                        required=True,
                        default='')
    parser.add_argument('-c',
                        '--csv',
                        type=str,
                        help='csv file, BLAST output',
                        dest='csvfile',
                        metavar='csvfile',
                        required=True,
                        default='')
    parser.add_argument('-d',
                        '--shift',
                        type=int,
                        help="""scan this lentgh at each end to create
                             stacks (default 10)""",
                        dest='minDec',
                        metavar='minDec',
                        required=False,
                        default='5')
    parser.add_argument('-l',
                        '--length',
                        type=int,
                        help="""length (default 50): i.e minimum length of an alignment
                                and uncovered zone to recycle""",
                        dest='minBpLength',
                        metavar='minBpLength',
                        required=False,
                        default='17')
    parser.add_argument('-I',
                        '--Id',
                        type=int,
                        help='loop number to identify files',
                        dest='recursId',
                        metavar='idrecursif',
                        required=True,
                        default=0)
    args = parser.parse_args()
    return(args)


def filter_algo():
    """
    0 Get args
    1 open files

        1.1 read fasta and csv files :
            1.1.1 if contig is in  csv :
                1.1.1.1 get contigSubjects (filtered):
                1.1.1.2 split contig uncovered zones and recycle

                next
            1.1.1 else :
                    contig -> neg.fasta

            next


    2 End exit
    """
    # Step 0 : Get arguments and create process id along start time
    args = opts_and_args()
    root = args.fastafile.split(".")

    p_id = os.getpid()
    print ("filtering started with proces ID: ", p_id)

    timestart = datetime.datetime.now()

    # Step 1 : opening files with  root string shared between all files
    with open(args.fastafile, 'r') as fasta, \
        open(args.csvfile, 'r') as csv, \
        open(root[0]+"."+str(args.recursId)+".bx_neg.fasta", 'w') \
        as neg_fasta,\
        open(root[0]+"."+str(args.recursId)+".bx.fas", 'w') as uncovered, \
        open(root[0]+"."+str(args.recursId)+".bx.aln", 'w') as alnfile, \
        open(root[0]+"."+str(args.recursId)+".bx.filtered.tsv", 'w') \
            as filteredcsv:

        seqio = SeqIO.parse(fasta, "fasta")
        csvio = CsvIO.CsvIO(csv, args.minDec, args.minBpLength)

        for s in seqio:
            contig = csvio.next(s.id)

            if contig:
                print ("filtering query sequence : ", contig.id)

                contig.filter_subjects()
                for sub in contig.selectedSubjects:
                    x = len(contig.selectedSubjects)
                    print ("found ", x, "subjects")
                    # print (sub.sacc,sub.pident,sub.qstart, sub.qend)

                    subjectdata = contig.id + "\t" + \
                        sub.staxid + "\t" + \
                        sub.sacc + "\t" + \
                        sub.pident + "\t" + \
                        sub.evalue + "\t" + \
                        sub.length + "\t" + \
                        sub.qstart + "\t" + \
                        sub.qend + "\t" + \
                        sub.sstart + "\t" + \
                        sub.send + "\n"

                    filteredcsv.write(subjectdata)

                    h = "aln : query=" + contig.id + \
                        " [ length =" + sub.qlen + "]" + \
                        " [frame =" + sub.qframe + "]" + \
                        " [" + sub.qstart + " " + \
                        sub.qend + "]" + \
                        " : subject=" + sub.sacc + \
                        " [" + sub.sstart + \
                        " " + sub.send + \
                        " " + str(sub.sens) + \
                        "]" + "\t" + sub.sskingdom + \
                        "\n" + sub.qseq + \
                        "\n" + sub.sseq + \
                        "----\n"

                    alnfile.write(h)

                contig.get_uncovered_zones()

                if contig.departs:

                    print ("recycling :")
                    rec_pool = contig.recycleseq(s.seq)
                    # print (rec_pool)
                    for r in rec_pool:
                        uncovered.write(r[0]+"\n")
                        uncovered.write(str(r[1])+"\n")

            else:
                print ("query sequence not in csv, negative: ", s.id)
                SeqIO.write(s, neg_fasta, "fasta")

    print ("fin de process, ID =", p_id)
    timeend = datetime.datetime.now()
    print ("start_time= ", timestart, "\n", "end_time=   ", timeend)


def main():
    filter_algo()

if __name__ == '__main__':
    main()
