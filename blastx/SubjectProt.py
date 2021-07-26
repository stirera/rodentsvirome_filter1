#!/usr/bin/python3
# -*-coding:Latin-1 -*

from Bio.Seq import Seq


class SubjectProt():
    def __init__(self, bl):
        #  format is:
        """qseqid sacc pident length mismatch gaps qstart"""
        #    0     1    2       3        4     5     6
        """qend sstart send evalue bitscore staxid sskingdom"""
        #    7    8      9    10      11      12        13
        """sframe qframe qlen qseq sseq"""
        #    14     15    16   17   18
        if (int(bl[15]) < 0):
            self.sacc, self.pident, self.length, self.mismatch, \
                self.gaps, self.qend, self.qstart, self.sstart, \
                self.send, self.evalue, self.bitscore, self.staxid, \
                self.sskingdom, self.sframe, self.qframe, self.qlen, \
                self.qseq, self.sseq = bl[1:19]
            self.sens = False
        else:
            self.sacc, self.pident, self.length, self.mismatch, \
                self.gaps, self.qstart, self.qend, self.send, \
                self.sstart, self.evalue, self.bitscore, self.staxid, \
                self.sskingdom, self.sframe, self.qframe, self.qlen, \
                self.qseq, self.sseq = bl[1:19]

            self.sens = True
            my_sseq = Seq(bl[16])
            # print "from : ", bl[16], "\nto : "
            rev_sseq = my_sseq.reverse_complement()
