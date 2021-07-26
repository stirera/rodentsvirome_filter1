#!/usr/bin/python3
# -*-coding:Latin-1 -*

from Bio.Seq import Seq


class Subject():

    def __init__(self, bl):

        self.sacc, self.pident, self.length, \
            self.mismatch, self.gaps, self.qstart, \
            self.qend, self.sstart, self.send, \
            self.evalue, self.bitscore, self.staxid, \
            self.sskingdom, self.qlen, self.qseq, self.sseq = bl[1:17]
        self.sens = True
        if (int(bl[9]) < int(bl[8])):
            self.sens = False
            my_sseq = Seq(bl[16])
            # print "from : ", bl[16], "\nto : "
            rev_sseq = my_sseq.reverse_complement()
