#!/usr/bin/python3
# -*-coding:Latin-1 -*
from pprint import pprint

from Bio.Seq import Seq

import Contig
import SubjectProt


class CsvIO():
    """
    This class takes input csv filehandle and returns contigs that
    matched in blast search.
    It returns "null" in case contig did not match or end of file.
    """

    def __init__(self, handle, minDec, minBpLength):
        self.handle = handle
        self.minDec = minDec
        self.minBpLength = minBpLength

    def next(self, s_id):
        pos = self.handle.tell()
        bl = self.handle.readline().split()

        if not bl:
            return  # end of file

        c = Contig.Contig(bl[0], self.minDec, self.minBpLength)  # init Contig
        # if Current contig (from fasta), does not appear in (output) csv file
        if c.id != s_id:
            self.handle.seek(pos)
            return

        c.subject.append(SubjectProt.SubjectProt(bl))

        # takes subsequent lines with the same query contig id
        while c.id == s_id:
            pos = self.handle.tell()
            bl = self.handle.readline().split()
            if not bl:
                return c   # end of file
            c.id = bl[0]
            c.subject.append(SubjectProt.SubjectProt(bl))

        else:  # end of this contig
            c.id = s_id
            c.subject = c.subject[:-1]
            self.handle.seek(pos)

        return c
