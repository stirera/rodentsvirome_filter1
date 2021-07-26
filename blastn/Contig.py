#!/usr/bin/python3
# -*-coding:Latin-1 -*

from pprint import pprint
from operator import attrgetter

from Bio.Seq import Seq


class Contig():

    def __init__(self, id, minDec, minBpLength):
        self.id = id
        self.subject = []
        self.minBpLength = minBpLength
        self.minDec = minDec
        self.selectedSubjects = []
        self.departs = []
        self.arrivees = []

    def find_subject_position(self, e, tab):
        """
            Find index before which this element (subject),
            called e, can be inserted
        """
        i = 0
        for t in range(len(tab)):
            # print (e, " compared to ", tab[t])

            if tab[t] > e:
                # print ("t= ", t, e, " compared to ", tab[t])
                i += t
                break
            else:
                i = 0

        return i

    def build_stacks(self, subject, matrix):
        """
        Given a subject and a matrix of subjects :
        Compare coordinates and determine if its alignment
        coordinates are covered by subject(s) in the matrix (+/- minDec)
        fusion is a flag variable set to 0, If covered flag = 1
        Returns fusion, coordinates in the matrix and subject as a table
        """
        fusion = 0
        tabfus = []

        for x in range(len(matrix)):
            for y in range(len(matrix[x])):
                dep = int(subject.qstart)
                arriv = int(subject.qend)
                # cellule => subject
                borders = [int(matrix[x][y].qstart), int(matrix[x][y].qend)]

                if (((dep - self.minDec) > borders[0]) and
                   ((arriv + self.minDec) < borders[1])):
                    # print("checking shift :", self.minDec)
                    # pprint(vars(subject))
                    fusion = 1
                    # print ("valeur de fusion if fusion :", fusion)
                    tabfus = [fusion, x, y, subject]
                    # print ("tabfus ALA2..// ", tabfus)
                    return tabfus
                    break
                m = 0
                while m <= self.minDec:
                    # from 0 to minDec (Shift parameter):
                    # search if border of (qstart/qend) fit in,
                    # in other words searching for overlaps <= minDec value
                    # print ("build_stacks loop number 3 while", y)
                    if ((int(dep-m) == borders[0]) or
                        (int(dep+m) == borders[0]) or
                       (int(arriv-m) == borders[1]) or
                       (int(arriv+m) == borders[1])):
                        fusion = 1
                        # print "valeur de fusion if fusion :", fusion
                        tabfus = [fusion, x, y, subject]
                        # print ("tabfus ALA1..// ", tabfus)
                        return tabfus
                        break
                    else:
                        # print ("fus else : ", fusion)
                        m += 1

        if (fusion == 0):
            tabfus = [fusion, x, y, subject]
            return tabfus

    def rearrange_items(self, tab1, tab2):
        tab1 = [x for x in tab1 if x not in tab2]

    def select_besthit(self, tab):
        """
        Selects the best bitscoring align in a stack
        """
        bitscores = []
        bestbitscore = 0.00
        for s in range(len(tab)):
            bitscores.append(float(tab[s].bitscore))
        bestbitscore = max(bitscores)
        i = bitscores.index(bestbitscore)

        return tab[i]

    def sort_subjects_on_starting(self, subjects):
        """
        Given a contig object with all hsps:
            filter on length (min = minBpLength)
            Order according to alignment starting coordinates
        """

        deps = []
        new_subjects = []
        for s in range(len(subjects)):
            if (int(subjects[s].length) >= int(self.minBpLength)):
                deps.append([int(subjects[s].qstart), subjects[s]])

        sorted_deps = sorted(deps, key=lambda d: int(d[0]))
        for d in sorted_deps:
            new_subjects.append(d[1])

        return new_subjects

    def filter_subjects(self):
        """
            Build a matrix of subjects with stacks, overlapping,
            or single subjects
        """
        new_subjects = self.sort_subjects_on_starting(self.subject)
        # print("new_subjects in filter")
        """
        for n in range(len(new_subjects)) :
                print (n,'*')
                pprint(vars(new_subjects[n]))
            print("####new_subjects in filter")
            classes = []
            smatrix = []
            c=0
            s=0
            new_subjects_seen=[]
            if (len(new_subjects)>1):
                for s in range(len(new_subjects)-1):
                    if (s not in new_subjects_seen):
                        overlappingSubjects = []
                        i=s+1
                        if (abs (int(new_subjects[i].qstart) -
                            int(new_subjects[s].qstart)
                            <=self.minDec)) or \
                            (abs(int(new_subjects[i].qend) -
                            int(new_subjects[s].qend)<=self.minDec)):
                            overlappingSubjects.append(new_subjects[s])
                            overlappingSubjects.append(new_subjects[i])
                            smatrix.append(overlappingSubjects)
                            # s=s-1
                            #s1=s1-1
                            new_subjects_seen.append(i)
                            c+=1
                        else :
                            smatrix.append([new_subjects[s]])
            else :
                smatrix=new_subjects
        """
        classes = []
        idtab = []
        smatrix = []
        for s in range(len(new_subjects)):
            tfus = []
            idtab = [new_subjects[s].sacc,
                     new_subjects[s].qstart,
                     new_subjects[s].qend]
            if idtab not in classes:
                classes.append(idtab)  # seen ?
                # pprint (vars(new_subjects[s]))
                if(len(smatrix) >= 1):
                    t = self.build_stacks(new_subjects[s], smatrix)
                    if t[0] == 1:
                        smatrix[t[1]].append(new_subjects[s])
                        # print ("taille de smatrix = ", len(smatrix))
                    else:
                        smatrix.append([new_subjects[s]])
                else:
                    smatrix.append([new_subjects[s]])
        # Finish with selecting best-hits
        for a in range(len(smatrix)):
            if (len(smatrix[a]) > 0):
                # print("best-hit...")
                # pprint(smatrix[a])
                s = self.select_besthit(smatrix[a])
                self.selectedSubjects.append(s)
            else:
                self.selectedSubjects.append(smatrix[a])
        # print ("Once filterd,...", len(self.selectedSubjects))

    def get_uncovered_zones(self):
        """
        Given final retained alignments:
        Get coordinates of retained alignments and
        cut where there is no alignments (min length applies)
        """

        new_selectedSubjects = \
            self.sort_subjects_on_starting(self.selectedSubjects)
        if len(new_selectedSubjects) == 1:
            s = new_selectedSubjects[0]
            if int(s.qstart) > self.minBpLength:
                splitdep = 1
                splitarr = int(s.qstart) - 1
                self.departs.append(splitdep)
                self.arrivees.append(splitarr)
                # print ("single case, front")
            if int(s.qend) + self.minBpLength < int(s.qlen):
                splitdep = int(s.qend)+1
                splitarr = int(s.qlen)
                self.departs.append(splitdep)
                self.arrivees.append(splitarr)
                # print ("single case, end")
        else:

            for s1 in (range(len(new_selectedSubjects) - 1)):
                s = new_selectedSubjects[s1]

                if s1 == 0 and int(s.qstart) > self.minBpLength:
                    splitdep = 1
                    splitarr = int(s.qstart)-1
                    self.departs.append(splitdep)
                    self.arrivees.append(splitarr)
                    # print("spliting between : ", splitdep, " and ", splitarr)
                    # print ("s1===0, multi")
                    if ((int(new_selectedSubjects[s1].qend) +
                       self.minBpLength) <
                       int(new_selectedSubjects[s1+1].qstart)):
                        splitarr = int(new_selectedSubjects[s1+1].qstart)-1
                        splitdep = int(s.qend)+1
                        # print("spliting between : ", splitdep,
                        #      " and ", splitarr)
                        self.departs.append(splitdep)
                        self.arrivees.append(splitarr)
                else:
                    if ((int(new_selectedSubjects[s1].qend) +
                       self.minBpLength) <
                       int(new_selectedSubjects[s1+1].qstart)):
                        splitarr = int(new_selectedSubjects[s1+1].qstart)-1
                        splitdep = int(s.qend)+1
                        # print("spliting between : ", splitdep,
                        #      " and ", splitarr)
                        self.departs.append(splitdep)
                        self.arrivees.append(splitarr)

            if (len(new_selectedSubjects) > 1):
                # s=new_selectedSubjects[-1]
                s = new_selectedSubjects.pop()
                if int(s.qend)+self.minBpLength < int(s.qlen):
                    splitdep = int(s.qend) + 1
                    splitarr = int(s.qlen)
                    # print("spliting between : ", splitdep, " and ", splitarr)
                    self.departs.append(splitdep)
                    self.arrivees.append(splitarr)
                    # print ("s1===n-1, multi, last")

            print("get_uncovered_zones method completed with :\n d=",
                  self.departs,
                  "\n a=",
                  self.arrivees)

    def recycleseq(self, sequence):
        """
        Given a sequence split, according to coordinates
        """

        pool = []
        for i in range(len(self.departs)):
            depart = int(self.departs[i]) - 1
            arrivee = int(self.arrivees[i]) - 1
            idseq = ">" + self.id + \
                    "_sect:" + str(self.departs[i]) + \
                    "_" + str(self.arrivees[i])
            section = sequence[depart:arrivee]
            pool.append([idseq, section])
        # section=sequence[1:10]

        return pool
