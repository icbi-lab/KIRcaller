#!/usr/bin/env python
#
#  callKIR3D.py
#
#  Copyright 2018 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from __future__ import division
import sys
import os
import getopt
import uuid
import tempfile
from subprocess import Popen, PIPE
from math import log
import pandas as pd

myPath = os.path.dirname(__file__)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:o:f:m:l:v", [
            "help", "bamDIR=", "output=", "featureCounts=", "maxReadsL0=",
            "minLogR="
        ])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    bamDIR = ""
    output = "./kir3D_genotypes.tsv"
    fcAlt = ""

    global verbose, maxReadsL0, minLogR
    verbose = False
    maxReadsL0 = 15
    minLogR = 2

    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-b", "--bamDIR"):
            bamDIR = os.path.realpath(a)
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-f", "--featureCounts"):
            fcAlt = a
        elif o in ("-m", "--maxReadsL0"):
            maxReadsL0 = a
        elif o in ("-l", "--minLogR"):
            minLogR = a
        else:
            assert False, "unhandled option"

    countRes = countKIR(bamDIR, output, fcAlt)
    zygRes = analyzeKIR(bamDIR, countRes)
    zygRes.to_csv(output, sep="\t", index=False)


def usage():
    print(
        sys.argv[0] +
        " --bamDIR <directory with BAM files> --output <KIR genotyping output> [--featureCounts <path to featureCounts>]"
    )


def mean(l):
    return sum(l) / float(len(l))


def find_nearest(array, value):
    n = [abs(i - value) for i in array]
    idx = n.index(min(n))
    return array[idx]


def most_common(lst):
    return max(set(lst), key=lst.count)


def checkBAMdir(bamDIR):

    if ((bamDIR == "") or (os.path.exists(bamDIR) == False)):
        print("Error: BAM dir: " + bamDIR + " not found!")
        usage()
        sys.exit(2)

    bamF = []
    for file in os.listdir(bamDIR + "/"):
        if file.endswith(".bam"):
            bamF.append(os.path.realpath(bamDIR + "/" + file))

    if (len(bamF) < 1):
        print("Error: no BAM files found in: " + bamDIR)
        sys.exit(2)
    return (bamF)


def countKIR(bamDIR, output, fcAlt):
    bamF = checkBAMdir(bamDIR)
    print("Analyzing " + str(len(bamF)) + " BAM files in: " + bamDIR +
          "\nSaving results in: " + output)

    fc = "featureCounts" if fcAlt == "" else fcAlt

    fcRes = tempfile.gettempdir() + "/KIR3Dfc_" + str(uuid.uuid4()) + ".txt"
    cmdLine = [
        fc, '-p', '-t', 'exon', '-g', 'Name', '-P', '-C', '--primary',
        '--fracOverlap', '0.90', '--maxMOp', '1', '-a',
        myPath + '/data/KIR3DS1_L1.gff', '-o', fcRes
    ]
    cmdLine.extend(bamF)

    if (verbose == True):
        print("Running featureCounts")
        print(" ".join(cmdLine))

    p = Popen(cmdLine, stdout=PIPE, stderr=PIPE)
    fcOut, fcErr = p.communicate()

    if (p.returncode != 0):
        print(fcErr)
        print("An error occured, exiting...")
        sys.exit(2)
    else:
        return (fcRes)


def analyzeKIR(bamDIR, countRes):
    countData = pd.read_table(
        countRes, sep="\t", index_col=0, skip_blank_lines=True, header=1)

    os.unlink(countRes)
    os.unlink(countRes + ".summary")

    countData.rename(
        columns=lambda x: x.replace(bamDIR + "/", "").replace(".bam", ""),
        inplace=True)
    subjects = countData.columns.tolist()[5:]

    zygData = []

    for subj in subjects:
        if (verbose == True):
            print(subj)

        l0 = countData.loc["KIR3DL1_uniq.0", subj]

        l1Orig = countData.loc["KIR3DL1_uniq.1", subj]
        l2Orig = countData.loc["KIR3DL1_uniq.2", subj]

        l1 = 1 if (countData.loc["KIR3DL1_uniq.1", subj] == 0
                  ) else countData.loc["KIR3DL1_uniq.1", subj]
        l2 = 1 if (countData.loc["KIR3DL1_uniq.2", subj] == 0
                  ) else countData.loc["KIR3DL1_uniq.2", subj]

        s1Orig = countData.loc["KIR3DS1_uniq.1", subj]
        s2Orig = countData.loc["KIR3DS1_uniq.2", subj]

        s1 = 1 if (countData.loc["KIR3DS1_uniq.1", subj] == 0
                  ) else countData.loc["KIR3DS1_uniq.1", subj]
        s2 = 1 if (countData.loc["KIR3DS1_uniq.2", subj] == 0
                  ) else countData.loc["KIR3DS1_uniq.2", subj]

        logR1 = round(log(l1 / s1, 2), 3)
        logR2 = round(log(l2 / s2, 2), 3)

        BAF1 = round((s1 / (l1 + s1)), 3)
        BAF2 = round((s2 / (l2 + s2)), 3)

        ZYG1 = "HOM_L" if (logR1 > minLogR) else ("HOM_S"
                                                  if ((logR1 < -minLogR) and
                                                      (l0 < maxReadsL0)) else
                                                  "HET")
        ZYG2 = "HOM_L" if (logR2 > minLogR) else ("HOM_S"
                                                  if ((logR2 < -minLogR) and
                                                      (l0 < maxReadsL0)) else
                                                  "HET")

        bafVal = find_nearest([0, 0.5, 1], BAF1)
        ZYG3 = "HOM_L" if (bafVal == 0) else ("HOM_S"
                                              if ((bafVal == 1) and
                                                  (l0 < maxReadsL0)) else "HET")

        bafVal = find_nearest([0, 0.5, 1], BAF2)
        ZYG4 = "HOM_L" if (bafVal == 0) else ("HOM_S"
                                              if ((bafVal == 1) and
                                                  (l0 < maxReadsL0)) else "HET")

        ZYGM = "HOM_L" if (mean([logR1, logR2]) >
                           1) else ("HOM_S" if ((mean([logR1, logR2]) < -1) and
                                                (l0 < maxReadsL0)) else "HET")
        ZYGC = most_common([ZYG1, ZYG2, ZYG3, ZYG4, ZYGM])

        zygData.append((subj, l0, l1Orig, s1Orig, l2Orig, s2Orig, logR1, logR2,
                        BAF1, BAF2, ZYG1, ZYG2, ZYG3, ZYG4, ZYGC))

    labels = [
        "Subject", "KIR3DL1_u0", "KIR3DL1_u1", "KIR3DS1_u1", "KIR3DL1_u2",
        "KIR3DS1_u2", "logR1", "logR2", "BAF1", "BAF2", "ZYG1", "ZYG2", "ZYG3",
        "ZYG4", "ZYGcons"
    ]
    zygRes = pd.DataFrame.from_records(zygData, columns=labels)

    return (zygRes)


if __name__ == "__main__":
    main()
