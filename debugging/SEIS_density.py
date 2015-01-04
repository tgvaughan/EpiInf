#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import argv, exit

class TreeEvent:
    def __init__(self, time, isLeaf):
        self.time = time
        self.isLeaf = isLeaf

    def __repr__(self):
        if self.isLeaf:
            return "Sample at t=" + str(self.time)
        else:
            return "Coalescence at t=" + str(self.time)

def loadTreeEvents(treeFile):
    treeEvents = []

    for line in treeFile.readlines():
        lineElements = line.strip().split(" ")
        time = float(lineElements[0])
        isLeaf = (lineElements[1] == "0")
        treeEvents.append(TreeEvent(time, isLeaf))

    return treeEvents


class ParticleState:
    def __init__(self, S, E, I, R, kE, kI):
        self.S = S
        self.I = I
        self.R = R
        self.kE = kE
        self.kI = kI


if __name__ == '__main__':

    parser = ArgumentParser(description="Compute parameter likelihood of SEIS transmission tree.")
    parser.add_argument("treefile", type=FileType('r'), help="File containing transmission tree (ExpoTree format)")
    parser.add_argument("-b", type=float, dest='beta', help="Infection rate")
    parser.add_argument("-a", type=float, dest='alpha', help="Activation rate")
    parser.add_argument("-g", type=float, dest='gamma', help="Recovery rate")
    parser.add_argument("-s", type=float, dest='psi', help="Psi sampling rate")

    args = parser.parse_args()

    print loadTreeEvents(args.treefile)
