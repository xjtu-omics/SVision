#!/usr/bin/env python3

class Sequence():
    def __init__(self, sequence):
        self.bases = sequence
        self.segSequence = []

    def addSequence(self, sequence):
        sequence = sequence.upper()
        self.bases += sequence

    def clear(self):
        self.bases = ""

    def length(self):
        return len(self.bases)

    def getBases(self):
        return self.bases

    def getReverseComplementBases(self):
        inv_seq = ""
        for i in range(len(self.bases) - 1, -1, -1):
            bp = self.bases[i]
            inv_bp = ''
            if bp == 'A':
                inv_bp = 'T'
            elif bp == 'T':
                inv_bp = 'A'
            elif bp == 'C':
                inv_bp = 'G'
            elif bp == 'G':
                inv_bp = 'C'
            else:
                inv_bp = 'N'

            inv_seq += inv_bp

        return inv_seq


class Segment():

    def __init__(self, x_start, y_start, length, forward, seg_id):
        self._xStart = x_start
        self._yStart = y_start
        self._length = length
        self._forward = forward
        self._segId = seg_id
        if forward:
            self._xEnd = self._xStart + (self._length - 1)
        else:
            self._xEnd = self._xStart - (self._length - 1)
        self._yEnd = self._yStart + (length - 1)

    def setxEnd(self, xEnd):
        self._xEnd = xEnd

    def setyEnd(self, yEnd):
        self._yEnd = yEnd

    def setLength(self, length):
        self._length = length

    def setxStart(self, xStart):
        self._xStart = xStart

    def setyStart(self, yStart):
        self._yStart = yStart

    def setForward(self, forward):
        self._forward = forward

    def setSegId(self, segId):
        self._segId = segId

    def toString(self):
        # coordS = "[" + str(self._xStart) + "," + str(self._yStart) + "]"
        # coordE = "[" + str(self._xEnd) + "," + str(self._yEnd) + "]"
        # return self._segId + "-" + coordS + "-" + coordE + "-" + self._forward + " " + self.length

        return str(self._xStart) + '\t' + str(self._xEnd) + '\t' + str(self._yStart) + '\t' + str(self._yEnd) + '\t' \
               + str(self._forward)

    def xStart(self):
        return self._xStart

    def yStart(self):
        return self._yStart

    def xEnd(self):
        return self._xEnd

    def yEnd(self):
        return self._yEnd

    def segId(self):
        return self._segId

    def length(self):
        return abs(self._xEnd - self._xStart)

    def Id(self):
        return self._segId

    def forward(self):
        return self._forward
