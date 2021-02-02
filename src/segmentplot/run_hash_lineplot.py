#!/usr/bin/env python3

from src.segmentplot.classes import Sequence
from src.segmentplot.classes import Segment
from src.segmentplot.hash_aligner import HashAligner


def select_longest(segments):
    longest_true = []
    longest_false = []
    for seg in segments:
        if seg.forward() == True:
            if len(longest_true) == 0:
                longest_true = [seg]
            else:
                if abs(seg.xEnd() - seg.xStart()) > abs(longest_true[0].xEnd() - longest_true[0].xStart()):
                    longest_true = [seg]
                elif abs(seg.xEnd() - seg.xStart()) == abs(longest_true[0].xEnd() - longest_true[0].xStart()):
                    longest_true.append(seg)

        else:
            if len(longest_false) == 0:
                longest_false = [seg]
            else:
                if abs(seg.xEnd() - seg.xStart()) > abs(longest_false[0].xEnd() - longest_false[0].xStart()):
                    longest_false = [seg]
                elif abs(seg.xEnd() - seg.xStart()) == abs(longest_false[0].xEnd() - longest_false[0].xStart()):
                    longest_false.append(seg)
    after_select = []
    after_select.extend(longest_true)
    after_select.extend(longest_false)

    return after_select

def cord_to_segments(cords):
    main_segs = []
    for cord in cords:
        if cord[2] == 0:
            forwoad = True
        else:
            forwoad = False
        x_start = cord[0][0]
        x_end = cord[0][1]
        y_start = cord[1][0]
        y_end = cord[1][1]
        length = int(y_end) - int(y_start) + 1

        main_segs.append(Segment(x_start, y_start, length, forwoad, 0))
    return main_segs


def hashplot_unmapped(ref, seq, k, min_accept):
    """
    Hashtable for realign unmapped seq
    :param ref: ref seq
    :param seq: unmapped seq
    :param k: k-mer size
    :param min_accept: min accept length for realigning
    :return: 
    """

    main_segs = None

    ref = Sequence(ref)
    read = Sequence(seq)

    repeat_thresh = 2
    aligner_ref = HashAligner(k, min_accept, 0, repeat_thresh)
    aligner_ref.run(ref, ref)
    diff_segs = aligner_ref.getSelfDiffSegs()

    y_hashvalues = aligner_ref.getHashValues()
    avoid_mers = aligner_ref.getAvoidKmer()

    aligner_merge = HashAligner(k, min_accept, 0, repeat_thresh)
    aligner_merge.run(read, ref, diff_segs, y_hashvalues, avoid_mers)
    segments_merge = aligner_merge.getMergeSegments()


    # select most long for both direction
    if len(segments_merge) >= 2:
        segments_merge = select_longest(segments_merge)


    return main_segs, segments_merge


