#!/usr/bin/env python3

# encoding: utf-8


from src.collection.classes import Signature
import pysam
import re
from src.segmentplot.run_hash_lineplot import hashplot_unmapped


def shift_left(ref_seq, ref_start, target_start, target_end):
    shift_len = 0
    # print('--', len(ref_seq), ref_start, target_start, target_end)
    relative_start = target_start - ref_start
    relative_end = target_end - ref_start
    ref_len = len(ref_seq)

    # out of range of ref seq
    if relative_start >= ref_len or relative_end >= ref_len:
        return target_start, target_end

    while True:

        # meet the leftmost of ref seq, then break
        if target_start - ref_start <= 0:
            break
        # previous base is the same with the last base, then we can shift
        # print(shift_len, ref_len, relative_start - shift_len - 1, relative_end - shift_len)

        if ref_seq[relative_start - shift_len - 1] == ref_seq[relative_end - shift_len]:
            shift_len += 1
            target_start -= 1
            target_end -= 1

        else:
            break

    return target_start, target_end


def fetch_ref_seq(ref_path, chr, start, end):
    ref = pysam.FastaFile(ref_path)

    ref_cutted = ref.fetch(chr, start, end)
    return ref_cutted


def cal_overlap_ratio(base_seg, target_seg, left_most, right_most):
    # print( base_seg['ref_start'], base_seg['ref_end'])

    # # SVision1.0.1 ADD. Fix bug
    if base_seg == target_seg:
        return 0
    if base_seg['ref_start'] < left_most:
        # print(-1, base_seg['ref_start'] , left_most)
        return 1.0
    if base_seg['ref_end'] > right_most:
        # print(-2, base_seg['ref_end'] , right_most)
        return 1.0
    # # End ADD

    base_len = base_seg['ref_end'] - base_seg['ref_start']
    # base is totaly coverd by target
    if base_seg['ref_start'] >= target_seg['ref_start'] and base_seg['ref_end'] <= target_seg['ref_end']:
        # print(1, target_seg['ref_start'], target_seg['ref_end'] )
        return 1.0
    #
    if base_seg['ref_end'] >= target_seg['ref_end'] > base_seg['ref_start'] and target_seg['ref_start'] < base_seg['ref_start']:
        covered_len = target_seg['ref_end'] - base_seg['ref_start']
        # print(2, target_seg['ref_start'], target_seg['ref_end'] )

        return covered_len / base_len
    if base_seg['ref_end'] < target_seg['ref_start'] < base_seg['ref_start'] and target_seg['ref_end'] > base_seg['ref_end']:
        covered_len = base_seg['ref_end'] - target_seg['ref_start']
        # print(3, target_seg['ref_start'], target_seg['ref_end'] )

        return covered_len / base_len

    return 0

def trim_segs(aligns_coverd, first_seg, last_seg):
    """
    Trim segs to balance them
    :param aligns_coverd:
    :param first_seg:
    :param last_seg:
    :return:
    """

    # # step 1: adjust the ratio between normal and gap
    distance_on_read = last_seg['q_start'] - first_seg['q_end']
    distance_on_ref = last_seg['ref_start'] - first_seg['ref_end']
    gap = max(distance_on_read, distance_on_ref)

    # find left and right-most cords determined by gap
    left_most = first_seg['ref_end'] - gap * 2
    right_most = last_seg['ref_start'] + gap * 2

    # print(left_most, right_most)
    for seg in aligns_coverd:
        if seg == first_seg:

            # if first seg is longer than left_most, then cut it
            if seg['ref_start'] < left_most:
                adjust_len = left_most - seg['ref_start']
                seg['ref_start'] = left_most
                seg['q_start'] = seg['q_start'] + adjust_len

            # if first seg is shorter than left_most, then lengthen it
            elif seg['ref_start'] > left_most:
                adjust_len = seg['ref_start'] - left_most
                seg['ref_start'] = left_most
                seg['q_end'] += adjust_len

                # if lengthered, then we need to adjust other segs' cord on read
                for tmp_seg in aligns_coverd:
                    if tmp_seg != first_seg:
                        tmp_seg['q_start'] += adjust_len
                        tmp_seg['q_end'] += adjust_len

            # same length, pass
            else:
                pass

        elif seg == last_seg:
            # if last seg is longer than right_most, then cut it
            if seg['ref_end'] > right_most:
                adjust_len = seg['ref_end'] - right_most
                seg['ref_end'] = right_most
                seg['q_end'] = seg['q_end'] - adjust_len

            # if last seg is shorter than right_most, the lengthen it
            elif seg['ref_end'] < right_most:
                adjust_len = right_most - seg['ref_end']
                seg['ref_end'] = right_most
                seg['q_end'] += adjust_len

                # if lengthened, other segs need do nothing
            # same length, pass
            else:
                pass

        # adjust other segs to make sure it locate in the range of (leftmost, rightmost)
        else:
            seg_len = seg['q_end'] - seg['q_start']
            if seg['ref_start'] < left_most:
                seg['ref_start'] = left_most
                seg['ref_end'] = seg['ref_start'] + seg_len
            if seg['ref_end'] > right_most:
                seg['ref_end'] = right_most
                seg['ref_start'] = seg['ref_end'] - seg_len


def analyze_gap(current_align, next_align, bam, options, help_aligns=[]):
    """
    Analyse gaps to find segment signatures
    :param current_align: the first major seg
    :param next_align: the second major seg
    :param bam: bamfile
    :param options: parameters
    :param help_aligns: all segs that between the two segs
    :return:
        Segment signature
    """


    aligns_coverd = []
    aligns_coverd.extend(help_aligns)

    qname = current_align['read_name']
    if current_align['ref_id'] == next_align['ref_id']:

        ref_chr = bam.getrname(current_align['ref_id'])
        # Same orientation
        if current_align['is_reverse'] == next_align['is_reverse']:

            # # DeepSV v1.1.0. ADD. shift left operation
            ref_cords = [current_align['ref_start'], current_align['ref_end'], next_align['ref_start'], next_align['ref_end']]
            ref_start = min(ref_cords)
            ref_end = max(ref_cords)
            ref_seq = fetch_ref_seq(options.genome, ref_chr, ref_start, ref_end)

            for seg in help_aligns:
                # reversed seg will not be shifted
                if seg['is_reverse'] == True:
                    continue
                else:
                    target_start = seg['ref_start']
                    target_end = seg['ref_end']
                    new_ref_start, new_ref_end = shift_left(ref_seq, ref_start, target_start, target_end)
                    seg['ref_start'] = new_ref_start
                    seg['ref_end'] = new_ref_end
            # End ADD


            distance_on_read = next_align['q_start'] - current_align['q_end']
            distance_on_ref = next_align['ref_start'] - current_align['ref_end']

            # no overlap on ref
            if distance_on_ref >= -options.min_sv_size:



                diff = distance_on_read - distance_on_ref
                # INS
                if diff >= options.min_sv_size:

                    aligns_coverd.append(current_align)
                    aligns_coverd.append(next_align)

                    sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                    # svision v1.2.1 modifiy
                    if next_align['ref_start'] <= current_align['ref_end']:
                        # modify
                        bkp_len = abs(distance_on_read) + abs(distance_on_ref)
                        all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]
                    else:
                        # modify
                        bkp_len = abs(distance_on_read)
                        all_bkps = [[current_align['ref_end'], next_align['ref_start'], bkp_len]]
                    # end

                    for align in sorted_aligns:
                        if align in help_aligns:
                            bkp_len = align['ref_end'] - align['ref_start']
                            all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                    all_left_cords = [all_bkps[0][0]]
                    all_right_cords = [all_bkps[0][1]]

                    for i in range(1, len(all_bkps)):
                        all_left_cords.append(all_bkps[i][0])
                        all_right_cords.append(all_bkps[i][1])

                    left_most_cord = min(all_left_cords)
                    right_most_cord = max(all_right_cords)

                    trim_segs(sorted_aligns, current_align, next_align)

                    if len(help_aligns) == 0:
                        return Signature(ref_chr, left_most_cord, right_most_cord + diff, "sigGap",
                                         qname, sorted_aligns, all_bkps, 'None')
                    else:
                        return Signature(ref_chr, left_most_cord, right_most_cord, "sigGap",
                                         qname, sorted_aligns, all_bkps, 'None')
                # DEL
                elif -options.max_sv_size <= diff <= -options.min_sv_size:

                    # # DeepSV v1.1.0. ADD. shift left operation
                    ref_cords = [current_align['ref_start'], current_align['ref_end'], next_align['ref_start'], next_align['ref_end']]
                    ref_start = min(ref_cords)
                    ref_end = max(ref_cords)
                    ref_seq = fetch_ref_seq(options.genome, ref_chr, ref_start, ref_end)

                    # shift to the left
                    target_start = current_align['ref_end']
                    target_end = next_align['ref_start']
                    new_target_start, new_target_end = shift_left(ref_seq, ref_start, target_start, target_end)

                    # update the cords
                    current_align['ref_end'] = new_target_start + 1
                    next_align['ref_start'] = new_target_end
                    # # End ADD


                    aligns_coverd.append(current_align)
                    aligns_coverd.append(next_align)

                    sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                    if next_align['ref_start'] <= current_align['ref_end']:
                        bkp_len = 1
                        all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]
                    else:
                        bkp_len = next_align['ref_start'] - current_align['ref_end']
                        all_bkps = [[current_align['ref_end'], next_align['ref_start'], bkp_len]]

                    for align in sorted_aligns:
                        if align in help_aligns:
                            bkp_len = align['ref_end'] - align['ref_start']
                            all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                    all_left_cords = [all_bkps[0][0]]
                    all_right_cords = [all_bkps[0][1]]

                    for i in range(1, len(all_bkps)):
                        all_left_cords.append(all_bkps[i][0])
                        all_right_cords.append(all_bkps[i][1])

                    left_most_cord = min(all_left_cords)
                    right_most_cord = max(all_right_cords)


                    # # Deepsv v1.1.4. Add. report repair mechanism
                    if len(help_aligns) is not 0:
                        mechanism = 'None'
                    else:
                        # breakpoint has insertion and insertion > 10bp
                        if distance_on_read > 10:
                            mechanism = 'MMBIR+{0}'.format(distance_on_read)
                        elif distance_on_read >= -2:
                            if distance_on_read >= 0:
                                mechanism = 'NHEJ+{0}'.format(distance_on_read)
                            else:
                                mechanism = 'NHEJ{0}'.format(distance_on_read)
                        elif distance_on_read >= -20:
                            mechanism = 'AltEJ{0}'.format(distance_on_read)
                        else:
                            mechanism = 'NAHR{0}'.format(distance_on_read)
                    # # End Add.

                    trim_segs(sorted_aligns, current_align, next_align)

                    return Signature(ref_chr, left_most_cord, right_most_cord, "sigGap", qname,
                                     sorted_aligns, all_bkps, mechanism)
                # INV
                else:

                    aligns_coverd.append(current_align)
                    aligns_coverd.append(next_align)

                    # cur and next is linear, indicates inv
                    if len(help_aligns) != 0:

                        sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                        if next_align['ref_start'] <= current_align['ref_end']:
                            bkp_len = 1
                            all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]
                        else:
                            bkp_len = next_align['ref_start'] - current_align['ref_end']
                            all_bkps = [[current_align['ref_end'], next_align['ref_start'], bkp_len]]

                        for align in sorted_aligns:
                            if align in help_aligns:
                                bkp_len = align['ref_end'] - align['ref_start']
                                all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                        all_left_cords = [all_bkps[0][0]]
                        all_right_cords = [all_bkps[0][1]]

                        for i in range(1, len(all_bkps)):
                            all_left_cords.append(all_bkps[i][0])
                            all_right_cords.append(all_bkps[i][1])

                        left_most_cord = min(all_left_cords)
                        right_most_cord = max(all_right_cords)

                        trim_segs(sorted_aligns, current_align, next_align)

                        if next_align['ref_start'] - current_align['ref_end'] > 0:
                            return Signature(ref_chr, left_most_cord, right_most_cord, "sigGap",
                                             qname, sorted_aligns, all_bkps, 'None')

            # overlap on ref, tDUP
            else:

                aligns_coverd.append(current_align)

                dup_len = abs(distance_on_ref)
                # extract the duo seg
                dup_seg = {
                    'q_start': next_align['q_start'],
                    'q_end': next_align['q_start'] + dup_len,

                    'qual': current_align['qual'],
                    'ref_id': current_align['ref_id'],

                    'ref_start': next_align['ref_start'],
                    'ref_end': next_align['ref_start'] + dup_len,
                    'is_reverse': current_align['is_reverse'],
                    'read_name': current_align['read_name']
                }
                aligns_coverd.append(dup_seg)
                # help_aligns.append(dup_seg)
                # next_seg is next_align - dup_seg
                new_next_align = {
                    'q_start': next_align['q_start'] + dup_len + 1,
                    'q_end': next_align['q_end'],

                    'qual': current_align['qual'],
                    'ref_id': current_align['ref_id'],

                    'ref_start': next_align['ref_start'] + dup_len + 1,
                    'ref_end': next_align['ref_end'],
                    'is_reverse': current_align['is_reverse'],
                    'read_name': current_align['read_name']
                }

                if new_next_align['q_end'] < new_next_align['q_start']:
                    new_next_align['q_end'] = dup_seg['q_end'] + dup_len
                    new_next_align['ref_end'] = dup_seg['ref_end'] + dup_len

                aligns_coverd.append(new_next_align)

                # svision v1.2.1 modifiy, pricious ins bkp
                sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                if new_next_align['ref_start'] <= current_align['ref_end']:
                    bkp_len = abs(distance_on_read) + abs(distance_on_ref)
                    all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]

                else:
                    bkp_len = abs(distance_on_read) + abs(distance_on_ref)
                    all_bkps = [[current_align['ref_end'], new_next_align['ref_start'], bkp_len]]
                # end

                for align in sorted_aligns:
                    if align in help_aligns or align == dup_seg:
                        bkp_len = align['ref_end'] - align['ref_start']
                        all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])


                all_left_cords = [all_bkps[0][0]]
                all_right_cords = [all_bkps[0][1]]

                for i in range(1, len(all_bkps)):
                    all_left_cords.append(all_bkps[i][0])
                    all_right_cords.append(all_bkps[i][1])

                left_most_cord = min(all_left_cords)
                right_most_cord = max(all_right_cords)

                trim_segs(sorted_aligns, current_align, new_next_align)

                return Signature(ref_chr, left_most_cord, right_most_cord, "sigDup",  qname,  sorted_aligns, all_bkps, 'None')

        # if cur and next are not same forward, we need to add a seg to make sure there are two forward segs
        else:

            aligns_coverd.append(current_align)
            aligns_coverd.append(next_align)
            # cur is forward and next is reverse
            if not current_align['is_reverse']:
                if len(help_aligns) == 0:
                    help_aligns = [next_align]
                    # new seg's len is the same as the for
                    new_len = current_align['q_end'] - current_align['q_start']
                    # overlap
                    if next_align['ref_end'] <= current_align['ref_end']:
                        added_next_align = {
                            'q_start': next_align['q_end'],
                            'q_end': next_align['q_end'] + new_len,
                            'qual': current_align['qual'],
                            'ref_id': current_align['ref_id'],
                            'ref_start': current_align['ref_end'],
                            'ref_end': current_align['ref_end'] + new_len,
                            'is_reverse': current_align['is_reverse'],
                            'read_name': current_align['read_name']
                        }

                        aligns_coverd.append(added_next_align)

                        sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                        if added_next_align['ref_start'] <= current_align['ref_end']:
                            bkp_len = 1
                            all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]
                        else:
                            bkp_len = added_next_align['ref_start'] - current_align['ref_end']
                            all_bkps = [[current_align['ref_end'], added_next_align['ref_start'], bkp_len]]

                        for align in sorted_aligns:
                            if align in help_aligns:
                                bkp_len = align['ref_end'] - align['ref_start']
                                all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                        all_left_cords = [all_bkps[0][0]]
                        all_right_cords = [all_bkps[0][1]]

                        for i in range(1, len(all_bkps)):
                            all_left_cords.append(all_bkps[i][0])
                            all_right_cords.append(all_bkps[i][1])

                        left_most_cord = min(all_left_cords)
                        right_most_cord = max(all_right_cords)

                        trim_segs(sorted_aligns, current_align, added_next_align)

                        return Signature(ref_chr, left_most_cord, right_most_cord, "sigUncovered", qname, sorted_aligns, all_bkps, 'None')

                    else:
                        fixed_inv_len = max((next_align['ref_end'] - current_align['ref_end']), (next_align['q_end'] - current_align['q_end']))
                        added_next_align = {
                            'q_start': current_align['q_end'] + fixed_inv_len,
                            'q_end': current_align['q_end'] + fixed_inv_len + new_len,
                            'qual': current_align['qual'],
                            'ref_id': current_align['ref_id'],
                            'ref_start': next_align['ref_end'],
                            'ref_end': next_align['ref_end'] + new_len,
                            'is_reverse': current_align['is_reverse'],
                            'read_name': current_align['read_name']
                        }
                        # next_align['ref_start'] = next_align['ref_end'] - fixed_inv_len
                        # next_align['q_end'] = next_align['q_start'] + fixed_inv_len
                        aligns_coverd.append(added_next_align)

                        sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                        if added_next_align['ref_start'] <= current_align['ref_end']:
                            bkp_len = 1
                            all_bkps = [[current_align['ref_end'], current_align['ref_end'] + 1, bkp_len]]
                        else:
                            bkp_len = added_next_align['ref_start'] - current_align['ref_end']
                            all_bkps = [[current_align['ref_end'], added_next_align['ref_start'], bkp_len]]

                        for align in sorted_aligns:
                            if align in help_aligns:
                                bkp_len = align['ref_end'] - align['ref_start']
                                all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                        all_left_cords = [all_bkps[0][0]]
                        all_right_cords = [all_bkps[0][1]]

                        for i in range(1, len(all_bkps)):
                            all_left_cords.append(all_bkps[i][0])
                            all_right_cords.append(all_bkps[i][1])

                        left_most_cord = min(all_left_cords)
                        right_most_cord = max(all_right_cords)


                        trim_segs(sorted_aligns, current_align, added_next_align)
                        return Signature(ref_chr, left_most_cord, right_most_cord, "sigUncovered", qname, sorted_aligns, all_bkps, 'None')
            else:
                if len(help_aligns) == 0:
                    help_aligns = [current_align]
                    new_len = next_align['q_end'] - next_align['q_start']
                    if current_align['ref_start'] >= next_align['ref_start']:

                        added_cur_align = {
                            'q_start': 0,
                            'q_end': 0 + new_len,
                            'qual': current_align['qual'],
                            'ref_id': current_align['ref_id'],
                            'ref_start': next_align['ref_start'] - new_len,
                            'ref_end': next_align['ref_start'] - 1,
                            'is_reverse': next_align['is_reverse'],
                            'read_name': current_align['read_name']
                        }
                        for align in aligns_coverd:
                            align['q_start'] += new_len
                            align['q_end'] += new_len
                        aligns_coverd.append(added_cur_align)

                        sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                        if next_align['ref_start'] <= added_cur_align['ref_end']:
                            bkp_len = 1
                            all_bkps = [[added_cur_align['ref_end'], added_cur_align['ref_end'] + 1], bkp_len]
                        else:
                            bkp_len = next_align['ref_start'] - added_cur_align['ref_end']
                            all_bkps = [[added_cur_align['ref_end'], next_align['ref_start'], bkp_len]]

                        for align in sorted_aligns:
                            if align in help_aligns:
                                bkp_len = align['ref_end'] - align['ref_start']
                                all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                        all_left_cords = [all_bkps[0][0]]
                        all_right_cords = [all_bkps[0][1]]

                        for i in range(1, len(all_bkps)):
                            all_left_cords.append(all_bkps[i][0])
                            all_right_cords.append(all_bkps[i][1])

                        left_most_cord = min(all_left_cords)
                        right_most_cord = max(all_right_cords)

                        trim_segs(sorted_aligns, added_cur_align, next_align)

                        return Signature(ref_chr, left_most_cord, right_most_cord,  "sigUncovered", qname, sorted_aligns, all_bkps, 'None')

                    else:
                        fixed_inv_len = max((next_align['ref_start'] - current_align['ref_start']), (next_align['q_start'] - current_align['q_start']))

                        added_cur_align = {
                            'q_start': 0,
                            'q_end': 0 + new_len,
                            'qual': current_align['qual'],
                            'ref_id': current_align['ref_id'],
                            'ref_start': next_align['ref_start'] - fixed_inv_len - new_len,
                            'ref_end': next_align['ref_start'] - fixed_inv_len - 1,
                            'is_reverse': next_align['is_reverse'],
                            'read_name': current_align['read_name']
                        }

                        new_len += abs((next_align['ref_start'] - current_align['ref_start']) - (next_align['q_start'] - current_align['q_start']))
                        for align in aligns_coverd:
                            align['q_start'] += new_len
                            align['q_end'] += new_len
                        aligns_coverd.append(added_cur_align)


                        sorted_aligns = sorted(aligns_coverd, key=lambda aln: (aln['q_start'], aln['q_end']))
                        if next_align['ref_start'] <= added_cur_align['ref_end']:
                            bkp_len = 1
                            all_bkps = [[added_cur_align['ref_end'], added_cur_align['ref_end'] + 1], bkp_len]
                        else:
                            bkp_len = next_align['ref_start'] - added_cur_align['ref_end']
                            all_bkps = [[added_cur_align['ref_end'], next_align['ref_start'], bkp_len]]

                        for align in sorted_aligns:
                            if align in help_aligns:
                                bkp_len = align['ref_end'] - align['ref_start']
                                all_bkps.append([align['ref_start'], align['ref_end'], bkp_len])

                        all_left_cords = [all_bkps[0][0]]
                        all_right_cords = [all_bkps[0][1]]

                        for i in range(1, len(all_bkps)):
                            all_left_cords.append(all_bkps[i][0])
                            all_right_cords.append(all_bkps[i][1])

                        left_most_cord = min(all_left_cords)
                        right_most_cord = max(all_right_cords)

                        trim_segs(sorted_aligns, added_cur_align, next_align)

                        return Signature(ref_chr, left_most_cord, right_most_cord, "sigUncovered", qname,  sorted_aligns, all_bkps, 'None')



def analyze_between_aligns(primary, supplmentary, bam, options):
    """
    Analyse between primary and supplementary aligns
    :param primary: primary align
    :param supplmentary: supplementary align
    :param bam: bam file
    :param options: parameters
    :return:
    """
    if options.contig == True:
        pass
    else:
        # Too much supplmentaries
        if len(supplmentary) > 4:
            return [], []
    read_name = primary.query_name

    alignments = [primary] + supplmentary
    primary_forward = primary.is_reverse

    major_segs = []
    minor_segs = []
    all_forward_segs = []

    # # traverse all aligns
    for alignment in alignments:
        ref_chr = bam.getrname(alignment.reference_id)
        # print(alignment.reference_start, alignment.reference_end, alignment.query_alignment_start, alignment.query_alignment_end,
        #       alignment.is_reverse, alignment.is_supplementary)

        # # focus on primary's forward, if other reads' forward is not same as primary, then adjust cords
        if alignment.is_reverse != primary_forward:
            q_start = alignment.query_length - alignment.query_alignment_end
            q_end = alignment.query_length - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        # # store the read's info in a dict
        alignment_dict = {  'q_start': q_start,
                            'q_end': q_end,
                            'qual': alignment.mapping_quality,
                            'ref_id': alignment.reference_id,
                            'ref_chr': ref_chr,
                            'ref_start': alignment.reference_start,
                            'ref_end': alignment.reference_end,
                            'read_name': read_name,
                            'cigarstring': alignment.cigarstring,
                            'read_seq': alignment.query_sequence[q_start: q_end]
                            }

        # set supplementary
        if alignment.is_supplementary:
            alignment_dict['is_supplementary'] = True
        else:
            alignment_dict['is_supplementary'] = False

        # reset seg's forward to the primary's forward
        # and those whose forward is the same as the primary, set reverse tag to False
        if alignment.is_reverse == primary_forward:
            alignment_dict['is_reverse'] = False
            all_forward_segs.append(alignment_dict)
        else:
            alignment_dict['is_reverse'] = True
            o_seg = alignment_dict
            o_seg['type'] = 'other'
            minor_segs.append(o_seg)
        # print(1, alignment_dict['ref_start'], alignment_dict['ref_end'], alignment_dict['q_start'], alignment_dict['q_end'], alignment_dict['is_reverse'])

    # # only one forward, then return
    if len(all_forward_segs) == 1:
        major_segs.append(all_forward_segs[0])
        for seg in major_segs:
            seg['type'] = 'main'
        for seg in minor_segs:
            seg['type'] = 'other'
        return major_segs, minor_segs

    sorted_forward_segs = sorted(all_forward_segs, key=lambda aln: (aln['q_start'], aln['q_end']))

    # # find segs that covered by other, if that, it is a other seg
    # left_most = sorted_forward_segs[0]['ref_start']
    # right_most = sorted_forward_segs[len(sorted_forward_segs) - 1]['ref_end']
    left_most = min([i['ref_start'] for i in sorted_forward_segs])
    right_most = max([i['ref_end'] for i in sorted_forward_segs])
    for i in range(len(sorted_forward_segs)):
        # the first and last seg are major segs
        if i == 0 or i == len(sorted_forward_segs) - 1:
            base_seg = sorted_forward_segs[i]
            base_seg['type'] = 'main'
            major_segs.append(base_seg)
            continue

        flag = 0
        base_seg = sorted_forward_segs[i]

        # other segs need to test if convered by others within a given ratio
        for j in range(len(sorted_forward_segs)):
            target_seg = sorted_forward_segs[j]

            overlap_ratio = cal_overlap_ratio(base_seg, target_seg, left_most, right_most)
            if overlap_ratio >= 0.8 and base_seg not in minor_segs:
                base_seg['type'] = 'other'
                minor_segs.append(base_seg)
                flag = 1
                break
        # not coverd by others, then it is a main seg
        if flag == 0:
            base_seg['type'] = 'main'
            major_segs.append(base_seg)


    if options.hash_table == True:
        # first find all the main segs' index
        main_segs_index = []
        all_segs = []
        all_segs.extend(major_segs)
        all_segs.extend(minor_segs)
        sorted_all_segs = sorted(all_segs, key=lambda aln: (aln['q_start'], aln['q_end']))
        for i in range(len(sorted_all_segs)):
            if sorted_all_segs[i]['type'] == 'main':
                main_segs_index.append(i)

        # # hashtable segs between main_segs when there is no other segs between them
        for i in range(len(main_segs_index) - 1):

            # there is no others between then
            if main_segs_index[i + 1] - main_segs_index[i] == 1:
                cur_main = sorted_all_segs[i].copy()
                next_main = sorted_all_segs[i + 1].copy()
                # and there is a gap on read
                if next_main['q_start'] - cur_main['q_end'] >= options.min_sv_size:
                    distance_on_read = next_main['q_start'] - cur_main['q_end']
                    distance_on_ref = next_main['ref_start'] - cur_main['ref_end']
                    diff = abs(distance_on_read - distance_on_ref)

                    if distance_on_ref >= -options.min_sv_size and diff >= options.min_sv_size:

                        ref_chr = cur_main['ref_chr']

                        # fetch read and ref seq
                        read_start = cur_main['q_end']
                        read_end = next_main['q_start']
                        read_seq = cur_main['read_seq'][read_start: read_end ]

                        ref_start = min(cur_main['ref_start'], next_main['ref_start'])
                        ref_end = max(cur_main['ref_end'], next_main['ref_end'])
                        ref_seq = fetch_ref_seq(options.genome, ref_chr, ref_start, ref_end)

                        # hashplot ins segment
                        # print(len(read_seq))
                        if len(read_seq) < options.max_hash_len:
                            # print('---------', len(read_seq), len(ref_seq))
                            m_segs, o_segs = hashplot_unmapped(ref_seq, read_seq, options.k_size, options.min_accept)
                            # print('done')
                            for seg in o_segs:
                                tmp_dict = {'q_start': seg.xStart() + read_start if seg.forward() else seg.xEnd() + read_start,
                                            'q_end': seg.xEnd() + read_start if seg.forward() else seg.xStart() + read_start,
                                            'qual': cur_main['qual'],
                                            'ref_id': cur_main['ref_id'],
                                            'ref_chr': ref_chr,
                                            'ref_start': seg.yStart() + ref_start,
                                            'ref_end': seg.yEnd() + ref_start,
                                            'read_name': cur_main['read_name'],
                                            'cigarstring': '',
                                            'type': 'other',
                                            'read_seq': read_seq,
                                            'is_reverse': False if seg.forward() else True,
                                            'is_supplementary': cur_main['is_supplementary']
                                            }

                                minor_segs.append(tmp_dict)

    for seg in major_segs:
        seg['type'] = 'main'
    for seg in minor_segs:
        seg['type'] = 'other'

    # for i in ((major_segs)):
    #     print('m', i['ref_start'], i['ref_end'], i['q_start'], i['q_end'], i['is_reverse'], i['is_supplementary'])
    # for i in ((minor_segs)):
    #     print('o', i['ref_start'], i['ref_end'], i['q_start'], i['q_end'], i['is_reverse'], i['is_supplementary'])
    return major_segs, minor_segs

    #
def analyze_inside_align(seg_dict, cigar_ops, cigar_lengths, options):
    """
    Analyse inside an align by processing its cigar
    :param seg_dict: a dict saving seg's info
    :param cigar_ops: cigar operations
    :param cigar_lengths: length for each cigar operations
    :param options: parameters
    :return:
    """

    readPos = 0
    refPos = seg_dict['ref_start']

    read_seq = seg_dict['read_seq']
    ref_start = seg_dict['ref_start']
    ref_end = seg_dict['ref_end']
    read_start = seg_dict['q_start']
    read_end = seg_dict['q_end']
    ref_chr = seg_dict['ref_chr']

    all_long_gaps = []
    all_ins_seqs = []

    # # traverse cigar to find gaps(I or D) longer than min_sv_size
    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        opLen = cigar_lengths[i]

        if op == "N" or op == "S":
            readPos += opLen

        elif op == "I":
            if opLen >= options.min_sv_size:
                all_long_gaps.append([[readPos, readPos + opLen], [refPos, refPos], 'I'])
                all_ins_seqs.append([readPos, readPos + opLen, refPos, refPos + 1, read_seq[readPos - read_start: readPos + opLen - read_start]])
            readPos += opLen

        elif op == "D":
            if opLen >= options.min_sv_size:
                all_long_gaps.append([[readPos, readPos], [refPos, refPos + opLen], 'D'])
            refPos += opLen

        elif op in ["M", "X", "E", '=']:
            refPos += opLen
            readPos += opLen

        elif op == "H":
            pass
        else:
            pass

    # no obvious gap, then return None
    if len(all_long_gaps) == 0:
        return None, None

    major_segs_cords = []
    minor_segs_cords = []
    major_segs_dicts = []
    minor_segs_dicts = []

    # # the followings: find major segs between long gaps
    # first gap, then seg before it is a main seg
    virtual_read_pos = read_start
    gap = all_long_gaps[0]
    m_length = gap[1][0] - ref_start
    major_segs_cords.append([virtual_read_pos, virtual_read_pos + m_length, ref_start, gap[1][0] - 1])
    virtual_read_pos += m_length

    # other gaps: there is a main seg between two gaps
    for i in range(len(all_long_gaps) - 1):

        cur_gap = all_long_gaps[i]
        # if curgap is a Ins, then add the length to virtual_read_pos
        if cur_gap[2] == 'I':
            m_length = cur_gap[0][1] - cur_gap[0][0]    # ins's length
            virtual_read_pos += m_length

        next_gap = all_long_gaps[i + 1]

        m_length = next_gap[1][0] - cur_gap[1][1]   # length between two gaps
        major_segs_cords.append([virtual_read_pos + 1, virtual_read_pos + m_length + 1, cur_gap[1][1], next_gap[1][0]])
        virtual_read_pos += m_length

    # last gap, then seg after this gap is a main seg
    gap = all_long_gaps[-1]
    if gap[2] == 'I':
        m_length = gap[0][1] - gap[0][0]
        virtual_read_pos += m_length

    m_length = ref_end - gap[1][1]
    major_segs_cords.append([virtual_read_pos + 1, virtual_read_pos + m_length + 1, gap[1][1], ref_end])
    # # end

    # # if there is a ins, then re-align to find if a dup exists
    if options.hash_table == True:

        for ins in all_ins_seqs:
            seg_read_start = ins[0]

            read_seq = ins[4]
            extend_num = 0
            ref_start = ref_start - extend_num
            ref_end = ref_end + extend_num
            ref_seq = fetch_ref_seq(options.genome, ref_chr, ref_start, ref_end)

            # hashplot ins segment
            if len(read_seq) < options.max_hash_len:
                main_segs, other_segs = hashplot_unmapped(ref_seq, read_seq, options.k_size, options.min_accept)
                for seg in other_segs:

                    tmp_dict = {'q_start': seg.xStart() + seg_read_start if seg.forward() else seg.xEnd() + seg_read_start,
                                'q_end': seg.xEnd() + seg_read_start if seg.forward() else seg.xStart() + seg_read_start,
                                'qual': seg_dict['qual'],
                                'ref_id': seg_dict['ref_id'],
                                'ref_chr': ref_chr,
                                'ref_start': seg.yStart() + ref_start,
                                'ref_end': seg.yEnd() + ref_start,
                                'read_name': seg_dict['read_name'],
                                'cigarstring': '',
                                'type': 'other',
                                'read_seq': '',
                                'is_reverse': False if seg.forward() else True,
                                'is_supplementary': seg_dict['is_supplementary']
                                }
                    minor_segs_dicts.append(tmp_dict)
    # print('-------------------------------------')
    # # convert seg cords to a dict format to save info
    for seg in major_segs_cords:

        new_seg_dict = {        'q_start': seg[0],
                            'q_end': seg[1],
                            'qual': seg_dict['qual'],
                            'ref_id': seg_dict['ref_id'],
                            'ref_chr': ref_chr,
                            'ref_start': seg[2],
                            'ref_end': seg[3],
                            'read_name': seg_dict['read_name'],
                            'cigarstring': '',
                            'type': 'main',
                            'read_seq': seg_dict['read_seq'][seg[0] - read_start: seg[1] - read_start],
                            'is_reverse': False,
                            'is_supplementary': seg_dict['is_supplementary']
                            }
        major_segs_dicts.append(new_seg_dict)
        # print(read_start, seg[0], seg[1], seg[1] - seg[0], len(new_seg_dict['read_seq']), 1)
    for seg in minor_segs_cords:

        new_seg_dict = {    'q_start': seg.xStart(),
                        'q_end': seg.xEnd(),
                        'qual': seg_dict['qual'],
                        'ref_id': seg_dict['ref_id'],
                        'ref_chr': ref_chr,
                        'ref_start': seg.yStart(),
                        'ref_end': seg.yEnd(),
                        'read_name': seg_dict['read_name'],
                        'cigarstring': '',
                         'type': 'other',
                         'read_seq': seg_dict[seg.xStart() - read_start: seg.xEnd() - read_start],
                        'is_reverse': False if seg.forward() == True else True,
                     'is_supplementary': seg_dict['is_supplementary']

                     }
        minor_segs_dicts.append(new_seg_dict)

    return major_segs_dicts, minor_segs_dicts
