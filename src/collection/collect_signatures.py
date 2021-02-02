#!/usr/bin/env python3

# encoding: utf-8

import logging
from src.collection.analyze_reads import analyze_between_aligns, analyze_inside_align
import pysam
import os

from src.collection.analyze_reads import analyze_gap
import re
from src.collection.graph import generate_graph

# Declarations from pysam doc
M = "M" # M
I = "I" # I
D = "D" # D
N = "N" # N
S = "S" # S
H = "H" # H
P = "P" # P
E = "E" # =
X = "X" # X

MATCH_SET = {M, X, E}

def cigar_to_list(cigar):

    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    return ops, lengths

def process_cigars(align_lengths, align_ops):
    '''
    Process cigars of an alignment, Find indels and remove corresponding operations
    :return:
    '''

    for i in range(1, len(align_lengths) - 1):
        if align_ops[i-1] != 'M' and align_ops[i+1] != 'M' \
                and align_ops[i-1] != align_ops[i + 1] \
                and align_ops[i] == 'M' and align_lengths[i - 1] == align_lengths[i + 1] \
                and align_lengths[i] < 4:

            align_lengths[i-1] = 0
            align_lengths[i+1] = 0

    new_lengths = []
    new_ops = []

    for i in range(0, len(align_lengths)):
        if align_lengths[i] != 0:
            new_lengths.append(align_lengths[i])
            new_ops.append(align_ops[i])

    packed_ops = align_ops
    packed_lengths = align_lengths

    return packed_ops, packed_lengths

def create_align(ref_id, old_align,):
    # nm = old_align.get_tag('NM')


    new_align = pysam.AlignedSegment()

    new_align.reference_id = ref_id
    new_align.reference_start = old_align.reference_start

    new_align.query_name = old_align.query_name
    #

    # if not old_align.is_reverse:
    #     new_align.flag = 2048
    # else:
    #     new_align.flag = 2064

    new_align.is_supplementary = old_align.is_supplementary
    new_align.is_reverse = old_align.is_reverse

    if not new_align.is_supplementary:
        new_align.query_sequence = old_align.query_sequence

    try:
        new_align.mapping_quality = old_align.mapq
    except OverflowError:
        new_align.mapping_quality = 0

    new_align.cigarstring = str(old_align.cigarstring).replace('H', "S")
    new_align.next_reference_id = -1
    new_align.next_reference_start = -1
    new_align.template_length = 0
    # new_align.query_qualities = old_align.query_qualities
    # new_align.set_tags([("NM", nm, "i")])


    return new_align



def debug_plot(sorted_segs_list, qname, options):
    # # DEBUG code2: output dotplots if required
    from src.segmentplot.plot_segment import PlotSingleImg2
    dotplots_out_path = os.path.join(options.out_path, 'dotplots')
    if not os.path.exists(dotplots_out_path):
        os.mkdir(dotplots_out_path)
    # plot
    ploter = PlotSingleImg2(sorted_segs_list, str(qname[: min(100, len(qname))]).replace('/', '_'),
                            dotplots_out_path, )
    ploter.plot()

def analyze_alignments(aligns, bam, options, part_num):
    """
    Analyse between and inside primary and supplementary segments, and then combining segs to find segment signatures
    :param aligns: Aligns need to be processed
    :param bam: original bam file
    :param options: parameter
    :param part_num:
    :return:
        seg_signatures: All found segment signatures
    """

    min_mapq = 0 if options.contig is True else options.min_mapq

    all_possible_chrs = pysam.FastaFile(options.genome).references

    # # collect reads's pm and sa
    reads_dict = {}
    for align in aligns:

        # if align.qname not in ['m54329U_190629_180018/78446969/ccs']:
        #     continue

        # # no cigar, then pass this align
        if align.cigarstring == None:
            continue
        # # unmapped or secondary or low mapping quality, then pass this align
        if align.is_unmapped or align.is_secondary or align.mapq < min_mapq:
            continue

        # if not align.is_supplementary and align.mapq < options.min_mapq:
        #     continue

        # # align to a ref that not in genome reference
        align_chr = align.reference_name
        if align_chr not in all_possible_chrs:
            print("[Warning]: '{0}' not in reference's .fa file, skip this read".format(align_chr))
            continue

        align_ref_id = bam.get_tid(align.reference_name)

        # create new align
        new_align = create_align(align_ref_id, align)

        if align.qname not in reads_dict.keys():
            reads_dict[align.qname] = [new_align]
        else:
            reads_dict[align.qname].append(new_align)


    # traverse reads
    read_num = 0
    seg_signatures = []
    for qname in reads_dict.keys():
        # print('---------------------------', qname)
        # for align in reads_dict[qname]:
        #     print(align.reference_start, align.reference_end, align.query_alignment_start, align.query_alignment_end, align.is_supplementary, align.is_reverse)

        # separate pm and sa
        pm_align = None
        supp_aligns = []
        for align in reads_dict[qname]:
            if not align.is_supplementary:
                pm_align = align
            else:
                supp_aligns.append(align)

        # no primary align, we cannot adjust aligns' forward, then we skip
        if pm_align == None:
            continue

        # set sa's query sequence to pm align's query sequence
        for sa in supp_aligns:
            sa.query_sequence = pm_align.query_sequence

        whole_read_seq = pm_align.query_sequence

        read_num += 1
        all_segs_dicts = []


        # # analyze between aligns: primary and its supplementary
        major_segs_dicts_from_intra, minor_segs_dicts_from_intra = analyze_between_aligns(pm_align, supp_aligns, bam, options)

        # # minor segments from 'between aligns' are saved directly and will not be processed any more
        all_segs_dicts.extend(minor_segs_dicts_from_intra)

        # # for all major segments from 'between aligns', apply inside analyze
        for seg_dict in major_segs_dicts_from_intra:

            # # cigar process, Find indels and remove corresponding operations
            align_cigar_ops, align_cigar_lengths = cigar_to_list(seg_dict['cigarstring'])

            # # analyze inside align, then we can get new major and minor segments
            major_segs_dicts, minor_segs_dicts = analyze_inside_align(seg_dict, align_cigar_ops, align_cigar_lengths, options)


            # # there is no SV inside seg
            if major_segs_dicts == None and minor_segs_dicts == None:
                all_segs_dicts.append(seg_dict)
            # # there is a SV inside seg
            else:
                all_segs_dicts.extend(major_segs_dicts)
                all_segs_dicts.extend(minor_segs_dicts)
        sorted_segs_list = sorted(all_segs_dicts, key=lambda aln: (aln['q_start'], aln['q_end']))

        # print(pm_align.query_sequence)
        # print('---------------------------', qname)
        # for align in sorted_segs_list:
        #     print(align['ref_start'], align['ref_end'], align['q_start'], align['q_end'], align['is_reverse'])

        # only one seg or no segs, means no signature
        if len(sorted_segs_list) == 1 or len(sorted_segs_list) == 0:
            continue

        elif len(sorted_segs_list) == 2:

            # debug_plot(sorted_segs_list, qname, options)


            current_align = sorted_segs_list[0].copy()
            next_align = sorted_segs_list[1].copy()

            graph = None
            if options.report_graph is True:
                graph = generate_graph(current_align, next_align, [], options.min_sv_size, whole_read_seq, options.genome)
            sig = analyze_gap(current_align, next_align, bam, options)
            if sig is not None:

                sig.set_graph(graph)
                seg_signatures.append(sig)

        else:

            # debug_plot(sorted_segs_list, qname, options)

            # in case the fist seg is reverse
            if sorted_segs_list[0]['is_reverse']:
                current_align = sorted_segs_list[0].copy()
                next_align = sorted_segs_list[1].copy()

                graph = None
                if options.report_graph is True:
                    graph = generate_graph(current_align, next_align, [], options.min_sv_size, whole_read_seq, options.genome)
                sig = analyze_gap(current_align, next_align, bam, options)
                if sig is not None:
                    sig.set_graph(graph)
                    seg_signatures.append(sig)

            # in case the last seg is reverse
            if sorted_segs_list[-1]['is_reverse']:

                current_align = sorted_segs_list[-2].copy()
                next_align = sorted_segs_list[-1].copy()

                graph = None
                if options.report_graph is True:
                    graph = generate_graph(current_align, next_align, [], options.min_sv_size, whole_read_seq, options.genome)
                sig = analyze_gap(current_align, next_align, bam, options)
                if sig is not None:
                    sig.set_graph(graph)
                    seg_signatures.append(sig)

            all_main_align = []
            all_main_index = []
            for i in range(len(sorted_segs_list)):
                if sorted_segs_list[i]['type'] == 'main':
                    all_main_align.append(sorted_segs_list[i])
                    all_main_index.append(i)

            # print('1111111111111111111111111', qname)
            # for align in all_main_align:
            #     print(align['ref_start'], align['ref_end'])

            for i in range(len(all_main_align) - 1):

                current_index = all_main_index[i]
                current_align = all_main_align[i].copy()
                next_index = all_main_index[i + 1]
                next_align = all_main_align[i + 1].copy()
                distance_on_read = next_align['q_start'] - current_align['q_end']
                if distance_on_read >= -25:
                    # segs between cur and next is help aligns
                    help_aligns = sorted_segs_list[current_index + 1: next_index]
                    graph = None
                    if options.report_graph is True:
                        if i != len(all_main_align) - 1 - 1:
                            graph = generate_graph(current_align, next_align, help_aligns, options.min_sv_size, whole_read_seq, options.genome, False)
                        else:
                            graph = generate_graph(current_align, next_align, help_aligns, options.min_sv_size, whole_read_seq, options.genome)

                    sig = analyze_gap(current_align, next_align, bam, options, help_aligns)
                    if sig is not None:
                        sig.set_graph(graph)

                        seg_signatures.append(sig)

    return seg_signatures

# def analyze_multi_alignments():
