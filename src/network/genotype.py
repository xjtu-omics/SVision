#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/7/2

'''

import pysam


def genotyper(candidate, support_reads, options):
    gt = './.'
    homo_thresh = options.homo_thresh
    hete_thresh = options.hete_thresh

    bam = pysam.AlignmentFile(options.bam_path, 'r')
    contig, start, end, svtype = candidate[0], candidate[1], candidate[2], candidate[3]

    contig_length = bam.get_reference_length(contig)
    aligns = bam.fetch(contig=contig, start=max(0, start-1000), stop=min(contig_length, end+1000))

    support_alt_reads = set([read for read in support_reads])

    support_ref_reads = set()

    aln_no = 0
    while aln_no < 500:
        try:
            current_alignment = next(aligns)
        except StopIteration:
            break

        if current_alignment.query_name in support_alt_reads:
            continue
        if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
            continue
        aln_no += 1

        if len(svtype) == 1:
            if svtype[0] == "DEL" or svtype[0] == "INV":
                minimum_overlap = min((end - start) / 2, 2000)
                if (current_alignment.reference_start < (end - minimum_overlap) and current_alignment.reference_end > (end + 100) or
                        current_alignment.reference_start < (start - 100) and current_alignment.reference_end > (start + minimum_overlap)):
                    support_ref_reads.add(current_alignment.query_name)

            if svtype[0] == "INS" or svtype[0] == "DUP":
                if current_alignment.reference_start < (start - 100) and current_alignment.reference_end > (end + 100):
                    support_ref_reads.add(current_alignment.query_name)
        else:
            support_ref_reads.add(current_alignment.query_name)

    alt_no = len(support_alt_reads)
    ref_no = len(support_ref_reads)

    if len(svtype) != 1:
        return gt, ref_no, alt_no

    if alt_no + ref_no >= options.min_gt_depth:
        ratio = alt_no / (alt_no + ref_no)
        if ratio >= homo_thresh:
            gt = '1/1'
        elif ratio >= hete_thresh and ratio < homo_thresh:
            gt = '0/1'
        elif ratio < hete_thresh:
            gt = '0/0'

    return gt, ref_no, alt_no