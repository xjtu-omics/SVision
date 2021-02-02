#!/usr/bin/env python3

# encoding: utf-8



from src.collection.collect_signatures import analyze_alignments
from src.collection.cluster_signatures import partition_and_cluster
import pysam
from src.collection.output_clusters import writer_cluster_to_file
import sys
import traceback

def run_detect(options, sample_path, chrom, part_num, window_size):
    if part_num == 0:
        print('[Processing]: Collecting ' + chrom)

    try:
        fai_file = options.genome + ".fai"
        genome_file = open(options.genome, "r")

        # # fctch from align_file to get the reads
        aln_file = pysam.AlignmentFile(sample_path)
        aligns = aln_file.fetch(chrom, part_num * window_size, (part_num + 1) * window_size)

        # # collect signatures
        sv_signatures = analyze_alignments(aligns, aln_file, options, part_num)

        # for sig in sv_signatures:
        #     print(sig.tstart, sig.tend, sig.qname)

        # # partition and clusters signatures to get cluster info
        clusters = partition_and_cluster(sv_signatures, fai_file, genome_file, chrom, sample_path, options)

        # # write cluster info to file
        writer_cluster_to_file(clusters, chrom, part_num, options)

        return None
    except:
        error_type, error_value, error_trace = sys.exc_info()
        # print(error_value, str(traceback.extract_tb(error_trace)))
        return "[ERROR]: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))



def run_refine(options, sample_path, chrom, start, end):
    part_num = 0
    try:
        fai_file = options.genome + ".fai"
        genome_file = open(options.genome, "r")

        # # fctch from align_file to get the reads
        aln_file = pysam.AlignmentFile(sample_path)
        aligns = aln_file.fetch(chrom, start, end)

        # # collect signatures
        sv_signatures = analyze_alignments(aligns, aln_file, options, part_num)

        # # partition and clusters signatures to get cluster info
        clusters = partition_and_cluster(sv_signatures, fai_file, genome_file, chrom, sample_path, options)

        # # write cluster info to file
        writer_cluster_to_file(clusters, chrom, part_num, options)

        return None
    except:
        error_type, error_value, error_trace = sys.exc_info()
        return "[ERROR]: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))