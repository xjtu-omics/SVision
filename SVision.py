#!/usr/bin/env python3

import sys
import os
# for file in os.listdir(os.getcwd()):
#     if os.path.isdir(file):
#         sys.path.append(os.path.join(os.getcwd(), file))

from src.collection import run_collection
from src.network.predict import Predict
import datetime
from src.network.output import merge_split_vcfs, cluster_original_callset, cal_scores_max_min
import shutil
import pysam
from src.collection.graph import collect_csv_same_format
from src.version import __version__

import argparse
import multiprocessing
import traceback


def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVision {0} \n \nShort Usage: Python SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>""".format(__version__))


    required_params = parser.add_argument_group("Input/Output parameters")
    required_params.add_argument('-o', dest="out_path", type=os.path.abspath, required=True, help='Absolute path to output ')
    required_params.add_argument('-b', dest='bam_path', type=os.path.abspath, required=True, help='Absolute path to bam file')
    required_params.add_argument('-m', dest="model_path", type=os.path.abspath, required=True, help='Absolute path to CNN predict model')
    required_params.add_argument('-g', dest='genome', type=os.path.abspath, required=True, help='Absolute path to your reference genome (.fai required in the directory)')
    required_params.add_argument('-n', dest='sample', type=str, required=True, help='Name of the BAM sample name')


    general_params = parser.add_argument_group("General parameters")
    general_params.add_argument('-t', dest="thread_num", type=int, default=1, help='Thread numbers [1]')
    general_params.add_argument('-s', dest="min_support", type=int, default=1, help='Min support read number for an SV [1]')
    general_params.add_argument('-c', dest="chrom", type=str, default=None, help='Specific region to detect, format: chr1:xxx-xxx or 1:xxx-xxx')
    general_params.add_argument('--hash_table', action="store_true", default=False,
                                 help='Activate hash table to align unmapped sequences')
    general_params.add_argument('--cluster_callset', action="store_true", default=False,
                                 help='Cluster original callset to merge uncovered event')
    general_params.add_argument('--report_mechanism', action="store_true", default=False,
                                 help='Report mechanisms for DEL event')
    general_params.add_argument('--report_graph', action="store_true", default=False,
                                 help='Report graph for events')

    general_params.add_argument('--contig', action="store_true", default=False,
                                 help='Activate contig mode')

    # general_params.add_argument('--rpmask', type=os.path.abspath, default='repeatmasker',
    #                              help='Path to RepeatMasker')
    # general_params.add_argument('--trf', type=os.path.abspath, default='trf',
    #                              help='Path to TRF')

    collect_params = parser.add_argument_group("Collect parameters")

    collect_params.add_argument("--min_mapq", type=int, default=10, help='Minimum mapping quality of reads to consider [10]')
    collect_params.add_argument("--min_sv_size", type=int, default=50, help='Minimum SV size to detect [50]')
    collect_params.add_argument("--max_sv_size", type=int, default=1000000, help='Maximum SV size to detect [1000000])')

    cluster_params = parser.add_argument_group("Cluster parameters")
    cluster_params.add_argument("--patition_max_distance", type=int, default=5000,
                                help='Maximum distance to partition signatures [5000]')
    cluster_params.add_argument("--cluster_max_distance", type=float, default=0.3,
                                help='Clustering maximum distance for a partition [0.3]]')


    hash_params = parser.add_argument_group("Hash table parameters")
    hash_params.add_argument("--k_size", type=int, default=10, help='Size of kmer [10]')
    hash_params.add_argument("--min_accept", type=int, default=50, help='Min match length to accept [50]')
    hash_params.add_argument("--max_hash_len", type=int, default=256, help='Max unmapped length to hashtable [256]')


    predict_params = parser.add_argument_group("Predict parameters")
    predict_params.add_argument("--batch_size", type=int, default=128, help='Batch size for the CNN prediction model [128]')

    options = parser.parse_args(arguments)



    return options

def main():
    options = parse_arguments()

    work_dir = options.out_path
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    # # SVision v1.0.3. ADD. log file
    log_file = open(os.path.join(work_dir, 'log.txt'), 'w')
    # # End ADD

    sample_path = options.bam_path
    aln_file = pysam.AlignmentFile(sample_path)

    print('[Processing]: Bam file at ', sample_path)

    window_size = 10000000  # slip on the chrom

    # chroms need to be detected
    refine_flag = 0
    start = 0
    end = 0
    if options.chrom == None:
        all_possible_chrs = pysam.FastaFile(options.genome).references

        chroms = []
        for line in aln_file.get_index_statistics():
            mapped_num = line[1]
            chrom = line[0]
            if mapped_num >= options.min_support and chrom in all_possible_chrs:
                chroms.append(chrom)
    else:
        refine_flag = 1

        region = options.chrom
        chroms = [str(region).split(':')[0]]
        cords = str(region).split(':')[1]
        start = int(cords.split('-')[0])
        end = int(cords.split('-')[1])

    if len(chroms) == 0:
        print('[ERROR]: No mapped reads in this bam. Exit!')
        exit()

    start_time = datetime.datetime.now()

    # clusters_out_path = os.path.join(options.out_path, 'clusters')
    segments_out_path = os.path.join(options.out_path, 'segments')
    predict_results_dir = os.path.join(work_dir, 'predict_results')
    if not os.path.exists(predict_results_dir):
        os.mkdir(predict_results_dir)
    if not os.path.exists(segments_out_path):
        os.mkdir(segments_out_path)
    if options.report_graph is True:
        graph_out_path = os.path.join(options.out_path, 'graphs')
        if os.path.exists(graph_out_path):
            shutil.rmtree(graph_out_path)
        os.mkdir(graph_out_path)

    # # DEBUG code1: create out dir
    # if os.path.exists(clusters_out_path):
    #     shutil.rmtree(clusters_out_path)
    # os.mkdir(clusters_out_path)

    process_pool = multiprocessing.Pool(processes=options.thread_num)
    pool_rets = []

    if refine_flag == 0:

        # traverse chroms to collect signatures and clusters
        for chrom in chroms:
            ref_len = aln_file.get_reference_length(chrom)
            part_num = 0

            if options.contig is True:
                window_size = ref_len
                run_collection.run_detect(options, sample_path, chrom, part_num, ref_len)
            else:

                while True:
                    if part_num * window_size > ref_len:
                        break

                    # SVision v1.0.3. MODIFY. more info to log out
                    pool_rets.append([process_pool.apply_async(run_collection.run_detect,
                                                               (options, sample_path, chrom, part_num, window_size)),
                                      chrom, part_num * window_size, (part_num + 1) * window_size])
                    # run_collection.run_detect(options, sample_path, chrom, part_num, window_size)
                    # End MODIFY

                    part_num += 1
    else:
        chrom = chroms[0]
        pool_rets.append(
            [process_pool.apply_async(run_collection.run_refine, (options, sample_path, chrom, start, end)),
             chrom, start, end])

    process_pool.close()
    process_pool.join()

    # exit()
    # # SVision v1.0.3. MODIFY. more info to log out
    # Error Exception
    error_num = 0
    for ret in pool_rets:
        ret_get = ret[0].get()
        if ret_get is not None:
            error_num += 1
            log_file.write('[Collecting Error]: at ' + ret[1] + '\t' + str(ret[2]) + '\t' + str(ret[3]) + '\n')
            log_file.write(ret_get + '\n')

    if error_num > 0:
        print('\n[Warning]: {0} threads failed when collecting. See log.txt for More!'.format(error_num))
    # # End MODIFY

    # merge splited files
    for chrom in chroms:
        # merge splited beds
        chr_segments_file = os.path.join(segments_out_path, chrom + ".segments.all.bed")
        cmd_str = 'cat {0}/{1}.segments.*.bed > {2}'.format(segments_out_path, chrom, chr_segments_file)
        os.system(cmd_str)

    end_time1 = datetime.datetime.now()
    cost_time = (end_time1 - start_time).seconds
    print("[Finished]: Collect segment signatures, Cost time: " + str(cost_time))

    # # begin to predict types
    def predict_one_chrom(chrom, predict_results_dir, options):
        try:
            segments_out_file = os.path.join(segments_out_path, chrom + ".segments.all.bed")
            chrom_predict_path = os.path.join(predict_results_dir, chrom + '.predict.' + 's' + str(options.min_support))

            predict = Predict(chrom, segments_out_file)
            predict.run(chrom_predict_path, options)
            return None
        except:
            error_type, error_value, error_trace = sys.exc_info()
            # print("[ERROR]: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace)))
            return "[ERROR]: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))

    process_pool = multiprocessing.Pool(processes=max(1, int(options.thread_num / 3)))
    pool_rets = []

    for chrom in chroms:
        # predict_one_chrom (chrom, predict_results_dir, options)
        pool_rets.append([process_pool.apply_async(predict_one_chrom, (chrom, predict_results_dir, options)), chrom])
    process_pool.close()
    process_pool.join()

    # # Modify. Deepsv v1.1.2
    error_num = 0
    for ret in pool_rets:
        ret_get = ret[0].get()
        if ret_get is not None:
            error_num += 1
            log_file.write('[Predicting Error]: at ' + ret[1] + '\n')
            log_file.write(ret_get + '\n')

    if error_num > 0:
        print('\n[Warning]: {0} threads failed when predicting. See log.txt for More!'.format(error_num))

    end_time2 = datetime.datetime.now()
    cost_time = (end_time2 - end_time1).seconds
    print("[Finished]: Predicting types, Cost time: " + str(cost_time))

    # # begin to score SV and merge chrom's results together
    max_score, min_score = cal_scores_max_min(predict_results_dir)
    merged_vcf_path = os.path.join(options.out_path, "{0}.svision.s{1}.vcf".format(options.sample, options.min_support))
    merge_split_vcfs(predict_results_dir, merged_vcf_path, options.genome, max_score, min_score)

    # # the following: run repeatmasker and trf to annotate TE and VNTR, non use for current version, but in v2.0
    # # # analyse mechanism
    # if options.report_mechanism == True:
    #     print("[Additional Func: Start]: Starting ananlyze mechanim......" )
    #     mechanism_dir = os.path.join(work_dir, 'mechanism')
    #     if not os.path.exists(mechanism_dir):
    #         os.mkdir(mechanism_dir)
    #     analyze_mechanism(merged_vcf_path, mechanism_dir, options)

    # # cluster original callset if required
    if options.cluster_callset is True:
        print("[Additional Func: Start]: Starting cluster original callset......")
        cluster_out_file = os.path.join(work_dir,
                                        "{0}.svision.s{1}.clusterd.vcf".format(options.sample, options.min_support))
        cluster_original_callset(merged_vcf_path, work_dir, sample_path, cluster_out_file)

    # # stats graphs file.
    if options.report_graph is True:
        graph_out_path = os.path.join(options.out_path, 'graphs')
        collect_csv_same_format(graph_out_path, merged_vcf_path, options.out_path, options.sample, options.min_support)

    # # rm tmp folder
    # shutil.rmtree(clusters_out_path)
    shutil.rmtree(segments_out_path)
    shutil.rmtree(predict_results_dir)

    log_file.close()

    end_time3 = datetime.datetime.now()
    cost_time = (end_time3 - start_time).seconds
    print("[Finished]: All. Total Cost time: " + str(cost_time) + "s")


if __name__ == '__main__':
    main()



 # python SVision.py -o /mnt/e/data_more/CCS/HG00733/ngmlr/sv_call/test/ -b /mnt/e/data_more/CCS/HG00733/ngmlr/chr22.bam -m /mnt/d/Workspace/Projects/DeepSV/src/model_ft/svision-cnn-model.ckpt -g /mnt/d/Data/Bams/REF/GRCh38/GRCh38_chr1-x.fa -t 4 -n test --report_graph -s 5