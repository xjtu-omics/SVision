#!/usr/bin/env python3

import os
import collections
import pysam
import numpy as np
from collections import Counter
from src.network.annotation import parse_trf, parse_rpmask
import multiprocessing
def cal_new_cluster(item_list):
    """
    Given a list of items in vcf, this function merge them together
    :param item_list: a list of items from vcf
    :return:
        new cluster : [chrom, start_cord, end_cord, SVLEN, most_type, most_type_support, sorted_types, all_sv_reads,
                    clusterd_id]
    """

    start_cord = 0
    end_cord = 0
    length = 0
    chrom = ""
    clusterd_id = ''
    sv_type = ''
    vaf = 0
    qual = 0

    all_supports = 0
    all_sv_reads = ''
    all_bkps = {}

    item_num = len(item_list)
    for item in item_list:

        chrom = item[0]
        start_cord += int(item[1])
        end_cord += int(item[2])
        length += int(item[3])

        sv_type = item[5]
        support = item[6]
        reads = item[7]
        bkps = item[8]
        id = item[9]
        vaf = item[10]

        qual += int(item[11])

        if clusterd_id == '':
            clusterd_id += str(id)
        else:
            clusterd_id += '_' + str(id)

        # collect support reads
        for read in reads:
            all_sv_reads += read + ','
        all_sv_reads = all_sv_reads[0: -1]

        # sum up supports
        all_supports += int(support)

        # collect breakpoints
        for bkp in bkps:
            bkp_split = bkp.split(':')
            sub_type = bkp_split[0]
            start = int(bkp_split[1].split('-')[0])
            end = int(bkp_split[1].split('-')[1])

            if sub_type not in all_bkps.keys():
                all_bkps[sub_type] = [[], []]
            all_bkps[sub_type][0].append(start)
            all_bkps[sub_type][1].append(end)


    start_cord = int(start_cord / len(item_list))
    end_cord = int(end_cord / len(item_list))
    length = int(length / len(item_list))
    qual = int(qual / len(item_list))
    new_cluster = [chrom, start_cord, end_cord, length, sv_type, all_supports, all_sv_reads, all_bkps, clusterd_id, vaf, item_num, qual]

    return new_cluster


def convert_to_vcf_format(new_cluster, sample_path, filter_type):
    """
    Convert new cluster's info from cal_new_cluster to vcf format
    :param new_cluster: new_cluster from cal_new_cluster
    :param sample_path:
    :return:
        a line in vcf format
    """

    chrom = new_cluster[0]
    start = new_cluster[1]
    end = new_cluster[2]
    length = new_cluster[3]
    sv_type = new_cluster[4]
    sv_supports = new_cluster[5]

    all_reads = new_cluster[6]
    clusterd_id = new_cluster[8]
    all_bkps = new_cluster[7]
    item_num = new_cluster[10]
    qual = new_cluster[11]

    if item_num == 1:
        coverage = 0
    else:
        # stats coverage
        coverage = 0
        aln_file = pysam.AlignmentFile(sample_path)
        reads = aln_file.fetch(chrom, start, end)
        for i in reads:
            coverage += 1

    # create INFO
    svtype_info = "SVTYPE=" + sv_type
    svsupp_info = "SUPPORT=" + str(sv_supports)
    svvaf_info = "VAF="
    svreads_info = "READS=" + all_reads
    svbkp_info = 'BKPS='

    sub_types = sv_type.split('+')
    for sub_type in sub_types:
        start_mean_cord = int(np.mean(all_bkps[sub_type][0]))
        end_mean_cord = int(np.mean(all_bkps[sub_type][1]))

        svbkp_info += sub_type + ':' + str(start_mean_cord) + '-' + str(end_mean_cord) + ','

    svbkp_info = svbkp_info[0: -1]

    if item_num == 1:
        vaf = new_cluster[9]
    else:
        if coverage == 0:
            vaf = 1.0
        else:
            vaf = round(sv_supports / coverage, 2)

    svvaf_info += str(vaf)


    info = "END={0};SVLEN={1};{2};{3};{4};{5};{6}".format(end, length, svtype_info, svsupp_info, svbkp_info, svvaf_info, svreads_info)
    line = chrom + '\t' + str(start) + '\t' + clusterd_id + '\t' + 'N' + '\t' + sv_type + '\t' + str(qual) + '\t' + filter_type + '\t' + info

    return line



def cluster_original_callset(callset_path, out_path, sample_path, cluster_out_file):
    """
    Custer unconvered region in callset together

    :param callset_path:  VCF path
    :param out_path:    out path
    :param sample_path:  bam path
    :param cluster_out_file:  out file path
    :return:
        write to results to cluster_out_file
    """
    uncovered_list = []

    thresh = 100

    # two tmp files, one for normal (covered) and one for cluster
    normal_file = open(os.path.join(out_path, 'normal_tmp.txt'), 'w')
    clusterd_file = open(os.path.join(out_path, 'cluster_tmp.txt'), 'w')

    fin = pysam.VariantFile(callset_path, 'r')
    normal_file.write(str(fin.header))  # write vcf's header

    # collect all unconvered SV, write other SV to file
    for record in fin:
        record_str = str(record)
        line_split = record_str.strip().split('\t')
        sig_type = line_split[6]

        # unconvered and convered are split
        if sig_type == "Uncovered":
            sv_type = record.info['SVTYPE'].replace('t', '')
            vaf = record.info['VAF']
            uncovered_list.append([record.contig, record.pos, record.stop, record.info['SVLEN'], sv_type, record.info['SVTYPE'], record.info['SUPPORT'], record.info['READS'], record.info['BKPS'], record.id, vaf, record.qual])
        else:
            normal_file.write(record_str)

    # sort
    sorted_uncovered_list = sorted(uncovered_list, key=lambda evi: evi[1])

    accessed_flag = [0 for i in range(len(sorted_uncovered_list))]  # a list to store if already accessed

    # traverse to merge them together
    for i in range(len(sorted_uncovered_list)):
        # already accessed
        if accessed_flag[i] == 1:
            continue

        # for each item, itself form a new cluster first
        item_list = [sorted_uncovered_list[i]]
        new_cluster = cal_new_cluster(item_list)

        # traverse others to merge
        for j in range(i + 1, len(sorted_uncovered_list)):
            # different chrom
            if sorted_uncovered_list[j][0] is not new_cluster[0]:
                continue

            # different svtype
            if sorted_uncovered_list[j][4].replace('t', '') != new_cluster[4].replace('t', ''):
                continue

            # already accessed
            if accessed_flag[j] == 1:
                continue

            # similar cords
            if abs(int(sorted_uncovered_list[j][1]) - new_cluster[1]) <= thresh or abs(int(sorted_uncovered_list[j][2]) - new_cluster[2]) < thresh:
                # print('add: ', sorted_uncovered_list[j])
                accessed_flag[j] = 1
                item_list.append(sorted_uncovered_list[j])

        new_cluster = cal_new_cluster(item_list)

        # change filter type
        if len(item_list) == 1:
            filter_type = 'Uncovered'
        else:
            filter_type = 'Clustered'

        # convert to vcf format
        line = convert_to_vcf_format(new_cluster, sample_path, filter_type)
        # write to file
        clusterd_file.write(line + '\n')


    clusterd_file.close()
    normal_file.close()

    # merge normal file and clusterd file
    cmd_str = 'cat {0} {1} > {2}'.format(os.path.join(out_path, 'normal_tmp.txt'), os.path.join(out_path, 'cluster_tmp.txt'), cluster_out_file)
    os.system(cmd_str)
    # remove tmp file
    os.remove(os.path.join(out_path, 'normal_tmp.txt'))
    os.remove(os.path.join(out_path, 'cluster_tmp.txt'))




def merge_split_vcfs(in_dir, merged_vcf_path, ref_path, max_score, min_score):
    files = os.listdir(in_dir)

    merged_vcf = open(merged_vcf_path, 'w')
    from src.version import __version__
    # Write header lines
    print("##fileformat=VCFv4.3", file=merged_vcf)
    print("##source=DeepSV v{0}".format(__version__) , file=merged_vcf)

    ref = pysam.FastaFile(ref_path)
    chroms = ref.references
    for chr in chroms:
        chr_length = ref.get_reference_length(chr)
        print('##contig=<ID={0},length={1}>'.format(chr, chr_length), file=merged_vcf)

    print("##CHROM=<CHROM=XXX,Description=\"Chromosome ID\">", file=merged_vcf)
    print("##POS=<POS=XXX,Description=\"Start position of the SV described in this region\">", file=merged_vcf)
    print("##ID=<ID=XXX,Description=\"ID of the SV described in this region\">", file=merged_vcf)
    print("##REF=<REF=N,Description=\"Ref's sequence in that region, default=N\">", file=merged_vcf)
    print("##QUAL=<QUAL=XXX,Description=\"The SV quality of the SV described in this region\">", file=merged_vcf)

    print("##ALT=<ID=SV,Description=\"Simple SVs\">", file=merged_vcf)
    print("##ALT=<ID=CSV,Description=\"Complex or nested SVs\">", file=merged_vcf)

    print("##FILTER=<ID=Covered,Description=\"Covered mean the SV is spanned by reads\">", file=merged_vcf)
    print("##FILTER=<ID=Uncovered,Description=\"UnCovered mean the SV is not spanned by reads\">", file=merged_vcf)
    print("##FILTER=<ID=Clustered,Description=\"Clustered mean the SV is not spanned by reads, but can be cluster together with others\">", file=merged_vcf)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV described in this region\">", file=merged_vcf)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=merged_vcf)
    print("##INFO=<ID=BKPS,Number=.,Type=String,Description=\"Possible SV breakpoints in this region\">", file=merged_vcf)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Possible SV types in this region\">", file=merged_vcf)
    print("##INFO=<ID=SUPPORT,Number=1,Type=String,Description=\"SV support number in this region\">", file=merged_vcf)
    print("##INFO=<ID=VAF,Number=1,Type=String,Description=\"SV allele frequency in this region\">", file=merged_vcf)
    print("##INFO=<ID=READS,Number=.,Type=String,Description=\"SV support read names in this region\">", file=merged_vcf)
    print("##INFO=<ID=MECHANISM,Number=1,Type=String,Description=\"DNA repair mechanism (MMBIR, NAHR, NHEJ, altEJ)\">", file=merged_vcf)
    print("##INFO=<ID=GraphID,Number=1,Type=String,Description=\"The corresponding graph id\">", file=merged_vcf)

    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=merged_vcf)
    print("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"high-quality reference reads\">", file=merged_vcf)
    print("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"high-quality variant reads\">", file=merged_vcf)

    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=merged_vcf)

    id_num = -1
    for vcf_file in files:
        if 'vcf' not in vcf_file:
            continue

        vcf_file = open(os.path.join(in_dir, vcf_file), 'r')

        previous_start = 0
        previous_end = 1
        sub_num = 1
        for record in vcf_file.readlines():

            line_split = str(record).split('\t')

            # set SV ID
            start = line_split[1]
            end = line_split[7].split(';')[0][4: ]
            if start == previous_start and end == previous_end:
                id_str = str(id_num) + '_' + str(sub_num)
                sub_num += 1
            else:
                previous_start = start
                previous_end = end
                id_num += 1
                sub_num = 1
                id_str = str(id_num)

            line_split[2] = id_str

            # set SV score
            old_score = float(line_split[5])

            new_score = int(100 - (round((old_score - min_score) / (max_score - min_score), 2) * 100))
            # print(new_score)
            line_split[5] = str(new_score)


            merged_vcf.write('\t'.join(line_split))

        vcf_file.close()

    merged_vcf.close()



def refine_type(original_type, original_bkps, options):
    """
    Refine original sv type by detailed bkps, primarily about:
    remove ins if ins_len is equal to dup_len; refind tdup to dup
    :param original_type:
    :param original_bkps:
    :param options:
    :return: new type
    """

    refined_type = []

    if 'INS' in original_type and 'tDUP' in original_type and 'DUP' not in original_type:
        ins_len = 0
        dup_len = 0
        for i in range(len(original_type)):
            if original_type[i] == 'INS':
                ins_len += int(original_bkps[i][2])
            elif original_type[i] == 'tDUP':
                dup_len += int(original_bkps[i][2])
        if ins_len - dup_len > options.min_sv_size:
            refined_type = original_type
        else:
            refined_type = [i for i in original_type if i != 'INS']
        return refined_type

    elif 'INS' in original_type and 'DUP' in original_type and 'tDUP' not in original_type:
        ins_len = 0
        dup_len = 0
        ins_pos = -1
        for i in range(len(original_type)):
            if original_type[i] == 'INS':
                ins_pos = int(original_bkps[i][0])
                ins_len += int(original_bkps[i][2])
            elif original_type[i] == 'DUP':
                dup_len += int(original_bkps[i][2])

                # refine dup and tdup by comparing with ins pos
                if ins_pos != -1:
                    dup_end = int(original_bkps[i][1])
                    if abs(ins_pos - dup_end) < 10:
                        original_type[i] = 'tDUP'

        if ins_len - dup_len > options.min_sv_size:
            refined_type = original_type
        else:
            refined_type = [i for i in original_type if i != 'INS']
        return refined_type

    elif 'INS' in original_type and 'DUP' in original_type and 'tDUP' in original_type:
        ins_len = 0
        dup_len = 0
        ins_pos = -1

        for i in range(len(original_type)):
            if original_type[i] == 'INS':
                ins_pos = int(original_bkps[i][0])
                ins_len += int(original_bkps[i][2])
            elif original_type[i] == 'DUP' or original_type[i] == 'tDUP':
                dup_len += int(original_bkps[i][2])

                # refine dup and tdup by comparing with ins pos
                if ins_pos != -1 and original_type[i] == 'DUP':
                    dup_end = int(original_bkps[i][1])
                    if abs(ins_pos - dup_end) < 10:
                        original_type[i] = 'tDUP'

        if ins_len - dup_len > options.min_sv_size:
            refined_type = original_type
        else:
            refined_type = [i for i in original_type if i != 'INS']
        return refined_type
    else:
        return original_type

def write_results_to_vcf(vcf_out, score_out, sv_stats, region, read_num_name_pair, sig_types, sig_score_pair, predict_scores, sig_mechanisms_pair, options):

    if len(sv_stats) > 0:


        avg_predict_score = (1 - round(np.mean(predict_scores), 2)) * 100
        all_support_reads = []
        all_mechanisms = []
        all_support_num = []
        all_vaf = []
        all_sv_types = []
        all_sv_bkps = []
        all_sig_scores = []

        region_split = region.split('+')
        chr = region_split[0]
        start = int(region_split[1])
        end = int(region_split[2])
        coverage = int(region_split[3])
        length = end - start

        for sv in sv_stats:
            sv_type = sv[0]
            sv_num = len(sv[1])
            sv_bkps = sv[2]
            if coverage == 0:
                vaf = 1.0
            else:
                vaf = round(len(sv[1]) / coverage, 2)

            # save this sv's info
            all_sv_types.append(sv_type)
            all_support_num.append(str(sv_num))
            all_vaf.append(str(vaf))
            all_sv_bkps.append(sv_bkps)

            cur_support_reads = []
            cur_sig_scores = []
            cur_mechanism = []

            for read_num in sv[1]:
                cur_support_reads.append(read_num_name_pair[read_num])
                cur_sig_scores.append(sig_score_pair[read_num])
                cur_mechanism.append(sig_mechanisms_pair[read_num])

            all_support_reads.append(cur_support_reads)
            all_sig_scores.append(cur_sig_scores)
            all_mechanisms.append(cur_mechanism)

        # calculate sig types
        sig_type_stat = collections.Counter(sig_types)
        if 'sigUncovered' in sig_type_stat.keys() and sig_type_stat['sigUncovered'] >= 0.75 * len(sig_types):
            filter_type = 'Uncovered'
        else:
            filter_type = 'Covered'



        for i in range(len(all_sv_types)):
            # create INFO
            svtype_info = "SVTYPE=" + all_sv_types[i]
            svsupp_info = "SUPPORT=" + all_support_num[i]
            svbkps_info = "BKPS="

            # Deepsv v1.1.6. Add
            if float(all_vaf[i]) > 1:
                all_vaf[i] = '1'
            svvaf_info = "VAF=" + all_vaf[i]

            svreads_info = "READS="
            mechanism_info = 'MECHANISM='

            svsig_scores = all_sig_scores[i]

            mechanisms = all_mechanisms[i]
            most_mechanism = Counter(mechanisms).most_common(1)[0][0]
            mechanism_info += most_mechanism


            for j in range(len(all_support_reads[i])):
                read_name = all_support_reads[i][j]
                if j == 0:
                    svreads_info += read_name
                else:
                    svreads_info += ',' + read_name

            split_types = all_sv_types[i].split('+')
            for j in range(len(all_sv_bkps[i])):
                if j == 0:
                    svbkps_info += split_types[j] + ':' + str(all_sv_bkps[i][j][2]) + '-' + str(all_sv_bkps[i][j][0]) + '-' + str(all_sv_bkps[i][j][1])
                else:
                    svbkps_info += ',' + split_types[j] + ':' + str(all_sv_bkps[i][j][2]) + '-' + str(all_sv_bkps[i][j][0]) + '-' + str(all_sv_bkps[i][j][1])


            # calculate sv score
            sv_score_std = np.std([int(score) for score in svsig_scores]) / int(all_support_num[i])
            sum_score = min(100, (sv_score_std + avg_predict_score))

            # classify SV to SV or CSV
            original_type = all_sv_types[i].split('+')

            # # svision v1.2.1 Add.
            # refine sv type
            refined_type = refine_type(original_type, all_sv_bkps[i], options)
            svtype_info = "SVTYPE=" + '+'.join(refined_type)
            # End Add

            # SVision v1.1.1 Modify.
            if len(refined_type) >= 2:
                new_type = '<CSV>'
            else:
                new_type = '<SV>'


            if options.report_mechanism == True:
                info = "END={0};SVLEN={1};{2};{3};{4};{5};{6};{7}".format(end, length, svtype_info, svsupp_info, svbkps_info, svvaf_info, svreads_info, mechanism_info)
            else:
                info = "END={0};SVLEN={1};{2};{3};{4};{5};{6}".format(end, length, svtype_info, svsupp_info, svbkps_info, svvaf_info, svreads_info)

            # SVision v1.1.6, ADD. report GT
            vaf = float(all_vaf[i])
            if vaf > 0.8:
                GT='1/1'
            elif vaf > 0.3:
                GT = '0/1'
            elif vaf >= 0:
                GT = '0/0'
            else:
                GT = './.'
            DV = int(all_support_num[i])
            DR = round(DV / (vaf + 0.001)) - DV

            gt_format = 'GT:DR:DV\t{0}:{1}:{2}'.format(GT, DR, DV)
            # Add, End

            line = chr + '\t' + str(start) + '\t' + '0' + '\t' + 'N' + '\t' + new_type + '\t' + str(sum_score) + '\t' + filter_type + '\t' + info + '\t' + gt_format


            print(sum_score, file=score_out)
            print(line, file=vcf_out)


def cal_scores_max_min(predict_path):
    files = os.listdir(predict_path)

    all_scores = []
    for file in files:
        if 'score.txt' in file:
            with open(os.path.join(predict_path, file)) as fin:
                for line in fin.readlines():
                    if line.strip() == '0':
                        continue
                    all_scores.append(float(line.strip()))



    return np.max(all_scores), np.min(all_scores)


def fetch_ref_seq(ref_path, chr, start, end):
    ref = pysam.FastaFile(ref_path)

    ref_cutted = ref.fetch(chr, start, end)
    return ref_cutted


# # the following: run repeatmasker and trf to annotate TE and VNTR, non use for current version, but in v2.0

# def analyze_record_mechanism(record, chrom, start, end, original_mechanism, out_dir, sub_process, rpmask_dir, trf_dir, out_split_path, options):
#     rpmask = options.rpmask
#     trf = options.trf
#
#     sv_region = chrom + '-' + str(start) + '-' + str(end)
#     fa_out_path = os.path.join(out_dir, '{0}.fa'.format(sv_region))
#     with open(fa_out_path, 'w') as fa_out:
#         ref_seq = fetch_ref_seq(options.genome, chrom, start, end)
#         fa_out.write('>' + str(sv_region) + '\n')
#         fa_out.write(ref_seq)
#
#     # run RepeatMaker
#     cmd_str = '{0} -parallel {1} -species human -gff -dir {2} {3} | grep nothing'.format(rpmask, sub_process, rpmask_dir,
#                                                                                 fa_out_path)
#     os.system(cmd_str)
#
#     # run TRF
#     cmd_str = '{0} {1} 2 7 7 80 10 50 500 -f -d -m | grep nothing '.format(trf, fa_out_path)
#     os.system(cmd_str)
#     cmd_str = 'mv {0}.fa.* {1}'.format(sv_region, trf_dir)
#     os.system(cmd_str)
#
#     # parse rp and trf mechanism
#     te_mechanism = []
#     rp_tbl_file = os.path.join(rpmask_dir, '{0}.fa.tbl'.format(sv_region))
#     if not os.path.exists(rp_tbl_file):
#         pass
#     else:
#         rp_mask_type = parse_rpmask(rp_tbl_file)
#         if rp_mask_type == -1:
#             pass
#         else:
#             te_mechanism = rp_mask_type
#
#     trf_file = os.path.join(trf_dir, '{0}.fa.2.7.7.80.10.50.500.1.txt.html'.format(sv_region))
#     if not os.path.exists(trf_file):
#         pass
#     else:
#         trf_type = parse_trf(trf_file)
#         if trf_type == -1:
#             pass
#         else:
#             te_mechanism.append(trf_type)
#
#     #
#     new_mechanism = ''
#     if len(te_mechanism) == 0:
#         new_mechanism = original_mechanism
#     else:
#         new_mechanism = ','.join(te_mechanism)
#
#     new_record = '='.join(str(record).strip().split('=')[0: -1]) + '=' + new_mechanism
#     if os.path.exists(out_split_path):
#         out_vcf = open(out_split_path)
#     else:
#         out_vcf = open(out_split_path, 'w')
#     out_vcf.write(new_record + '\n')
#     out_vcf.close()
#
# def analyze_mechanism(merged_vcf_path, out_dir, options):
#     """
#
#     :param merged_vcf_path:
#     :param out_dir:
#     :param options:
#     :return:
#     """
#
#     trf_dir = os.path.join(out_dir, 'trf')
#     if not os.path.exists(trf_dir):
#         os.mkdir(trf_dir)
#     trf_log = os.path.join(trf_dir, 'trf.log')
#
#     rpmask_dir = os.path.join(out_dir, 'rpmask')
#     if not os.path.exists(rpmask_dir):
#         os.mkdir(rpmask_dir)
#     rpmask_log = os.path.join(rpmask_dir, 'rpmask.log')
#
#     merged_vcf =  pysam.VariantFile(merged_vcf_path, 'r')
#     header = merged_vcf.header
#
#     out_vcf_path = os.path.join(options.out_path, "SVision.s{0}.mechanism.vcf".format(options.min_support))
#     out_vcf = open(out_vcf_path, 'w')
#     out_vcf.write(str(header))
#     out_vcf.close()
#
#     if options.thread_num > 2:
#         sub_process = 2
#     else:
#         sub_process = 1
#
#     process_num = max(1, int(options.thread_num / sub_process))
#     process_pool = multiprocessing.Pool(processes=process_num)
#     out_split_files = []
#     for i in range(process_num):
#         out_split_files.append(os.path.join(options.out_path, "SVision.mechanism.{0}.vcf".format(i)))
#
#     i = 0
#     for record in merged_vcf.fetch():
#         i += 1
#         out_split_path = out_split_files[i % process_num]
#
#         chrom = record.contig
#         start = record.pos
#         end = record.stop
#         original_mechanism = record.info['MECHANISM']
#
#         process_pool.apply_async(analyze_record_mechanism, (str(record),  chrom, start, end, original_mechanism, out_dir, sub_process, rpmask_dir, trf_dir, out_split_path, options))
#         # analyze_record_mechanism(record, out_dir, sub_process, rpmask_dir, trf_dir, out_split_file, options)
#     process_pool.close()
#     process_pool.join()
#
#     # merge split files, then delete split files
#     cmd_str = "cat "
#     for i in range(process_num):
#         cmd_str += os.path.join(options.out_path, "SVision.mechanism.{0}.vcf".format(i)) + ' '
#     cmd_str += '>> {0}'.format(out_vcf_path)
#     os.system(cmd_str)
#
#     # delete
#     for i in range(process_num):
#         os.remove(os.path.join(options.out_path, "SVision.mechanism.{0}.vcf".format(i)))

