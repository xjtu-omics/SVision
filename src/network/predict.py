#!/usr/bin/env python3


import numpy as np
import tensorflow as tf
from src.network.alexnet import AlexNet
from src.network.create_batch import BatchGenerator
from src.network.output import write_results_to_vcf



class Predict:
    def __init__(self, chrom, segments_out_file):
        # Path to the textfiles for the trainings and validation set
        self.segments_out_file = segments_out_file
        self.chrom = chrom
        # Learning

        # Network params
        self.dropout_rate = 1.
        self.num_classes = 5
        self.train_layers = ['fc8', 'fc7']

        # How often we want to write the tf.summary data to disk
        self.display_step = 1

    def get_region_potential_svtypes(self, reads_dict):
        """
        stats sv types from predict results

        :param reads_dict: format (allow same predicted types): {read_id: [(svtype, (bkp_start, bkp_end, bkp_len))]}
                                # e.g {'1': [(1, (859987, 859988, 2389)), (4, (857598, 859987, 2389))], '2': [(1, (859987, 859988, 2388)), (4, (857623, 859987, 2364))]}
                         format (not allow same predicted types): {read_id: [svtype, [bkp_start, bkp_end, bkp_len]}
                                # e.g {'1': {1: [906269, 906270, 200], 4: [905711, 906269, 200]}, '2': {1: [906281, 906282, 200]}}
        :return:
            sv_stats, format: [(type, support_read_id, avg_bkps), ()...]
                e.g. [('INS', ['1', '5'], [[13340320, 13340322]]), ('INS+DUP', ['2', '3', '4', '6',], [[13340323, 13340326], [13338229, 13338722]])]

        """

        # # allow same predicted types
        # stats = {}
        #
        # for read_id, sv_type_infos in reads_dict.items():
        #     sv_type_infos_srt_by_svtype = sorted(sv_type_infos, key=lambda x: x[0])
        #     read_final_svtype = ''.join([str(sv_info[0]) for sv_info in sv_type_infos_srt_by_svtype])
        #
        #     # meet a new final_svtype, collect each subtype's bkps and push to stats
        #     if read_final_svtype not in stats.keys():
        #
        #         new_bkps = []
        #         for i in range(len(read_final_svtype)):
        #             new_bkps.append(sv_type_infos_srt_by_svtype[i][1])
        #         stats[read_final_svtype] = [[read_id], new_bkps]
        #
        #     # meet a already_in final_svtype, calculate avg bkps and
        #     else:
        #
        #         new_bkps = []
        #         for i in range(len(read_final_svtype)):
        #             new_bkps.append(sv_type_infos_srt_by_svtype[i][1])
        #
        #         old_bkps = stats[read_final_svtype][1]
        #         old_read_num = len(stats[read_final_svtype][0])
        #
        #         # # calculate avg bkps
        #         avg_bkps = []
        #
        #         for i in range(len(new_bkps)):
        #
        #             # when meet a same svtype, calculate the avg bkps
        #             start_avg = int((new_bkps[i][0] + old_bkps[i][0] * old_read_num) / (old_read_num + 1))
        #             end_avg = int((new_bkps[i][1] + old_bkps[i][1] * old_read_num) / (old_read_num + 1))
        #             len_avg = int((new_bkps[i][2] + old_bkps[i][2] * old_read_num) / (old_read_num + 1))
        #
        #             avg_bkps.append([start_avg, end_avg, len_avg])
        #
        #
        #         stats[read_final_svtype][0].append(read_id)
        #         stats[read_final_svtype][1] = avg_bkps

        # # not allow same predicted types
        stats = {}

        for read_id, sv_type_infos in reads_dict.items():
            sv_type_str = ''.join(str(i) for i in sorted(sv_type_infos.keys()))
            if sv_type_str not in stats.keys():

                new_bkps = []
                for i in range(len(sv_type_str)):
                    new_bkps.append(sv_type_infos[int(sv_type_str[i])])
                stats[sv_type_str] = [[read_id], new_bkps]
            else:

                new_bkps = []
                for i in range(len(sv_type_str)):
                    new_bkps.append(sv_type_infos[int(sv_type_str[i])])

                old_bkps = stats[sv_type_str][1]
                old_read_num = len(stats[sv_type_str][0])

                avg_bkps = []

                for i in range(len(new_bkps)):

                    # when meet a same svtype, calculate the avg bkps
                    start_avg = int((new_bkps[i][0] + old_bkps[i][0] * old_read_num) / (old_read_num + 1))
                    end_avg = int((new_bkps[i][1] + old_bkps[i][1] * old_read_num) / (old_read_num + 1))
                    len_avg = int((new_bkps[i][2] + old_bkps[i][2] * old_read_num) / (old_read_num + 1))

                    avg_bkps.append([start_avg, end_avg, len_avg])

                stats[sv_type_str][0].append(read_id)
                stats[sv_type_str][1] = avg_bkps


        stats = sorted(stats.items(), key=lambda x: len(x[1][0]), reverse=True)
        sv_stats = []

        # convert type id to type string: 123 -> DEL+INS+INV
        for j in range(len(stats)):
            sv_type_infos, sv_info = stats[j]
            support_read_ids = sv_info[0]
            bkps = sv_info[1]

            read_final_svtype = ""
            for i in range(len(sv_type_infos)):
                if i is not 0:
                    read_final_svtype += "+"

                if sv_type_infos[i] == '0':
                    read_final_svtype += "DEL"
                elif sv_type_infos[i] == '1':
                    read_final_svtype += "INS"
                elif sv_type_infos[i] == '2':
                    read_final_svtype += "INV"
                elif sv_type_infos[i] == '3':
                    read_final_svtype += "DUP"
                elif sv_type_infos[i] == '4':
                    read_final_svtype += "tDUP"
            sv_stats.append((read_final_svtype, support_read_ids, bkps))

        return sv_stats


    def run(self, out_path_prefix, options):
        """
        Begin to predict types
        :param out_path_prefix:
        :param options:
        :return:
        """
        tf.compat.v1.disable_eager_execution()

        out_score_file = open(out_path_prefix + '.score.txt', 'w')
        out_vcf_file = open(out_path_prefix + '.vcf', 'w')

        batch_size = options.batch_size

        # # convert original data to batch format data
        batch_generator = BatchGenerator(self.segments_out_file, shuffle=False, nb_classes=self.num_classes, batch_size=batch_size)

        tf.compat.v1.reset_default_graph()

        x = tf.compat.v1.placeholder(tf.float32, [batch_size, 227, 227, 3])

        keep_prob = tf.compat.v1.placeholder(tf.float32)

        # Initialize model_ft
        model = AlexNet(x, keep_prob, self.num_classes, self.train_layers)

        # how many batches do we have
        val_batches_per_epoch = np.floor(batch_generator.data_size / batch_size).astype(np.int16)

        # Link variable to model_ft output
        score = model.fc8
        with tf.compat.v1.Session() as sess:
            # Initialize all variables
            model_path = options.model_path
            sess.run(tf.compat.v1.global_variables_initializer())
            saver = tf.compat.v1.train.Saver()
            saver.restore(sess, (model_path))

            # new_saver = tf.train.import_meta_graph(model_path)
            # new_saver.restore(sess, model_path)

            sess.run(tf.compat.v1.local_variables_initializer())

            reads_dict = {}     # save predict results
                                # format: {read_id: [svtype, [bkp_start, bkp_end]}
                                # e.g {'1': {1: [906269, 906270], 4: [905711, 906269]}, '2': {1: [906281, 906282]}}

            read_num_name_pair = {}     # save read_name:read_id pair
            sig_score_pair = {}         # save read_name:sig score pair
            sig_mechanisms_pair = {}         # save read_name:mechanism pair

            last_region = ""
            region = ""

            sig_types = []
            predict_scores = []
            print('[Processing]: Predicting ' + self.chrom)

            for _ in range(val_batches_per_epoch):
                batch_px, batch_py = batch_generator.next_batch(batch_size)
                # # predict
                score_value, predict_value, softmax_value = sess.run([score, tf.argmax(score, 1), tf.nn.softmax(score)],
                                                                     feed_dict={x: batch_px, keep_prob: self.dropout_rate})

                # # parse predict results
                for i in range(len(batch_py)):
                    if 'complement' in batch_py[i]:
                        continue

                    # parse batch info
                    split_item = batch_py[i].split('svision')
                    read_num = split_item[0]
                    region = split_item[1]
                    read_name = split_item[2]

                    bkp_start = int(split_item[4])
                    bkp_end = int(split_item[5])
                    bkp_len = int(split_item[9])
                    mechanism = split_item[8]

                    # # SVision v1.0.1. ADD, Fix bug: remove wrong 'INV' prediction
                    forward = split_item[7]
                    if forward == 'True' and predict_value[i] == 2:
                        continue
                    # # End ADD

                    # meet a new region, then write region to file
                    if region != last_region:
                        if last_region != "":
                            region_potential_svtypes = self.get_region_potential_svtypes(reads_dict)
                            write_results_to_vcf(out_vcf_file, out_score_file, region_potential_svtypes, last_region, read_num_name_pair, sig_types, sig_score_pair, predict_scores, sig_mechanisms_pair, options)

                        last_region = region
                        reads_dict = {}
                        read_num_name_pair = {}
                        sig_score_pair = {}
                        sig_types = []
                        predict_scores = []
                        sig_mechanisms_pair = {}
                    read_num_name_pair[read_num.replace('m', '')] = read_name
                    sig_types.append(split_item[3])

                    predict_scores.append(round(softmax_value[i][predict_value[i]], 2))
                    sig_score_pair[read_num.replace('m', '')] = split_item[6]

                    sig_mechanisms_pair[read_num.replace('m', '')] = mechanism

                    # only main segs can be predicted to ins or del
                    # # allow same predicted types
                    # if "m" not in read_num:
                    #     if predict_value[i] == 0 or predict_value[i] == 1:
                    #         continue
                    #     else:
                    #         # if meet a new read
                    #         if read_num not in reads_dict.keys():
                    #             reads_dict[read_num] = [[predict_value[i], [bkp_start, bkp_end, bkp_len]]]
                    #         else:
                    #             reads_dict[read_num].append([predict_value[i], [bkp_start, bkp_end, bkp_len]])
                    #
                    # else:
                    #     read_num = read_num.replace('m', '')
                    #     # read_num = batch_py[i]
                    #     if read_num not in reads_dict.keys():
                    #         reads_dict[read_num] = [[predict_value[i], [bkp_start, bkp_end, bkp_len]]]
                    #     else:
                    #         reads_dict[read_num].append([predict_value[i], [bkp_start, bkp_end, bkp_len]])

                    # # not allow same predicted types
                    if "m" not in read_num:
                        if predict_value[i] == 0 or predict_value[i] == 1:
                            continue
                        else:
                            # if meet a new read
                            if read_num not in reads_dict.keys():
                                reads_dict[read_num] = {predict_value[i]: [bkp_start, bkp_end, bkp_len]}
                            else:
                                reads_dict[read_num][predict_value[i]] = [bkp_start, bkp_end, bkp_len]

                    else:
                        read_num = read_num.replace('m', '')
                        # read_num = batch_py[i]
                        if read_num not in reads_dict.keys():
                            reads_dict[read_num] = {predict_value[i]: [bkp_start, bkp_end, bkp_len]}
                        else:
                            reads_dict[read_num][predict_value[i]] = [bkp_start, bkp_end, bkp_len]


                # break
            # meet the end
            region_potential_svtypes = self.get_region_potential_svtypes(reads_dict)
            write_results_to_vcf(out_vcf_file, out_score_file, region_potential_svtypes, last_region, read_num_name_pair, sig_types, sig_score_pair, predict_scores, sig_mechanisms_pair, options)

        out_vcf_file.close()
