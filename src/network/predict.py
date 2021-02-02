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

    def stats_sv_predicts(self, reads_dict):
        """
        stats sv predict results

        :param reads_dict: format: {read_id: [svtype, [bkp_start, bkp_end, bkp_len]}
                                # e.g {'1': {1: [906269, 906270, 200], 4: [905711, 906269, 200]}, '2': {1: [906281, 906282, 200]}}
        :return:
            sv_stats, format: [(type, support_read_id, avg_bkps), ()...]
                e.g. [('INS', ['1', '5'], [[13340320, 13340322]]), ('INS+DUP', ['2', '3', '4', '6',], [[13340323, 13340326], [13338229, 13338722]])]

        """

        stats = {}

        for read_id, sv_type_info in reads_dict.items():
            sv_type_str = ''.join(str(i) for i in sorted(sv_type_info.keys()))
            if sv_type_str not in stats.keys():

                new_bkps = []
                for i in range(len(sv_type_str)):
                    new_bkps.append(sv_type_info[int(sv_type_str[i])])
                stats[sv_type_str] = [[read_id], new_bkps]
            else:

                new_bkps = []
                for i in range(len(sv_type_str)):
                    new_bkps.append(sv_type_info[int(sv_type_str[i])])

                old_bkps = stats[sv_type_str][1]
                old_read_num = len(stats[sv_type_str][0])

                avg_bkps = []

                for i in range(len(new_bkps)):

                    # when meet a same svtype, calculate the avg bkps
                    start_avg = int((new_bkps[i][0] + old_bkps[i][0] * old_read_num) / (old_read_num + 1))
                    end_avg = int((new_bkps[i][1] + old_bkps[i][1] * old_read_num) / (old_read_num + 1))
                    len_avg = int((new_bkps[i][2] + old_bkps[i][2] * old_read_num) / (old_read_num + 1))

                    avg_bkps.append([start_avg, end_avg, len_avg])

                    # print('-', new_bkps[i][0], old_bkps[i][0], old_read_num, start_avg)
                    # print('-', new_bkps[i][1], old_bkps[i][1], old_read_num, end_avg)

                # print(old_read_num, new_bkps, old_bkps, avg_bkps)

                stats[sv_type_str][0].append(read_id)
                stats[sv_type_str][1] = avg_bkps

        # print(stats.items())
        stats = sorted(stats.items(), key=lambda x: len(x[1][0]), reverse=True)
        # print(stats)
        # exit()
        sv_stats = []

        # convert type id to type string
        for j in range(len(stats)):
            sv_type_info, sv_info = stats[j]
            support_read_ids = sv_info[0]
            bkps = sv_info[1]

            sv_type_str = ""
            for i in range(len(sv_type_info)):
                if i is not 0:
                    sv_type_str += "+"

                if sv_type_info[i] == '0':
                    sv_type_str += "DEL"
                elif sv_type_info[i] == '1':
                    sv_type_str += "INS"
                elif sv_type_info[i] == '2':
                    sv_type_str += "INV"
                elif sv_type_info[i] == '3':
                    sv_type_str += "DUP"
                elif sv_type_info[i] == '4':
                    sv_type_str += "tDUP"
            sv_stats.append((sv_type_str, support_read_ids, bkps))

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
                # exit()
                score_value, predict_value, softmax_value = sess.run([score, tf.argmax(score, 1), tf.nn.softmax(score)],
                                                                     feed_dict={x: batch_px, keep_prob: self.dropout_rate})


                for i in range(len(batch_py)):
                    if 'complement' in batch_py[i]:
                        continue

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


                    # write region to file
                    if region != last_region:
                        if last_region != "":

                            sv_stats = self.stats_sv_predicts(reads_dict)
                            write_results_to_vcf(out_vcf_file, out_score_file, sv_stats, last_region, read_num_name_pair, sig_types, sig_score_pair, predict_scores, sig_mechanisms_pair, options)

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
            sv_stats = self.stats_sv_predicts(reads_dict)
            write_results_to_vcf(out_vcf_file, out_score_file, sv_stats, last_region, read_num_name_pair, sig_types, sig_score_pair, predict_scores, sig_mechanisms_pair, options)

        out_vcf_file.close()
