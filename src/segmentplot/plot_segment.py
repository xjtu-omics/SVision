#!/usr/bin/env python3

import numpy as np
import cv2 as cv
import os


class PlotSingleImg():
    def __init__(self, segments_orignal, refLength, readLength,
                 type, outDir, cnt, mode="save"):

        self.ratio = (float(max(readLength, refLength) / 227.0))
        # print(self.ratio)
        if self.ratio < 1:
            self.ratio = 1

        # self.ratio = 1

        self.segments_orignal = segments_orignal

        self.readLength = int(readLength / self.ratio)
        self.refLength = int(refLength / self.ratio)

        # self.region = region
        self.type = type
        self.outDir = outDir
        self.cnt = cnt

        self.line = 1
        self.mode = mode


    def plot(self):
        # plt.subplots(figsize=(2.24, 2.24))
        # plt.subplot(111)
        # fig, ax = plt.subplots()

        img_length = 227
        img = np.zeros((img_length, img_length, 3))

        first_channel = np.zeros((img_length, img_length))
        third_channel = np.zeros((img_length, img_length))
        for seg in self.segments_orignal:

            if seg.forward():
                cv.line(first_channel, (int(seg.yStart() / self.ratio), int(seg.xStart() / self.ratio)),
                        (int(seg.yEnd() / self.ratio), int(seg.xEnd() / self.ratio)), 255, self.line)
            else:
                cv.line(first_channel, (int(seg.yEnd() / self.ratio), int(seg.xEnd() / self.ratio)),
                        (int(seg.yStart() / self.ratio), int(seg.xStart() / self.ratio)), 255, self.line)
                cv.line(third_channel, (int(seg.yEnd() / self.ratio), int(seg.xEnd() / self.ratio)),
                        (int(seg.yStart() / self.ratio), int(seg.xStart() / self.ratio)), 255, self.line)


        img[:, :, 0] = first_channel

        second_channel = np.zeros((img_length, img_length))

        for i in range(img_length):

            dot_pos = np.where(first_channel[:, i] != 0)[0]
            dot_num = len(dot_pos)
            if dot_num >= 2:

                second_channel[dot_pos, i] = 255

        img[:, :, 1] = second_channel
        img[:, :, 2] = third_channel

        if self.mode == 'save':
            cv.imwrite(os.path.join(self.outDir, self.type + "-" + str(self.cnt) + ".png"), img[:, :, :])
        else:
            return img[:, :, :]



class PlotSingleImg2():
    def __init__(self, segments_orignal, name, outDir, mode="save"):

        ref_length = max([align['ref_end'] for align in segments_orignal]) - min([align['ref_start'] for align in segments_orignal])
        read_length = max([align['q_end'] for align in segments_orignal]) - min([align['q_start'] for align in segments_orignal])

        if ref_length < read_length:
            if ref_length < 1000:
                self.ratio = 1
            else:
                self.ratio = 10
                while ref_length / self.ratio > 1000:
                    self.ratio = self.ratio * 10
        else:
            if read_length < 1000:
                self.ratio = 1
            else:
                self.ratio = 10
                while read_length / self.ratio > 1000:
                    self.ratio = self.ratio * 10
        # print(self.ratio)
        if self.ratio < 1:
            self.ratio = 1

        # self.ratio = 1

        self.segments_orignal = segments_orignal

        self.read_length = int(read_length / self.ratio)
        self.ref_length = int(ref_length / self.ratio)

        # print(self.read_length, self.ref_length, self.ratio)

        # self.region = region
        self.name = name
        self.outDir = outDir

        self.line = 1
        self.mode = mode


    def plot(self):

        ref_left_most = min([align['ref_start'] for align in self.segments_orignal])
        for align in self.segments_orignal:
            align['ref_start'] -= ref_left_most
            align['ref_end'] -= ref_left_most




        img = np.ones((self.read_length, self.ref_length, 1)) * 255
        for seg in self.segments_orignal:
            if not seg['is_reverse']:
                cv.line(img, (int(seg['ref_start'] / self.ratio), int(seg['q_start'] / self.ratio)),
                        (int(seg['ref_end'] / self.ratio), int(seg['q_end'] / self.ratio)), 0, self.line)
            else:
                cv.line(img, (int(seg['ref_end'] / self.ratio), int(seg['q_end'] / self.ratio)),
                        (int(seg['ref_start'] / self.ratio), int(seg['q_start'] / self.ratio)), 0, self.line)
        cv.imwrite(os.path.join(self.outDir, str(self.name) + ".png"), img[:, :, :])

        #



        # img = np.zeros((self.read_length, self.ref_length, 3))
        #
        # first_channel = np.zeros((self.read_length, self.ref_length))
        # third_channel = np.zeros((self.read_length, self.ref_length))
        # for seg in self.segments_orignal:
        #
        #     if not seg['is_reverse']:
        #         cv.line(first_channel, (int(seg['ref_start'] / self.ratio), int(seg['q_start'] / self.ratio)),
        #                 (int(seg['ref_end']  / self.ratio), int(seg['q_end'] / self.ratio)), 255, self.line)
        #     else:
        #         cv.line(first_channel, (int(seg['ref_end'] / self.ratio), int(seg['q_end'] / self.ratio)),
        #                 (int(seg['ref_start'] / self.ratio), int(seg['q_start'] / self.ratio)), 255, self.line)
        #         cv.line(third_channel, (int(seg['ref_end'] / self.ratio), int(seg['q_end'] / self.ratio)),
        #                 (int(seg['ref_start'] / self.ratio), int(seg['q_start'] / self.ratio)), 255, self.line)
        #
        #
        # img[:, :, 0] = first_channel
        #
        # second_channel = np.zeros((self.read_length, self.ref_length))
        #
        # for i in range(self.ref_length):
        #
        #     dot_pos = np.where(first_channel[:, i] != 0)[0]
        #     dot_num = len(dot_pos)
        #     if dot_num >= 2:
        #
        #         second_channel[dot_pos, i] = 255
        #
        # img[:, :, 1] = second_channel
        # img[:, :, 2] = third_channel
        #
        #
        # cv.imwrite(os.path.join(self.outDir, str(self.name) + ".png"), img[:, :, :])





