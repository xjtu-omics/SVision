#!/usr/bin/env python3

import numpy as np
import cv2
import math

from src.segmentplot.classes import Segment
from src.segmentplot.plot_segment import PlotSingleImg


class BatchGenerator:
    def __init__(self, segments_file, horizontal_flip=False, shuffle=False,
                 mean=np.array([104., 117., 124.]), scale_size=(227, 227),
                 nb_classes=2, batch_size=128):

        # Init params
        self.horizontal_flip = horizontal_flip
        self.n_classes = nb_classes
        self.shuffle = shuffle
        self.mean = mean
        self.scale_size = scale_size
        self.pointer = 0
        self.batch_size = batch_size
        self.read_class_list(segments_file)

        if self.shuffle:
            self.shuffle_data()

    def read_class_list(self, segments_file):
        """
        Scan the image file and get the image paths and labels
        """

        with open(segments_file) as f:

            self.images = []
            self.labels = []

            lines = f.readlines()
            for line in lines:

                line_split = line.strip('\n').split('\t')

                # rubost
                data = '_'.join(line_split[1: 13])

                self.images.append(data)
                label = line_split[13] + 'svision' + line_split[0] + 'svision' + line_split[15] + 'svision' + line_split[16] + 'svision' + line_split[17] + 'svision' + line_split[18] + 'svision' + line_split[19] + 'svision' + line_split[20] + 'svision' + line_split[21] + 'svision' + line_split[22]
                self.labels.append(label)

            # store total number of data
            self.data_size = len(self.labels)

            complemented_num = self.batch_size * math.ceil(self.data_size / self.batch_size) - self.data_size
            complemented_data = '0_1_0_1_True_1_1_1_1_True_2_2'
            complemented_label = 'complement-complement'
            for i in range(complemented_num):
                self.images.append(complemented_data)
                self.labels.append(complemented_label)

            self.data_size = len(self.labels)


    def shuffle_data(self):
        """
        Random shuffle the images and labels
        """
        images = self.images.copy()
        labels = self.labels.copy()
        self.images = []
        self.labels = []

        # create list of permutated index and shuffle data accoding to list
        idx = np.random.permutation(len(labels))
        for i in idx:
            self.images.append(images[i])
            self.labels.append(labels[i])

    def reset_pointer(self):
        """
        reset pointer to begin of the list
        """
        self.pointer = 0

        if self.shuffle:
            self.shuffle_data()

    def next_batch(self, batch_size):
        """
        This function gets the next n ( = batch_size) images from the path list
        and labels and loads the images into them into memory
        """
        # Get next batch of image (path) and labels
        paths = self.images[self.pointer:self.pointer + batch_size]
        labels = self.labels[self.pointer:self.pointer + batch_size]

        # update pointer
        self.pointer += batch_size

        # Read images
        images = np.ndarray([batch_size, self.scale_size[0], self.scale_size[1], 3])

        for i in range(len(paths)):
            items = paths[i].split('_')
            # print(paths[i])
            seg1_x_start = int(items[0])
            seg1_x_end = int(items[1])
            seg1_y_start = int(items[2])
            seg1_y_end = int(items[3])

            if items[4] == 'True':
                seg1_forward = True
            elif items[4] == 'False':
                seg1_forward = False
            else:
                seg1_forward = None

            seg1 = Segment(seg1_x_start, seg1_y_start, (seg1_y_end - seg1_y_start), seg1_forward, 0)

            seg2_x_start = int(items[5])
            seg2_x_end = int(items[6])
            seg2_y_start = int(items[7])
            seg2_y_end = int(items[8])

            if items[9] == "True":
                seg2_forward = True
            elif items[9] == 'False':
                seg2_forward = False
            else:
                seg2_forward = None

            seg2 = Segment(seg2_x_start, seg2_y_start, (seg2_y_end - seg2_y_start), seg2_forward, 0)

            segments_pair = [seg1, seg2]

            read_len = int(items[10])
            ref_len = int(items[11])

            ploter = PlotSingleImg(segments_pair, read_len, ref_len, '', '', '', mode='return')
            img = ploter.plot()

            if self.horizontal_flip and np.random.random() < 0.5:
                img = cv2.flip(img, 1)

            # rescale image
            img = cv2.resize(img, (self.scale_size[0], self.scale_size[1]))
            img = img.astype(np.float32)

            # subtract mean
            img -= self.mean

            images[i] = img

        # return array of images and labels
        return images, labels
