#!/usr/bin/env python3

# encoding: utf-8

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
# from sklearn.cluster import DBSCAN
from src.collection.classes import Cluster, read_fai_file


def partition_and_cluster(signatures, fai_file, genome, chr, this_sample_path, options):

    partitions = signature_partition(signatures, options.patition_max_distance)

    clusters = cluster_partitions(partitions, fai_file, genome, this_sample_path, options)
    return clusters


def cluster_positions_simple(positions, max_delta):
    """Form partitions of (contig, position) pairs using simple distance."""
    sorted_positions = sorted(positions, key=lambda position: (position[0], position[1]))
    partitions = []
    current_partition = []
    for position in sorted_positions:
        if len(current_partition) < 1:
            current_partition.append(position)
            continue
        if distance_positions(current_partition[0], position) > max_delta:
            partitions.append(current_partition[:])
            while len(current_partition) > 0 and distance_positions(current_partition[0], position) > max_delta:
                current_partition.pop(0)
        current_partition.append(position)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions

def distance_positions(position1, position2):
    return float("inf") if position1[0] != position2[0] else abs(position1[1] - position2[1])

def signature_partition(signatures, max_distance):

    sorted_signatures = sorted(signatures, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for signature in sorted_signatures:
        if len(current_partition) > 0 and current_partition[-1].position_distance_to(signature) > max_distance:
            partitions.append(current_partition[:])
            current_partition = []
        current_partition.append(signature)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])

    return partitions

def cluster_partitions(partitions, fai_file, genome, this_sample_path, options):
    '''
    Using hierarchical clustering to cluster signatures in each partition
    :param partitions:
    :param options:
    :return:
    '''
    fai_dict = read_fai_file(fai_file)

    clusters = []
    for partition in partitions:
        if len(partition) == 1:
            this_cluster = Cluster(partition, this_sample_path)

            # # SVision v1.0.2. ADD. remove cluster with wrong cord
            if this_cluster.abandon == 0:
                clusters.append(this_cluster)
            continue

        data = np.array([[signature.get_source()[1], signature.get_source()[2], 1000] for signature in partition])

        # DBSCAN
        # db = DBSCAN(eps=0.3, min_samples=1, metric=span_position_distance).fit(data)
        # cluster_indices = db.labels_
        # new_clusters = [[] for i in range(max(cluster_indices) + 1)]
        #
        # for signature_index, cluster_index in enumerate(cluster_indices):
        #     if cluster_index == -1:
        #         continue
        #     cur_sig = partition[signature_index]
        #     new_clusters[cluster_index].append(cur_sig)
        #
        # for cluster_list in new_clusters:
        #     this_cluster = Cluster(cluster_list, this_sample_path)
        #     this_cluster.annotate_cluster(cluster_list, genome, fai_dict)
        #     clusters.append(this_cluster)

        # # he cluster
        Z = linkage(data, method="average", metric=span_position_distance)
        cluster_indices = list(fcluster(Z, options.cluster_max_distance, criterion='distance'))

        new_clusters = [[] for i in range(max(cluster_indices))]

        for signature_index, cluster_index in enumerate(cluster_indices):
            cur_sig = partition[signature_index]
            new_clusters[cluster_index-1].append(cur_sig)

        for cluster_list in new_clusters:
            this_cluster = Cluster(cluster_list, this_sample_path)

            # # SVision v1.0.2. ADD. remove cluster with wrong cord
            if this_cluster.abandon == 0:
                clusters.append(this_cluster)

    return sorted(clusters, key=lambda cluster:(cluster.contig, (cluster.cstart + cluster.cend) / 2))

def span_position_distance(signature1, signature2):
    distance_normalizer = signature1[2]
    span1 = signature1[1] - signature1[0]
    span2 = signature2[1] - signature2[0]
    center1 = (signature1[0] + signature1[1]) // 2
    center2 = (signature2[0] + signature2[1]) // 2
    position_distance = min(abs(signature1[0] - signature2[0]), abs(signature1[1] - signature2[1]), abs(center1 - center2)) / distance_normalizer
    span_distance = abs(span1 - span2) / max(span1, span2)

    return position_distance + span_distance





