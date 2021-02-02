#!/usr/bin/env python3

import os
from src.segmentplot.run_hash_lineplot import cord_to_segments
from src.collection.graph import graph_is_same_as
def linearOrNot(i, j):

    distance_on_ref = j.yStart() - i.yEnd()
    distance_on_read = j.xStart() - i.xEnd()
    # diff = abs(distance_on_read - distance_on_ref)
    if distance_on_read == 0:
        distance_on_read = 1
    diff = distance_on_ref / distance_on_read

    if i.forward() != j.forward():
        return False
    # if diff >= 50:
    if diff >= 1.5 or diff <= 0.7:
        return False
    else:
        return True


def write_graph_to_file(file_path, graphs):

    out_graphs = [graphs[0]]

    for i in range(1, len(graphs)):
        graph = graphs[i]

        same_flags = []
        for target_graph in out_graphs:
            if graph_is_same_as(graph, target_graph, strict=True) == False:
                same_flags.append(0)
            else:
                same_flags.append(1)

        if 1 not in same_flags:
            out_graphs.append(graph)
        else:
            for i in range(len(same_flags)):
                if same_flags[i] == 1:
                    out_graphs[i].appear_time += 1

    for graph in sorted(out_graphs, key=lambda g: g.appear_time, reverse=True):

        nodes = graph.nodes
        edges = graph.edges

        # #in case: for a same region, there are more than one SV
        graph_file_path = file_path
        cnt = 1
        while True:
            if os.path.exists(graph_file_path):
                graph_file_path = '.'.join(file_path.split('.')[0: -1]) + "_{0}.gfa".format(cnt)
                cnt += 1
            else:
                break

        # # output
        with open(graph_file_path, 'w') as fout:
            for node in nodes:
                node_seq = node.seq if node.seq is not "" else 'N'
                if 'I' in node.id:
                    if node.node_is_dup == True:
                        fout.write("S\t{0}\t{1}\tSN:Z:{2}\tSO:i:{3}\tSR:i:0\tLN:i:{4}\tDP:S:{5}:{6}\n".format(node.id, node_seq, node.host, node.read_start, len(node_seq), node.dup_from, node.ref_start))
                    else:
                        fout.write("S\t{0}\t{1}\tSN:Z:{2}\tSO:i:{3}\tSR:i:0\tLN:i:{4}\n".format(node.id, node_seq, node.host, node.read_start, len(node_seq)))
                else:
                    fout.write("S\t{0}\t{1}\tSN:Z:{2}\tSO:i:{3}\tSR:i:0\tLN:i:{4}\n".format(node.id, node_seq, node.host, node.ref_start, len(node_seq)))

            for edge in edges:
                fout.write("L\t{0}\t{1}\t{2}\t{3}\t0M\tSR:i:0\n".format(edge.node1, '-' if edge.node1_is_reverse else '+', edge.node2, '-' if edge.node2_is_reverse else '+'))

def writer_cluster_to_file(clusters, chr, part_num, options):
    """
    Traverse each cluster, collect it's segments's info, and finally write them to file
    :param clusters: all cluster form he-cluster
    :param chr: current chromosome
    :param part_num: current part num on this chromosome
    :param options: parameters
    :return: No return
    """

    cluster_num = 0

    all_segs = []   # save segments's info in a list

    # res = []      # DEBUG code1

    for cluster in clusters:

        if int(cluster.cend) - int(cluster.cstart) > options.max_sv_size:
            continue

        if cluster.read_num >= options.min_support:
            cluster_num += 1
            cluster, segs_of_cluster = proc_one_cluster(cluster, options)

            # res.append(cluster)       # DEBUG code1

            # # Add. Deepsv v1.2. write graph to file
            if options.report_graph is True:
                graph_out_path = os.path.join(options.out_path, 'graphs')

                graphs = [sig.graph for sig in cluster.get_signatures()]
                cluster_graph_file = os.path.join(graph_out_path, "{0}-{1}-{2}.gfa".format(cluster.contig, int(cluster.cstart), int(cluster.cend)))
                write_graph_to_file(cluster_graph_file, graphs)
            # # End Add

            all_segs.extend(segs_of_cluster)


    # # DEBUG code1 : write cluster info to file
    # clusters_out_path = os.path.join(options.out_path, 'clusters')
    # cluster_info_file = os.path.join(clusters_out_path, chr + '.clusters.' + str(part_num) + '.bed')
    # cluster_writer = open(cluster_info_file, 'w')
    # for cluster in res:
    #     cluster_writer.write(cluster.to_string() + '\n')
    # cluster_writer.close()
    # #


    # write segments' info to file
    segments_out_path = os.path.join(options.out_path, 'segments')
    chr_segments_out_file = os.path.join(os.path.join(segments_out_path, chr + '.segments.' + str(part_num) + '.bed'))
    chr_segments_out_file = open(chr_segments_out_file, 'w')
    for seg in all_segs:
        chr_segments_out_file.write(seg)
    chr_segments_out_file.close()


def proc_one_cluster(cluster, options):
    """
    for each cluster, process sigs that contained in that cluster, then process segments in each sig
    :param cluster: input cluster which contains sig's info
    :param options: parameters
    :return:
        cluster: cluster's info
        segs_of_cluster: all segments's info in that cluster
    """

    signatures = cluster.get_signatures()

    cluster_region = cluster.contig + '+' + str(int(cluster.cstart)) + '+' + str(int(cluster.cend)) \
                     + '+' + str(cluster.coverage)

    # traverse sigs
    segs_of_cluster = []
    sig_cnt = 0
    for sig in signatures:
        sig_cnt += 1

        # collect segment's info within this sig
        segs_of_sig = proc_one_sig(cluster_region, sig, sig_cnt, options)

        # pass this sig if it has only one major segment
        if segs_of_sig == -1:
            continue

        segs_of_cluster.extend(segs_of_sig)

    return cluster, segs_of_cluster

def proc_one_sig(cluster_region, sig, sig_cnt, options):
    """
    Process each sig, split it's segment-plot
    :param cluster_region: cluster's info which this sig belongs to
    :param sig: sig need to process
    :param sig_cnt: number id of this sig
    :param options: parameters
    :return:
        segs_of_sig: all segments's info in that sig
    """
    refLength, readLenght, main_segs_cord, other_segs_cord = sig.get_segs_cords()
    # pass this sig if it has only one major segment
    if main_segs_cord == -1:
        return -1

    main_segs = cord_to_segments(main_segs_cord)
    other_segs = cord_to_segments(other_segs_cord)

    qname = sig.qname
    sig_type = sig.type
    sig_bkps = sig.bkps

    sig_mechanism = sig.mechanism

    segs_of_sig = []
    all_segs = main_segs + other_segs

    """
    # # DEBUG code2: output dotplots if required
    from src.segmentplot.plot_segment import PlotSingleImg
    dotplots_out_path = os.path.join(options.out_path, 'dotplots')
    if not os.path.exists(dotplots_out_path):
        os.mkdir(dotplots_out_path)
    cluster_dotplot_out_path = os.path.join(dotplots_out_path, cluster_region)
    if not os.path.exists(cluster_dotplot_out_path):
        os.mkdir(cluster_dotplot_out_path)
    # plot
    ploter = PlotSingleImg(all_segs, refLength, readLenght, str(sig_cnt), cluster_dotplot_out_path, str('0-' + str(qname[: min(100, len(qname))]).replace('/', '_')))
    ploter.plot()
    """

    #
    # calculate non-linear score for each dotplot
    non_linear_score = cal_non_linear(all_segs)

    sub_cnt = 0

    # conbining major segments
    for i in range(len(main_segs) - 1):
        sub_cnt += 1
        if not linearOrNot(main_segs[i], main_segs[i + 1]):
            segments_pair = [main_segs[i], main_segs[i + 1]]
            segs_of_sig.append(cluster_region + '\t' + segments_pair[0].toString() + '\t' + segments_pair[1].toString()
                                        + '\t' + str(readLenght) + '\t' + str(refLength) + '\t' + str(sig_cnt) + "m"
                               +'\t' + str(sub_cnt) + '\t' + qname + '\t' + sig_type + '\t' + str(sig_bkps[0][0]) + '\t' + str(sig_bkps[0][1]) + '\t' + str(non_linear_score) + '\t' + 'True' + '\t' + sig_mechanism + '\t' + str(sig_bkps[0][2]) + '\n')

            # # # DEBUG code2: output dotplots if required
            # ploter = PlotSingleImg(segments_pair, refLength, readLenght, str(sig_cnt), cluster_dotplot_out_path,
            #                        str(str(sub_cnt) + '-' + str(qname[: min(100, len(qname))]).replace('/', '_')))
            # ploter.plot()
    # combining major segment with one other segment
    for seg1 in main_segs:
        for i in range(len(other_segs)):
            seg2 = other_segs[i]

            sub_cnt += 1
            if seg2 in main_segs and seg1 in main_segs:
                pass
            else:
                if not linearOrNot(seg1, seg2):
                    segments_pair = [seg1, seg2]

                    # # SVision v1.0.1. ADD, Fix bug: removewrong 'INV' prediction
                    if seg1.forward() == False or seg2.forward() == False:
                        forward = 'False'
                    else:
                        forward = 'True'
                    # # End ADD

                    segs_of_sig.append(cluster_region + '\t' + segments_pair[0].toString() + '\t' + segments_pair[1].toString()
                        + '\t' + str(readLenght) + '\t' + str(refLength) + '\t' + str(sig_cnt) + '\t' + str(
                            sub_cnt) + '\t' + qname + '\t' + sig_type + '\t' + str(sig_bkps[i + 1][0]) + '\t' + str(sig_bkps[i + 1][1]) + '\t' + str(non_linear_score) + '\t' + forward + '\t' + sig_mechanism + '\t' + str(sig_bkps[i + 1][2]) + '\n')

                    # # DEBUG code2: output dotplots if required
                    # ploter = PlotSingleImg(segments_pair, refLength, readLenght, str(sig_cnt), cluster_dotplot_out_path,
                    #                        str(str(sub_cnt) + '-' + str(qname[: min(100, len(qname))]).replace('/',
                    #                                                                                            '_')))
                    # ploter.plot()
    return segs_of_sig

def cal_non_linear(all_segs):
    """
    Calculate non-linear score
    :param all_segs: all seg's info
    :return:
        non-linear score
    """
    ref_cords = []

    sumup_non_linear = 0
    for seg in all_segs:
        ref_cords.append(seg.yStart())
        ref_cords.append(seg.yEnd())

        ref_mid = (seg.xStart() + seg.xEnd()) / 2
        read_mid = (seg.yStart() + seg.yEnd()) / 2
        diff = abs(ref_mid - read_mid)

        length = seg.length()

        non_linear_num = diff * length
        sumup_non_linear += non_linear_num


    # normalize by span distance on ref
    ref_span = max(ref_cords) - min(ref_cords)
    final_score = int(sumup_non_linear / ref_span)

    return final_score