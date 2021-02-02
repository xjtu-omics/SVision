#!/usr/bin/env python3

# encoding: utf-8

import pysam

class Signature:
    """
    Signature class for abnormal alignment observed from BAM file
    """
    def __init__(self, contig, tstart, tend, type, qname, sorted_aligns, all_bkps, mechanism):
        self.contig = contig
        self.tstart = tstart
        self.tend = tend

        self.qname = qname
        self.type = type

        self.bkps = all_bkps
        self.sorted_aligns = sorted_aligns

        self.mechanism = mechanism

        if self.tend < self.tstart:
            print("[WARNING]: " + "Signature with invalid coordinates (end < start): " + self.to_string())


    def get_source(self):
        return (self.contig, self.tstart, self.tend)


    def get_key(self):
        return (self.contig, (self.tstart + self.tend) // 2)

    def position_distance_to(self, signature2):
        """Return position distance between two signatures."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = signature2.get_source()
        this_center = (this_start + this_end) // 2
        other_center = (other_start + other_end) // 2
        if this_contig == other_contig:
            return min(abs(this_start - other_start), abs(this_end - other_end), abs(this_center - other_center))
        else:
            return float("inf")

    def set_graph(self, graph):
        self.graph = graph
    def to_string(self):
        return "\t".join(["{0}","{1}","{2}","{3}"]).format(self.contig, self.tstart, self.tend, "{0};{1}".format(self.type, self.qname))

    # def get_D_seq(self, genome_file, fai_dict, flank):
    #
    #     if self.tend > fai_dict[self.contig][0]:
    #         sys.stderr.write('ERROR! Signature past the genome end. ' + str(self.tend))
    #
    #     del_start = max(self.tstart - flank, 0)
    #     del_end = min(self.tend + flank, fai_dict[self.contig][0])
    #
    #     del_seq = extract_sequence([self.contig, del_start, del_end], genome_file, fai_dict)
    #
    #     return del_seq
    #
    # def get_I_seq(self, fai_dict, flank):
    #
    #     ins_start = max(self.qstart - flank, 0)
    #     ins_end = min(self.qend + flank, fai_dict[self.contig][0])
    #
    #     ins_seq = self.seq[ins_start, ins_end]
    #
    #     return ins_seq

    def get_segs_cords(self):
        # # sort by read cords
        # sorted_alignment_list = sorted(self.sorted_aligns, key=lambda aln: (aln['q_start'], aln['q_end']))
        # print('out-----------------')
        # for first_seg in sorted_alignment_list:
        #     print(first_seg['ref_start'], first_seg['ref_end'], first_seg['q_start'], first_seg['q_end'],
        #           first_seg['is_reverse'], first_seg['is_supplementary'])

        sorted_alignment_list = self.sorted_aligns
        read_start = sorted_alignment_list[0]['q_start']
        read_end = sorted_alignment_list[-1]['q_end']
        ref_start = sorted_alignment_list[0]['ref_start']
        ref_end = sorted_alignment_list[-1]['ref_end']

        main_segs_cord = []
        other_segs_cord = []
        # print('start cords: ', ref_start, read_start)

        # # adjust cords again to get cords start at zero
        for i in range(len(sorted_alignment_list)):
            cur_ref_start = sorted_alignment_list[i]['ref_start']
            cur_ref_end = sorted_alignment_list[i]['ref_end']
            cur_read_start = sorted_alignment_list[i]['q_start']
            cur_read_end = sorted_alignment_list[i]['q_end']
            sorted_alignment_list[i]['ref_start'] = cur_ref_start - ref_start
            sorted_alignment_list[i]['ref_end'] = cur_ref_end - ref_start
            sorted_alignment_list[i]['q_start'] = cur_read_start - read_start
            sorted_alignment_list[i]['q_end'] = cur_read_end - read_start

            # # find main_segs, the first and the last seg are main_segs
            if i == 0 or i == (len(sorted_alignment_list) - 1):
                main_segs_cord.append([[sorted_alignment_list[i]['q_start'], sorted_alignment_list[i]['q_end']],
                                       [sorted_alignment_list[i]['ref_start'], sorted_alignment_list[i]['ref_end']], 0])
            else:
                if sorted_alignment_list[i]['is_reverse']:
                    other_segs_cord.append([[sorted_alignment_list[i]['q_end'], sorted_alignment_list[i]['q_start'] ],
                                       [sorted_alignment_list[i]['ref_start'], sorted_alignment_list[i]['ref_end']], 1])
                else:
                    other_segs_cord.append([[sorted_alignment_list[i]['q_start'], sorted_alignment_list[i]['q_end']],
                                        [sorted_alignment_list[i]['ref_start'], sorted_alignment_list[i]['ref_end']], 0])

        read_len = sorted_alignment_list[-1]['q_end']
        ref_len = sorted_alignment_list[-1]['ref_end']

        # print(main_segs_cord)
        return ref_len, read_len, main_segs_cord, other_segs_cord




class Cluster:

    def __init__(self, cluster, this_sample_path, cstart_end=None):

        self.sample_path = this_sample_path

        self.sigs = cluster
        self.contig = cluster[0].contig
        self.read_num = len(cluster)

        # Debug usage
        self.query_info = ""
        self.signatures = cluster
        self.coverage = 0

        starts = []
        ends = []

        for sig in cluster:
            starts.append(sig.tstart)
            ends.append(sig.tend)

            self.query_info += "{0},{1}:{2},{3};".format(str(sig.bkps), sig.tstart, sig.tend, sig.type)

        if cstart_end == None:
            self.cstart = sum(starts) / len(cluster)
            self.cend = sum(ends) / len(cluster)
        else:
            self.cstart = cstart_end[0]
            self.cend = cstart_end[1]

        # # deesv v 1.0.2. ADD. filter clusters with wrong cord
        self.abandon = 0
        if self.cstart < 0 or self.cend < 0 or self.cstart > self.cend:
            # print(self.sigs, self.cstart, self.cend)
            # for i in self.sigs:
            #     print(i.qname, i.tstart, i.tend)
            self.abandon = 1

        if self.abandon == 0:
            self.stats_coverd_reads()
        # # End ADD

    def stats_coverd_reads(self):
        aln_file = pysam.AlignmentFile(self.sample_path)
        reads = aln_file.fetch(self.contig, self.cstart, self.cend)

        for i in reads:
            self.coverage += 1


    def get_signatures(self):
        return self.signatures

    def to_string(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.contig, int(self.cstart), int(self.cend), int(self.cend) - int(self.cstart), 'None', self.read_num, self.query_info)




def read_fai_file(fai_file_name):
    """
    Get a dictionary of an FAI file (reference index).

    """

    fai_record_dict = {}

    with open(fai_file_name, 'r') as fai_file:

        for line in fai_file:

            line = line.strip()

            if not line:
                continue

            fai_record = line.split()

            fai_record_dict[fai_record[0]] = [int(x) for x in fai_record[1:]]

    return fai_record_dict
