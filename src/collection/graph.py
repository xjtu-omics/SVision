import pysam
import os


class Node:
    def __init__(self, chr, ref_start, ref_end, read_start, read_end, seq, is_reverse, id, host):

        self.chr = chr
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.is_reverse = is_reverse
        self.id = id
        self.seq = seq
        self.host = host
        self.depth = 0


        self.node_is_dup = False
        self.dup_from = -1

    def add_depth(self, num):
        self.depth += num

    def set_dup_node(self, dup_from):
        if dup_from != -1:

            self.node_is_dup = True
            self.dup_from = dup_from

    def to_string(self):
        if self.node_is_dup:
            return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(self.id, self.chr, self.ref_start, self.ref_end, self.read_start, self.read_end, '-' if self.is_reverse else '+', self.dup_from)
        else:
            return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(self.id, self.chr, self.ref_start, self.ref_end, self.read_start, self.read_end, '-' if self.is_reverse else '+')


class Edge:
    def __init__(self, node1, node1_is_reverse, node2, node2_is_reverse, id):

        self.node1 = node1
        self.node1_is_reverse = node1_is_reverse

        self.node2 = node2
        self.node2_is_reverse = node2_is_reverse

        self.id = id

        self.edge_is_dup = False

    def set_dup_edge(self):
        self.edge_is_dup = True

    def to_string(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.id, self.node1, '-' if self.node1_is_reverse else '+', self.node2, '-' if self.node2_is_reverse else '+')


class Graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges

        self.appear_time = 1



def parse_gfa_file(gfa_path):
    nodes = []
    edges = []
    with open(gfa_path) as fin:
        for line in fin.readlines():
            line_split = line.strip().split('\t')
            if line_split[0] == 'S':
                node_id = line_split[1]
                node_seq = line_split[2]
                node_host = line_split[3].split(':')[-1]
                node_start = line_split[4].split(':')[-1]

                if len(line_split) == 8:
                    node_is_dup = True
                    dup_from = line_split[7].split(':')[2]
                else:
                    node_is_dup = False
                    dup_from = -1

                node = Node(-1, -1, -1, -1, -1, -1, False, node_id, node_host)
                node.set_dup_node(dup_from)
                nodes.append(node)
            elif line_split[0] == 'L':
                edge_node1 = line_split[1]

                edge_node1_is_reverse = True if line_split[2] == '-' else False
                edge_node2 = line_split[3]
                edge_node2_is_reverse = True if line_split[4] == '-' else False
                edges.append(Edge(edge_node1, edge_node1_is_reverse, edge_node2, edge_node2_is_reverse, 0))
            else:
                pass

    return Graph(nodes, edges)

def graph_is_same_as(graph1, graph2, strict=False, symmetry=False):
    # # collect nodes, edges number and count nodes to dict
    nodes1 = graph1.nodes
    nodes1_num = len(nodes1)
    nodes1_dict = {}
    for node in nodes1:
        node_type = node.id[0]
        if node_type not in nodes1_dict:
            nodes1_dict[node_type] = 0
        nodes1_dict[node_type] += 1

        # count dup node's number
        if node.node_is_dup is True:
            node_type2 = 'D'
            if node_type2 not in nodes1_dict:
                nodes1_dict[node_type2] = 0
            nodes1_dict[node_type2] += 1
    # edge feature
    edges1 = graph1.edges
    edges1_num = len(edges1)
    edges1_path = ''
    for edge in edges1:
        edges1_path += edge.node1
        edges1_path += '-' if edge.node1_is_reverse else '+'
        edges1_path += edge.node2
        edges1_path += '-' if edge.node2_is_reverse else '+'

    # # collect nodes, edges number and count nodes to dict
    nodes2 = graph2.nodes
    nodes2_num = len(nodes2)

    nodes2_dict = {}
    for node in nodes2:
        node_type = node.id[0]
        if node_type not in nodes2_dict:
            nodes2_dict[node_type] = 0
        nodes2_dict[node_type] += 1

        # count dup node's number
        if node.node_is_dup is True:
            node_type2 = 'D'
            if node_type2 not in nodes2_dict:
                nodes2_dict[node_type2] = 0
            nodes2_dict[node_type2] += 1
    # edge feature
    edges2 = graph2.edges
    edges2_num = len(edges2)
    edges2_path = ''
    for edge in edges2:
        edges2_path += edge.node1
        edges2_path += '-' if edge.node1_is_reverse else '+'
        edges2_path += edge.node2
        edges2_path += '-' if edge.node2_is_reverse else '+'

    # nodes number or edge number are not the same, return False
    if nodes1_num != nodes2_num or edges1_num != edges2_num:
        return False
    else:
        # compare node type count
        for node_t in nodes1_dict.keys():
            # unique node type
            if node_t not in nodes2_dict.keys():
                return False
            # same node type but different node number
            if nodes1_dict[node_t] != nodes2_dict[node_t]:
                return False

        if symmetry == True:
            # update graph2's node id to its reverse
            update_dict = {}
            for node in nodes2:
                node_type = node.id[0]
                node_num = int(node.id[1: ])
                new_node_id = "{0}{1}".format(node_type, nodes2_dict[node_type] - node_num - 1)
                update_dict[node.id] = new_node_id

            # update edges's node id
            new_edges2_path = ""
            for i in range(len(edges2) - 1, -1, -1):
                new_edges2_path += update_dict[edges2[i].node2]
                new_edges2_path += '-' if edges2[i].node2_is_reverse else '+'
                new_edges2_path += update_dict[edges2[i].node1]
                new_edges2_path += '-' if edges2[i].node1_is_reverse else '+'

            if edges1_path != new_edges2_path:
                return False

        if strict == True:
            if edges1_path != edges2_path:
                return False
        return True


def cal_overlap_ratio(base_node, target_node, left_most, right_most):

    if target_node == None:
        return 0
    if base_node == target_node:
        return 0
    if base_node.ref_start < left_most:
        return 1.0
    if base_node.ref_end > right_most:
        return 1.0

    base_len = base_node.ref_end - base_node.ref_start
    # base is totaly coverd by target
    if base_node.ref_start >= target_node.ref_start and base_node.ref_end <= target_node.ref_end:
        return 1.0
    #
    if base_node.ref_end >= target_node.ref_end > base_node.ref_start and target_node.ref_start < base_node.ref_start:
        covered_len = target_node.ref_end - base_node.ref_start
        return covered_len / base_len
    if base_node.ref_end < target_node.ref_start < base_node.ref_start and target_node.ref_end > base_node.ref_end:
        covered_len = base_node.ref_end - target_node.ref_start
        return covered_len / base_len

    return 0




def generate_graph(cur_align, next_align, help_aligns, min_sv_size, whole_read_seq, ref_path, next_is_last=True):
    skeleton_nodes_num = 0
    insert_node_num = 0
    skeleton_nodes = []
    insert_nodes = []

    # # add cur and next align to skeleton nodes
    cur_node = Node(cur_align['ref_chr'], cur_align['ref_start'], cur_align['ref_end'], cur_align['q_start'], cur_align['q_end'], cur_align['read_seq'], cur_align['is_reverse'], 'S{0}'.format(skeleton_nodes_num), cur_align['ref_chr'])
    skeleton_nodes.append(cur_node)
    skeleton_nodes_num += 1

    #  overlap on ref
    distance_on_ref = next_align['ref_start'] - cur_align['ref_end']
    dup_len = abs(distance_on_ref)

    if distance_on_ref <= -min_sv_size:


        dup_seg = {
            'ref_chr': next_align['ref_chr'],

            'q_start': next_align['q_start'],
            'q_end': next_align['q_start'] + dup_len,
            'qual': cur_align['qual'],
            'ref_id': cur_align['ref_id'],
            'read_seq': next_align['read_seq'][0: dup_len],
            'ref_start': next_align['ref_start'],
            'ref_end': next_align['ref_start'] + dup_len,
            'is_reverse': cur_align['is_reverse'],
            'read_name': cur_align['read_name']
        }
        help_aligns.append(dup_seg)
        new_next_align = {
            'ref_chr': next_align['ref_chr'],
            'q_start': next_align['q_start'] + dup_len + 1,
            'q_end': next_align['q_end'],

            'qual': cur_align['qual'],
            'ref_id': cur_align['ref_id'],
            'read_seq': next_align['read_seq'][dup_len: ],
            'ref_start': next_align['ref_start'] + dup_len + 1,
            'ref_end': next_align['ref_end'],
            'is_reverse': cur_align['is_reverse'],
            'read_name': cur_align['read_name']
        }

        if new_next_align['ref_start'] < new_next_align['ref_end']:
            next_node = Node(new_next_align['ref_chr'], new_next_align['ref_start'], new_next_align['ref_end'], new_next_align['q_start'], new_next_align['q_end'], new_next_align['read_seq'], new_next_align['is_reverse'],'S{0}'.format(skeleton_nodes_num), new_next_align['ref_chr'])
            skeleton_nodes.append(next_node)
            skeleton_nodes_num += 1
        else:
            if next_is_last is True:
                next_node = None
            else:
                next_node = Node(new_next_align['ref_chr'], new_next_align['ref_start'], new_next_align['ref_start'] + 500,new_next_align['q_start'], new_next_align['q_start'] + 500, new_next_align['read_seq'], new_next_align['is_reverse'], 'S{0}'.format(skeleton_nodes_num), new_next_align['ref_chr'])
                skeleton_nodes.append(next_node)
                skeleton_nodes_num += 1

    else:

        next_node = Node(next_align['ref_chr'], next_align['ref_start'], next_align['ref_end'], next_align['q_start'], next_align['q_end'], next_align['read_seq'], next_align['is_reverse'], 'S{0}'.format(skeleton_nodes_num), next_align['ref_chr'])
        skeleton_nodes.append(next_node)
        skeleton_nodes_num += 1



    # # check help aligns whether they are insert_nodes
    left_most = cur_align['ref_start']
    right_most = next_align['ref_end']
    for align in help_aligns:
        tmp_node = Node(align['ref_chr'], align['ref_start'], align['ref_end'], align['q_start'], align['q_end'], align['read_seq'], align['is_reverse'], 'None', align['read_name'])

        overlap_with_cur = cal_overlap_ratio(tmp_node, cur_node, left_most, right_most)
        overlap_with_next = cal_overlap_ratio(tmp_node, next_node, left_most, right_most)

        # reverse help align, check whether it is overlaped with cur and next align
        if align['is_reverse'] == True:

            if overlap_with_cur > 0.8:
                tmp_node.id = 'I{0}'.format(insert_node_num)
                # set dup flag
                tmp_node.node_is_dup = True
                tmp_node.dup_from = cur_node.id

                insert_nodes.append(tmp_node)
                insert_node_num += 1

            elif overlap_with_next > 0.8:
                tmp_node.id = 'I{0}'.format(insert_node_num)
                # set dup flag
                tmp_node.node_is_dup = True
                tmp_node.dup_from = next_node.id

                insert_nodes.append(tmp_node)
                insert_node_num += 1
            else:
                tmp_node.id = 'S{0}'.format(skeleton_nodes_num)
                tmp_node.host = align['ref_chr']
                skeleton_nodes.append(tmp_node)
                skeleton_nodes_num += 1

        # non-reverse help align is directly assigned to insert nodes
        else:
            if overlap_with_cur > 0.8:
                # set dup flag
                tmp_node.node_is_dup = True
                tmp_node.dup_from = cur_node.id
            elif overlap_with_next > 0.8:
                # set dup flag
                tmp_node.node_is_dup = True
                tmp_node.dup_from = next_node.id

            tmp_node.id = 'I{0}'.format(insert_node_num)
            insert_nodes.append(tmp_node)
            insert_node_num += 1

    # generate edges
    nodes_srt_by_read = sorted(skeleton_nodes + insert_nodes, key=lambda node: node.read_start)
    edges = []
    edge_num = 0
    for i in range(1, len(nodes_srt_by_read)):
        # check whether there is a ins-gap on read between two nodes
        if nodes_srt_by_read[i].read_start - nodes_srt_by_read[i - 1].read_end > min_sv_size:
            tmp_node = Node(nodes_srt_by_read[i].chr, nodes_srt_by_read[i].ref_start, nodes_srt_by_read[i].ref_start, nodes_srt_by_read[i - 1].read_end + 1, nodes_srt_by_read[i].read_start - 1, whole_read_seq[nodes_srt_by_read[i - 1].read_end + 1: nodes_srt_by_read[i].read_start - 1], False, 'I{0}'.format(insert_node_num), cur_align['read_name'])
            insert_nodes.append(tmp_node)
            insert_node_num += 1

            edges.append(Edge(nodes_srt_by_read[i - 1].id, nodes_srt_by_read[i - 1].is_reverse, tmp_node.id, tmp_node.is_reverse, 'E{0}'.format(edge_num)))
            edge_num += 1
            edges.append(Edge(tmp_node.id, tmp_node.is_reverse, nodes_srt_by_read[i].id, nodes_srt_by_read[i].is_reverse, 'E{0}'.format(edge_num)))
            edge_num += 1


        else:
            edges.append(Edge(nodes_srt_by_read[i - 1].id, nodes_srt_by_read[i - 1].is_reverse, nodes_srt_by_read[i].id, nodes_srt_by_read[i].is_reverse, 'E{0}'.format(edge_num)))
            edge_num += 1


    # add nodes if there is a del-gap between every two skeleton nodes
    skeleton_srt_by_ref = sorted(skeleton_nodes, key=lambda node: node.ref_start)
    for i in range(1, len(skeleton_srt_by_ref)):
        node_chr = skeleton_srt_by_ref[i].chr
        gap_on_ref = skeleton_srt_by_ref[i].ref_start - skeleton_srt_by_ref[i - 1].ref_end
        # there is gap on ref, then we add a node
        if gap_on_ref > min_sv_size:
            node_ref_start = skeleton_srt_by_ref[i - 1].ref_end + 1
            node_ref_end = skeleton_srt_by_ref[i].ref_start - 1
            ref_seq = pysam.FastaFile(ref_path).fetch(node_chr, node_ref_start, node_ref_end)

            skeleton_nodes.append(Node(node_chr, node_ref_start, node_ref_end, -1, -1, ref_seq, False, 'S{0}'.format(skeleton_nodes_num), skeleton_srt_by_ref[i].host))
            skeleton_nodes_num += 1

    # # update the ids of nodes
    # update skeleton nodes
    update_id = {}
    skeleton_srt_by_ref = sorted(skeleton_nodes, key=lambda node: node.ref_start)
    i = 0
    for node in skeleton_srt_by_ref:
        update_id[node.id] = "S{0}".format(i)
        node.id = "S{0}".format(i)
        i += 1

    # # update insert nodes
    insert_srt_by_read = sorted(insert_nodes, key=lambda node: node.read_start)
    i = 0
    for node in insert_srt_by_read:
        update_id[node.id] = "I{0}".format(i)
        node.id = "I{0}".format(i)
        i += 1

    # # update ids of edge
    for edge in edges:
        edge.node1 = update_id[edge.node1]
        edge.node2 = update_id[edge.node2]

    return Graph(skeleton_srt_by_ref + insert_srt_by_read, edges)


def parse_graph_features(graph):
    # # collect nodes, edges number and count nodes to dict
    nodes = graph.nodes
    nodes_dict = {}

    for node in nodes:
        node_type = node.id[0]
        if node_type not in nodes_dict:
            nodes_dict[node_type] = 0
        nodes_dict[node_type] += 1

        # count dup node's number
        if node.node_is_dup is True:
            node_type2 = 'D'
            if node_type2 not in nodes_dict:
                nodes_dict[node_type2] = 0
            nodes_dict[node_type2] += 1

    nodes_feature = ""
    for key in nodes_dict:
        nodes_feature += "{0}:{1},".format(key, nodes_dict[key])
    nodes_feature = nodes_feature[: -1]

    # edge feature
    edges = graph.edges
    edges_num = len(edges)
    edges_path = ''
    for edge in edges:


        edges_path += edge.node1
        edges_path += '-' if edge.node1_is_reverse else '+'
        edges_path += edge.node2
        edges_path += '-' if edge.node2_is_reverse else '+'

    return nodes_feature, edges_num, edges_path

def collect_csv_same_format(gfa_path, vcf_path, out_path, sample, min_support):
    # svision v1.2.1 add
    graph_vcf = open(os.path.join(out_path, "{0}.svision.s{1}.graph.vcf".format(sample, min_support)), 'w')
    # end

    exactly_matching = {}
    symmetry_matching = {}
    symmetry_matching_pair = {}

    cnt = 0
    vcf_file = pysam.VariantFile(vcf_path)
    graph_vcf.write(str(vcf_file.header))

    for record in vcf_file:

        cnt += 1

        # if not satisfy, output directly
        # if 'CSV' not in str(record) or 'Uncovered' in str(record) or int(record.info['SUPPORT']) < 5:
        if 'CSV' not in str(record):
            # SVision v1.2.1, ADD. r
            # # add graph id to vcf info column, we need recalculate GT info because pysam does not output GT
            vaf = float(record.info['VAF'])
            if vaf > 0.8:
                GT = '1/1'
            elif vaf > 0.3:
                GT = '0/1'
            elif vaf >= 0:
                GT = '0/0'
            else:
                GT = './.'
            DV = int(record.info['SUPPORT'])
            DR = round(DV / (vaf + 0.001)) - DV

            gt_format = 'GT:DR:DV\t{0}:{1}:{2}'.format(GT, DR, DV)
            graph_vcf.write(str(record).strip() + ';GraphID=-1'+ '\t' + gt_format + '\n')
            # Add, End

            continue

        chr = record.contig
        start = record.start + 1
        end = record.stop

        sv_type = record.info['SVTYPE']

        if '_' in record.id:
            sub_id = str(record.id).split('_')[1]
            target_gfa = '{0}-{1}-{2}_{3}_{4}'.format(chr, start, end, sub_id, sv_type)
        else:
            target_gfa = '{0}-{1}-{2}_{3}'.format(chr, start, end, sv_type)

        target_gfa_path = '_'.join(target_gfa.split('_')[0: -1])
        if not os.path.exists(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path))):
            # print(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa)))
            continue

        # # DEBUG: check path
        # target_graph = parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path)))
        # # edge feature
        # edges1 = target_graph.edges
        # edges1_num = len(edges1)
        # edges1_path = ''
        # for edge in edges1:
        #     edges1_path += edge.node1
        #     edges1_path += '-' if edge.node1_is_reverse else '+'
        #     edges1_path += edge.node2
        #     edges1_path += '-' if edge.node2_is_reverse else '+'
        # if edges1_path == 'S0+I0+I0+S1+':
        #     continue
        # # # end

        # # exactly matching
        exactly_flag = 0
        exactly_base_gfas = ''
        cnt = -1
        for base_gfa in exactly_matching.keys():
            cnt += 1
            base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
            target_gfa_path = '_'.join(target_gfa.split('_')[0: -1])
            if graph_is_same_as(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path))),
                                parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(base_gfa_path))), strict=True):
                exactly_flag = 1
                exactly_base_gfas = base_gfa
                break

        if exactly_flag == 0:
            exactly_matching[target_gfa] = [target_gfa]
            graph_id = len(exactly_matching.keys()) - 1
        elif exactly_flag == 1:
            exactly_matching[exactly_base_gfas].append(target_gfa)
            graph_id = cnt
        else:
            graph_id = -1
            pass

        # SVision v1.2.1, ADD. r
        # # add graph id to vcf info column, we need recalculate GT info because pysam does not output GT
        vaf = float(record.info['VAF'])
        if vaf > 0.8:
            GT = '1/1'
        elif vaf > 0.3:
            GT = '0/1'
        elif vaf >= 0:
            GT = '0/0'
        else:
            GT = './.'
        DV = int(record.info['SUPPORT'])
        DR = round(DV / (vaf + 0.001)) - DV

        gt_format = 'GT:DR:DV\t{0}:{1}:{2}'.format(GT, DR, DV)
        graph_vcf.write(str(record).strip() + ';GraphID={0}'.format(graph_id) + '\t' + gt_format + '\n')
        # Add, End

    graph_vcf.close()


    # deal symmetry
    gfas = list(exactly_matching.keys())
    for i in range(len(gfas)):
        for j in range(i + 1, len(gfas)):
            base_gfa =gfas[i]
            target_gfa = gfas[j]

            base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
            target_gfa_path = '_'.join(target_gfa.split('_')[0: -1])
            if not graph_is_same_as(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path))),
                                    parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(base_gfa_path))), strict=True):

                if graph_is_same_as(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path))),
                                    parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(base_gfa_path))), strict=False, symmetry=True):
                    symmetry_matching[base_gfa] = exactly_matching[base_gfa] + exactly_matching[target_gfa]
                    symmetry_matching_pair[base_gfa] = ['{0},{1}'.format(i, j), '{0},{1}'.format(len(exactly_matching[base_gfa]), len(exactly_matching[target_gfa]))]
                else:
                    pass

    # write summary
    with open(os.path.join(out_path, '{0}.graph_exactly_match.txt'.format(sample)), 'w') as fout:
        cnt = 0

        for base_gfa in exactly_matching.keys():

            base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])

            nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(base_gfa_path))))
            fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}\n'.format(cnt, len(exactly_matching[base_gfa]), nodes_feat, edges_feat, path_feat))
            fout.write('\t'.join(exactly_matching[base_gfa]))
            fout.write('\n')

            cnt += 1
    with open(os.path.join(out_path, '{0}.graph_symmetry_match.txt'.format(sample)), 'w') as fout:

        for base_gfa in symmetry_matching.keys():
            if len(symmetry_matching[base_gfa]) != 0:
                base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])

                # write to file
                nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(base_gfa_path))))
                fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}'.format(symmetry_matching_pair[base_gfa][0], symmetry_matching_pair[base_gfa][1], nodes_feat, edges_feat, path_feat))

                # for target_gfa in symmetry_matching[base_gfa]:
                target_gfa_path = '_'.join(symmetry_matching[base_gfa][-1].split('_')[0: -1])
                _, _, path_feat = parse_graph_features(parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format(target_gfa_path))))
                fout.write(',{0}\n'.format(path_feat))

                fout.write('\t'.join(symmetry_matching[base_gfa]))
                fout.write('\n')
    return exactly_matching, symmetry_matching



def merge_samples_results():
    out_path = '/home/DATA/CCS/graph2/'

    tot_exactly_matching = {}
    tot_symmetry_matching = {}
    tot_symmetry_matching_pair = {}

    for sample in ['HG00731', 'HG00732', 'HG00733']:
    # for sample in ['NA19240', 'HG00512']:
        print('-------------------' + sample + '--------------------')

        gfa_path = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(sample)
        sample_exact_file = os.path.join(out_path, '{0}_graph_exactly_match.txt'.format(sample))
        for line in open(sample_exact_file).readlines():
            if '>' in line:
                continue

            target_gfa = line.split('\t')[0]
            found_flag = 0
            for base_gfa in tot_exactly_matching.keys():
                base_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
                if graph_is_same_as(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format('_'.join(base_gfa.split('_')[0: -1])))),
                                    parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format('_'.join(target_gfa.split('_')[0: -1])))), strict=True):
                    for gfa in line.split('\t'):
                        tot_exactly_matching[base_gfa].append(gfa.strip() + '_' + sample)
                    found_flag = 1
                    break
            if found_flag == 0:
                tot_exactly_matching[target_gfa] = []
                for gfa in line.split('\t'):
                    tot_exactly_matching[target_gfa].append(gfa.strip() + '_' + sample)

    # # deal tot symmetry match
    tot_exactly_matching = dict(sorted(tot_exactly_matching.items(), key=lambda d: len(d[1]), reverse=True))

    gfas = list(tot_exactly_matching.keys())
    for i in range(len(gfas)):
        for j in range(i + 1, len(gfas)):
            base_gfa = gfas[i]
            target_gfa = gfas[j]

            base_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
            target_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_exactly_matching[target_gfa][0].split('_')[-1])

            base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
            target_gfa_path = '_'.join(target_gfa.split('_')[0: -1])
            if not graph_is_same_as(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))),
                                    parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))), strict=True):

                if graph_is_same_as(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))),
                                    parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))), strict=False, symmetry=True):
                    tot_symmetry_matching[base_gfa] = tot_exactly_matching[base_gfa] + tot_exactly_matching[target_gfa]
                    tot_symmetry_matching_pair[base_gfa] = ['{0},{1}'.format(i, j), '{0},{1}'.format(len(tot_exactly_matching[base_gfa]), len(tot_exactly_matching[target_gfa]))]
                else:
                    pass

    # write summary
    with open(os.path.join(out_path, 'Tot_graph_exactly_match.txt'), 'w') as fout:
        cnt = 0

        for base_gfa in tot_exactly_matching.keys():

            base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
            base_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
            nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))))
            fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}\n'.format(cnt, len(tot_exactly_matching[base_gfa]), nodes_feat, edges_feat, path_feat))
            fout.write('\t'.join(tot_exactly_matching[base_gfa]))
            fout.write('\n')

            cnt += 1
    with open(os.path.join(out_path, 'Tot_graph_symmetry_match.txt'), 'w') as fout:

        for base_gfa in tot_symmetry_matching.keys():
            if len(tot_symmetry_matching[base_gfa]) != 0:
                base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
                base_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])

                # write to file
                nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))))
                fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}'.format(tot_symmetry_matching_pair[base_gfa][0], tot_symmetry_matching_pair[base_gfa][1], nodes_feat, edges_feat, path_feat))

                # for target_gfa in symmetry_matching[base_gfa]:
                target_gfa_prefix = '/home/DATA/CCS/{0}/ngmlr/sv_call_graph2/graphs/'.format(tot_symmetry_matching[base_gfa][-1].split('_')[-1])
                target_gfa_path = '_'.join(tot_symmetry_matching[base_gfa][-1].split('_')[0: -2])
                _, _, path_feat = parse_graph_features(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))))

                fout.write(',{0}\n'.format(path_feat))
                fout.write('\t'.join(tot_symmetry_matching[base_gfa]))
                fout.write('\n')



# if __name__ == '__main__':
    # gfa1 = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/SVision/results/graph/graph_test_data/chr1-209761662-209762731-ins-dup.gfa'
    # gfa2 = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/SVision/results/graph/graph_test_data/chr1-239952716-239953381-ins-inv-dup.gfa'
    #
    # graph1 = parse_gfa_file(gfa1)
    # graph2 = parse_gfa_file(gfa2)
    #
    # print(graph_is_same_as(graph1, graph2))
    #

    # sample = 'NA19240'
    # #
    # vcf_path = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/SVision.s5.vcf'.format(sample)
    # gfa_path = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(sample)
    # out_path = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/'.format(sample)
    # sample_exactly_matching, sample_symmetry_matching = collect_csv_same_format(gfa_path, vcf_path, out_path, sample)
    #
    # exit()

    # merge_samples()

    # tot_exactly_matching = {}
    # tot_symmetry_matching = {}
    # tot_symmetry_matching_pair = {}
    #
    # out_path = '/mnt/d/Data/Bams/TGS/CCS/NA19240/ngmlr/sv_call_graph/'
    #
    # # for sample in ['HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'HG02818', 'HG03125', 'HG03486', 'NA12878', 'NA19238', 'NA19239', 'NA19240']:
    # for sample in ['NA19240', 'HG00512']:
    #     print('-------------------' + sample + '--------------------')
    #     vcf_path = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/SVision.s5.vcf'.format(sample)
    #     gfa_path = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(sample)
    #     sample_exactly_matching, sample_symmetry_matching = collect_csv_same_format(gfa_path, vcf_path, out_path, sample)
    #
    #
    #     # # deal tot_exactly_match
    #     for target_gfa in sample_exactly_matching.keys():
    #         found_flag = 0
    #         for base_gfa in tot_exactly_matching.keys():
    #             base_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
    #             if graph_is_same_as(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format('_'.join(base_gfa.split('_')[0: -1])))),
    #                                 parse_gfa_file(os.path.join(gfa_path, '{0}.gfa'.format('_'.join(target_gfa.split('_')[0: -1])))), strict=True):
    #                 for gfa in sample_exactly_matching[target_gfa]:
    #                     tot_exactly_matching[base_gfa].append(gfa + '_' + sample)
    #                 found_flag = 1
    #                 break
    #         if found_flag == 0:
    #             tot_exactly_matching[target_gfa] = []
    #             for gfa in sample_exactly_matching[target_gfa]:
    #                 tot_exactly_matching[target_gfa].append(gfa + '_' + sample)
    #
    # # # deal tot symmetry match
    # tot_exactly_matching = dict(sorted(tot_exactly_matching.items(), key=lambda d: len(d[1]), reverse=True))
    #
    #
    # gfas = list(tot_exactly_matching.keys())
    # for i in range(len(gfas)):
    #     for j in range(i + 1, len(gfas)):
    #         base_gfa = gfas[i]
    #         target_gfa = gfas[j]
    #
    #         base_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
    #         target_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_exactly_matching[target_gfa][0].split('_')[-1])
    #
    #         base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
    #         target_gfa_path = '_'.join(target_gfa.split('_')[0: -1])
    #         if not graph_is_same_as(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))),
    #                                 parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))), strict=True):
    #
    #             if graph_is_same_as(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))),
    #                                 parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))), strict=False, symmetry=True):
    #                 tot_symmetry_matching[base_gfa] = tot_exactly_matching[base_gfa] + tot_exactly_matching[target_gfa]
    #                 tot_symmetry_matching_pair[base_gfa] = ['{0},{1}'.format(i, j), '{0},{1}'.format(len(tot_exactly_matching[base_gfa]), len(tot_exactly_matching[target_gfa]))]
    #             else:
    #                 pass
    #
    # # write summary
    # with open(os.path.join(out_path, 'Tot_graph_exactly_match.txt'), 'w') as fout:
    #     cnt = 0
    #
    #     for base_gfa in tot_exactly_matching.keys():
    #
    #         base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
    #         base_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
    #         nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))))
    #         fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}\n'.format(cnt, len(tot_exactly_matching[base_gfa]), nodes_feat, edges_feat, path_feat))
    #         fout.write('\t'.join(tot_exactly_matching[base_gfa]))
    #         fout.write('\n')
    #
    #         cnt += 1
    # with open(os.path.join(out_path, 'Tot_graph_symmetry_match.txt'), 'w') as fout:
    #
    #     for base_gfa in tot_symmetry_matching.keys():
    #         if len(tot_symmetry_matching[base_gfa]) != 0:
    #             base_gfa_path = '_'.join(base_gfa.split('_')[0: -1])
    #             base_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_exactly_matching[base_gfa][0].split('_')[-1])
    #
    #             # write to file
    #             nodes_feat, edges_feat, path_feat = parse_graph_features(parse_gfa_file(os.path.join(base_gfa_prefix, '{0}.gfa'.format(base_gfa_path))))
    #             fout.write('> GraphId={0}\tNumber={1}\tNodes={2}\tEdges={3}\tPath={4}'.format(tot_symmetry_matching_pair[base_gfa][0], tot_symmetry_matching_pair[base_gfa][1], nodes_feat, edges_feat, path_feat))
    #
    #             # for target_gfa in symmetry_matching[base_gfa]:
    #             target_gfa_prefix = '/mnt/d/Data/Bams/TGS/CCS/{0}/ngmlr/sv_call_graph/graphs/'.format(tot_symmetry_matching[base_gfa][-1].split('_')[-1])
    #             target_gfa_path = '_'.join(tot_symmetry_matching[base_gfa][-1].split('_')[0: -2])
    #             _, _, path_feat = parse_graph_features(parse_gfa_file(os.path.join(target_gfa_prefix, '{0}.gfa'.format(target_gfa_path))))
    #
    #             fout.write(',{0}\n'.format(path_feat))
    #             fout.write('\t'.join(tot_symmetry_matching[base_gfa]))
    #             fout.write('\n')