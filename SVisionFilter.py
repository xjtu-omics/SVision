#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/6/9

'''
import sys
import argparse
import pysam
import os
import pandas as pd
from intervaltree import IntervalTree


AUTOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

def parse_exclude_regions(exclude):
    exclude_dict = {}
    with open(exclude, 'r') as f:
        for line in f:
            entries = line.strip().split("\t")
            start = int(entries[1])
            end = int(entries[2])
            if entries[0] not in exclude_dict:
                exclude_dict[entries[0]] = IntervalTree()
                exclude_dict[entries[0]][start:end] = (start, end)
            else:
                exclude_dict[entries[0]][start:end] = (start, end)
    return exclude_dict

def contains_gaps(chrom, start, end, ref):
    seq = ref.fetch(chrom, start - 2000, end + 2000)
    in_gap = False
    for base in seq:
        if base == 'N':
            in_gap = True
            break
    return in_gap

def get_incomplete_graphID(graph):
    incomplete_graphs = []
    complete_graphs = []
    with open(graph, 'r') as f:
        for line in f:
            if ">" in line:
                entries = line.strip().split("\t")
                id = entries[0].split("=")[1]
                path = entries[4].split("=")[1]

                if path[-3] != "S":
                    incomplete_graphs.append(id)
                else:
                    complete_graphs.append(id)

    return incomplete_graphs, complete_graphs

def overlap(a,b,c,d):
    r = 0 if a==c and b==d else min(b,d)-max(a,c)
    if r>=0: return r


def run_filter(svision_vcf, svision_exact_graph, ref_fasta, exclude_file, exclude_graphid, min_sr, max_sv_size, outdir):

    filtered_prefix = '.'.join(os.path.basename(svision_vcf).split('.')[0:-1])

    exclude_dict = parse_exclude_regions(exclude_file)
    ref_file = pysam.FastaFile(ref_fasta)
    incomplete_graphs, complete_graphs = get_incomplete_graphID(svision_exact_graph)

    csv_by_id = {}
    all_sv_num = 0
    filtered_vcf_writer = open(outdir + '/{0}.filtered.vcf'.format(filtered_prefix), 'w')

    with open(svision_vcf, 'r') as f:
        for line in f:
            if "#" in line:
                continue
            entries = line.strip().split("\t")
            chrom = entries[0]
            id = entries[2]
            if "_" in id:
                id = id.split("_")[0]

            if chrom in AUTOSOMES:
                start = int(entries[1])

                info_tokens = entries[7].split(";")
                info_dict = {}

                for token in info_tokens:
                    info_dict[token.split("=")[0]] = token.split("=")[1]

                this_chrom_exclude_tree = exclude_dict[chrom]
                if this_chrom_exclude_tree.overlaps(start, int(info_dict['END'])):
                    continue

                if contains_gaps(chrom, start, int(info_dict['END']), ref_file):
                    continue

                if int(info_dict["SVLEN"]) < max_sv_size and int(info_dict['SUPPORT']) >= min_sr:
                    all_sv_num += 1
                    if info_dict['GraphID'] != "-1":
                        if info_dict['GraphID'] not in incomplete_graphs and entries[6] != 'Uncovered':
                            if id in csv_by_id:
                                csv_by_id[id].append((chrom, start, info_dict['END'], entries[2], info_dict['GraphID'], int(info_dict['SUPPORT'])))
                            else:
                                csv_by_id[id] = [(chrom, start, info_dict['END'], entries[2], info_dict['GraphID'], int(info_dict['SUPPORT']))]

                            filtered_vcf_writer.write(line)
                    else:
                        filtered_vcf_writer.write(line)

    high_conf_csv_list = list()
    for id, csvs in csv_by_id.items():

        sorted_csvs = sorted(csvs, key=lambda x: x[-1], reverse=True)
        selected_csv = sorted_csvs[0]
        if selected_csv[4] in exclude_graphid:
           continue

        high_conf_csv_list.append((selected_csv[0], int(selected_csv[1]), int(selected_csv[2]), selected_csv[3], selected_csv[4]))

    df_high_conf_csvs = pd.DataFrame(high_conf_csv_list, columns=['chrom', 'start', 'end', 'id', 'graphid'])
    sorter_index = dict(zip(AUTOSOMES, range(len(AUTOSOMES))))
    df_high_conf_csvs['chrom_rank'] = df_high_conf_csvs['chrom'].map(sorter_index)

    # Write high confident csvs
    df_high_conf_csvs.sort_values(by=['chrom_rank', 'start'], inplace=True)
    df_high_conf_csvs.drop('chrom_rank', 1, inplace=True)
    df_high_conf_csvs.to_csv(outdir + "/{0}.HQ-CSVs.tsv".format(filtered_prefix), sep="\t", header=True, index=False)

    print('SV after filtering {0}, containing {1} high-quality CSVs'.format(all_sv_num, len(high_conf_csv_list)))


def main():

    ## Filtering for the results in the paper
    # svision_vcf = '/Users/apple/SVision/HG00733/svision/HG00733.svision.s5.graph.vcf'
    # svision_exact_graph = '/Users/apple/SVision/HG00733/svision/HG00733.graph_exactly_match.txt'
    # ref_fasta = "/Users/apple/Data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # exclude_file = "/Users/apple/Data/genome/grch38.exclude_regions_cen.bed"
    # outdir = '/Users/apple/SVision/HG00733/svision'
    # exclude_graphid = ['0', '4']

    arguments = sys.argv[1:]

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""This is used to filter raw SVision calls and collects high-quality CSVs, 
which is consisted of (default settings used to produce results in the publication):
 1. SV length between 50bp and 100Kbp.
 2. SV supported by at least 5 reads.                                    
 3. SV detected in gapped regions.
 4. SV detected in low mapping quality regions, such as centromere specified by --region parameter.
 5. SV of partial graph representation because of read length limitation.
 6. CSV graph indicate insertion associated duplication events, indicating by -i parameter.
""")

    required = parser.add_argument_group("Required parameters")
    required.add_argument('-v', dest='vcf', type=str, help='Raw SVision callset in VCF format')
    required.add_argument('-g', dest='graph', type=str, help='Exactly matched graphs obtained from SVision output')
    required.add_argument('-r', dest='ref', type=str, help='Path to the reference file in .fa or .fasta format')
    required.add_argument('-o', dest='outdir', type=os.path.abspath, help='Output directory of filtered calls and CSVs')
    required.add_argument('-i', dest='id', type=str, help='Comma separated graph id to exclude. (e.g 0,4)')


    optional = parser.add_argument_group("Optional parameters")
    optional.add_argument('--region', type=os.path.abspath, help='Path to the BED file of exclude regions')
    optional.add_argument('--min_sr', type=int, default=5, help='Filtering events by supporting reads (default: %(default)s)')
    optional.add_argument('--max_sv_size', type=int, default=100000, help='Maximum event size to include (default: %(default)s)')

    options = parser.parse_args(arguments)

    svision_vcf = options.vcf
    svision_exact_graph = options.graph
    ref_fasta = options.ref
    exclude_file = options.region
    outdir = options.outdir
    exclude_graphid = options.id.split(',')
    min_sr = options.min_sr
    max_sv_size = options.max_sv_size

    run_filter(svision_vcf, svision_exact_graph, ref_fasta, exclude_file, exclude_graphid, min_sr, max_sv_size, outdir)

if __name__ == '__main__':
    main()