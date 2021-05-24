import numpy as np
import os
from bs4 import BeautifulSoup as bs
from bs4.element import NavigableString

def process_tbl(tbl_file):

    sv_rep_info_list = list()
    line_counter = 1

    for line in open(tbl_file, 'r'):

        if line_counter == 6:
            value = line.strip().split(':')[1].strip()
            per_str = 'masked:{0}'.format(value.split(' ')[3])
            sv_rep_info_list.append(per_str)

        if line_counter == 11:
            value = line.strip().split(':')[1].strip()
            sine_str = 'SINE:{0}'.format(value.split(' ')[-2])
            sv_rep_info_list.append(sine_str)

        # if line_counter == 12:
        #     alu_str = 'ALU:{0}'.format(line.strip().split(' ')[-2])
        #     sv_rep_info_list.append(alu_str)

        if line_counter == 15:
            value = line.strip().split(':')[1].strip()
            line_str = 'LINE:{0}'.format(value.split(' ')[-2])
            sv_rep_info_list.append(line_str)

        if line_counter == 20:
            value = line.strip().split(':')[1].strip()
            ltr_str = 'LTR:{0}'.format(value.split(' ')[-2])
            sv_rep_info_list.append(ltr_str)

        # if line_counter == 30:
        #     value = line.strip().split(':')[1].strip()
        #     other_str = 'Others:{0}'.format(value.split(' ')[-2])
        #     sv_rep_info_list.append(other_str)

        # if line_counter == 32:
        #     value = line.strip().split(':')[1].strip()
        #     spersed_str = 'total_interspersed:{0}'.format(value.split(' ')[-2])
        #     sv_rep_info_list.append(spersed_str)

        # if line_counter == 35:
        #     value = line.strip().split(':')[1].strip()
        #     rna_str = 'RNA:{0}%'.format(value.split(' ')[-2])
        #     sv_rep_info_list.append('\t{0}'.format(rna_str))

        if line_counter == 37:
            value = line.strip().split(':')[1].strip()
            satellite_str = 'satellite:{0}'.format(value.split(' ')[-2])
            sv_rep_info_list.append(satellite_str)

        # if line_counter == 38:
        #     value = line.strip().split(':')[1].strip()
        #     simRep_str = 'simRep:{0}'.format(value.split(' ')[-2])
        #     sv_rep_info_list.append(simRep_str)

        if line_counter == 39:
            value = line.strip().split(':')[1].strip()
            lowCom_str = 'lowCom:{0}'.format(value.split(' ')[-2])
            sv_rep_info_list.append(lowCom_str)

        line_counter += 1

    return sv_rep_info_list


def parse_rpmask(tbl_file):

    rp_list = process_tbl(tbl_file)
    if float(rp_list[0].split(":")[1]) == 0.0:
        return -1
    else:

        te_stats = {}
        for i in range(1, len(rp_list)):
            te_type = rp_list[i].split(':')[0]
            te_val = float(rp_list[i].split(':')[1])
            te_stats[te_type] = te_val

        non_0_te = []
        for te_type in te_stats.keys():
            te_val = float(te_stats[te_type])
            if te_val != 0:
                non_0_te.append('{0}:{1}'.format(te_type, te_val))

        return non_0_te


def seperate_seq_from_string(input):
    nucleotide = ['A','T','C','G']
    seqStart = 0
    for i in range(len(input)):
        if input[i] in nucleotide:
            seqStart = i
            break

    return input[0:seqStart], input[seqStart:len(input)]

def process_navigable_string(infoTokens):
    # string contains sequence info
    nucleotideContentInfo = infoTokens[17].replace('ACGTcount:', '').strip()

    return nucleotideContentInfo.replace(' ','')

def process_trf_navigable_string(inputStr):
    # trf infos
    trfList = []
    trfInfo = inputStr.strip().split('\n')

    motif = ''
    motifNecleotideInfoStr = ''
    startIdx = 0
    endIdx = 0
    score = 0
    copy_num = 0
    match = 0
    for i in range(len(trfInfo)):
        ele = trfInfo[i]
        if 'Consensus pattern' in ele:
            motif = trfInfo[i + 1]
        if 'Period size' in ele:
            startIdx = i + 2
            ele = ele.replace(":", '').split(' ')
            copy_num = float(ele[ele.index("Copynumber") + 1])
        if 'Statistics' in ele:

            stats_info = trfInfo[i + 1].replace(" ","")
            stats_info_tokens = stats_info.split(",")
            match = int(stats_info_tokens[0].split(":")[1])
            endIdx = i - 3

        if 'Score' in ele:
            score = int(ele.split(' ')[-1])
        if 'ACGTcount' in ele:
            motifNecleotideInfoStr = ele.replace('ACGTcount:', '').strip()

    for j in range(startIdx, endIdx):
        tmp = trfInfo[j].replace(' ', '').split(' ')
        if '' in tmp:
            continue
        nextTmp = trfInfo[j + 1].replace(' ', '').split(' ')
        if '' not in nextTmp:
            pos, seq = seperate_seq_from_string(tmp[0])
            posNext, seqNext = seperate_seq_from_string(nextTmp[0])
            if len(seq) == len(seqNext):
                str = '{0}:{1}'.format(pos, seq)
                trfList.append(str)

    return motif, trfList, score, copy_num, match


def parse_trf(trf_file):


    soup = bs(open(trf_file), 'html.parser')

    contentList = soup.find('pre').contents
    preInfo = contentList[0].replace(' ','').strip().split('\n')

    sv_region_str = ""
    for ele in preInfo:
        if "Sequence" in ele:
            sv_region_str = ele.split(":")[1]
    sv_size = int(sv_region_str.split("-")[2]) - int(sv_region_str.split("-")[1]) + 1

    motif = ''
    max_matches = 0

    seqInfo = contentList[0].strip().split('\n')


    # SV has TRF annotation
    if len(contentList) > 5:
        for i in range(1, len(contentList)):
            ele = contentList[i]
            if isinstance(ele, NavigableString):
                curMotif, curTrfsList, curScore, copies, matches = process_trf_navigable_string(ele)

                if matches > max_matches:
                    motif = curMotif
                    # score = curScore
                    max_matches = matches


        masked_perc = max_matches / float(sv_size)
        masked_perc_str = round(masked_perc, 4) * 100

        if masked_perc_str > 0:
            new_trf_annt = "STRs:{0}".format(round(masked_perc_str, 4))
            if len(motif) >= 7:
                new_trf_annt = "VNTRs:{0}".format(round(masked_perc_str, 4))
            return new_trf_annt
        else:
            return -1

    # # No TRF annotation
    else:
        return -1








if __name__ == '__main__':
    a = parse_rpmask('D:\\Workspace\\test\\rpmask\\chr20-5015546-5015860.fa.tbl')
    print(a)

    sv_trf_annt_dict = parse_trf('D:\\Workspace\\test\\trf\\chr20-5453209-5453361.fa.2.7.7.80.10.50.500.1.txt.html')
    print(sv_trf_annt_dict)
