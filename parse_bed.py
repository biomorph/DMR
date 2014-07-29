__author__ = 'ravi'
import pybedtools
from pybedtools import BedTool
import sys
import argparse
import os

#methyl_cov = BedTool('../../bismark_cov/SRR1036970_1_val_1.fq.gz_bismark_pe.deduplicated.bismark.cov')


def make_chrom_bins_bed(chromSizes_file, window_size):
    chrom_bin_bed_file = str(window_size) +  'windowGenome.bed'
    if os.path.exists(chrom_bin_bed_file): return BedTool(chrom_bin_bed_file)
    else:
        a = pybedtools.example_bedtool('a.bed')
        chrom_windows_bed = a.window_maker(genome='hg19',w=window_size)
        chrom_windows_bed.saveas(chrom_bin_bed_file)
        return chrom_windows_bed

def methyl_stats_for_bin(bin,window):
    methylated_Cs = 0
    unmethylated_Cs = 0
    return methylated_Cs,unmethylated_Cs


def get_bin_intervals(feature,bin,window):
    return feature.chrom == bin[0] and feature.start >= bin[1] and feature.start <= bin[1] + window


def get_windowed_methylation_bed(replicates, chrom_bins_bed, minCov, window):
    window_methylation_bed = chrom_bins_bed
    windowed_bed_file_name = str(window) + 'window' + str(minCov) + 'minCov'
    for replicate in replicates:
        replicate_bed = BedTool(replicate)
        sorted_bed = replicate_bed.sort()
        window_methylation_bed = window_methylation_bed.map(b=sorted_bed,c='5,6')
        windowed_bed_file_name += '-' + replicate.split('.')[0]
    window_methylation_bed.filter(lambda d: d[3] != '.' and d[4] != '.' and int(d[3])+int(d[4])>10).saveas(windowed_bed_file_name)

def group_input_files(files, groups):
    file_group_dict = {}
    file_groups = files.split(':')
    groups = groups.split(':')
    if len(file_groups) != len(groups):
        sys.exit('Number of file groups and number of conditions do not match')
    for i, group in enumerate(groups):
        if group in file_group_dict:
            file_group_dict[group] += [file_groups[i].split(',')]
        else:
            file_group_dict[group] = [file_groups[i].split(',')]
    return file_group_dict


parser = argparse.ArgumentParser()

parser.add_argument('-data', type = str, help = 'Input matched replicate coverage files separated by comma and diff condition files separated by semicolon, ex. 1A,1B:2A,2B etc')

parser.add_argument('-groups', type = str, help = 'Specify the grouping of the samples your entered separated by colon , ctrl:treatment etc')

parser.add_argument('-window', type = int, help = 'Specify window size for binning methylation coverage')

parser.add_argument('-chromSizes', type = str, help = 'chromSizes file with chr# size')

parser.add_argument('-minCov', type = int, help = 'minimum number of cytosines covered, defaults to 10')

args = parser.parse_args()


condition_dict = group_input_files(args.data, args.groups)
chrom_bins_bed = make_chrom_bins_bed(args.chromSizes,args.window)

minCov = 10
if args.minCov: minCov = args.minCov

for condition in condition_dict:
    biological_replicates = condition_dict[condition]
    for tech_replicate in biological_replicates:
        get_windowed_methylation_bed(tech_replicate, chrom_bins_bed, minCov, args.window)







