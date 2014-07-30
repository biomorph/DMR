__author__ = 'ravi'
import pybedtools
from pybedtools import BedTool
import sys
import argparse
import os


# This function makes window_sized genome windows bed
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


# This function takes the bismark coverage files and adds up the number of methylated and unmethylated C's over each window and filters windows with missing data and with coverage < minCov
def get_filtered_window_methylation(cov_bed, chrom_bins_bed, minCov, window):
    window_methylation_bed = chrom_bins_bed
    windowed_bed_file_name = str(window) + 'window-' + str(minCov) + 'minCov-' + cov_bed.split('.')[0] + '.bed'
    bed = BedTool(cov_bed)
    sorted_bed = bed.sort()
    window_methylation_bed = window_methylation_bed.map(b=sorted_bed,c='5,6')
    window_methylation_bed.filter(lambda d: d[3] != '.' and d[4] != '.').saveas(windowed_bed_file_name)
    return windowed_bed_file_name

# Filters the windowed beds to common set of windows among all samples. This will make plotting and statistics easier
def get_common_windows (bed_list):
    a = BedTool()
    common_bed = a.multi_intersect(i=bed_list).filter(lambda z:int(z[3]) == len(bed_list)).saveas('common.bed')
    for bed in bed_list:
        BedTool(bed).intersect(b=common_bed).saveas(bed+'-common')



'''
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
'''

parser = argparse.ArgumentParser()

parser.add_argument('-data', type = str, help = 'Input matched replicate coverage files separated by comma and diff condition files separated by semicolon, ex. 1A,1B:2A,2B etc')

parser.add_argument('-groups', type = str, help = 'Specify the grouping of the samples your entered separated by colon , ctrl:treatment etc')

parser.add_argument('-window', type = int, help = 'Specify window size for binning methylation coverage')

parser.add_argument('-chromSizes', type = str, help = 'chromSizes file with chr# size')

parser.add_argument('-minCov', type = int, help = 'minimum number of cytosines covered, defaults to 10')

args = parser.parse_args()


file_array = args.data.strip().split(',')
chrom_bins_bed = make_chrom_bins_bed(args.chromSizes,args.window)

minCov = 10
if args.minCov: minCov = args.minCov
filtered_bed_array = []

for file in file_array:
    filtered_bed = get_filtered_window_methylation (file, chrom_bins_bed, minCov, args.window)
    filtered_bed_array.append(filtered_bed)

get_common_windows(filtered_bed_array)









