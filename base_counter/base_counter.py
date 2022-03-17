'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 19 Sep 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Count frequency of DNA bases at genomic positions.
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pysam
import pathlib
import csv
from umi_tools import UMIClusterer
from collections import defaultdict, namedtuple
import statistics


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "base_counter"
DEFAULT_WINDOW = 1
MAX_READ_DEPTH = 1000000000


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Count frequency of DNA bases at genomic positions'
    parser = ArgumentParser(description=description)
    parser.add_argument('--targets', type=str, required=True,
        help='Bed file coordinates of genomic targets to consider')
    parser.add_argument('--sample', type=str, required=True,
        help='Sample ID')
    parser.add_argument('--window', type=int, required=False, default=DEFAULT_WINDOW, help='Window of loci centered on target for base counting (default: %(default)s)')
    parser.add_argument('--mapqual', type=float, required=False, help='Minimum mapping quality threshold for reads')
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('bam',
                        metavar='BAM_FILE',
                        type=str,
                        help='Input BAM file')
    return parser.parse_args()

# format: CHR10_10023321_A_CMHS234_T_0.1286

def parse_target_label(label):
    fields = label.split('_')
    _chrom, _start, ref, carrier_sample, alt, vaf = fields[:6]
    vaf = float(vaf)
    return (carrier_sample, ref, alt, vaf)


def get_targets(options):
    result = {} 
    with open(options.targets) as targets_file:
        for line in targets_file:
            fields = line.split()
            chrom, start, _end, label = fields[:4]
            carrier_sample, ref, alt, vaf = parse_target_label(label)
            start = int(start)
            #end = int(end)
            result[(chrom, start)] = (carrier_sample, ref, alt, vaf)
    return result

class Counts(object):
    def __init__(self):
        self.A = 0
        self.A_strand = {'fwd': 0, 'rev': 0} 
        self.T = 0
        self.T_strand = {'fwd': 0, 'rev': 0} 
        self.G = 0
        self.G_strand = {'fwd': 0, 'rev': 0} 
        self.C = 0
        self.C_strand = {'fwd': 0, 'rev': 0} 
        self.N = 0
        self.N_strand = {'fwd': 0, 'rev': 0} 
        self.D = 0    # this counts absent bases, such as deletions and refskips
        self.D_strand = {'fwd': 0, 'rev': 0} 
        

    def increment_base_count(self, base, strand):
        if base == 'A':
            self.A += 1
            self.A_strand[strand] += 1
        elif base == 'T':
            self.T += 1
            self.T_strand[strand] += 1
        elif base == 'G':
            self.G += 1
            self.G_strand[strand] += 1
        elif base == 'C':
            self.C += 1
            self.C_strand[strand] += 1
        elif base == 'N':
            self.N += 1
            self.N_strand[strand] += 1
        elif base == 'D':
            self.D += 1
            self.D_strand[strand] += 1
        else:
            exit(f"Unrecognised base: {base}")

    def __str__(self):
        return f"A:{self.A}, T:{self.T}, G:{self.G}, C:{self.C}"

def is_variant(alt, counts):
    if alt == 'A':
        return counts.A > 0
    elif alt == 'T':
        return counts.T > 0
    elif alt == 'G':
        return counts.G > 0
    elif alt == 'C':
        return counts.C > 0
    else:
        return False


VALID_DNA_BASES = "ATGCN"
valid_dna_bases_set = set(VALID_DNA_BASES)

def chrom_name_remove_chr(chrom):
    chrom_no_chr = chrom
    if chrom_no_chr.startswith('chr'):
        chrom_no_chr = chrom_no_chr[3:] 
    return chrom_no_chr


class ReadStats:
    def __init__(self):
        self.total_reads = 0
        self.retained_reads = 0

    def percent_retained_reads(self):
        if self.total_reads > 0:
            return (self.retained_reads / self.total_reads) * 100
        else:
            return 0

BaseQualStats = namedtuple("BaseQualStats", ["mean", "stdev", "quartiles"])


def get_pileup_stats(mapqual_threshold, read_stats, samfile, chrom, pos):
    base_read_lens = {'A': [], 'T': [], 'G': [], 'C': [], 'N': [], 'D': [] } 
    base_counts = Counts()
    # ignore_orphans (bool) – ignore orphans (paired reads that are not in a proper pair). 
    # ignore_overlaps (bool) – If set to True, detect if read pairs overlap and only take the higher quality base. 
    # NOTE: this loop should only happen once for each variant, since they are SNVs and occupy one genomics position
    base_qualities = []
    strand = None
    for pileupcolumn in samfile.pileup(chrom, pos, pos+1, truncate=True, stepper='samtools',
                                       ignore_overlaps=False, ignore_orphans=False,
                                       max_depth=MAX_READ_DEPTH):
        pileupcolumn.set_min_base_quality(0)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                this_base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                base_qualities.append(this_base_qual)
            read_stats.total_reads += 1
            this_map_qual = pileupread.alignment.mapping_quality
            if mapqual_threshold is None or this_map_qual >= mapqual_threshold:
                read_stats.retained_reads += 1
                # thus_cluster_read_lens.append(read.alignment.query_alignment_length)
                fragment_length = abs(pileupread.alignment.template_length)
                if pileupread.alignment.is_reverse:
                    strand = 'rev'
                else:
                    strand = 'fwd' 
                if not pileupread.is_del and not pileupread.is_refskip:
                    this_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                else:
                    # deletion and refskip get special 'D' (deletion) base
                    this_base = 'D'
                base_counts.increment_base_count(this_base, strand)
                base_read_lens[this_base].append(fragment_length)
    try:
        base_qual_mean = statistics.mean(base_qualities)
        base_qual_stdev = statistics.stdev(base_qualities)
        base_qual_quartiles = statistics.quantiles(base_qualities, n=4)
    except statistics.StatisticsError:
        base_qual_mean = None
        base_qual_stdev = None
        base_qual_quartiles = [None, None, None] 
    base_qual_stats = BaseQualStats(mean=base_qual_mean, stdev=base_qual_stdev, quartiles=base_qual_quartiles)
    base_read_len_means = {} 
    base_read_len_means['A'] = statistics.mean(base_read_lens['A']) if base_read_lens['A'] else ''
    base_read_len_means['T'] = statistics.mean(base_read_lens['T']) if base_read_lens['T'] else ''
    base_read_len_means['G'] = statistics.mean(base_read_lens['G']) if base_read_lens['G'] else ''
    base_read_len_means['C'] = statistics.mean(base_read_lens['C']) if base_read_lens['C'] else ''
    base_read_len_means['N'] = statistics.mean(base_read_lens['N']) if base_read_lens['N'] else ''
    return base_qual_stats, base_read_len_means, base_counts 


def compute_vaf(counts, allele):
    if allele == 'A':
        num = counts.A
    elif allele == 'T':
        num = counts.T
    elif allele == 'G':
        num = counts.G
    elif allele == 'C':
        num = counts.C
    elif allele == 'N':
        num = counts.N
    else:
        num = 0
    total = counts.A + counts.T + counts.G + counts.C + counts.N
    if total > 0:
       return num / total
    else:
       return ''

def get_noise(counts, ref, alt):
    noise_read_count = 0
    ref_alt = [ref, alt]
    for base in ['A', 'T', 'G', 'C', 'N']:
        if base not in ref_alt:
            noise_read_count += getattr(counts, base)
    return noise_read_count
    
fieldnames = ['chrom', 'pos', 'ref', 'alt', 'is_target', 'sample', 'is_carrier', 'tumour_vaf', 'A', 'T', 'G', 'C', 'N', 'DP', 'vaf', 'base_qual_mean', 'base_qual_stdev', 'base_qual_quartile_1', 'base_qual_quartile_2', 'base_qual_quartile_3', 'A_read_len_mean', 'T_read_len_mean', 'G_read_len_mean', 'C_read_len_mean', 'N_read_len_mean', 'A_fwd', 'A_rev', 'T_fwd', 'T_rev', 'C_fwd', 'C_rev', 'G_fwd', 'G_rev', 'N_fwd', 'N_rev', 'maybe_variant']

def process_bam_file(options, targets):
    sample = options.sample
    samfile = pysam.AlignmentFile(options.bam, "rb" )
    logging.info(f"Processing BAM file from {options.bam} for sample {sample}")
    writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=fieldnames) 
    read_stats = ReadStats()
    writer.writeheader()

    seen_coords = set()

    for coord, meta in targets.items():
        chrom, start = coord
        chrom_no_chr = chrom_name_remove_chr(chrom) 
        half_window = options.window // 2
        window_start = start - half_window
        window_end = start + (options.window - half_window)

        assert(options.window == (window_end - window_start))

        for pos in range(window_start, window_end):

            is_target = pos == start

            if is_target or ((chrom, pos) not in targets and (chrom, pos) not in seen_coords):

                seen_coords.add((chrom, pos))

                # get all the reads that align to this position
                base_qual_stats, base_read_lens, base_counts = get_pileup_stats(options.mapqual, read_stats, samfile, chrom_no_chr, pos)
                locus_depth = base_counts.A + base_counts.T + base_counts.G + base_counts.C + base_counts.N

                if is_target:
                    carrier_sample, ref, alt, tumour_vaf = meta
                    vaf = compute_vaf(base_counts, alt)
                    is_carrier = sample == carrier_sample
                    is_target = True
                    maybe_variant = is_carrier and is_variant(alt, base_counts)
                else:
                    maybe_variant = ''
                    vaf = ''
                    vaf_raw = ''
                    ref = ''
                    alt = ''
                    tumour_vaf = ''
                    is_carrier = ''
                    is_target = False
                # bed file input coordinates are zero based, but output format is 1-based to be similar with VCF and
                # standard variant nomenclature 
                this_row = {'chrom': chrom, 'pos': pos+1, 'ref': ref, 'alt': alt, 'is_target': is_target, 'sample': sample, 'is_carrier': is_carrier, 'tumour_vaf': tumour_vaf,
                            'A': base_counts.A, 'T': base_counts.T, 'G': base_counts.G, 'C': base_counts.C, 'N': base_counts.N, 'DP': locus_depth, 'vaf': vaf,
                            'maybe_variant': maybe_variant, 'base_qual_mean': base_qual_stats.mean, 'base_qual_stdev': base_qual_stats.stdev, 
                            'base_qual_quartile_1': base_qual_stats.quartiles[0], 'base_qual_quartile_2': base_qual_stats.quartiles[1], 'base_qual_quartile_3': base_qual_stats.quartiles[2],
                            'A_read_len_mean': base_read_lens['A'], 'T_read_len_mean': base_read_lens['T'], 'G_read_len_mean': base_read_lens['G'], 'C_read_len_mean': base_read_lens['C'], 'N_read_len_mean': base_read_lens['N'],
                            'A_fwd': base_counts.A_strand['fwd'], 'A_rev': base_counts.A_strand['rev'],
                            'T_fwd': base_counts.T_strand['fwd'], 'T_rev': base_counts.T_strand['rev'],
                            'G_fwd': base_counts.G_strand['fwd'], 'G_rev': base_counts.G_strand['rev'],
                            'C_fwd': base_counts.C_strand['fwd'], 'C_rev': base_counts.C_strand['rev'],
                            'N_fwd': base_counts.N_strand['fwd'], 'N_rev': base_counts.N_strand['rev'],
                            }
                writer.writerow(this_row)
            #else:
            #    print(f"skipping duplicate pos {chrom_no_chr} {pos}")
    logging.info(f"Total number of reads in input: {read_stats.total_reads}")
    logging.info(f"Number of reads retained after quality filtering: {read_stats.retained_reads}, {read_stats.percent_retained_reads():.2f}%")
    samfile.close()


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))



def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    targets = get_targets(options) 
    process_bam_file(options, targets)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
