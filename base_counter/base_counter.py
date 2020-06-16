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
EXIT_INVALID_UMI = 3
EXIT_REPEAT_READ_UMI = 4
PROGRAM_NAME = "base_counter"
UMI_LENGTH = 10


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
    parser.add_argument('--umi', type=str, required=True, help='Path of FASTQ file containing UMIs')
    #parser.add_argument('--umistats', type=str, required=False, help='Path of CSV file to output UMI stats for this sample')
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
            chrom, start, end, label = fields[:4]
            carrier_sample, ref, alt, vaf = parse_target_label(label)
            start = int(start)
            end = int(end)
            result[(chrom, start, end)] = (carrier_sample, ref, alt, vaf)
    return result

class Counts(object):
    def __init__(self):
        self.A = 0
        self.T = 0
        self.G = 0
        self.C = 0
        self.N = 0

    def increment_base_count(self, base):
        if base == 'A':
            self.A += 1
        elif base == 'T':
            self.T += 1
        elif base == 'G':
            self.G += 1
        elif base == 'C':
            self.C += 1
        elif base == 'N':
            self.N += 1
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


def group_runs_by(xs, project, comp):
    result = []
    current = 0
    this_group = []
    while current < len(xs):
        this_item = xs[current]
        this_key = project(this_item)
        if current == 0:
            this_group.append(this_item)
        else:
            previous_key = project(xs[current - 1])
            if comp(this_key, previous_key):
                this_group.append(this_item)
            else:
                result.append(this_group)
                this_group = [this_item]
        current += 1
    if this_group:
        result.append(this_group)
    return result

def adjacent(x, y):
    return x - y <= 1

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

    
def get_pileup_reads(mapqual_threshold, read_stats, samfile, chrom, start, end):
    pileup_reads = []
    # ignore_orphans (bool) – ignore orphans (paired reads that are not in a proper pair). 
    # ignore_overlaps (bool) – If set to True, detect if read pairs overlap and only take the higher quality base. 
    # NOTE: this loop should only happen once for each variant, since they are SNVs and occupy one genomics position
    base_qualities = []
    for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True, stepper='samtools',
                                       ignore_overlaps=False, ignore_orphans=False,
                                       max_depth=1000000000):
        pileupcolumn.set_min_base_quality(0)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                this_base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                base_qualities.append(this_base_qual)
            read_stats.total_reads += 1
            this_map_qual = pileupread.alignment.mapping_quality
            if mapqual_threshold is not None and this_map_qual >= mapqual_threshold:
                read_stats.retained_reads += 1
                pileup_reads.append(pileupread)
    base_qual_mean = statistics.mean(base_qualities)
    base_qual_stdev = statistics.stdev(base_qualities)
    #base_qual_quartiles = statistics.quantiles(base_qualities, n=4)
    # XXX fixme
    base_qual_quartiles = [0, 0, 0]
    base_qual_stats = BaseQualStats(mean=base_qual_mean, stdev=base_qual_stdev, quartiles=base_qual_quartiles)
    return base_qual_stats, pileup_reads


def validate_clusters(counts, clusters):
    counts_umis = sorted(counts.keys())
    clusters_umis = sorted([u for c in clusters for u in c])
    return counts_umis == clusters_umis

# umis_to_reads is a dictionary mapping umi to [reads]
UmiCluster = namedtuple("UmiCluster", ["umis_to_reads"])

'''
read_groups is a list of read groups, where each group is aligned
to the same location (or similar location).

This function clusters the umis for each group.

Result is [UmiCluster]
'''
def cluster_reads_by_umi(umis, read_groups):
    umi_clusterer = UMIClusterer(cluster_method="directional")
    result = []
    for group in read_groups:
        # number of reads for each umi
        umi_counts = defaultdict(int) 
        # reads associated with each umi
        umi_reads = defaultdict(list)
        for read in group:
            this_alignment = read.alignment
            query_name = this_alignment.query_name
            if query_name in umis:
                # UMIClusterer requires a bytestring
               this_umi = umis[query_name].encode()
               umi_counts[this_umi] += 1 
               umi_reads[this_umi].append(read)
        umi_clusters = umi_clusterer(umi_counts, threshold=1)
        for this_cluster in umi_clusters:
            this_umis_to_reads = {}
            for this_umi in this_cluster:
                this_umis_to_reads[this_umi] = umi_reads[this_umi]
            result_cluster = UmiCluster(umis_to_reads = this_umis_to_reads)
            result.append(result_cluster)
    return result 
    

def display_umi_clusters(clusters):
    for this_cluster in clusters:
        umis = this_cluster.umis_to_reads.keys()
        cluster_size = len(umis) 
        if cluster_size >= 8:
            print(30 * '-')
            print(f"size: {cluster_size}")
            for this_umi, this_umi_reads in this_cluster.umis_to_reads.items():
                print(f"\t{this_umi}")
                for read in this_umi_reads:
                    pos = read.alignment.reference_start
                    print(f"{pos} {read.alignment.query_sequence}") 


def get_cluster_dominant_base(bases):
    counts = Counts()
    for base in bases:
        counts.increment_base_count(base)
    sorted_counts = sorted([(counts.A, 'A'), (counts.T, 'T'), (counts.G, 'G'), (counts.C, 'C')])
    max_count, max_base = sorted_counts[-1]
    if max_count > 0:
        return max_base
    else:
        return None

def count_bases(clusters):
    # base counts for collapsed UMI clusters
    cluster_counts = Counts()
    # base counts for all reads at this locus
    all_counts = Counts()
    # this_cluster is a collection of UMIs that cluster together and should be
    # considered the same. We should flatten them into a single base call
    for this_cluster in clusters: 
        this_cluster_bases = []
        for umi, this_reads in this_cluster.umis_to_reads.items():
            for read in this_reads:
                if not read.is_del and not read.is_refskip:
                    this_base = read.alignment.query_sequence[read.query_position].upper()
                    this_cluster_bases.append(this_base)
                    all_counts.increment_base_count(this_base)
        dominant_base = get_cluster_dominant_base(this_cluster_bases) 
        if dominant_base is not None:
            cluster_counts.increment_base_count(dominant_base)
    return cluster_counts, all_counts


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
    
fieldnames = ['chrom', 'pos', 'ref', 'alt', 'sample', 'is_carrier', 'tumour_vaf', 'A', 'T', 'G', 'C', 'N', 'DP', 'vaf', 'A_raw', 'T_raw', 'G_raw', 'C_raw', 'N_raw', 'DP_raw', 'vaf_raw', 'base_qual_mean', 'base_qual_stdev', 'base_qual_quartile_1', 'base_qual_quartile_2', 'base_qual_quartile_3', 'maybe_variant']

def process_bam_file(options, umis, targets):
    sample = options.sample
    samfile = pysam.AlignmentFile(options.bam, "rb" )
    logging.info(f"Processing BAM file from {options.bam} for sample {sample}")
    writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=fieldnames) 
    read_stats = ReadStats()
    writer.writeheader()
    for coord, meta in targets.items():
        chrom, start, end = coord
        #print(30 * '#')
        #print(f"{chrom} {start} {end}")
        # we are only considering SNVs, so the start and end coord should always be exactly 1 bp apart
        assert end - start == 1
        chrom_no_chr = chrom_name_remove_chr(chrom) 
        carrier_sample, ref, alt, tumour_vaf = meta
        is_carrier = sample == carrier_sample
        # get all the reads that align to this position
        base_qual_stats, pileup_reads = get_pileup_reads(options.mapqual, read_stats, samfile, chrom_no_chr, start, end)
        # sort the pileup reads based on their leftmost alignment coordinate
        sorted_pileup_reads = sorted(pileup_reads, key=lambda r: r.alignment.reference_start)
        # group the pileup reads into runs where consecutive reads are at most 1bp apart
        # grouped_pileup_reads_pos = group_runs_by(pileup_reads, lambda r: r.alignment.reference_start, adjacent)
        grouped_pileup_reads_pos = group_runs_by(sorted_pileup_reads, lambda r: r.alignment.reference_start, adjacent)
        umi_clusters = cluster_reads_by_umi(umis, grouped_pileup_reads_pos)
        #display_umi_clusters(umi_clusters)
        cluster_counts, raw_counts = count_bases(umi_clusters)
        depth_uncorrected = len(pileup_reads)
        depth_corrected = cluster_counts.A + cluster_counts.T + cluster_counts.G + cluster_counts.C + cluster_counts.N
        maybe_variant = is_carrier and is_variant(alt, cluster_counts)
        # bed file input coordinates are zero based, but output format is 1-based to be similar with VCF and
        # standard variant nomenclature 
        vaf = compute_vaf(cluster_counts, alt)
        vaf_raw = compute_vaf(raw_counts, alt)
        this_row = {'chrom': chrom, 'pos': start+1, 'ref': ref, 'alt': alt, 'sample': sample, 'is_carrier': is_carrier, 'tumour_vaf': tumour_vaf,
                    'A': cluster_counts.A, 'T': cluster_counts.T, 'G': cluster_counts.G, 'C': cluster_counts.C, 'N': cluster_counts.N, 'DP': depth_corrected, 'vaf': vaf,
                    'A_raw': raw_counts.A, 'T_raw': raw_counts.T, 'G_raw': raw_counts.G, 'C_raw': raw_counts.C, 'N_raw': raw_counts.N, 'DP_raw': depth_uncorrected, 'vaf_raw': vaf_raw,
                    'maybe_variant': maybe_variant, 'base_qual_mean': base_qual_stats.mean, 'base_qual_stdev': base_qual_stats.stdev, 
                    'base_qual_quartile_1': base_qual_stats.quartiles[0], 'base_qual_quartile_2': base_qual_stats.quartiles[1], 'base_qual_quartile_3': base_qual_stats.quartiles[2]}
        writer.writerow(this_row)
    logging.info(f"Total number of reads in input: {read_stats.total_reads}")
    logging.info(f"Number of reads retained after quality filtering: {read_stats.retained_reads}, {read_stats.percent_retained_reads():.2f}%")
    samfile.close()


def is_valid_umi(umi):
    return len(umi) == UMI_LENGTH and valid_dna_bases_set.issuperset(set(umi)) 

def get_umis(options):
    # mapping from read IDs to UMIs
    result = {}
    with pysam.FastxFile(options.umi) as file:
        for entry in file:
            if entry.name not in result:
                if is_valid_umi(entry.sequence):
                   result[entry.name] = entry.sequence
                else:
                   exit_with_error(f"Invalid UMI: {entry.sequence}", EXIT_INVALID_UMI) 
            else:
                exit_with_error(f"Read ID repeated in UMI file: {entry.name}", EXIT_REPEAT_READ_UMI)
    return result

def umi_stats(options, umis):
    total_umis = 0
    unique_umis = set()
    for read, umi in umis.items():
        total_umis += 1
        unique_umis.add(umi)
    num_unqiue_umis = len(unique_umis)
    num_umis_with_n = sum([1 for u in unique_umis if u.count("N") > 0])
    if options.umistats:
        with open(options.umistats, "w") as file:
            print(f"{total_umis},{num_unqiue_umis},{num_umis_with_n}", file=file)
        


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
    # umis is a dictionary mapping read ID to UMI
    umis = get_umis(options)
    #umi_stats(options, umis)
    targets = get_targets(options) 
    process_bam_file(options, umis, targets)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
