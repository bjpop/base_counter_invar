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
    parser.add_argument('--umistats', type=str, required=False, help='Path of CSV file to output UMI stats for this sample')
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
        self.UMIs = set()

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


VALID_DNA_BASES = "ATGCN"
valid_dna_bases_set = set(VALID_DNA_BASES)

def process_bam_file(options, targets):
    total_reads = 0 
    fieldnames = ['chrom', 'pos', 'ref', 'alt', 'sample', 'is carrier', 'A', 'T', 'G', 'C', 'N', 'DP', 'UMI', 'read/UMI', 'maybe variant']
    writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=fieldnames) 
    writer.writeheader()
    sample = options.sample
    samfile = pysam.AlignmentFile(options.bam, "rb" )
    logging.info(f"Processing BAM file from {options.bam} for sample {sample}")
    for coord, meta in targets.items():
        chrom, start, end = coord
        chrom_no_chr = chrom
        if chrom_no_chr.startswith('chr'):
            chrom_no_chr = chrom_no_chr[3:]
        carrier_sample, ref, alt, vaf = meta
        is_carrier = sample == carrier_sample
      
        # ignore_orphans (bool) – ignore orphans (paired reads that are not in a proper pair). 
        # ignore_overlaps (bool) – If set to True, detect if read pairs overlap and only take the higher quality base. 
        for pileupcolumn in samfile.pileup(chrom_no_chr, start, end, truncate=True, stepper='samtools',
                                           ignore_overlaps=False, ignore_orphans=False,
                                           max_depth=1000000000):
            this_pos_zero_based = pileupcolumn.pos
            this_pos_one_based = this_pos_zero_based + 1
            counts = Counts()
            # the maximum number of reads aligning to this position before filtering
            unfiltered_coverage = pileupcolumn.nsegments
            # collect all bases in the current column, and assign them to either read1s or read2s
            for pileupread in pileupcolumn.pileups:
                this_alignment = pileupread.alignment
                query_name = this_alignment.query_name
                mapping_quality = this_alignment.mapping_quality
                alignment_length = this_alignment.query_alignment_length
                is_read1 = this_alignment.is_read1
                is_read2 = this_alignment.is_read2
                if not pileupread.is_del and not pileupread.is_refskip:
                    this_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    counts.increment_base_count(this_base)
            total_counts = counts.A + counts.T + counts.G + counts.C + counts.N
            maybe_variant = is_carrier and is_variant(alt, counts)
            output_row = {'chrom': chrom, 'pos': this_pos_one_based, 'ref': ref, 'alt': alt, 'sample': sample, 'is carrier': is_carrier, 'A': counts.A, 'T': counts.T, 'G': counts.G, 'C': counts.C, 'N': counts.N, 'DP': total_counts, 'UMI': 0, 'read/UMI': 0, 'maybe variant': maybe_variant }
            writer.writerow(output_row) 
    samfile.close()


def is_valid_umi(umi):
    return len(umi) == UMI_LENGTH and valid_dna_bases_set.issuperset(set(umi)) 

def get_umis(options):
    # mapping from read IDs to UMIs
    result = {}
    with pysam.FastxFile(options.umi) as fh:
        for entry in fh:
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
    for read,umi in umis.items():
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
    umis = get_umis(options)
    umi_stats(options, umis)
    exit(0)
    targets = get_targets(options) 
    process_bam_file(options, targets)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
