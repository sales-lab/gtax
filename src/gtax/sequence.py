import argparse
import gzip
import os
import random
import uuid
from functools import partial
from multiprocessing import Pool

import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gtax.taxonomy import Taxonomy


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def load_reference_transcripts(reference_file, total_sequences, short_seq_len):
    transcripts = []
    sum_length = 0
    with gzip.open(reference_file, 'rt') as trans_stream:
        for r in SeqIO.parse(trans_stream, "fasta"):
            transcripts.append(r)
            sum_length += len(r.seq)
    print('{} transcripts loaded'.format(len(transcripts)))
    number_seq_factor = total_sequences * short_seq_len / sum_length
    print('Using factor {:.2f} to generate {} number of short sequences'.
          format(number_seq_factor, total_sequences))
    return transcripts, number_seq_factor


def fasta_paired_end_generator(records, number_seq_factor, short_seq_len, min_seq_size):
    pre = uuid.uuid4()
    print("Opening tmp files: temp_{0}_1.fa.gz and temp_{0}_1.fa.gz".format(pre))
    with gzip.open('temp_{}_1.fa.gz'.format(pre), 'wt') as fout_1, \
            gzip.open('temp_{}_2.fa.gz'.format(pre), 'wt') as fout_2:
        for r in records:
            length = len(r.seq)
            no_reads = int(number_seq_factor * length / short_seq_len)
            i = 0
            iterations = no_reads * 100
            while i < no_reads:
                iterations -= 1
                insert_size = random.randrange(200, 350)
                pos = random.randrange(length)
                f_start = pos - short_seq_len
                if f_start < 0:
                    f_start = 0
                f_read = r.seq[f_start:pos]
                r_start = pos + insert_size
                r_read = r.seq[r_start:r_start + short_seq_len].reverse_complement()
                if 'N' not in f_read and 'N' not in r_read and \
                        len(f_read) > min_seq_size and \
                        len(r_read) > min_seq_size:
                    f_rec = SeqRecord(Seq(f_read),
                                      id='{}_{}/1'.format(r.id, f_start), description='')
                    fout_1.write(f_rec.format("fasta"))

                    r_rec = SeqRecord(Seq(r_read),
                                      id='{}_{}/2'.format(r.id, f_start), description='')
                    fout_2.write(r_rec.format("fasta"))
                    i += 1
                if iterations <= 0:
                    break
    return pre


def fasta_single_end_generator(records, number_seq_factor, short_seq_len, min_seq_size):
    pre = uuid.uuid4()
    with gzip.open('temp_{}.fa.gz'.format(pre), 'wt') as fout:
        for r in records:
            length = len(r.seq)
            no_reads = int(number_seq_factor * length / short_seq_len)
            i = 0
            iterations = no_reads * 100
            while i < no_reads:
                iterations -= 1
                pos = random.randrange(length)
                f_start = pos - short_seq_len
                if f_start < 0:
                    f_start = 0
                f_read = r.seq[f_start:pos]
                if 'N' not in f_read and len(f_read) > min_seq_size:
                    f_rec = SeqRecord(Seq(f_read), id='{}_{}'.format(r.id, f_start), description='')
                    fout.write(f_rec.format("fasta"))
                    i += 1
                if iterations <= 0:
                    break
    return pre


def create_random_paired_end_files(reference_file, output_prefix,
                                   number_sequences, short_seq_len, min_seq_size,
                                   number_files, threads, split_file_by):
    transcripts, number_seq_factor = load_reference_transcripts(reference_file,
                                                                number_sequences, short_seq_len)

    for i in range(1, number_files + 1):
        with Pool(processes=threads) as p:
            results = p.map(partial(fasta_paired_end_generator,
                                    number_seq_factor=number_seq_factor,
                                    short_seq_len=short_seq_len,
                                    min_seq_size=min_seq_size),
                            [d for d in list(chunks(transcripts, 1000))])

        # p = Pool(processes=8)
        # results_cont = p.map(contamination_read_generator, taxonomy_groups_reads.keys())
        # p.close()

        print('Creating file {}_{}'.format(output_prefix, i))
        with gzip.open('{}_{}_1.fa.gz'.format(output_prefix, i), 'wt') as fout_1, \
                gzip.open('{}_{}_2.fa.gz'.format(output_prefix, i), 'wt') as fout_2, \
                gzip.open('{}_{}_cont_1.fa.gz'.format(output_prefix, i), 'wt') as fout2_1, \
                gzip.open('{}_{}_cont_2.fa.gz'.format(output_prefix, i), 'wt') as fout2_2:
            for pre in results:
                with gzip.open('temp_{}_1.fa.gz'.format(pre), 'rt') as fin:
                    for r in SeqIO.parse(fin, "fasta"):
                        fout_1.write(r.format('fasta'))
                        fout2_1.write(r.format('fasta'))
                os.remove('temp_{}_1.fa.gz'.format(pre))
                with gzip.open('temp_{}_2.fa.gz'.format(pre), 'rt') as fin:
                    for r in SeqIO.parse(fin, "fasta"):
                        fout_2.write(r.format('fasta'))
                        fout2_2.write(r.format('fasta'))
                os.remove('temp_{}_2.fa.gz'.format(pre))
            # for k in taxonomy_groups_reads:
            #     with gzip.open('{}_1.fa.gz'.format(k), 'rt') as fin:
            #         for r in SeqIO.parse(fin, "fasta"):
            #             fout2_1.write(r.format('fasta'))
            #     os.remove('{}_1.fa.gz'.format(k))
            #     with gzip.open('{}_2.fa.gz'.format(k), 'rt') as fin:
            #         for r in SeqIO.parse(fin, "fasta"):
            #             fout2_2.write(r.format('fasta'))
            #     os.remove('{}_2.fa.gz'.format(k))


def contamination_fasta_single_end_generator(fasta_file, total_read_seq, short_seq_len, min_seq_size):
    pre = uuid.uuid4()
    with open(fasta_file) as fin, gzip.open('temp_{}.fa.gz'.format(pre), 'wt') as fout:
        for r in SeqIO.parse(fin, "fasta"):
            length = len(r.seq)
            i = 0
            iterations = total_read_seq * 100
            while i < total_read_seq:
                iterations -= 1
                pos = random.randrange(length)
                f_start = pos - short_seq_len
                if f_start < 0:
                    f_start = 0
                f_read = r.seq[f_start:pos]
                if 'N' not in f_read and len(f_read) > min_seq_size:
                    f_rec = SeqRecord(Seq(f_read), id='{}_{}'.format(r.id, f_start), description='')
                    fout.write(f_rec.format("fasta"))
                    i += 1
                if iterations <= 0:
                    break
    return pre


def process_by_taxonomy_group(tax_group, gtax_fasta_dir, total_reads, short_seq_len, min_seq_size):
    print('Processing tax group: {}'.format(tax_group))
    fasta_file = os.path.join(gtax_fasta_dir, '{}.fsa'.format(tax_group))
    idx_file = os.path.join(gtax_fasta_dir, '{}.idx'.format(tax_group))
    df = pandas.read_csv(idx_file, sep='\t', header=None)
    total_seq = len(df)
    total_read_seq = int(total_reads / total_seq)
    if total_read_seq < 1:
        total_read_seq = 1
    return contamination_fasta_single_end_generator(fasta_file, total_read_seq * total_seq, short_seq_len, min_seq_size)


def create_random_single_end_files(reference_file, output_prefix,
                                   number_sequences, short_seq_len, min_seq_size,
                                   number_files, threads, split_file_by,
                                   taxonomy, foreign_contamination_percent,
                                   gtax_fasta_dir, taxonomy_group):
    print("Creating single-end short sequences")
    transcripts, number_seq_factor = load_reference_transcripts(reference_file,
                                                                number_sequences, short_seq_len)

    for i in range(1, number_files + 1):
        print('Creating reference sequences')
        with Pool(processes=threads) as p:
            results = p.map(partial(fasta_single_end_generator,
                                    number_seq_factor=number_seq_factor,
                                    short_seq_len=short_seq_len,
                                    min_seq_size=min_seq_size),
                            [d for d in list(chunks(transcripts, 1000))])
        print('Creating foreign sequences')
        with Pool(processes=threads) as p:
            results_cont = p.map(partial(process_by_taxonomy_group,
                                         gtax_fasta_dir=gtax_fasta_dir,
                                         total_reads=number_sequences * foreign_contamination_percent,
                                         short_seq_len=short_seq_len,
                                         min_seq_size=min_seq_size),
                                 [d for d in taxonomy.taxonomy_groups if d != taxonomy_group])
        print('Printing final files')
        if split_file_by == 0:
            print('Creating file {}_{}'.format(output_prefix, i))
            with gzip.open('{}_{}.fa.gz'.format(output_prefix, i), 'wt') as fout, \
                    gzip.open('{}_{}_cont.fa.gz'.format(output_prefix, i), 'wt') as fout2:
                for pre in results:
                    with gzip.open('temp_{}.fa.gz'.format(pre), 'rt') as fin:
                        for r in SeqIO.parse(fin, "fasta"):
                            fout.write(r.format('fasta'))
                            fout2.write(r.format('fasta'))
                    os.remove('temp_{}.fa.gz'.format(pre))
                for pre in results_cont:
                    with gzip.open('temp_{}.fa.gz'.format(pre), 'rt') as fin:
                        for r in SeqIO.parse(fin, "fasta"):
                            fout2.write(r.format('fasta'))
                os.remove('temp_{}.fa.gz'.format(pre))
        else:
            count = 0
            suffix = 1
            fout = gzip.open('{}_{}_{}.fa.gz'.format(output_prefix, i, suffix), 'wt')
            fout2 = gzip.open('{}_{}_{}_cont.fa.gz'.format(output_prefix, i, suffix), 'wt')
            for pre in results:
                with gzip.open('temp_{}.fa.gz'.format(pre), 'rt') as fin:
                    for r in SeqIO.parse(fin, "fasta"):
                        if count == split_file_by:
                            fout.close()
                            fout2.close()
                            count = 0
                            suffix += 1
                            fout = gzip.open('{}_{}_{}.fa.gz'.format(output_prefix, i, suffix), 'wt')
                            fout2 = gzip.open('{}_{}_{}_cont.fa.gz'.format(output_prefix, i, suffix), 'wt')
                        count += 1
                        fout.write(r.format('fasta'))
                        fout2.write(r.format('fasta'))
                os.remove('temp_{}.fa.gz'.format(pre))
            fout.close()
            for pre in results_cont:
                with gzip.open('temp_{}.fa.gz'.format(pre), 'rt') as fin:
                    for r in SeqIO.parse(fin, "fasta"):
                        if count == split_file_by:
                            fout2.close()
                            count = 0
                            suffix += 1
                            fout2 = gzip.open('{}_{}_{}_cont.fa.gz'.format(output_prefix, i, suffix), 'wt')
                        count += 1
                        fout2.write(r.format('fasta'))
            os.remove('temp_{}.fa.gz'.format(pre))
            fout2.close()


def create_random_short_sequences():
    parser = argparse.ArgumentParser()

    parser.add_argument('--reference_file', help='Reference Transcriptome gzip file',
                        required=True)
    parser.add_argument('--taxonomy_group', help='Taxonomy group name for the target organism',
                        required=True)
    parser.add_argument('--gtax_fasta_dir', help='Gtax database directory contaning FASTA and IDX files',
                        required=True)
    parser.add_argument('--output_prefix', help='Output file prefix', required=True)
    parser.add_argument('--number_files', help='Number of random file to create',
                        default=4, type=int, required=False)
    parser.add_argument('--threads', help='Number of threads',
                        type=int, required=True)
    parser.add_argument('--foreign_contamination_percent',
                        help='Foreign Contamination percent by taxonomy group',
                        type=float, required=True)
    parser.add_argument('--short_seq_len', help='Length of the short read',
                        default=100, type=int, required=False)
    parser.add_argument('--number_sequences',
                        help='The number of short reads sequences to generate ',
                        type=int, required=True)
    parser.add_argument('--split_file_by',
                        help='Generate output files with this number of sequences ',
                        default=0,
                        type=int, required=False)
    parser.add_argument('--paired_end',
                        help='Generate paired-end mate sequence',
                        action="store_true", default=False, required=False)
    args = parser.parse_args()
    min_seq_size = args.short_seq_len * 75 / 100

    taxonomy = Taxonomy(tax_pickle_file=os.path.join(args.gtax_fasta_dir,
                                                     'taxonomy.pickle'),
                        group_pickle_file=os.path.join(args.gtax_fasta_dir,
                                                       'taxonomy_groups.pickle'))

    if args.paired_end:
        create_random_paired_end_files(args.reference_file,
                                       args.output_prefix,
                                       args.number_sequences,
                                       args.short_seq_len,
                                       min_seq_size,
                                       args.number_files,
                                       args.threads,
                                       args.split_file_by)
    else:
        create_random_single_end_files(args.reference_file,
                                       args.output_prefix,
                                       args.number_sequences,
                                       args.short_seq_len,
                                       min_seq_size,
                                       args.number_files,
                                       args.threads,
                                       args.split_file_by,
                                       taxonomy,
                                       args.foreign_contamination_percent,
                                       args.gtax_fasta_dir,
                                       args.taxonomy_group)


if __name__ == "__main__":
    create_random_short_sequences()
