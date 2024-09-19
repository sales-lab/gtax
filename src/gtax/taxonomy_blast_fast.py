import os
import gzip
import loky
import pandas
import argparse
from Bio import SeqIO
from functools import partial
from gtax.taxonomy import Taxonomy
from tqdm import tqdm


def transcript_contamination(filename, blast_columns, tax_ids):
    blast_df = pandas.read_csv(filename, sep='\t', header=None,
                               names=blast_columns.split(' '),
                               low_memory=False)
    groupby_qseqid = blast_df.groupby('qseqid')
    data = []
    for g in groupby_qseqid.groups.keys():
        df = groupby_qseqid.get_group(g)
        df = df[df['evalue'] == df['evalue'].min()]
        if not all(elem in tax_ids for elem in df['staxid'].unique()):
            data.append(g)

    return data


def taxonomy_blast_fast():
    parser = argparse.ArgumentParser(prog='taxonomy_blast_fast',
                                     description='This tools process BLAST output to find contamination.')

    parser.add_argument('--threads', help='No. of threads',
                        required=True)
    parser.add_argument('--prefix', help='Prefix for output files',
                        required=True)
    parser.add_argument('--fastq', help='Reference Transcriptome FASTQ file',
                        required=True)
    parser.add_argument('--taxid', help='Parent TaxID to use as filter',
                        required=True)
    parser.add_argument('--blastdir', help='Directory with BLAST gzip results *.out.gz',
                        required=True)
    parser.add_argument('--blast_columns', help='BLAST -outfmt columns. Must include '
                                                'qseqid, saccver, evalue, qcovs, staxid',
                        required=True)
    args = parser.parse_args()

    tax_ids = [int(i) for i in Taxonomy().successors(args.taxid)]
    print('{} taxonomies IDs in the list'.format(len(tax_ids)))

    records = {}
    with gzip.open(args.fastq, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            records[record.id] = record
    print(f"{len(records)} sequences loaded")

    files = [os.path.join(args.blastdir, f)
             for f in os.listdir(args.blastdir) if f.endswith('.out.gz')]

    contamination = 0
    cont_ids = set()

    executor = loky.get_reusable_executor(max_workers=int(args.threads))
    results = executor.map(
        partial(transcript_contamination,
                blast_columns=args.blast_columns,
                tax_ids=tax_ids),
        files
    )
    results = tqdm(results, total=len(files))
    for idx, data in enumerate(results):
        if (idx % 10) == 0:
            results.refresh()

        contamination += len(data)
        cont_ids |= set(data)

    with open('{}_clean.fastq'.format(args.prefix), 'w') as f_fsa:
        for r in records:
            if r.id not in cont_ids:
                f_fsa.write(records[r].format('fastq'))
    print(f'Input Transcripts: {len(records)}\n'
          f'Clean Transcripts: {len(records) - contamination}\n'
          f'Contaminated transcripts: {contamination}')
