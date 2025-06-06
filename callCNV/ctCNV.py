#!/usr/bin/env python
# -*- coding: utf-8 -*-
# updata 190428
#    1): autocall  Edd parameter -e  default Exon_counts > 5.0 == > Exon_counts > 2.0
#    2): scatter   add parameter -d / --output--dir
# updata 190816
#    1): scatter annotate gene
# updata 211230
#    1): trans to python3
""" callCNV : an novel algorithm to detect SCNA in cfDNA samples. """

import argparse
import collections
import functools
import io
import logging
import math
import os
import random
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from itertools import takewhile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import pyfaidx

__version__ = "0.2.0"


def print_version(_args):
    """Display version"""
    print("callCNV version: %r" % __version__)


def report_bad_line(line_parser):
    """Raise value error when reading bad line in bed file."""

    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError:
            raise ValueError("Bad line: %r" % line)

    return wrapper


def read_bed(infile):
    """Read standard multi-column bed file."""

    @report_bad_line
    def _parse_line(line):
        fields = line.split('\t', 7)
        chrom, start, end = fields[:3]
        gene = (fields[3].rstrip() if len(fields) >= 4 else '-')
        multiple = (fields[4].rstrip() if len(fields) >= 5 else 'N')
        enID = (fields[5].rstrip() if len(fields) >= 6 else '-')
        strand = (fields[6].rstrip() if len(fields) >= 7 else '.')
        return chrom, int(start), int(end), gene, multiple, enID, strand

    with open(infile, 'rb') as handle:
        rows = list(map(_parse_line, (line for line in handle)))
        return pd.DataFrame.from_records(
            rows,
            columns=[
                "Chr",
                "start",
                "end",
                "gene",
                "multiple",
                "refseqID",
                "strand"])


def read_bed4(infile):
    """Read four columns from bed file."""
    table = read_bed(infile)
    table2 = table.loc[:, ['Chr', 'start', 'end', 'gene', 'multiple']]
    table3 = table2.sort_values(['Chr', 'start'])
    return table3


def ensure_path(fname):
    """Create dirs and move an existing file to avoid overwriting, if necessary.
    If a file already exists at the given path, it is renamed with an integer
    suffix to clear the way.
    """
    if '/' in os.path.normpath(fname):
        # Ensure the output directory exists
        dname = os.path.dirname(os.path.abspath(fname))
        if dname and not os.path.isdir(dname):
            try:
                os.makedirs(dname)
            except OSError as exc:
                raise OSError("Output path " + fname +
                              " contains a directory " + dname +
                              " that cannot be created: %s" % exc)
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = "%s.%d" % (fname, cnt)
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = "%s.%d" % (fname, cnt)
        os.rename(fname, bak_fname)
        logging.info("Moved existing file %s -> %s", fname, bak_fname)
    return True


def join_path(output_dir, fname):
    """ Join the output dir name and the file name."""
    if output_dir and isinstance(fname, str):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        outfname = os.path.join(output_dir, fname)
        return outfname
    elif isinstance(fname, str):
        return fname
    else:
        raise ValueError("Cannot access the dir: %r" % output_dir)


@contextmanager
def safe_write(outfile):
    """Check existance of path before creat the file."""
    if isinstance(outfile, str):
        dirname = os.path.dirname(outfile)
        if dirname and not os.path.isdir(dirname):
            os.mkdir(dirname)
        with open(outfile, 'w') as handle:
            yield handle
    else:
        yield outfile


def tab_read(infile):
    """Read tab seperated file into GenomicArray object."""
    try:
        # dframe = pd.read_table(infile, dtype={'Chr': 'str'})
        dframe = pd.read_csv(infile, dtype={'Chr': 'str'}, sep='\t')
    except pd.io.common.EmptyDataError:
        dframe = []
    result = sorted(GenomicArray(dframe))
    return result


def tab_write(GA, outfile=None):
    """Write GenomicArray object into tab seperated file."""
    if isinstance(GA, GenomicArray):
        dframe = GA.data
    else:
        dframe = GA
    with safe_write(outfile or sys.outfilestdout) as handle:
        dframe.to_csv(
            handle,
            header=True,
            index=False,
            sep='\t',
            float_format='%.8g')


def sorter_chrom(label):
    """Create a sorting key from chromosome label.
    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.
    E.g. chr1 < chr2 < chr10 < chrX < chrY < chrM
    """
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key


class GenomicArray(object):
    """Array class for genomic intervals."""
    _required_columns = ("Chr", "start", "end")
    _required_dtypes = (str, int, int)

    def __init__(self, data_table, meta_dict=None):
        if not isinstance(data_table, pd.DataFrame):
            data_table = pd.DataFrame(data_table)
        if not all(c in data_table.columns for c in self._required_columns):
            raise ValueError(
                "data table must have at least columns %r; got %r" % (
                    self._required_columns,
                    tuple(data_table.columns)))
        self.data = data_table
        self.meta = (
            dict(meta_dict) if meta_dict is not None and len(meta_dict) else {}
        )

    @property
    def Chr(self):
        return self.data['Chr']

    @property
    def start(self):
        return self.data['start']

    @property
    def end(self):
        return self.data['end']

    @property
    def _chr_x_label(self):
        if 'chr_x' in self.meta:
            return self.meta['chr_x']
        chr_x = ('chrX' if self.Chr.iat[0].startswith('chr') else 'X')
        self.meta['chr_x'] = chr_x
        return chr_x

    @property
    def _chr_y_label(self):
        if 'chr_y' in self.meta:
            return self.meta['chr_y']
        chr_y = ('chrY' if self.Chr.iat[0].startswith('chr') else 'Y')
        self.meta['chr_y'] = chr_y
        return chr_y

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.data.iloc[index]
        elif isinstance(index, str):
            return self.data[index]
        elif (isinstance(index, tuple) and len(index) == 2
              and index[1] in self.data.columns):
            return self.data.loc[index]
        elif isinstance(index, slice):
            return self.as_dataframe(self.data[index])
        else:
            try:
                if isinstance(index, type(None) or len(index) == 0):
                    empty = pd.DataFrame(columns=self.data.columns)
                    return self.as_dataframe(empty)
            except TypeError:
                raise TypeError("object of type %r" % type(index) +
                                "cannot be used as an index into a " +
                                self.__class__.__name__)
            return self.as_dataframe(self.data[index])

    def __setitem__(self, index, value):
        if isinstance(index, int):
            self.data.iloc[index] = value
        elif isinstance(index, str):
            self.data[index] = value
        elif isinstance(index, tuple) and len(index) == 2 and index[1] in self.data.columns:
            self.data.loc[index] = value
        else:
            assert isinstance(index, slice) or len(index) > 0
            self.data[index] = value

    def __iter__(self):
        return self.data.itertuples(index=False)

    def as_dataframe(self, dframe):
        return self.__class__(dframe.reset_index(drop=True), self.meta.copy())

    def by_chromosome(self):
        """Iterate over regions grouped by chromosome."""
        self.data.Chr = self.data.Chr.astype(str)
        for chrom, subtable in self.data.groupby("Chr", sort=False):
            yield chrom, self.as_dataframe(subtable)

    def autosomes(self):
        """Filter autosomes."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            is_auto = self.Chr.str.contains(r"\d")
        if not is_auto.any():  # 这句有啥用？
            return self
        return self[is_auto]

    def cutnames(self):
        """Remove multiple gene names"""
        # result = self.copy()
        if not self.data[self.data.gene.str.contains(';')].empty:
            self.data.gene = self.data.gene.astype(str).str.replace(';.*$', '')

    def coords(self):
        coordframe = self.data.loc[:, list(self._required_columns)]
        return coordframe.itertuples(index=False)

    def sort(self):
        """Sort this array's bins in-place, with smart chromosome ordering."""
        sort_key = self.data.Chr.apply(sorter_chrom)
        self.data = (self.data.assign(_sort_key_=sort_key)
                     .sort_values(by=['_sort_key_', 'start', 'end'],
                                  kind='mergesort',
                                  ignore_index=True)
                     .drop('_sort_key_', axis=1))

    def sort_columns(self):
        """Sort this array's columns in-place, per class definition."""
        extra_cols = []
        for col in self.data.columns:
            if col not in self._required_columns:
                extra_cols.append(col)
        sorted_colnames = list(self._required_columns) + sorted(extra_cols)
        assert len(sorted_colnames) == len(self.data.columns)
        self.data = self.data.reindex(columns=sorted_colnames)

    def data_print(self):
        print(self.data)


def Pool_map(Command, Parameter_List, Thread):
    with ThreadPoolExecutor(max_workers=Thread) as executor:
        future_tasks = executor.map(Command, Parameter_List)
        return future_tasks


def multi_bam_coverage(para):
    bed_f, min_mapq, bam = para
    result = bed_f.multi_bam_coverage(q=min_mapq, bams=bam)
    # data_table = pd.read_table(StringIO.StringIO(result), header=None)
    data_table = pd.read_csv(io.StringIO(result), header=None, sep='\t')
    index1 = data_table[0].astype(str).str.cat(
        data_table[1].astype(str), sep='\t').str.cat(
        data_table[2].astype(str), sep='\t').str.cat(
        data_table[3].astype(str), sep='\t').str.cat(
        data_table[4].astype(str), sep='\t')
    series = pd.Series(list(data_table.iloc[:, -1]), index=index1)
    return series


def multicov(bed="/work-a/user/guoh/Data/panel6.bed",
             bamfs=None,
             min_mapq=0, from_string=False, thread=1,
             # fa_fname="/lustre/rdi/user/guoh/Data/gatk_bundle/2.8/hg19/ucsc.hg19.fasta"):
             fa_fname="/mnt/share02/renhd/980_CNV/original_info/reference/hg19_reference_with_NC.fasta"):
    """Calculte coverage for multiple bam files or single bam file.
    using min_mapq to filter reads with low mapping qulity.
    using StringIO.StringIO(), if bed is not a file path but a data_frame.
    """

    if bamfs is None:
        bamfs = ['./mapped/P2198_blood_1a.bam']
    if not from_string:
        bed_f = pybedtools.BedTool(bed)
    else:
        bed_f = pybedtools.BedTool(bed, from_string=True)
    """
    result1 = bed_f.multi_bam_coverage(q=min_mapq, bams=bamfs)
    # bed_col_count = bed_f.field_count()
    data_table = pd.read_table(StringIO.StringIO(result1), header=None)
    """
    bam_count = len(bamfs)
    paralist = [[bed_f, min_mapq, bamfs[i]] for i in range(bam_count)]
    data_table_tmp = Pool_map(multi_bam_coverage, paralist, thread)
    data_table = pd.concat(data_table_tmp, axis=1)
    # tmplist = [tmp[i].split('\t') for i in xrange(len(list(data_table_tmp2.index)))]
    # tmparray = np.array(tmplist)
    data_table.insert(0, 'e', data_table.index.map(lambda x: x.split('\t')[4]))
    data_table.insert(0, 'd', data_table.index.map(lambda x: x.split('\t')[3]))
    data_table.insert(
        0, 'c', data_table.index.map(
            lambda x: int(
                x.split('\t')[2])))
    data_table.insert(
        0, 'b', data_table.index.map(
            lambda x: int(
                x.split('\t')[1])))
    data_table.insert(0, 'a', data_table.index.map(lambda x: x.split('\t')[0]))
    data_table.index = list(range(len(data_table)))

    # pybedtools.cleanup(remove_all=True)
    col_names = detect_multicov_columns(data_table, bamfs)
    data_table.columns = col_names
    if 'gene' in data_table:
        data_table['gene'] = data_table['gene'].fillna('-')
    else:
        data_table['gene'] = '-'
    if 'multiple' in data_table:
        data_table['multiple'] = data_table['multiple'].fillna('N')
    else:
        data_table['multiple'] = 'N'
    if 'depth' in data_table:
        sample_name = os.path.basename(bamfs[0]).split('.')[0]
        data_table.loc[data_table['depth'] == 0, 'depth'] = 0.001
        # ok_idx = (data_table['depth'] > 0)
        """
        ok_idx1 = (data_table['multiple'] == 'Y') & (data_table['depth'] > 0)
        ok_idx2 = (data_table['multiple'] == 'N') & (data_table['depth'] > 0)
        data_table.loc[ok_idx1, sample_name+'-log2_ND'] = (
            np.log2(data_table.loc[ok_idx1, 'depth']) -
            np.median(np.log2(data_table.loc[ok_idx1, 'depth'])))
        data_table.loc[ok_idx2, sample_name+'-log2_ND'] = (
            np.log2(data_table.loc[ok_idx2, 'depth']) -
            np.median(np.log2(data_table.loc[ok_idx2, 'depth'])))
        """
        for mul in set(data_table['multiple']):
            ok_idx_i = (
                               data_table['multiple'] == mul) & (
                               data_table['depth'] > 0)
            data_table.loc[ok_idx_i, sample_name + '-log2_ND'] = (
                    np.log2(data_table.loc[ok_idx_i, 'depth']) -
                    np.median(np.log2(data_table.loc[ok_idx_i, 'depth'])))

        tmp_GA = sorted(GenomicArray(data_table))
        gc, rmask = get_fasta_stats(tmp_GA, fa_fname)
        exon_len = np.array(tmp_GA.data['end'] - tmp_GA.data['start'])
        fix_b1 = center_by_window(
            tmp_GA, 0.1, gc, depth_column=sample_name + '-log2_ND')
        fix_b1 = center_by_window(
            fix_b1, 0.1, exon_len, depth_column=sample_name + '-log2_ND')
        data_table = fix_b1.data
        del data_table['multiple']
    else:
        data_table = data_table.fillna(0)
        for sample in col_names[5:]:
            data_table.loc[data_table[sample] == 0, sample] = 0.001
            """
            #ok_idx = (data_table[sample] > 0)
            ok_idx1 = (data_table['multiple'] == 'Y') & (data_table[sample] > 0)
            ok_idx2 = (data_table['multiple'] == 'N') & (data_table[sample] > 0)
            data_table.loc[ok_idx1, sample+'-log2_ND'] = (
                np.log2(data_table.loc[ok_idx1, sample]) -
                np.median(np.log2(data_table.loc[ok_idx1, sample])))
            data_table.loc[ok_idx2, sample+'-log2_ND'] = (
                np.log2(data_table.loc[ok_idx2, sample]) -
                np.median(np.log2(data_table.loc[ok_idx2, sample])))
            """
            for mul in set(data_table['multiple']):
                ok_idx_i = (
                                   data_table['multiple'] == mul) & (
                                   data_table[sample] > 0)
                data_table.loc[ok_idx_i, sample + '-log2_ND'] = (
                        np.log2(data_table.loc[ok_idx_i, sample]) -
                        np.median(np.log2(data_table.loc[ok_idx_i, sample])))
            tmp_GA = sorted(GenomicArray(data_table))
            gc, rmask = get_fasta_stats(tmp_GA, fa_fname)
            exon_len = np.array(tmp_GA.data['end'] - tmp_GA.data['start'])
            fix_b1 = center_by_window(
                tmp_GA, 0.1, gc, depth_column=sample + '-log2_ND')
            fix_b1 = center_by_window(
                fix_b1, 0.1, exon_len, depth_column=sample + '-log2_ND')
            fix_b1.sort()
            data_table = fix_b1.data
        data_table['mean_log2ND'] = (data_table
                                     .iloc[:, -int(bam_count):]
                                     .apply(np.mean, axis=1))
        data_table['std_log2ND'] = (data_table
                                    .iloc[:, -int(bam_count) - 1:-1]
                                    .apply(np.std, axis=1))
        data_table['median_log2ND'] = (data_table
                                       .iloc[:, -int(bam_count) - 2:-2]
                                       .apply(np.median, axis=1))
        del data_table['multiple']
    return data_table


def detect_multicov_columns(data_frame, bamfs):
    """Detect and rename the columns of the coverage data frame."""
    col_count = len(data_frame.columns)
    if col_count == 6:
        return ['Chr', 'start', 'end', 'gene', 'multiple', 'depth']
    elif col_count == 4:
        return ['Chr', 'start', 'end', 'depth']
    elif col_count == 5:
        return ['Chr', 'start', 'end', 'gene', 'depth']
    # elif col_count >= 6 and 'multiple' not in data_frame.columns:
    #    return ['Chr', 'start', 'end', 'gene'] + \
    #            [os.path.basename(a).split('.')[0] for a in bamfs]
    elif col_count >= 7:
        return ['Chr', 'start', 'end', 'gene', 'multiple'] + \
               [os.path.basename(a).split('.')[0] for a in bamfs]
    raise RuntimeError("No columns from multicov: \n%r" % data_frame.columns)


def ref_region_Z(ref_data, bamfs):
    """Caculate region Z score for reference."""
    bam_count = len(bamfs)
    h1 = list(ref_data.columns[0:4])
    h2 = list(ref_data.columns[-3:])
    ref_data2 = ref_data.iloc[:, -3:]
    # ref_data.to_csv('ref_data.tsv', sep='\t')
    # ref_data2.to_csv('ref_data2.tsv', sep='\t')
    sl = []
    for s in ref_data.iloc[:, -3 - bam_count:-3].columns:
        s_short = s.split("-log2")[0]
        print(s, s_short)
        sl.append(s_short)
        ref_data[s_short] = (ref_data[s] - ref_data2.iloc[:, -3]) / ref_data2.iloc[:, -2]
    return ref_data[h1 + h2 + sl]


def ref_GCS(ref_Z, bamfs):
    """Caculate gene specific score for reference."""
    bamcount = len(bamfs)
    GCS = ref_Z.loc[ref_Z['gene'].apply(len) > 1]
    ref_gcs = GCS.groupby('gene').apply(
        lambda subg: subg.iloc[:, -bamcount:].sum() / math.sqrt(subg['start'].count()))
    ref_gcs['Median_GCS'] = ref_gcs.apply(
        lambda x: np.median(x), axis=1)
    ref_gcs['cut_off'] = ref_gcs.apply(lambda x: abs(
        x).sort_values(ascending=True)[int(0.95 * len(x))], axis=1)
    return ref_gcs


def sample_region_Z(sample_data, ref_data):
    """Caculate region Z score for tumor sample."""
    name = sample_data.columns[-1].split('-')[0]
    sample_data.index = list(
        sample_data.Chr.astype(str).str
            .cat(sample_data.start.astype(str), sep=':')
            .astype(str).str.cat(sample_data.end.astype(str), sep='-'))
    ref_data.index = list(
        ref_data.Chr.astype(str).str.cat(ref_data.start.astype(str), sep=':')
            .astype(str).str.cat(ref_data.end.astype(str), sep='-'))
    result = pd.concat(
        [ref_data['gene'], ref_data.iloc[:, -3:], sample_data.iloc[:, -1]],
        axis=1,
        join_axes=[ref_data.index])
    result[name + '-RZ'] = (
                                   result.iloc[:, -1] - result.iloc[:, -4]) / result.iloc[:, -3]
    result = result.dropna()
    result[name + '-CopyRatio'] = (result.iloc[:, -2] - result.iloc[:, -3]) - np.median(
        (result.iloc[:, -2] - result.iloc[:, -3]))
    return result[result[name + '-RZ'].notnull()]


def sample_GCS(sample_Z):
    """Caculate gene specific score for tumor sample."""
    GCS = sample_Z.loc[sample_Z['gene'].apply(len) > 1]
    try:
        sample_gcs = GCS.groupby('gene').apply(
            lambda subg: subg.iloc[:, -2].sum() / math.sqrt(subg.iloc[:, -2].count()))
    except BaseException:
        print(sample_Z)
    sample_ratio = GCS.groupby('gene').apply(
        lambda subg: np.median(subg.iloc[:, -1]))
    sample_exons = GCS.groupby('gene').apply(
        lambda subg: subg.iloc[:, 1].count())
    sample_result = pd.concat(
        [sample_gcs, sample_ratio, sample_exons],
        axis=1,
        join_axes=[sample_gcs.index])
    return sample_result


def CNV_raw(sample_result, ref_gcs, bams2):
    """Call copy number and gene specific score."""
    col_name = (
            list(ref_gcs.columns) + [bams2 + '-GCS'] +
            [bams2 + '-CopyRatio'] + ['Exon_counts'])
    result = pd.concat(
        [ref_gcs, sample_result], axis=1, join_axes=[ref_gcs.index])
    result.columns = col_name
    result = result.assign(
        Filter=lambda x: (
                np.abs(x.iloc[:, -3]) > x.iloc[:, -4]
        ).map(pd.Series({True: "PASS", False: "."})))
    result = result.assign(
        CopyNumber=lambda x: 2 * np.power(2, x.iloc[:, -3]))
    return result


def center_by_window(cnarr, fraction, sort_key, depth_column='log2'):
    """Smooth out biases according to the trait specified by sort_key.
    E.g. correct GC-biased bins by windowed averaging across similar-GC
    bins; or for similar interval sizes.
    """
    depth_column = str(depth_column)
    df = cnarr.data.reset_index(drop=True)
    np.random.seed(0xA5EED)
    shuffle_order = np.random.permutation(df.index)
    df = df.reindex(shuffle_order)
    # Apply the same shuffling to the key array as to the target probe set
    assert isinstance(sort_key, (np.ndarray, pd.Series))
    sort_key = sort_key[shuffle_order]
    # Sort the data according to the specified parameter
    order = np.argsort(sort_key, kind='mergesort')
    df = df.iloc[order]
    biases = rolling_median(df[depth_column], fraction)
    df[depth_column] -= biases
    fixarr = sorted(cnarr.as_dataframe(df))
    return fixarr


def get_fasta_stats(cnarr, fa_fname):
    """Calculate GC and RepeatMasker content of each bin in the
    FASTA genome."""
    gc_rm_vals = [calculate_gc_lo(subseq)
                  for subseq in fasta_extract_regions(fa_fname, cnarr)]
    gc_vals, rm_vals = list(zip(*gc_rm_vals))
    return np.asfarray(gc_vals), np.asfarray(rm_vals)


def calculate_gc_lo(subseq):
    """Calculate the GC and lowercase (RepeatMasked) content of a string."""
    cnt_at_lo = subseq.count('a') + subseq.count('t')
    cnt_at_up = subseq.count('A') + subseq.count('T')
    cnt_gc_lo = subseq.count('g') + subseq.count('c')
    cnt_gc_up = subseq.count('G') + subseq.count('C')
    tot = float(cnt_gc_up + cnt_gc_lo + cnt_at_up + cnt_at_lo)
    if not tot:
        return 0.0, 0.0
    frac_gc = (cnt_gc_lo + cnt_gc_up) / tot
    frac_lo = (cnt_at_lo + cnt_gc_lo) / tot
    return frac_gc, frac_lo


def fasta_extract_regions(fa_fname, intervals):
    """Extract an iterable of regions from an indexed FASTA file.
    Input: FASTA file name; iterable of (seq_id, start, end) (1-based)
    Output: iterable of string sequences.
    """
    with pyfaidx.Fasta(fa_fname, as_raw=True) as fa_file:
        for chrom, subarr in intervals.by_chromosome():
            for _chrom, start, end in subarr.coords():
                # print(_chrom, start, end)
                # pyfaidx fetch sequence chrN: (start+1)-end from fasta file
                # bed file don't include end base
                # yield fa_file[_chrom][start.item()-1:end.item()-1]
                yield fa_file[_chrom][start - 1:end - 1]


def check_inputs(x, width):
    """Transform width into a half-window size.
    `width` is either a fraction of the length of `x` or an integer size of the
    whole window. The output half-window size is truncated to the length of `x`
    if needed.
    """
    x = np.asfarray(x)
    if 0 < width < 1:
        wing = int(math.ceil(len(x) * width * 0.5))
    elif 2 <= width == int(width):
        wing = width // 2
    else:
        raise ValueError("width fraction must be between 0 and 1 (got %s)"
                         % width)
    wing = max(wing, 5)
    wing = min(wing, len(x) - 1)
    assert wing > 0, "Wing must be greater than 0 (got %s)" % wing
    # Pad the edges of the original array with mirror copies
    signal = pd.Series(np.concatenate((x[wing - 1::-1], x, x[:-wing - 1:-1])))
    return x, wing, signal


def rolling_median(x, width):
    """Rolling median with mirrored edges."""
    x, wing, signal = check_inputs(x, width)
    rolled = signal.rolling(2 * wing + 1, 1, center=True).median()
    # if rolled.hasnans:
    #     rolled = rolled.interpolate()
    return np.asfarray(rolled[wing:-wing])


def chromosome_sizes(probes):
    """Caculate target region count for each autosome."""
    chrom_sizes = collections.OrderedDict()
    for chrs, rows in probes.by_chromosome():
        chrom_sizes[chrs] = rows['id'].max() - rows['id'].min() + 1
    return chrom_sizes


def plot_x_split(axis, chrom_sizes):
    """Plot the vline to seperate different Chromosomes."""
    x_dividers = []
    x_centers = []
    x_starts = collections.OrderedDict()
    offset = 0
    for lable, size in list(chrom_sizes.items()):
        x_starts[lable] = offset
        x_centers.append(offset + 0.5 * size)
        x_dividers.append(offset + size - 1)
        offset += size
    axis.set_xlim(0, offset + 1)
    for xposn in x_dividers[:-1]:
        axis.axvline(x=xposn, color='k', linestyle=':')
    # Use chromosome names as x-axis labels (instead of base positions)
    axis.set_xticks(x_centers)
    axis.set_xticklabels(list(chrom_sizes.keys()), rotation=60)
    axis.tick_params(labelsize='small')
    axis.tick_params(axis='x', length=0)
    axis.get_yaxis().tick_left()
    return x_starts


def run(refbams, tbam, bed, fasta=None):
    # bedfile = bed
    a = read_bed4(bed)
    cnarr = GenomicArray(a)
    cnarr2 = cnarr.autosomes()
    cnarr2.cutnames()
    data_string = cnarr2.data.to_string(header=False, index=False)
    bams = [x for x in refbams.strip().split(',')]
    # bams=['./mapped/P2920_blood_1a.bam','./mapped/P3040_blood_1a.bam','./mapped/P3043_blood_1a.bam']
    b = multicov(data_string, bamfs=bams, from_string=1)
    b.to_csv('call63-reference_depth.txt', sep='\t')
    c = ref_region_Z(b, bams)
    d = ref_GCS(c, bams)
    bams2 = [x for x in tbam.strip().split(',')]
    sample_name = os.path.basename(bams2[0]).split('.')[0]
    # bams2 = ['./mapped/P3039_blood_1a.bam']
    sample = multicov(data_string, bamfs=bams2, from_string=1)
    sample = sample[sample.depth > 0]
    sample.to_csv(sample_name + "-Coverage.txt", sep='\t')
    result = sample_region_Z(sample, b)
    result.to_csv(sample_name + "-RZ.txt", sep='\t')
    sample_gcs = sample_GCS(result)
    f = CNV_raw(sample_gcs, d, bams2)
    f.to_csv(sample_name + "-SCNA_results.txt", sep='\t')
    SCNA_passFilter = f.query("(Exon_counts>2.0) & (Filter=='PASS')")
    SCNA_passFilter.to_csv(
        sample_name + "-SCNA_results_passed_Filter.txt", sep='\t')


def _cmd_reference(args):
    """Make normal-pool reference for downstream CNV calling."""
    bed_raw = read_bed4(args.bed)
    cnarr = GenomicArray(bed_raw)
    cnarr2 = cnarr.autosomes()
    cnarr2.cutnames()
    data_string = cnarr2.data.to_string(header=False, index=False)
    fasta_name = args.fasta
    bams = []
    output_dir = args.output_dir
    thread_num = args.thread_num
    for path in args.inBams:
        if os.path.isdir(path):
            bams.extend(
                os.path.join(path, f) for f in os.listdir(
                    path) if f.endswith('.bam'))
        else:
            bams.append(path)
    ref_data = multicov(
        data_string, bamfs=bams, min_mapq=args.min_mapq,
        from_string=1, thread=thread_num, fa_fname=fasta_name)
    ref_fname = (
                        args.output and args.output + "_ref_COV.txt") or "cnv_ref_COV.txt"
    ref_fname = join_path(output_dir, ref_fname)
    ensure_path(ref_fname)
    tab_write(ref_data, outfile=ref_fname)
    ref_Z = ref_region_Z(ref_data, bams)
    ref_gcs = ref_GCS(ref_Z, bams)
    ref_gcs_fname = (
                            args.output and args.output + "_ref_GCS.txt") or "cnv_ref_GCS.txt"
    ref_gcs_fname = join_path(output_dir, ref_gcs_fname)
    ensure_path(ref_gcs_fname)
    ref_gcs.to_csv(ref_gcs_fname, sep='\t')


def _cmd_coverage(args):
    """Calculate coverage in the given regions from BAM file."""
    bed_raw = read_bed4(args.bed)
    cnarr = GenomicArray(bed_raw)
    cnarr2 = cnarr.autosomes()
    cnarr2.cutnames()
    data_string = cnarr2.data.to_string(header=False, index=False)
    fasta_name = args.fasta
    bams2 = args.inBam
    sample_name = os.path.basename(bams2).split('.')[0]
    sample_data = multicov(
        data_string, bamfs=[bams2], min_mapq=args.min_mapq,
        from_string=1, fa_fname=fasta_name)
    sample_fname = (
                           args.output and args.output + "_COV.txt") or sample_name + "_COV.txt"
    ensure_path(sample_fname)
    tab_write(sample_data, outfile=sample_fname)


def _cmd_sex(args):
    """Detect sample sex."""
    bed_raw = read_bed4(args.bed)
    cnarr2 = GenomicArray(bed_raw)
    cnarr2.cutnames()
    data_string = cnarr2.data.to_string(header=False, index=False)
    fasta_name = args.fasta
    bams2 = args.inBam
    sample_name = os.path.basename(bams2).split('.')[0]
    sample_data = multicov(
        data_string, bamfs=[bams2], min_mapq=args.min_mapq,
        from_string=1, fa_fname=fasta_name)
    sex_dic = {}
    sex_dic.setdefault('chrX_RegionCount', str(
        sample_data.loc[sample_data['Chr'] == 'chrX'].iloc[:, -1].count()))
    sex_dic.setdefault('chrY_RegionCount', str(
        sample_data.loc[sample_data['Chr'] == 'chrY'].iloc[:, -1].count()))
    sex_dic.setdefault('chrX_log2Ratio', str(
        sample_data.loc[sample_data['Chr'] == 'chrX'].iloc[:, -1].median()))
    sex_dic.setdefault('chrY_log2Ratio', str(
        sample_data.loc[sample_data['Chr'] == 'chrY'].iloc[:, -1].median()))
    x_cp = round(
        2 * 2 ** sample_data.loc[sample_data['Chr'] == 'chrX'].iloc[:, -1].median()
    )
    y_cp = round(
        2 * 2 ** sample_data.loc[sample_data['Chr'] == 'chrY'].iloc[:, -1].median()
    )
    sex_dic.setdefault('chrX_CopyNumber', str(x_cp))
    sex_dic.setdefault('chrY_CopyNumber', str(y_cp))
    sample_sex = (
                         x_cp > 1.8 and y_cp < 0.2 and "Female") or (
                         x_cp > 0.6 and y_cp > 0.5 and "Male") or "Und"
    sex_dic.setdefault('sample_sex', sample_sex)
    df = pd.DataFrame(sex_dic, index=[0])
    sample_fname = (
                           args.output and args.output + "_sex.txt") or sample_name + "_sex.txt"
    ensure_path(sample_fname)
    tab_write(df, outfile=sample_fname)


def _cmd_sigtest(args):
    """Estimate the significance of gene specific CNV."""
    ref_GA = tab_read(args.refCov)
    tumor_GA = tab_read(args.tumorCov)
    ref_data = ref_GA.data
    tumor_data = tumor_GA.data
    sample_Z = sample_region_Z(tumor_data, ref_data)
    sample_name = os.path.basename(args.tumorCov).split('_COV')[0]
    fname1 = (
                     args.output and args.output + "_RZ.txt") or sample_name + "_RZ.txt"
    ensure_path(fname1)
    tab_write(sample_Z, outfile=fname1)
    sample_gcs = sample_GCS(sample_Z)
    fname2 = (
                     args.output and args.output + "_GCS.txt") or sample_name + "_GCS.txt"
    ensure_path(fname2)
    sample_gcs.columns = ['GCS', 'log2_ratio', 'exons']
    sample_gcs.to_csv(fname2, sep='\t')


def _cmd_call(args):
    """Call somatic gene specific copy number variants."""
    bams2 = args.output
    refGCS_data = pd.read_table(args.refGCS, sep='\t', index_col=0)
    tumorGCS_data = pd.read_table(args.tumorGCS, sep='\t', index_col=0)
    CNV_result = CNV_raw(tumorGCS_data, refGCS_data, bams2)
    fname1 = args.output + "_SCNA_results.txt"
    ensure_path(fname1)
    CNV_result.to_csv(fname1, header=True, index=True, sep='\t')
    fname2 = args.output + "-SCNA_results_passed_Filter.txt"
    ensure_path(fname2)
    SCNA_passFilter = CNV_result.query("(Exon_counts>4.0) & (Filter=='PASS')")
    SCNA_passFilter.to_csv(fname2, header=True, index=True, sep='\t')


def _cmd_autocall(args):
    """A batch commands to run callCNV pipeline based on existing
    reference files (COV and GCS) and tumor sample bam.
    """
    output_dir = args.output_dir
    # calculate tumor sample coverage
    bed_raw = read_bed4(args.bed)
    bams2 = args.inBam
    sample_name = os.path.basename(bams2).split('.')[0]
    cnarr = GenomicArray(bed_raw)
    cnarr.cutnames()
    test = cnarr.data
    gene_dic = {x[3]: x[0] for x in test.itertuples(index=False)}
    cnarr2 = cnarr.autosomes()
    cnarr2.cutnames()
    data_string = cnarr2.data.to_string(header=False, index=False)
    fasta_name = args.fasta
    sample_data = multicov(
        data_string, bamfs=[bams2], min_mapq=args.min_mapq,
        from_string=1, fa_fname=fasta_name)
    sample_fname = (
                           args.output and args.output + "_COV.txt") or sample_name + "_COV.txt"
    sample_fname = join_path(output_dir, sample_fname)
    ensure_path(sample_fname)
    tab_write(sample_data, outfile=sample_fname)

    # calculate CNV significance
    ref_GA = tab_read(args.refCov)
    ref_data = ref_GA.data
    tumor_data = sample_data
    sample_Z = sample_region_Z(tumor_data, ref_data)
    fname1 = (
                     args.output and args.output + "_RZ.txt") or sample_name + "_RZ.txt"
    fname1 = join_path(output_dir, fname1)
    ensure_path(fname1)
    tab_write(sample_Z, outfile=fname1)
    sample_gcs = sample_GCS(sample_Z)
    fname2 = (
                     args.output and args.output + "_GCS.txt") or sample_name + "_GCS.txt"
    fname2 = join_path(output_dir, fname2)
    ensure_path(fname2)
    sample_gcs.columns = ['GCS', 'log2_ratio', 'exons']
    sample_gcs.to_csv(fname2, sep='\t')

    # call somatic gene specific copy number variants
    # refGCS_data = pd.read_table(args.refGCS, sep='\t', index_col=0)
    refGCS_data = pd.read_csv(args.refGCS, sep='\t', index_col=0)
    tumorGCS_data = sample_gcs
    CNV_result = CNV_raw(tumorGCS_data, refGCS_data, args.output)
    fname3 = args.output + "_SCNA_results.txt"
    fname3 = join_path(output_dir, fname3)
    ensure_path(fname3)
    CNV_result.to_csv(fname3, header=True, index=True, sep='\t')
    fname4 = args.output + "-SCNA_results_passed_Filter.txt"
    fname4 = join_path(output_dir, fname4)
    ensure_path(fname4)
    Exon_counts = args.min_exonCount
    SCNA_passFilter = CNV_result[(CNV_result["Exon_counts"] > Exon_counts) & (
            CNV_result["Filter"] == "PASS")]
    # SCNA_passFilter = CNV_result.query("(Exon_counts>2.0) & (Filter=='PASS')")
    SCNA_passFilter.to_csv(fname4, header=True, index=True, sep='\t')

    # generate reformated tsv fiel for GeneticTest report
    fname5 = args.output + "-SCNA_results_gtest.txt"
    fname5 = join_path(output_dir, fname5)
    ensure_path(fname5)
    f = SCNA_passFilter.copy()
    ff = f.reset_index()
    ff['log2_ratio'] = ff.iloc[:, -4]
    ff['Chr'] = ff['gene'].map(gene_dic)
    reformated_results = ff[['Chr', 'gene', 'log2_ratio', 'CopyNumber']]
    reformated_results.to_csv(fname5, header=True, index=False, sep='\t')


def _cmd_scatter(args):
    """Plot the copy number variations."""
    ref_COV = args.refCov
    sample_COV = args.tumorCov
    SCNA_pass = args.SCNAtxt
    SampleType = args.SampleType

    ref_GA = tab_read(ref_COV)
    sample_GA = tab_read(sample_COV)
    sample_name = os.path.basename(sample_COV).split("_COV")[0]
    sample_data = sample_GA.data
    ref_data = ref_GA.data
    sample_data.index = list(
        sample_data.Chr.astype(str).str.cat(
            sample_data.start.astype(str),
            sep=':').astype(str).str.cat(
            sample_data.end.astype(str),
            sep='-'))
    ref_data.index = list(
        ref_data.Chr.astype(str).str.cat(
            ref_data.start.astype(str),
            sep=':').astype(str).str.cat(
            ref_data.end.astype(str),
            sep='-'))
    SN = int((len(ref_data.columns) - 7) / 2)
    result = pd.concat([ref_data.iloc[:, 0:4], ref_data.iloc[:, SN + 4:],
                        sample_data.iloc[:, -1]], axis=1, join_axes=[ref_data.index])
    rr = result.reset_index()
    del rr['index']
    rr['id'] = rr.index
    dic = {}
    for x in rr.columns[4:4 + SN]:
        dic.setdefault(x, x.split("-log2")[0])
    dic.setdefault(rr.columns[-2], 'tumor')
    cc = rr.rename(columns=dic)
    for x in cc.columns[4:4 + SN]:
        cc[x] = cc[x] - cc['median_log2ND']
    cc['tumor'] = cc['tumor'] - cc['median_log2ND']
    # ref_list = list(cc.columns[4:SN+4]) + [str(cc.columns[-1])]
    # sample_list =  list(cc.columns[-2:])+['gene']
    sample_list2 = ['Chr', 'start', 'end', 'gene', 'tumor', 'id']
    pdata = cc[sample_list2]
    pdata2 = pdata.copy()
    pdata2['median'] = pdata2['gene'].map(
        dict(
            pdata2.groupby("gene").apply(
                lambda subg: subg['tumor'].median())))

    pass_dic = {}
    with open(SCNA_pass, 'rb') as f:
        for line in f:
            pass_dic.setdefault(line.strip().split()[
                                    0], line.strip().split()[-2])
    pdata2['color'] = pdata2['gene'].map(pass_dic).map(
        {".": "blue", "PASS": "red"}).fillna("gray")
    pdata2['tumor'] = 2 * 2 ** pdata2['tumor'].astype(float)
    # pdata2 = pdata2.fillna("gray")
    cnarr1 = GenomicArray(pdata2)
    cnarr1 = cnarr1[cnarr1['Chr'].str.len() < 6]
    # matplotlib parameters
    plt.switch_backend('agg')
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'DejaVuSerif'
    plt.rcParams['font.monospace'] = 'DroidSansMono'
    plt.style.use('classic')
    _fig, axis = plt.subplots(figsize=(15, 6))
    chrom_sizes = chromosome_sizes(cnarr1)
    axis.set_title(sample_name)
    x_starts = plot_x_split(axis, chrom_sizes)
    axis.axhline(y=2, color='k')
    axis.set_ylabel("Copy Number")
    ymax = 2 * 2 ** cnarr1.data[cnarr1.data['color']
                                == 'red']['median'].max() + 1
    axis.set_ylim(-0.1, 2 * 2 **
                  cnarr1.data[cnarr1.data['color'] == 'red']['median'].max() + 1)
    # Plot target region log2 ratio
    axis.scatter(
        cnarr1['id'],
        cnarr1['tumor'],
        color='#808080',
        edgecolor='none',
        alpha=0.8,
        marker='o')
    cnarr1 = cnarr1[cnarr1['gene'].str.len() > 1]
    seg_lines = []
    for gene, subdata in cnarr1.data.groupby("gene", sort=False):
        line_color = (
                             subdata['median'].median() > 0 and list(
                         subdata['color'])[0] == 'red' and '#ff4500') or (
                             subdata['median'].median() < 0 and list(
                         subdata['color'])[0] == 'red' and '#7cfc00') or "#555555"
        seg_lines.append(
            (subdata['median'].median(),
             subdata['id'].min(),
             subdata['id'].max(),
             line_color,
             gene))
    # Plot gene log2 ratio
    cutoff = 2.5 if SampleType == 1 else 3
    geneList = [
        'ALK',
        'AR',
        'CDK4',
        'EGFR',
        'ERBB2',
        'FGFR1',
        'FGFR2',
        'FGFR3',
        'FLT1',
        'KIT',
        'KRAS',
        'MET',
        'PIGF',
        'PIK3CA',
        'PTEN',
        'RICTOR',
        'VEGFA',
        'CDK6',
        'MDM2',
        'MDM4',
        'STK11',
        'B2M',
        'POLE',
        'CCND1',
        'FGF3',
        'FGF4',
        'FGF19',
        'MYC',
        'TP53',
        'RB1',
        'BRAF',
        'SMAD4',
        'APC',
        'ATM',
        'TSC1',
        'NF1',
        'PBRM1',
        'FBXW7']
    x_marker = [0]
    y_marker = [0]
    sepList = list(set(range(1, 20)) - {5})
    for line in seg_lines:
        y1, x1, x2, segment_color, gene = line
        y1 = 2 * 2 ** y1
        axis.plot((x1, x2), (y1, y1),
                  color=segment_color, linewidth=4, solid_capstyle='round')
        if ((gene in geneList) and (y1 > cutoff) and segment_color != "#555555") or (
                (gene in geneList) and (y1 < 1.5) and segment_color != "#555555"):
            if y1 > 2:
                ytest = y1 + \
                        random.uniform(min(0.5, 0.6 * (ymax - y1)), max(0.5, 0.6 * (ymax - y1)))
                y_clust = 0
                count = 0
                for i in range(-1, max(-6, -1 * len(y_marker)), -1):
                    if abs(
                            ytest - y_marker[i]) < 0.5 and y_marker[i] > y_clust:
                        y_clust = y_marker[i]
                        count += 1
                if count > 0:
                    ytest = min(
                        ymax,
                        y_clust +
                        np.sign(
                            ytest -
                            y_clust) *
                        sepList[count] /
                        5.00)
                xtest = max(x_marker[-1], x1) + 100
            else:
                ytest = y1 - \
                        random.uniform(min(0.5, 0.6 * y1), max(0.5, 0.6 * y1))
                y_clust = 100000
                count = 0
                for i in range(-1, max(-6, -1 * len(y_marker)), -1):
                    if abs(
                            ytest - y_marker[i]) < 0.5 and y_marker[i] < y_clust:
                        y_clust = y_marker[i]
                        count += 1
                if count > 0:
                    ytest = max(-1, y_clust + np.sign(ytest -
                                                      y_clust) * sepList[count] / 5.00)
                xtest = max(x_marker[-1], x1) + 100
            x_marker.append(xtest)
            y_marker.append(ytest)
            plt.annotate(gene,
                         xy=(x1 + 0.5 * (x2 - x1),
                             y1),
                         xytext=(xtest,
                                 ytest),
                         arrowprops=dict(arrowstyle='-|>',
                                         connectionstyle='arc3,rad=0.2',
                                         color='darkslategrey'))

    plt.tight_layout()
    output_dir = args.output_dir
    fname = sample_name + '_CopyNumber_scatter.png'
    sample_fname = join_path(output_dir, fname)
    plt.savefig(sample_fname, dpi=300)


def main():
    """ Main command line parser"""
    parser = argparse.ArgumentParser(
        description="""callCNV: A novel algorithm for identification of Somatic
                    Copy Number Variations in tumor plasma cfDNA sample.""",
        epilog='Contact Hao Guo <guo.hao@genecast.com.cn> for help.')
    subparsers = parser.add_subparsers(
        help='Sub-commands, using -h for more info.')

    P_ref = subparsers.add_parser('reference', help=_cmd_reference.__doc__)
    P_ref.add_argument(
        '-i', '--inBams', nargs='*',
        help="Input normal sample bam files.")
    P_ref.add_argument(
        '-o', '--output',
        help="Output basename followed by ext like XXX_ref_COV.txt")
    P_ref.add_argument('-b', '--bed', help="Input target region BED file.")
    P_ref.add_argument(
        '-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
    P_ref.add_argument(
        '-q', '--min-mapq', type=int, default=0,
        help="""Minimum mapping quality score (MAPQ) allowed.
                                [Default: %(default)s]""")
    P_ref.add_argument(
        '-p', '--thread_num', type=int, default=1,
        help="number of threads[1]")
    P_ref.add_argument(
        '-d', '--output-dir', default='.', help="Output directory.")
    P_ref.set_defaults(func=_cmd_reference)

    P_cov = subparsers.add_parser('coverage', help=_cmd_coverage.__doc__)
    P_cov.add_argument('-i', '--inBam', help='Input tumor sample bam file.')
    P_cov.add_argument(
        '-o', '--output', help='Output basename of sample coverage file.')
    P_cov.add_argument(
        '-b', '--bed', help="Input target region BED file.")
    P_cov.add_argument(
        '-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
    P_cov.add_argument(
        '-q', '--min-mapq', type=int, default=0,
        help="""Minimum mapping quality score (MAPQ).
                                [Default: %(default)s]""")
    P_cov.set_defaults(func=_cmd_coverage)

    P_sig = subparsers.add_parser('sigtest', help=_cmd_sigtest.__doc__)
    P_sig.add_argument(
        '-r', '--refCov', help='Input baseline reference coverage file.')
    P_sig.add_argument(
        '-t', '--tumorCov', help='Input tumor sample coverage file.')
    P_sig.add_argument(
        '-o', '--output',
        help='Output basename of region score and gene copy-number score.')
    P_sig.set_defaults(func=_cmd_sigtest)

    P_call = subparsers.add_parser('call', help=_cmd_call.__doc__)
    P_call.add_argument(
        '-r', '--refGCS', help='Input baseline reference GCS file.')
    P_call.add_argument(
        '-t', '--tumorGCS', help='Input tumor sample GCS file.')
    P_call.add_argument(
        '-o', '--output', help='Output basename of SCNV result file')
    # P_call.add_argument('', '', help='')
    P_call.set_defaults(func=_cmd_call)

    P_autocall = subparsers.add_parser('autocall', help=_cmd_autocall.__doc__)
    P_autocall.add_argument(
        '-i', '--inBam', help='Input tumor sample bam file.')
    P_autocall.add_argument(
        '-c', '--refCov', help='Input baseline reference coverage file.')
    P_autocall.add_argument(
        '-g', '--refGCS', help='Input baseline reference GCS file.')
    P_autocall.add_argument(
        '-o', '--output', help='Output basename of SCNV result file')
    P_autocall.add_argument(
        '-d', '--output-dir', default='.', help="Output directory.")
    P_autocall.add_argument(
        '-b', '--bed', help="Input target region BED file.")
    P_autocall.add_argument(
        '-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
    P_autocall.add_argument(
        '-e', '--min-exonCount', type=int, default=2,
        help="""Minimum Exon Count.
                                [Default: %(default)s]""")
    P_autocall.add_argument(
        '-q', '--min-mapq', type=int, default=0,
        help="""Minimum mapping quality score (MAPQ).
                                [Default: %(default)s]""")
    P_autocall.set_defaults(func=_cmd_autocall)

    P_scatter = subparsers.add_parser('scatter', help=_cmd_scatter.__doc__)
    P_scatter.add_argument(
        '-r', '--refCov', help='Input baseline reference coverage file.')
    P_scatter.add_argument(
        '-t', '--tumorCov', help='Input tumor sample coverage file.')
    P_scatter.add_argument(
        '-s', '--SCNAtxt', help='Input SCNA result txt file with PASS Filter.')
    P_scatter.add_argument(
        '-st',
        '--SampleType',
        type=int,
        default=2,
        help='Input Sample Type(CF:1, non-CF:2).')
    P_scatter.add_argument(
        '-d', '--output-dir', help='Output directory.')
    P_scatter.set_defaults(func=_cmd_scatter)

    P_sex = subparsers.add_parser('guessex', help=_cmd_sex.__doc__)
    P_sex.add_argument(
        '-i', '--inBam', help='Input tumor sample bam file.')
    P_sex.add_argument(
        '-o', '--output', help='Output basename of sample coverage file.')
    P_sex.add_argument(
        '-b', '--bed', help="Input target region BED file.")
    P_sex.add_argument(
        '-f', '--fasta',
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)")
    P_sex.add_argument(
        '-q', '--min-mapq', type=int, default=0,
        help="""Minimum mapping quality score (MAPQ).
                                [Default: %(default)s]""")
    P_sex.set_defaults(func=_cmd_sex)

    P_version = subparsers.add_parser(
        'version', help=print_version.__doc__)
    P_version.set_defaults(func=print_version)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
