#cython: language_level=3
# cython: profile=False

from __future__ import absolute_import
import click
import numpy as np
import logging

from dodi import io_funcs
from dodi.samclips cimport fixsam
from dodi.pairing cimport pairing_process
from dodi.io_funcs cimport Template, Params, Reader, Writer, make_template, sam_to_array, choose_supplementary
from libc.stdio cimport FILE, stdout, stdin


cdef extern from "stdio.h":  # manually wrapped, for some reason this needs to be done
    FILE *fopen(const char *, const char *)
    int fclose(FILE *)
    int fputs(const char *, FILE *)


def echo(*arg):
    click.echo(arg, err=True)


cdef to_output(Template template, Params params, Writer writer):
    if template.outstr:
        sam = list(template.outstr)
    else:
        sam = fixsam(template, params)

    if not sam:
        logging.error(f"Read dropped QNAME={template.name}, covert to sam failed")

    for i, item in enumerate(sam):
        seq_len = len(item[8])
        qual_len = len(item[9])
        if seq_len != qual_len and item[8] and item[9] != "*":
            logging.critical(
                f'SEQ and QUAL not same length, index={i}, qlen={seq_len}, quallen={qual_len} ' +
                f'{template.name} {sam}'
            )
    template_name_with_tab = template.name + "\t"
    formatted_lines = [template_name_with_tab + "\t".join(item) for item in sam]
    for item in formatted_lines:
        writer.write(item.encode('ascii'))


def median(L):
    # Stolen from https://github.com/arq5x/lumpy-sv/blob/master/scripts/pairend_distro.py
    if len(L) == 0:
        return 500
    elif len(L) == 1:
        return L[0]
    if len(L) % 2 == 1:
        return L[int(len(L)/2)]  # cast to int since divisions always return floats in python3
    mid = int(len(L) / 2) - 1
    return (L[mid] + L[mid+1]) / 2.0


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    xs.sort()
    med = median(xs)
    umad = median([x - med for x in xs if x > med])
    return med, umad


def mean_std(L):
    s = sum(L)
    mean = np.median(L)
    sq_sum = 0.0
    for v in L:
        sq_sum += (v - mean)**2.0
    var = sq_sum / float(len(L))
    return mean, var**0.5


def get_insert_params(L, mads=8):
    c = len(L)
    med, umad = unscaled_upper_mad(L)
    upper_cutoff = int(min(med + mads * umad, 10000))
    L = [v for v in L if v < upper_cutoff]
    new_len = len(L)
    removed = c - new_len
    mean, stdev = mean_std(L)
    mean = int(mean)
    stdev = int(stdev)
    if logging.DEBUG >= logging.root.level:
        logging.debug(f"dodi insert size {mean} +/- {stdev}")
    return mean, stdev


def insert_size(batch):

    required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted
    cdef int flag
    tlens = []
    for b in batch:
        for aln in b.inputdata:
            flag = int(aln[1])
            if not flag & 2:
                continue
            tlen = int(aln[8])
            rname = aln[2]
            rnext = aln[6]
            if rnext == '=':
                rnext = rname
            if rname == rnext and flag & flag_mask == required and tlen >= 0:
                tlens.append(tlen)
                break
        if len(tlens) > 1000:
            break

    if len(tlens) > 0:
        insert_m, insert_stdev = get_insert_params(tlens)
    else:
        logging.info("dodi insert size, not enough pairs, using 350 +/- 200")
        insert_m, insert_stdev = 350, 200
    return insert_m, insert_stdev


cdef void process_batch(batch, Params params, Writer writer, tree):
    cdef char * to_write
    cdef bytes byte_string
    if params.find_insert_size:
        insert, insert_std = insert_size(batch)
        max_d = insert + (4 * insert_std)
        params.mu = insert
        params.sigma = insert_std
        params.default_max_d = max_d

    cdef int done
    cdef Template temp
    c = 0
    for temp in batch:
        done = sam_to_array(temp, params, tree)
        if done < 0:
            logging.error(f"Read dropped QNAME={temp.name}, sam_to_array failed")
        elif done:
            to_output(temp, params, writer)
        else:
            pairing_process(temp, params)
            choose_supplementary(temp)
            if not temp.passed:
                logging.error(f"Read dropped QNAME={temp.name}, process failed")
            to_output(temp, params, writer)
            c += 1


def process_reads(args, batch_size):
    # Set Writer
    cdef FILE * outsam
    cdef char *fname
    if args["output"] in {"-", "stdout"} or args["output"] is None:
        outsam = stdout
    else:
        filename_byte_string = args["output"].encode("UTF-8")
        fname = filename_byte_string
        outsam = fopen(fname, "rb")
    if outsam == NULL:
        raise IOError("Unable to open output file")

    cdef Writer writer = Writer()
    writer.set_ptr(outsam)

    # Set Reader
    cdef FILE * insam
    if args['sam'] in '-stdin':
        insam = stdin
    else:
        filename_byte_string = args['sam'].encode("UTF-8")
        fname = filename_byte_string
        insam = fopen(fname, "rb")
    if insam == NULL:
        raise FileNotFoundError(f"No such file or directory: {args['sam']}")

    cdef Reader reader = Reader()
    reader.set_ptr(insam)

    # For checking if alignments fall in target regions
    tree = io_funcs.overlap_regions(args["include"])

    cdef Params params = Params(args)
    cdef int count = 0

    first_line = reader.header_to_file(writer, io_funcs.header_info_line(args))
    if not first_line:
        return 0

    batch = []
    rows = [first_line.split("\t", 9)]

    while True:
        line = reader.read_line()
        if not line:
            break

        l = line.split("\t", 9)
        if len(rows) == 0 or rows[-1][0] == l[0]:
            rows.append(l)

        elif rows:
            count += 1
            batch.append(make_template(rows))
            if len(batch) >= batch_size:
                process_batch(batch, params, writer, tree)
                batch = []
            rows = [l]

    if rows:
        count += 1
        batch.append(make_template(rows))
    if batch:
        process_batch(batch, params, writer, tree)

    return count
