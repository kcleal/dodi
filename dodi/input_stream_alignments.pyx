#cython: language_level=3
# cython: profile=False

from __future__ import absolute_import
import multiprocessing
import pkg_resources
import click
import time
import datetime
import numpy as np
import logging


from . import pairing, io_funcs, samclips
from dodi.io_funcs import make_template, sam_to_str
from dodi.io_funcs cimport Params
from libc.stdio cimport FILE, stdout


cdef extern from "stdio.h":
    FILE *fopen(const char *, const char *)
    int fclose(FILE *)
    int fputs(const char *, FILE *)


def echo(*arg):
    click.echo(arg, err=True)


cdef void process_template(read_template, params):
    done = io_funcs.sam_to_array(read_template, params)
    if done:
        return

    pairing.process(read_template, params)
    if not read_template.passed:
        return

    io_funcs.choose_supplementary(read_template)
    # if read_template.secondary:
    #     io_funcs.score_alignments(read_template, params)


cpdef list to_output(template, params):
    if template.outstr:
        return list(template.outstr)
    return samclips.fixsam(template, params)


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
    logging.info(f"dodi insert size {mean} +/- {stdev}")
    return mean, stdev


def insert_size(batch):

    required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted

    tlens = []
    for b in batch:
        for aln, _ in b.inputdata:
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

    if len(tlens) > 0:
        insert_m, insert_stdev = get_insert_params(tlens)
    else:
        logging.info("dodi insert size, not enough pairs, using 350 +/- 200")
        insert_m, insert_stdev = 350, 200
    return insert_m, insert_stdev


def process_reads(args):
    t0 = time.time()
    if args['template_size'] != 'auto':
        insert_std = args["template_size"].split(",")
        args["insert_median"] = float(insert_std[0])
        args["insert_stdev"] = float(insert_std[1])
    else:
        args["insert_median"] = 300.
        args["insert_stdev"] = 200.

    if not args["include"]:
        args["bias"] = 1.0
    else:
        logging.info("Elevating alignments in --include with --bias {}".format(args["bias"]))

    cdef FILE * outsam
    cdef char *fname
    if args["output"] in {"-", "stdout"} or args["output"] is None:
        outsam = stdout
    else:
        filename_byte_string = args["output"].encode("UTF-8")
        fname = filename_byte_string
        outsam = fopen(fname, "rb")

    #
    count = 0

    version = pkg_resources.require("dodi")[0].version
    itr = io_funcs.iterate_mappings(args, version)

    params = Params(args)

    n_jobs = args['procs']

    cdef char* to_write
    cdef bytes byte_string

    if params.paired_end:
        batch_size = 10_000
    else:
        batch_size = 10

    if n_jobs == 1:

        header_string = next(itr)

        byte_string = header_string.encode('ascii')
        to_write = byte_string
        fputs(to_write, outsam)

        jobs = []
        batch = []
        for rows, last_seen_chrom in itr:
            count += 1
            temp = make_template(rows, last_seen_chrom)
            if len(batch) < batch_size:
                batch.append(temp)
            else:
                # process one batch
                max_d = params.default_max_d
                insert = params.mu
                insert_std = params.sigma
                if params.find_insert_size:
                    insert, insert_std = insert_size(batch)
                    max_d = insert + (4 * insert_std)
                    params.mu = insert
                    params.sigma = insert_std
                    params.default_max_d = max_d
                for temp in batch:
                    process_template(temp, params)
                    if temp.passed:
                        sam = to_output(temp, params)
                        if sam:
                            byte_string = sam_to_str(temp.name, sam).encode('ascii')
                            to_write = byte_string
                            fputs(to_write, outsam)
                batch = []

        if len(batch) > 0:

            max_d = params.default_max_d
            insert = params.mu
            insert_std = params.sigma
            if params.find_insert_size:
                insert, insert_std = insert_size(batch)
                max_d = insert + (4 * insert_std)
                params.mu = insert
                params.sigma = insert_std
                params.default_max_d = max_d

            for temp in batch:
                process_template(temp, params)
                if temp.passed:
                    sam = to_output(temp, params)
                    if sam:
                        byte_string = sam_to_str(temp.name, sam).encode('ascii')
                        to_write = byte_string
                        fputs(to_write, outsam)

            batch = []

    fclose(outsam)

    logging.info("dodi completed in {} h:m:s".format(str(datetime.timedelta(seconds=int(time.time() - t0)))),
               )