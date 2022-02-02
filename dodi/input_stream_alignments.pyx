#cython: language_level=3
# cython: profile=True

from __future__ import absolute_import
import multiprocessing
import sys
import pkg_resources
import click
import time
import datetime
import pickle
import numpy as np
from sys import stderr
import logging

from . import pairing, io_funcs, samclips
from sys import stderr


def echo(*arg):
    click.echo(arg, err=True)


cdef void process_template(read_template):
    paired = io_funcs.sam_to_array(read_template)

    # if read_template['name'] == 'V300096939L3C004R0130394504':
    #     print(read_template, file=stderr)
    #     print("PAIRED", file=stderr)
    #     print(paired, file=stderr)

    if paired:
        return

    res = pairing.process(read_template)

    # if read_template['name'] == 'V300096939L3C004R0130394504':
    #     #     print(read_template)
    #     #     print(read_template['data'].astype(int))
    #     #     print(res)
    if res:
        read_template["passed"] = True
        io_funcs.add_scores(read_template, *res)
        io_funcs.choose_supplementary(read_template)
        if read_template['secondary']:
            io_funcs.score_alignments(read_template)


cpdef list to_output(dict template):

    if "outstr" in template:
        return list(template["outstr"])

    return samclips.fixsam(template)


def load_mq_model(pth):
    if pth:
        click.echo("Loading MapQ recalibrator {}".format(pth), err=True)
        return pickle.load(open(pth, "rb"))
    else:
        click.echo("No MapQ recalibration", err=True)
        return None


def phred_from_model(p):
    if p == 1:
        return 30
    if p < 0.5:
        return 0
    # return 40
    P = 1 - p
    v = int(round(-10 * np.log10(P)))
    return v if v <= 30 else 30


def predict_mapq(xtest, model):
    return list(map(phred_from_model, model.predict_proba(xtest)[:, 1]))


cdef write_records(sam, mq_model, outsam):
    # Model features are "AS", "DA", "DN", "DP", "NP", "PS", "XS", "kind_key", "mapq", "DS
    # if not mq_model:
    for name, record in sam:
        outsam.write(io_funcs.sam_to_str(name, record))

    # else:
    #     map_qs = []
    #     for name, alns in sam:
    #
    #         for a in alns:
    #             mapq = a[3]
    #
    #             t = {i[:2]: i[5:] for i in a[10:]}
    #             kind_key = 0 if int(a[0]) & 2048 else 1
    #             f = [t["AS"], t["DA"], t["DN"], t["DP"], t["NP"], t["PS"], t["XS"] if "XS" in t else 0, kind_key, mapq, t["DS"]]
    #             map_qs.append(f)
    #
    #     mqs = iter(predict_mapq(np.array(map_qs).astype(float), mq_model))
    #
    #     for name, alns in sam:
    #         for i in range(len(alns)):
    #             flag = int(alns[i][0])
    #             # cigar = alns[i][4]
    #             mq = next(mqs)
    #
    #             if flag & 2048:
    #                 if flag & 4:
    #                     alns[i][3] = "0"
    #                 else:
    #                     alns[i][3] = str(mq)
    #
    #         outsam.write(io_funcs.sam_to_str(name, alns))


cpdef list job(data_tuple):

    temp = io_funcs.make_template(*data_tuple)

    process_template(temp)
    sam_temp = []
    if temp['passed']:
        sam = to_output(temp)
        if sam:
            sam_temp.append((temp["name"], sam))
    return sam_temp


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

    cdef int required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted

    tlens = []
    for b in batch:
        for aln, _ in b['inputdata']:
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

    insert_std = args["template_size"].split(",")
    args["insert_median"] = float(insert_std[0])
    args["insert_stdev"] = float(insert_std[1])

    if not args["include"]:
        args["bias"] = 1.0
    else:
        logging.info("Elevating alignments in --include with --bias {}".format(args["bias"]))

    if args["output"] in {"-", "stdout"} or args["output"] is None:
        outsam = sys.stdout
    else:
        outsam = open(args["output"], "w")

    if (args["fq1"] or args["fq2"]) and args["procs"] > 1:
        raise ValueError("Cant use procs > 1 with fq input")

    map_q_recal_model = None  # load_mq_model(args["mq"])

    count = 0

    version = pkg_resources.require("dodi")[0].version
    itr = io_funcs.iterate_mappings(args, version)

    isize = (args["insert_median"], args["insert_stdev"])
    match_score = args["match_score"]
    pairing_params = (args["max_insertion"], args["min_aln"], args["max_overlap"], args["ins_cost"],
                      args["ol_cost"], args["inter_cost"], args["u"])
    paired_end = int(args["paired"] == "True")
    bias = args["bias"]
    replace_hard = int(args["replace_hardclips"] == "True")
    secondary = args['secondary'] == 'True'

    default_max_d = args["insert_median"] + 4*args["insert_stdev"]  # Separation distance threshold to call a pair discordant

    # if args["procs"] != 1:
    #
    #     header_string = next(itr)
    #     outsam.write(header_string)
    #
    #     temp = []
    #     for rows, last_seen_chrom, fq in itr:
    #         temp.append((rows, max_d, last_seen_chrom, fq, pairing_params, paired_end, isize,
    #                                      match_score, bias, replace_hard))
    #
    #         if len(temp) > 10000:
    #             with multiprocessing.Pool(args["procs"]) as p:
    #                 res = p.map(job, temp)
    #
    #             res = [item for sublist in res for item in sublist if item]  # Remove [] and flatten list
    #             write_records(res, map_q_recal_model, outsam)
    #             temp = []
    #
    #     if len(temp) > 0:
    #
    #         with multiprocessing.Pool(args["procs"]) as p:
    #             res = p.map(job, temp)
    #         res = [item for sublist in res for item in sublist if item]
    #         write_records(res, map_q_recal_model, outsam)

    # Use single process

    find_insert_size = True
    if not args['paired']:
        find_insert_size = False

    if True:

        header_string = next(itr)
        outsam.write(header_string)

        sam_temp = []

        batch = []
        for rows, last_seen_chrom, fq in itr:

            count += 1

            # rows, max_d, last_seen_chrom, fq
            temp = io_funcs.make_template(rows, default_max_d, last_seen_chrom, fq, pairing_params, paired_end, isize,
                                         match_score, bias, replace_hard, secondary)

            if len(batch) < 100e3:
                batch.append(temp)

            else:
                max_d = default_max_d
                if find_insert_size:
                    insert, insert_std = insert_size(batch)
                    max_d = insert + (4 * insert_std)

                for temp in batch:

                    if max_d != default_max_d:
                        temp['max_d'] = max_d

                    process_template(temp)

                    if temp['passed']:
                        sam = to_output(temp)
                        if sam:
                            sam_temp.append((temp["name"], sam))

                        if len(sam_temp) > 50000:
                            write_records(sam_temp, map_q_recal_model, outsam)
                            sam_temp = []
                batch = []


        if len(batch) > 0:

            max_d = default_max_d
            if find_insert_size:
                insert, insert_std = insert_size(batch)
                max_d = insert + (4 * insert_std)

            for temp in batch:

                if max_d != default_max_d:
                    temp['max_d'] = max_d

                process_template(temp)

                if temp['passed']:
                    sam = to_output(temp)
                    if sam:
                        sam_temp.append((temp["name"], sam))
            batch = []

        if len(sam_temp) > 0:
            write_records(sam_temp, map_q_recal_model, outsam)

    if args["output"] != "-" or args["output"] is not None:
        outsam.close()

    logging.info("dodi completed in {} h:m:s".format(str(datetime.timedelta(seconds=int(time.time() - t0)))),
               )
