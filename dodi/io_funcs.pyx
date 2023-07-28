#!python
#cython: language_level=2, boundscheck=False
#cython: profile=False
#distutils: language=c++

# Language level 2 is needed for char map

import numpy as np
cimport numpy as np
from collections import defaultdict
import os
import click
from sys import stderr, stdin
import logging

ctypedef np.float_t DTYPE_t

from libc.stdlib cimport malloc
from libcpp.vector cimport vector as cpp_vector
import re


def echo(*arg):
    click.echo(arg, err=True)


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))


def to_dict(self):
    return {v: self.__getattribute__(v) for v in dir(self) if "__" not in v and v not in ('to_dict', 'from_dict', 'data', 'rows')} # != "to_dict" and v != "from_dict"}


def from_dict(self, d):
    allowed = set(dir(self))
    for k, v in d.items():
        if k in allowed:
            self.__setattr__(k, v)
    return self


cdef class Params:
    def __init__(self, args):
        self.match_score = args["match_score"]
        self.mu = args["insert_median"]
        self.sigma = args["insert_stdev"]
        self.min_aln = args["min_aln"]
        self.max_homology = args["max_overlap"]
        self.inter_cost = args["inter_cost"]
        self.U = args["u"]
        # cdef float zero_cost_boundary = args["zero_cost_boundary"]
        # cdef float max_gap_cost = args["max_gap_cost"]
        self.ins_cost = args['ins_cost']
        self.ol_cost = args['ol_cost']

        self.paired_end = int(args["paired"])
        self.bias = args["bias"]
        self.secondary = args['secondary']

        self.default_max_d = self.mu + (4 * self.sigma)  # Separation distance threshold to call a pair discordant

        self.find_insert_size = False if (not args['paired'] or args['template_size'] == 'auto') else True
        self.modify_mapq = args['modify_mapq']
        self.add_tags = args['tags']
    def __repr__(self):
        return ', '.join([f'{k}={v}' for k, v in {'match_score': self.match_score, 'mu': self.mu, 'sigma': self.sigma, 'min_aln': self.min_aln,
                'max_homology': self.max_homology, 'inter_cost': self.inter_cost, 'u': self.U, 'ins_cost': self.ins_cost,
                'ol_cost': self.ol_cost, 'paired_end': self.paired_end, 'bias': self.bias, 'secondary': self.secondary,
                'default_max_d': self.default_max_d, 'find_insert_size': self.find_insert_size, 'modify_mapq': self.modify_mapq,
                'add_tags': self.add_tags
        }.items()])


cdef class Template:
    def __init__(self, rows, last_seen_chrom): # , paired_end, match_score, bias,
                         # secondary, min_aln, max_hom, inter_cost, U, zero_cost_boundary, max_gap_cost):

        # self.match_score = match_score
        #
        # self.min_aln = min_aln
        # self.max_homology = max_hom
        # self.inter_cost = inter_cost
        # self.U = U
        # self.zero_cost_bound = zero_cost_boundary
        # self.max_gap_cost = max_gap_cost
        #
        # self.paired_end = paired_end
        # self.bias = bias
        # self.secondary = secondary

        self.inputdata = rows

        self.read1_length = 0
        self.read2_length = 0
        self.score_mat = {}
        self.passed = 0
        self.name = rows[0][0][0]
        self.last_seen_chrom = last_seen_chrom
        self.read1_seq = ""  # Some sam records may have seq == '*' , need a record of full seq for adding back in
        self.read2_seq = ""
        self.read1_q = ""
        self.read2_q = "",
        self.outstr = ""
        self.read1_reverse = False  # Set to true if aligner has reverse complemented the sequence
        self.read2_reverse = False
        self.read1_unmapped = False
        self.read2_unmapped = False

        self.first_read2_index = 0

    def __repr__(self):
        return str(to_dict(self))


cpdef Template make_template(rows, last_seen_chrom):
                             # bint paired_end, float match_score, float bias,
                         # bint secondary, float min_aln, float max_hom, float inter_cost, float U,
                         #     float zero_cost_boundary, float max_gap_cost):

    return Template(rows, last_seen_chrom)  #, paired_end, match_score, bias,
                    # secondary, min_aln=min_aln, max_hom=max_hom, inter_cost=inter_cost,
                    # U=U, zero_cost_boundary=zero_cost_boundary, max_gap_cost=max_gap_cost)


def sam_to_str(template_name, sam):
    for i, item in enumerate(sam):
        if len(item[8]) != len(item[9]):
            logging.critical(f'SEQ and QUAL not same length, index={i}, qlen={len(item[8])}, quallen={len(item[9])} ' + template_name + " " + str(sam))
            quit()
    return "".join(template_name + "\t" + "\t".join(i) + "\n" for i in sam)


def get_include_reads(include_regions, bam):

    if not include_regions:
        for r in bam:
            yield r

    regions = [i.strip().split("\t")[:3] for i in open(include_regions, "r") if i[0] != "#"]
    for c, s, e in regions:
        logging.info("Reading {}:{}-{}".format(c, s, e))
        for r in bam.fetch(c, int(s), int(e)):
            yield r


from libc.stdio cimport FILE, stdin

cdef extern from "stdio.h":
    FILE *fopen(const char *, const char *)
    int fclose(FILE *)
    ssize_t getline(char **, size_t *, FILE *)


def file_iter(sam):

    filename_byte_string = sam.encode("ascii")
    cdef char * fname = filename_byte_string

    cdef FILE * cfile
    if sam in '-stdin':
        cfile = stdin
    else:
        cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(f"No such file or directory: {sam}")

    cdef char * buffer = NULL
    cdef size_t bufsize = 0
    cdef ssize_t read

    while True:
        read = getline(&buffer, &bufsize, cfile)
        if read == -1:
            break
        yield str(buffer.decode("ascii"))
    fclose(cfile);


def sam_itr(args):

    # itr = args["sam"]
    itr = file_iter(args['sam'])
    tree = overlap_regions(args["include"])

    # First get header
    header_string = ""
    last_seen_chrom = ""
    first_line = ""
    for t in itr:
        if t[0] == "@":
            header_string += t
            continue
        first_line = t.split("\t", 9)
        last_seen_chrom = first_line[2]
        yield header_string
        break

    try:
        pos = int(first_line[3])
    except IndexError:
        raise IOError("No sam lines detected")

    ol = intersecter(tree, first_line[2], pos, pos + 250)

    yield first_line, last_seen_chrom, ol

    for t in itr:
        line = t.split("\t", 9)
        if line[3] != last_seen_chrom:
            last_seen_chrom = line[2]

        pos = int(line[3])
        ol = intersecter(tree, line[2], pos, pos + 250)

        yield line, last_seen_chrom, ol


def iterate_mappings(args, version):

    params = {'paired',
              'template_size',
              'keep_mapq',
              'secondary',
              'tags',
              'min_aln',
              'max_overlap',
              'zero_cost_boundary',
              'max_gap_cost',
              'inter_cost',
              'pairing_cost',
              'match_score',
              'include',
              'bias',
              'u'}
    cp_args = {k: v for k, v in args.items() if k in params}

    arg_str = ", ".join(["{}={}".format(i, j) for i, j in args.items() if i in params])
    inputstream = sam_itr(args)

    total = 0
    name = ""
    rows = []
    header_string = next(inputstream)
    header_string += "@PG\tID:DODI\tPN:dodi\tVN:{}\tCL:{}\n".format(version, arg_str)

    yield header_string

    last_seen_chrom = ""
    for m, last_seen_chrom, ol in inputstream:  # Alignment
        nm = m[0]
        if name != nm:
            if len(rows) > 0:
                total += 1
                yield rows, last_seen_chrom
            rows = []
            name = nm
        rows.append((m, ol))  # String, ol states if alignment overlaps ROI

    # Deal with last record
    if len(rows) > 0:
        total += 1
        yield rows, last_seen_chrom


cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0',
                       '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0',  'a',  'a' ]

np.random.seed(0)


cpdef str reverse_complement(str seq, int seq_len):
    """https://bioinformatics.stackexchange.com/questions/3583/\
    what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho/3595#3595"""

    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = '\0'

    cdef bytes py_bytes = seq.encode('ascii')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('ascii')


cpdef get_start_end(cigar):
    c = re.split(r'(\d+)', cigar)
    cdef int end = 0
    cdef int start = 0
    cdef int template_length = 0
    cdef int i

    for i in range(1, len(c), 2):
        ci = int(c[i])
        if i == 1 and (c[i+1] == "S" or c[i+1] == "H"):
            start += ci
            end += ci
        elif c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += ci
        if c[i+1] != "D":
            template_length += ci
    return start, end, template_length


cpdef int get_align_end_offset(cigar):
    c = re.split(r'(\d+)', cigar)
    cdef int end = 0
    cdef int i
    for i in range(1, len(c), 2):
        if c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return end


cdef tuple check_for_good_pairing(template, add_tags, max_d):

    # Check to see if this pair of alignments needs pairing
    r1 = template.inputdata[0]
    r2 = template.inputdata[1]

    cdef int aflag = int(r1[0])
    cdef int bflag = int(r2[0])

    if aflag & 4:
        template.read1_unmapped = 1
    if bflag & 4:
        template.read2_unmapped = 1

    # Check if rnext and pnext set properly
    cdef int p1, p2
    cdef int proper_pair = 1 if aflag & 2 else 0

    if not proper_pair and not aflag & 12:  # Read unmapped/mate unmapped. Double check in case flag not set properly
        p1 = int(r1[2])  # Position
        p2 = int(r2[2])

        if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

            if abs(p1 - p2) < max_d:
                # Check for FR or RF orientation
                if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                    proper_pair = 1

    if add_tags:
        if proper_pair == 1:
            dn = '0.0'
        else:
            dn = '250.0'

        if aflag & 4 and bflag & 4:  # Read unmapped, mate unmapped
            r1 += f"ZM:f:0\tZA:i:0\tZP:f:0\tZN:f:{dn}\tZS:f:0\tZO:i:0"
            r2 += f"ZM:f:0\tZA:i:0\tZP:f:0\tZN:f:{dn}\tZS:f:0\tZO:i:0"

        elif aflag & 4:
            r1 += f"ZM:f:0\tZA:i:0\tZP:f:0\tZN:f:{dn}\tZS:f:0\tZO:i:0"
            r2 += f"ZM:f:0\tZA:i:100\tZP:f:125\tZN:f:{dn}\tZS:f:125\tZO:i:0"

        elif bflag & 4:
            r1 += f"ZM:f:0\tZA:i:100\tZP:f:125\tZN:f:{dn}\tZS:f:125\tZO:i:0"
            r2 += f"ZM:f:0\tZA:i:0\tZP:f:0\tZN:f:{dn}\tZS:f:0\tZO:i:0"

        else:
            r1 += f"ZM:f:1\tZA:i:100\tZP:f:250\tZN:f:{dn}\tZS:f:250\tZO:i:0"
            r1 += f"ZM:f:1\tZA:i:100\tZP:f:250\tZN:f:{dn}\tZS:f:250\tZO:i:0"

    return r1, r2


cpdef int sam_to_array(template, params) except -1:
    # Expect read1 and read2 alignments to be concatenated, not mixed together
    data, overlaps = list(zip(*template.inputdata))

    # split the rest of the columns
    template.inputdata = [i[1:-1] + i[-1].strip().split("\t") for i in data]

    # If only one alignment for read1 and read2, no need to try pairing, just send sam to output
    if params.paired_end and len(data) == 2:
        pair_str = check_for_good_pairing(template, params.add_tags, params.default_max_d)

        if pair_str:
            template.passed = 1
            template.outstr = pair_str
            return 1

    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_mis-matches, original_aln_score]
    cdef np.ndarray[np.float_t, ndim=2] arr = np.zeros((len(data), 10))

    chrom_ids = {}

    if template.inputdata[0][1] == "*":
        template.inputdata[0][1] = template.last_seen_chrom

    cdef int cc = 0
    cdef int idx, pos, flag, seq_len, query_start, query_end, start_temp
    cdef str chromname, cigar, k, t, v
    cdef float bias = params.bias

    cdef int read1_strand_set, read2_strand_set, current_l
    cdef int first_read2_index = len(template.inputdata) + 1
    read1_set = 0  # Occasionally multiple primaries, take the longest
    read2_set = 0

    for idx in range(len(template.inputdata)):

        l = template.inputdata[idx]

        flag = int(l[0])
        pos = int(l[2])
        if l[1] != "*":
            chromname = l[1]
            template.last_seen_chrom = chromname
        else:
            l[1] = template.last_seen_chrom
            chromname = l[1]
        if chromname not in chrom_ids:
            chrom_ids[chromname] = cc
            cc += 1

        arr[idx, 0] = chrom_ids[chromname]
        arr[idx, 1] = pos
        arr[idx, 5] = idx
        arr[idx, 6] = -1 if flag & 16 else 1

        if idx == 0 and flag & 4:
            template.read1_unmapped = 1

        # Sometimes both first and second are not flagged. Assume first
        if not flag & 64 and not flag & 128:
            arr[idx, 7] = 1
        else:
            if flag & 64:
                arr[idx, 7] = 1
            else:
                arr[idx, 7] = 2
                if idx < first_read2_index:
                    first_read2_index = idx
                    if flag & 4:
                        template.read2_unmapped = 1

        tags = [i.split(":") for i in l[11:]]
        seq_len = len(l[8])
        for k, t, v in tags:
            if k == "NM":
                arr[idx, 8] = float(v)
            elif k == "AS":
                arr[idx, 9] = float(v)  # Keep copy of original alignment score
                if overlaps[idx]:
                    arr[idx, 4] = float(v) * bias
                else:
                    arr[idx, 4] = float(v)

        current_l = len(l[8])

        if params.paired_end:

            if flag & 64 and read1_set < current_l and len(l[8]) > 1:  # First in pair
                template.read1_seq = l[8]
                template.read1_q = l[9]
                template.read1_length = seq_len
                read1_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template.read1_reverse = 1 if flag & 16 else 0

            elif flag & 128 and read2_set < current_l and len(l[8]) > 1:  # Second in pair
                template.read2_seq = l[8]
                template.read2_q = l[9]
                template.read2_length = seq_len
                read2_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template.read2_reverse = 1 if flag & 16 else 0

        else:
            if not flag & 256 and len(l[8]) > 1 and read1_set < current_l:
                template.read1_seq = l[8]
                template.read1_q = l[9]
                template.read1_length = len(l[8])
                read1_set = current_l

        cigar = l[4]
        if not cigar:
            query_start = 0  # Unmapped read? no cigar
            query_end = 0
        else:

            query_start, query_end, template_length = get_start_end(cigar)

            # If current alignment it not primary, and on different strand from primary, count from other end
            if params.paired_end:
                if flag & 64 and template.read1_reverse != bool(flag & 16):  # First in pair, read1_rev != read_rev
                    start_temp = template.read1_length - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

                elif flag & 128 and template.read2_reverse != bool(flag & 16):  # Second in pair
                    start_temp = template.read2_length - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

            else:
                # if template.read1_length == 0:
                if template_length == 0:
                    template_length = len(l[8])  # if length from cigar failed, use seq length
                    if template_length <= 1:
                        echo('Could not infer template length')
                        return -1
                        # raise ValueError('Could not infer template length')
                template.read1_length = template_length

                if flag & 16:  # Single end Reverse strand, count from end
                    start_temp = template.read1_length - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

        arr[idx, 2] = query_start
        arr[idx, 3] = query_end

    # if template.name == '01726306-45be-4337-bf31-c1b71083d2e9.21q1F_False':
    #     echo('Hiiii')
    #     echo(len(template.read1_seq))
    # todo is this needed?
    # if first_read2_index == len(arr) + 1:
    #     template.paired_end = 0

    cdef int j
    if params.paired_end:  # Increment the contig position of read 2
        for j in range(first_read2_index, len(arr)):
            if arr[j, 7] == 2:  # Second in pair
                arr[j, 2] += template.read1_length
                arr[j, 3] += template.read1_length
            if arr[j, 3] > template.read1_length + template.read2_length:
                echo('Infered template length greather than found read1 + read2 length. Is the file name sorted?', arr[j, 3], template.read1_length, template.read2_length)
                # raise ValueError
                return -1
    template.first_read2_index = first_read2_index
    template.data_ori = arr
    template.chrom_ids = chrom_ids

    return 0


cpdef void choose_supplementary(template):
    # Final alignments have been chosen, but need to decide which is supplementary
    cdef int j = 0

    cdef double read1_max = -1
    cdef double read2_max = -1
    cdef int i = 0
    cdef int row_idx

    primary_1 = -1
    primary_2 = -1

    cdef double[:, :] data_table = template.data_ori
    cdef long[:] rows = template.rows
    for row_idx in rows:
        if data_table[row_idx, 7] == 1 and data_table[row_idx, 9] > read1_max:  # Use original alignment score, not biased (idx==4)
            read1_max = data_table[row_idx, 9]
            primary_1 = row_idx

        elif data_table[row_idx, 7] == 2 and data_table[row_idx, 9] > read2_max:
            read2_max = data_table[row_idx, 9]
            primary_2 = row_idx

    template.primary1 = primary_1
    template.primary2 = primary_2


cpdef void score_alignments(template):
    return
    # Scans all alignments for each query, slow for long reads but ok for short read data
    # the goal is to try and keep extra supplementary alignments that are not on the main alignment path
    # Used for DN, similar to XS
    # all_xs = []
    # cdef int i, actual_row, idx
    # cdef float xs = -1
    # cdef float size = 0
    # cdef float qstart = 0
    # cdef float qend = 0
    # cdef float isize = 0
    # cdef float istart = 0
    # cdef float iend = 0
    # cdef float iscore = 0
    # cdef float ol = 0
    # cdef float ori_aln_score = 0
    # cdef int twoi = template["first_read2_index"]  # Use to skip read1/read2 alignments
    # cdef int row_indx
    #
    # idx = 0
    #
    # inputdata = template['inputdata']
    # tabledata = template['data_ori']
    #
    # extra_row_idxs = []
    # # for item in template_rows:
    # for row_indx in template['rows']:
    #
    #     item = tabledata[5]
    #     qstart = item[2]
    #     qend = item[3]
    #     readn = item[7]
    #     ori_aln_score = item[9]
    #     size = qend - qstart
    #
    #     if template["paired_end"]:
    #         if readn == 2:
    #             rr = range(twoi, len(tabledata))
    #         else:
    #             rr = range(0, twoi)
    #     else:
    #         rr = range(len(tabledata))
    #
    #     for i in rr:
    #         if i == row_indx:
    #             continue
    #         istart, iend, iscore = tabledata[i, [2, 3, 4]]  # Note use biased align score, otherwise gets confusing
    #         isize = (iend - istart) + 1e-6
    #         # Check for overlap
    #         ol = max(0, min(qend, iend) - max(qstart, istart))
    #         if ol and (ol / size > .85) and (ol / isize > .85):  # check for 85% reciprocal overlap with current
    #             if iscore > xs:
    #                 xs = iscore
    #
    #     if xs == -1:
    #         xs = ori_aln_score
    #     template["score_mat"][template["locs"][idx]][1] = xs
    #     idx += 1


cdef list get_bed_regions(str bed):
    b = [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r")
         if i[0] != "#" and len(i) > 0 and "\t" in i]
    if len(b) == 0:
        raise ValueError("Bed regions not formatted correctly")
    return b


def merge_intervals(intervals, srt=True, pad=0, add_indexes=False):
    """
    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]
    """
    if srt:
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end (index)
    else:
        sorted_by_lower_bound = intervals

    if pad:
        if not add_indexes:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad] for c, i, j in sorted_by_lower_bound]
        else:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad, k] for c, i, j, k in sorted_by_lower_bound]

    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                if not add_indexes:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
                else:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]), lower[3] + [higher[3]])
            else:
                if not add_indexes:
                    merged.append(higher)
                else:
                    merged.append(list(higher)[:3] + [[higher[3]]])
    return merged





cdef class Py_BasicIntervalTree:
    def __cinit__(self):
        self.thisptr = new BasicIntervalTree()
    def __dealloc__(self):
        del self.thisptr
    cpdef void add(self, int start, int end, int index):
        self.thisptr.add(start, end, index)
    cpdef bint searchInterval(self, int pos, int pos2):
        return self.thisptr.searchInterval(pos, pos2)
    cpdef overlappingInterval(self, int pos, int pos2):
        cdef Interval* res = self.thisptr.overlappingInterval(pos, pos2)
        if res[0] is None:  # [0] dereferences pointer
            return None
        else:
            return res[0].low, res[0].high
    cpdef void index(self):
        self.thisptr.index()
    cpdef allOverlappingIntervals(self, int start, int end):
        cdef cpp_vector[int] res
        self.thisptr.allOverlappingIntervals(start, end, res)
        return list(res)
    cpdef int countOverlappingIntervals(self, int pos, int pos2):
        return self.thisptr.countOverlappingIntervals(pos, pos2)


def iitree(a, add_value=False):
    # sorted input list a will not lead to a balanced binary-search-tree if added in sequential order
    # This function reorders the list to ensure a balanced BST when added using Py_BasicIntervalTree
    tree = Py_BasicIntervalTree()
    index = 0
    for h in a:
        if add_value:
            tree.add(h[0], h[1], h[2])
        else:
            tree.add(h[0], h[1], index)
        index += 1
    tree.index()
    return tree


cpdef dict overlap_regions(str bed):
    if not bed:
        return {}
    regions = merge_intervals(get_bed_regions(bed))
    chrom_intervals = defaultdict(list)
    for c, s, e in regions:
        chrom_intervals[c].append((s, e))

    chrom_intervals = {k: iitree(v) for k, v in chrom_intervals.items()}
    return chrom_intervals


cpdef int intersecter(tree, chrom, int start, int end):
    cdef bint found = 0
    if tree and chrom in tree:
        found = tree[chrom].searchInterval(start, end)
    return found
