#!python
# cython: language_level=3, boundscheck=False
# cython: profile=False
# distutils: language=c++


import numpy as np
cimport numpy as np
from collections import defaultdict
import os
import click
import logging

ctypedef np.float_t DTYPE_t

from libc.stdlib cimport malloc
from libcpp.vector cimport vector as cpp_vector

from cython.operator import dereference


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
        self.ins_cost = args['ins_cost']
        self.ol_cost = args['ol_cost']
        self.paired_end = int(args["paired"])
        self.bias = args["bias"]
        self.secondary = args['secondary']
        self.default_max_d = self.mu + (4 * self.sigma)  # Separation distance threshold to call a pair discordant
        self.find_insert_size = False if (not args['paired'] or args['template_size'] != 'auto') else True
        self.modify_mapq = args['modify_mapq']

    def __repr__(self):
        return ', '.join([f'{k}={v}' for k, v in {'match_score': self.match_score, 'mu': self.mu, 'sigma': self.sigma, 'min_aln': self.min_aln,
                'max_homology': self.max_homology, 'inter_cost': self.inter_cost, 'u': self.U, 'ins_cost': self.ins_cost,
                'ol_cost': self.ol_cost, 'paired_end': self.paired_end, 'bias': self.bias, 'secondary': self.secondary,
                'default_max_d': self.default_max_d, 'find_insert_size': self.find_insert_size, 'modify_mapq': self.modify_mapq,
        }.items()])


cdef class Template:
    # def __init__(self, rows, last_seen_chrom):
    def __cinit__(self, rows):
        self.inputdata = rows
        self.read1_length = 0
        self.read2_length = 0
        # self.score_mat = {}
        self.passed = 0
        # self.name = rows[0][0][0]
        self.name = rows[0][0]
        # self.last_seen_chrom = last_seen_chrom
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
        self.primary1 = -1
        self.primary2 = -1
    def __init__(self, rows):
        pass

    def __repr__(self):
        return str(to_dict(self))


cdef Template make_template(rows):
    return Template(rows)


def get_include_reads(include_regions, bam):
    if not include_regions:
        for r in bam:
            yield r
    regions = [i.strip().split("\t")[:3] for i in open(include_regions, "r") if i[0] != "#"]
    for c, s, e in regions:
        logging.info("Reading {}:{}-{}".format(c, s, e))
        for r in bam.fetch(c, int(s), int(e)):
            yield r


cdef class Writer:
    def __init__(self):
        pass
    cdef void set_ptr(self, FILE *out):
        self.outsam = out
    cdef void write(self, bytes byte_string):
        cdef char * to_write = byte_string
        fputs(to_write, self.outsam)
    def __dealloc__(self):
        if self.outsam:
            fclose(self.outsam)


cdef class Reader:
    def __cinit__(self):
        self.bufsize = 0
    def __init__(self):
        pass
    cdef void set_ptr(self, FILE *in_sam):
        self.insam = in_sam

    cdef header_to_file(self, Writer writer, extra_header_line):
        cdef int hl = 0
        cdef bytes byte_string
        cdef char * to_write
        while True:
            read = getline(&self.buffer, &self.bufsize, self.insam)
            if read == -1:
                return
            if self.buffer[0] == b'@':
                fputs(self.buffer, writer.outsam)
                hl += 1
            else:
                if hl > 0:
                    byte_string = extra_header_line.encode('ascii')
                    to_write = byte_string
                    fputs(to_write, writer.outsam)
                return self.buffer.decode("ascii")
    cdef read_line(self):
        read = getline(&self.buffer, &self.bufsize, self.insam)
        if read == -1:
            return
        return self.buffer.decode("ascii")

    def __dealloc__(self):
        if self.insam:
            fclose(self.insam)


def header_info_line(args):
    params = ['paired',
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
              'u']
    arg_str = ", ".join(["{}={}".format(i, j) for i, j in args.items() if i in params])
    return "@PG\tID:DODI\tPN:dodi\tVN:{}\tCL:{}\n".format(args['version'], arg_str)


cdef char *basemap = [ b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b'T', b'\0',  b'G', b'\0', b'\0', b'\0',  b'C', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',  b'N', b'\0',
                       b'\0', b'\0', b'\0', b'\0',  b'A',  b'A', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b't', b'\0',  b'g', b'\0', b'\0', b'\0',  b'c', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                       b'\0', b'\0', b'\0', b'\0',  b'a',  b'a' ]


cpdef str reverse_complement(str seq, int seq_len):
    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = b'\0'

    cdef bytes py_bytes = seq.encode('ascii')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('ascii')


cdef tuple get_start_end(cigar):
    cdef int start = 0
    cdef int end = 0
    cdef int template_length = 0
    cdef int num = 0  # To accumulate the number as we parse through the string
    cdef char op
    cdef bytes t = cigar.encode('utf8')
    cdef char *c_cigar = t
    cdef bint first_op = True
    while dereference(c_cigar):
        if b'0' <= dereference(c_cigar) <= b'9':
            # Convert digit char to int and accumulate
            num = num * 10 + (dereference(c_cigar) - 48)  # ord(b'0') = 48
        else:
            op = dereference(c_cigar)
            if first_op and (op == b'S' or op == b'H'):
                start += num
                end += num
            elif op != b'D' and op != b'S' and op != b'H':
                end += num
            if op != b'D':
                template_length += num
            num = 0
            first_op = False
        c_cigar += 1
    return start, end, template_length


cdef int get_align_end_offset(cigar):
    cdef int end = 0
    cdef int num = 0
    cdef char op
    cdef bytes t = cigar.encode('utf8')
    cdef char *c_cigar = t
    while dereference(c_cigar):
        if b'0' <= dereference(c_cigar) <= b'9':
            # Convert digit char to int and accumulate
            num = num * 10 + (dereference(c_cigar) - 48)  # ord(b'0') = 48
        else:
            op = dereference(c_cigar)
            if op not in b'DSH':  # Don't count deletions, or soft/hard clips at right-hand side
                end += num
            num = 0
        c_cigar += 1
    return end


DEF AS_FOUND = 1
DEF NM_FOUND = 2

cdef int sam_to_array(Template template, Params params, tree):
    # Expect read1 and read2 alignments to be concatenated, not mixed together
    data = template.inputdata

    # split the rest of the columns
    template.inputdata = [i[1:-1] + i[-1].split("\t") for i in data]

    # If only one alignment for read1 and read2, no need to try pairing, just send sam to output
    if params.paired_end and len(data) == 2:
        template.passed = 1
        template.outstr = (template.inputdata[0], template.inputdata[1])
        return 1

    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, mq, original_aln_score]
    cdef np.ndarray[np.float_t, ndim=2] arr = np.zeros((len(data), 10))

    chrom_ids = {}

    cdef int cc = 0
    cdef int idx, pos, flag, seq_len, query_start, query_end, start_temp, tc, j
    cdef str chromname, cigar, k, t, v
    cdef float bias = params.bias

    cdef int read1_strand_set, read2_strand_set, current_l
    cdef int first_read2_index = len(template.inputdata) + 1
    cdef int read1_set = 0  # Occasionally multiple primaries, take the longest
    cdef int read2_set = 0

    for idx in range(len(template.inputdata)):

        l = template.inputdata[idx]

        flag = int(l[0])
        pos = int(l[2])
        chromname = l[1]
        if chromname not in chrom_ids:
            chrom_ids[chromname] = cc
            cc += 1

        arr[idx, 0] = chrom_ids[chromname]
        arr[idx, 1] = pos
        arr[idx, 5] = idx
        arr[idx, 6] = -1 if flag & 16 else 1

        arr[idx, 8] = int(l[3])  # mapq

        seq_len = len(l[8])
        tc = 0
        for j in range(11, len(l)):
            tag = l[j]
            if not tc & AS_FOUND and tag.startswith("AS"):
                arr[idx, 9] = float(tag[5:])  # Keep copy of original alignment score
                if tree and intersecter(tree, chromname, pos, pos + 1):
                    arr[idx, 4] = arr[idx, 9] * bias
                else:
                    arr[idx, 4] = arr[idx, 9]
                tc &= AS_FOUND
            elif not tc & NM_FOUND and tag.startswith("NM"):
                arr[idx, 8] = float(tag[5:])
                tc &= NM_FOUND
            if tc == 3:
                break

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
                        logging.error('Could not infer template length')
                        return -1
                        # raise ValueError('Could not infer template length')
                template.read1_length = template_length

                if flag & 16:  # Single end Reverse strand, count from end
                    start_temp = template.read1_length - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

        arr[idx, 2] = query_start
        arr[idx, 3] = query_end

    if params.paired_end:  # Increment the contig position of read 2
        for idx in range(first_read2_index, len(arr)):
            if arr[idx, 7] == 2:  # Second in pair
                arr[idx, 2] += template.read1_length
                arr[idx, 3] += template.read1_length
            if arr[idx, 3] > template.read1_length + template.read2_length:
                logging.error('Inferred template length greater than found read1 + read2 length. Is the file name sorted?', arr[idx, 3], template.read1_length, template.read2_length)

    template.first_read2_index = first_read2_index
    template.data_ori = arr
    template.chrom_ids = chrom_ids

    return 0


cdef void choose_supplementary(Template template):
    # Final alignments have been chosen, but need to decide which is supplementary
    cdef int j = 0
    cdef double read1_max = -1
    cdef double read2_max = -1
    cdef int i = 0
    cdef size_t row_idx
    cdef double[:, :] data_table = template.data_ori
    for row_idx in template.path_result.path:
        if data_table[row_idx, 7] == 1 and data_table[row_idx, 9] > read1_max:  # Use original alignment score
            read1_max = data_table[row_idx, 9]
            template.primary1 = row_idx
        elif data_table[row_idx, 7] == 2 and data_table[row_idx, 9] > read2_max:
            read2_max = data_table[row_idx, 9]
            template.primary2 = row_idx


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
