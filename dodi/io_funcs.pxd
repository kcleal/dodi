#!python
#cython: language_level=3, boundscheck=False
#cython: profile=False
#distutils: language=c++

from libcpp.vector cimport vector as cpp_vector


# cdef class Template:
#     cdef public float max_d, match_score, bias
#     cdef public int read1_length, read2_length, first_read2_index
#
#     cdef public bint read1_unmapped, read2_unmapped, read1_reverse, read2_reverse, passed, paired_end, secondary
#     cdef public object pairing_params, score_mat, name, last_seen_chrom, inputdata, r1_len, r2_len, isize, outstr, \
#         chrom_ids, read1_q, read2_q, primary1, primary2
#     cdef public object read1_seq, read2_seq
#     cdef public double[:, :] data, data_ori
#     cdef public long[:] rows
#
#     # pairing params
#     cdef public float min_aln, max_homology, inter_cost, U, zero_cost_bound, max_gap_cost, mu, sigma


cdef class Params:
    cdef public float match_score, mu, sigma, min_aln, max_homology, inter_cost, U, ins_cost, ol_cost, bias, default_max_d
    cdef public bint paired_end, secondary, find_insert_size, modify_mapq, add_tags


cdef class Template:
    cdef public object inputdata, score_mat, name, last_seen_chrom, read1_seq, read2_seq, read1_q, read2_q, outstr,\
        data_ori, chrom_ids, rows
    # cdef public float max_d
    cdef public int read1_length, read2_length, passed, first_read2_index, primary1, primary2
    cdef public bint read1_reverse, read2_reverse, read1_unmapped, read2_unmapped


cdef extern from "itree_wrapper.h":
    ctypedef struct Interval:
        int low
        int high


cdef extern from "itree_wrapper.h":
    cdef cppclass BasicIntervalTree:
        BasicIntervalTree()
        void add(int, int, int)
        bint searchInterval(int, int)
        Interval* overlappingInterval(int, int)
        void index()
        void allOverlappingIntervals(int, int, cpp_vector[int])
        int countOverlappingIntervals(int, int)


cdef class Py_BasicIntervalTree:
    cdef BasicIntervalTree *thisptr

    cpdef void add(self, int start, int end, int index)
    cpdef bint searchInterval(self, int pos, int pos2)
    cpdef overlappingInterval(self, int pos, int pos2)
    cpdef void index(self)
    cpdef allOverlappingIntervals(self, int pos, int pos2)
    cpdef int countOverlappingIntervals(self, int pos, int pos2)
