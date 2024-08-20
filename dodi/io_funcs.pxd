#!python
#cython: language_level=3, boundscheck=False
#cython: profile=False
#distutils: language=c++

from libcpp.vector cimport vector as cpp_vector


cdef class Params:
    cdef:
        public float match_score, mu, sigma, min_aln, max_homology, inter_cost, U
        public float ins_cost, ol_cost, bias, default_max_d
        public bint paired_end, secondary, find_insert_size, modify_mapq, add_tags


cdef struct PathResult:
    cpp_vector[size_t] path
    double path_score, second_best


cdef class Template:
    cdef:
        public object inputdata, data_ori, chrom_ids
        public object name, read1_seq, read2_seq, read1_q, read2_q, outstr
        public int read1_length, read2_length, passed, first_read2_index, primary1, primary2
        public bint read1_reverse, read2_reverse, read1_unmapped, read2_unmapped
        PathResult path_result


from libc.stdio cimport FILE, stdin

cdef extern from "stdio.h":
    FILE *fopen(const char *, const char *)
    int fclose(FILE *)
    ssize_t getline(char **, size_t *, FILE *)
    int fputs(const char *, FILE *)
    size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);


cdef class Writer:
    cdef FILE *outsam
    cdef void set_ptr(self, FILE *out)
    cdef void write(self, bytes byte_string)


cdef class Reader:
    cdef FILE *insam
    cdef char * buffer
    cdef size_t bufsize
    cdef ssize_t read
    cdef void set_ptr(self, FILE *in_sam)
    cdef header_to_file(self, Writer writer, extra_header_line)
    cdef read_line(self)


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


cdef Template make_template(rows)


cdef tuple get_start_end(cigar)


cdef int get_align_end_offset(cigar)


cdef int sam_to_array(Template template, Params params, tree)


cdef void choose_supplementary(Template template)
