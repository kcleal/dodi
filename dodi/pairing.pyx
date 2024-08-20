#cython: language_level=3, boundscheck=False, wraparound=False
# cython: profile=False
#distutils: language=c++

"""
Find the best path through a collection of alignments
Works with paried-end reads or single contigs. Borrows the pairing heuristic from bwa.
"""
from __future__ import absolute_import
import numpy as np
cimport numpy as np
from libc.math cimport exp, log, sqrt, fabs
from libcpp.vector cimport vector
from dodi.io_funcs cimport Template, Params, PathResult
import time
from sys import stderr


ctypedef np.float_t DTYPE_t


cdef float erfcc(float x):
    """Complementary error function."""
    cdef float z, t, r, p1, p2, p3, p4, p5, p6, p7, p8, p9
    z = fabs(x)
    t = (1. / (1. + 0.5*z))
    p1 = -.82215223+t*.17087277
    p2 = 1.48851587+t*p1
    p3 = -1.13520398+t*p2
    p4 = .27886807+t*p3
    p5 = -.18628806+t*p4
    p6 = .09678418+t*p5
    p7 = .37409196+t*p6
    p8 = 1.00002368+t*p7
    p9 = 1.26551223+t*p8
    r = t * exp(-z*z-p9)
    if x >= 0.:
        return r
    else:
        return 2. - r


cdef float normcdf(float x, float mu, float sigma):
    cdef float t, y
    t = x-mu
    y = 0.5*erfcc(-t/(sigma*sqrt(2.0)))
    if y > 1.0:
        y = 1.0
    return y


cdef int is_proper_pair(float d, float pos1, float pos2, float strand1, float strand2, float mu, float sigma):
    if d <= (mu + 4*sigma):
        if pos1 < pos2 and strand1 == 1. and strand2 == -1.:
            return 1
        if pos2 < pos1 and strand2 == 1. and strand1 == -1.:
            return 1
    return 0


cdef float bwa_pair_score(float d, float mu, float sigma, float match_score, float u):
    """
    Calculate the bwa pairing cost
    :param d: seperation distance
    :param mu: insert size mean
    :param sigma: insert size std
    :param match_score: score gained from a matched base
    :param u: constant parameter, worst cost possible
    :return: The bwa pairing cost (float), whether the read is FR (int)
    """
    cdef float prob, c
    prob = (1 - normcdf(d, mu, sigma))
    c = -match_score * (log(prob)/log(4))
    if c < u:
        return c
    else:
        return u


cdef void optimal_path(double[:, :] segments,
                        float contig_length,
                        float mu,
                        float sigma,
                        float min_aln,
                        float max_homology,
                        float ins_cost,
                        float hom_cost,
                        float inter_cost,
                        float U,
                        float match_score,
                        int paired_end,
                        long[:] pred,
                        double[:] node_scores,
                        double[:] nb_node_scores,
                        size_t n_rows,
                        PathResult& result
                        ):
    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.
    # Use a special start and end node to score the first and last alignments

    cdef:
        double chr1, pos1, start1, end1, score1, row_index1, strand1, r1
        double chr2, pos2, start2, end2, score2, row_index2, strand2, r2
        double micro_h, ins, best_score, next_best_score, best_normal_orientation, current_score, total_cost
        double S, path_score, jump_cost, distance

        int proper_pair
        size_t i, j, p, normal_end_index

    # Deal with first score
    for i in range(n_rows):
        node_scores[i] = segments[i, 4] - (segments[i, 2] * ins_cost)
        pred[i] = -1
        nb_node_scores[i] = -1e6  # Must set to large negative, otherwise a value of zero can imply a path to that node

    best_score = 0  # Declare here in case only 1 alignment
    next_best_score = 0
    best_normal_orientation = 0  # Keep track of the best normal pairing score, F first R second
    normal_end_index = -1  # Keep track of the last normal-index for updating the normal-score later on

    for i in range(1, n_rows):

        chr1 = segments[i, 0]
        pos1 = segments[i, 1]
        start1 = segments[i, 2]
        end1 = segments[i, 3]
        score1 = segments[i, 4]
        row_index1 = segments[i, 5]
        strand1 = segments[i, 6]
        r1 = segments[i, 7]

        p = -1  # -1 means there is no predecessor
        best_score = score1 - (start1 * ins_cost)  # Implies all preceding alignments skipped!
        next_best_score = - (start1 * ins_cost)  # Worst case

        # Walking backwards
        for j in range(i-1, -1, -1):
            chr2 = segments[j, 0]
            pos2 = segments[j, 1]
            start2 = segments[j, 2]
            end2 = segments[j, 3]
            score2 = segments[j, 4]
            row_index2 = segments[j, 5]
            strand2 = segments[j, 6]
            r2 = segments[j, 7]

            # Allow alignments with minimum sequence and max overlap
            if start1 > end2 - max_homology and end1 - end2 > min_aln:

                # Microhomology and insertion lengths between alignments on same read only
                micro_h = 0
                ins = 0
                micro_h = end2 - start1
                if micro_h < 0:
                    ins = fabs(micro_h)
                    micro_h = 0

                # Define jump cos
                proper_pair = 0
                if chr1 == chr2:
                    if r1 == r2:  # If the jump occurs on the same read, means there is an SV
                        jump_cost = U

                    else:
                        distance = fabs(pos1 - pos2)
                        proper_pair = is_proper_pair(distance, pos1, pos2, strand1, strand2, mu, sigma)
                        if proper_pair:
                            jump_cost = bwa_pair_score(distance, mu, sigma, match_score, U)
                        else:
                            jump_cost = U

                else:
                    jump_cost = inter_cost + U

                total_cost = (micro_h * hom_cost) + (ins * ins_cost) + jump_cost
                S = score1 - total_cost
                current_score = node_scores[j] + S

                # Update best score and next best score
                if current_score > best_score:
                    next_best_score = best_score
                    best_score = current_score
                    p = j

                elif current_score > next_best_score:
                    next_best_score = current_score

                if proper_pair:
                    if current_score > best_normal_orientation:
                        best_normal_orientation = current_score
                        normal_end_index = i  # index i is towards the right-hand side of the array

        node_scores[i] = best_score
        pred[i] = p
        nb_node_scores[i] = next_best_score

    # Update the score for jumping to the end of the sequence
    # Basically calculates what the best and secondary scores are for the 'virtual' end node
    cdef float node_to_end_cost
    cdef float secondary = 0

    path_score = -(contig_length * ins_cost)  # Worst case

    cdef int end_i = -1

    for i in range(n_rows):
        node_to_end_cost = node_scores[i] - (ins_cost * (contig_length - segments[i, 3]))
        node_scores[i] = node_to_end_cost  # Need to know for secondary path
        if node_to_end_cost > path_score:
            path_score = node_to_end_cost
            end_i = i
        elif node_to_end_cost > secondary:
            secondary = node_to_end_cost
        if i == normal_end_index:  # Update the best normal score, deals with right-hand-side soft-clips
            best_normal_orientation = node_to_end_cost

    # Get the path using traceback
    cdef:
        float s, potential_secondary
        int last_i, next_i, path_length, current_i

    # do one pass to get path length, then another to fill path array
    path_length = 0
    last_i = end_i
    if end_i != -1:
        path_length += 1
        while True:
            next_i = pred[last_i]
            if next_i == -1:
                break
            path_length += 1
            last_i = next_i

        result.path.resize(path_length)

        current_i = 0
        last_i = end_i
        while True:
            next_i = pred[last_i]
            result.path[path_length - 1 - current_i] = <size_t> segments[last_i, 5]
            if next_i == -1:
                break
            last_i = next_i
            current_i += 1

    if secondary < 0:
        secondary = 0
    if path_score < 0:
        path_score = 0
    if best_normal_orientation < 0:
        best_normal_orientation = 0

    result.path_score = path_score
    result.second_best = secondary


def sort_func(row):
    return row[7], row[2], -row[4]


cdef void pairing_process(Template rt, Params params):
    """
    Assumes that the reads are ordered read1 then read2, in the FR direction
    :param rt: Read_template object, contains all parameters within the pairing_params array
    :returns 5-tuple: list path, float length, float second_best, float dis_to_normal, float norm_pairings
    """

    cdef bint paired_end = params.paired_end
    t = np.array(sorted(rt.data_ori, key=sort_func))

    cdef size_t n_rows = t.shape[0]

    cdef np.ndarray[np.int64_t, ndim=1] pred = np.empty(n_rows, dtype=np.int64)
    cdef np.ndarray[np.float64_t, ndim=1] node_scores = np.empty(n_rows, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] nb_node_scores = np.empty(n_rows, dtype=np.float64)  # Next best node score, for finding the secondary path

    rt.passed = True

    cdef PathResult path1_res, path2_res

    if not paired_end:
        if rt.read1_unmapped:
            rt.path_result.path.push_back(0)
            return
        optimal_path(t, rt.read1_length, params.mu, params.sigma, params.min_aln, params.max_homology,
                            params.ins_cost, params.ol_cost, params.inter_cost, params.U, params.match_score, 0,
                                                                                    pred, node_scores, nb_node_scores, n_rows, rt.path_result)

    else:
        if (rt.read1_length == 0 and rt.read1_length == 0) or (rt.read1_unmapped and rt.read2_unmapped):
            rt.path_result.path.push_back(0)
            rt.path_result.path.push_back(rt.first_read2_index)

        elif rt.read2_unmapped:
            read1_arr = t[t[:, 7] == 1]
            optimal_path(read1_arr, rt.read1_length, params.mu, params.sigma, params.min_aln, params.max_homology,
                            params.ins_cost, params.ol_cost, params.inter_cost, params.U, params.match_score, 1,
                                                                                        pred, node_scores, nb_node_scores, len(read1_arr), rt.path_result)

        elif rt.read1_unmapped:
            read2_arr = t[t[:, 7] == 2]
            optimal_path(read2_arr, rt.read2_length, params.mu, params.sigma, params.min_aln, params.max_homology,
                            params.ins_cost, params.ol_cost, params.inter_cost, params.U, params.match_score, 1,
                                                                                        pred, node_scores, nb_node_scores, len(read2_arr), rt.path_result)

        else:
            optimal_path(t, rt.read1_length + rt.read2_length, params.mu, params.sigma, params.min_aln, params.max_homology,
                            params.ins_cost, params.ol_cost, params.inter_cost, params.U, params.match_score, 1, pred, node_scores, nb_node_scores, n_rows, path1_res)
            if path1_res.path.size() == n_rows:
                rt.path_result = path1_res
                return

            # put read 2 first, and read 1 second
            t[:rt.first_read2_index, 2:4] += rt.read2_length
            t[rt.first_read2_index:, 2:4] -= rt.read1_length
            t = t[(-t[:, 7]).argsort(kind='stablesort')]

            optimal_path(t, rt.read1_length + rt.read2_length, params.mu, params.sigma, params.min_aln, params.max_homology,
                            params.ins_cost, params.ol_cost, params.inter_cost, params.U, params.match_score, 1,
                                                                                        pred, node_scores, nb_node_scores, n_rows, path2_res)

            if path2_res.path_score > path1_res.path_score:
                rt.path_result = path2_res
            else:
                rt.path_result = path1_res
