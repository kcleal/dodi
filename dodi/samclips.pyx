#cython: language_level=3
# cython: profile=False

"""
Utils to generate proper sam output and flag information
"""
from __future__ import absolute_import
import re
import click
from dodi import io_funcs
from dodi.io_funcs cimport get_align_end_offset, Template, Params
import logging


def echo(*arg):
    click.echo(arg, err=True)


cdef int set_bit(int v, int index, int x):
    """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
    cdef int mask
    mask = 1 << index  # Compute mask, an integer with just bit 'index' set.
    v &= ~mask  # Clear the bit indicated by the mask (if x is False)
    if x:
        v |= mask  # If x was True, set the bit indicated by the mask.
    return v


cdef list set_tlen(out):

    pri_1 = out[0][1]
    pri_2 = out[1][1]

    cdef int flg1 = pri_1[0]
    cdef int flg2 = pri_2[0]

    cdef int tlen1, tlen2, t1, t2, p1_pos, p2_pos
    if flg1 & 12 or flg2 & 12 or pri_1[1] != pri_2[1]:  # Read or mate unmapped, translocation
        tlen1 = 0
        tlen2 = 0
        t1 = 0
        t2 = 0

    else:
        p1_pos = int(pri_1[2])
        p2_pos = int(pri_2[2])

        # Use the end position of the alignment if read is on the reverse strand, or start pos if on the forward
        if flg1 & 16:  # Read rev strand
            aln_end1 = get_align_end_offset(pri_1[4])
            t1 = p1_pos + aln_end1
        else:
            t1 = p1_pos

        if flg2 & 16:
            aln_end2 = get_align_end_offset(pri_2[4])
            t2 = p2_pos + aln_end2
        else:
            t2 = p2_pos

        if t1 <= t2:
            tlen1 = t2 - t1  # Positive
            tlen2 = t1 - t2  # Negative
        else:
            tlen1 = t2 - t1  # Negative
            tlen2 = t1 - t2  # Positive

    pri_1[7] = str(tlen1)
    pri_2[7] = str(tlen2)

    out2 = [(out[0][0], pri_1, out[0][2]), (out[1][0], pri_2, out[1][2])]

    # Set tlen's of supplementary
    cdef int sup_flg, sup_pos
    for sup_tuple in out[2:]:
        sup_tuple = list(sup_tuple)
        sup_flg = sup_tuple[1][0]
        sup_chrom = sup_tuple[1][1]
        sup_pos = int(sup_tuple[1][2])

        sup_end = get_align_end_offset(sup_tuple[1][4])
        if sup_flg & 16:  # If reverse strand, count to end
            sup_pos += sup_end

        if sup_flg & 64:  # First in pair, mate is second
            other_end = t2
            other_chrom = pri_2[1]
            other_flag = pri_2[0]
        else:
            other_end = t1
            other_chrom = pri_1[1]
            other_flag = pri_1[0]
        # This is a bit of a hack to make the TLEN identical to bwa
        # Make sure they are both on same chromosome
        if sup_chrom == other_chrom:
            if sup_pos < other_end:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end
            else:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end

            sup_tuple[1][7] = str(tlen)
        out2.append(tuple(sup_tuple))

    return out2


cdef set_mate_flag(a, b, max_d, read1_rev, read2_rev):

    if not a or not b:  # No alignment, mate unmapped?
        return False, False

    # Make sure chromosome of mate is properly set not "*"
    chrom_a, mate_a = a[2], a[5]
    chrom_b, mate_b = b[2], b[5]
    if chrom_a != mate_b:
        b[5] = chrom_a
    if chrom_b != mate_a:
        a[5] = chrom_b

    cdef int aflag = a[0]
    cdef int bflag = b[0]

    cdef bint reverse_A = False
    cdef bint reverse_B = False

    # If set as not primary, and has been aligned to reverse strand, and primary is mapped on forward
    # the sequence needs to be rev complement
    if aflag & 256:
        if (aflag & 16) and (not read1_rev):
            reverse_A = True
        elif (not aflag & 16) and read1_rev:
            reverse_A = True

    if bflag & 256:
        if (bflag & 16) and (not read2_rev):
            reverse_B = True
        elif (not bflag & 16) and read2_rev:
            reverse_B = True

    # Turn off proper pair flag, might be erroneously set
    aflag = set_bit(aflag, 1, 0)  # Bit index from 0
    bflag = set_bit(bflag, 1, 0)

    # Turn off supplementary pair flag
    aflag = set_bit(aflag, 11, 0)
    bflag = set_bit(bflag, 11, 0)

    # Set paired
    aflag = set_bit(aflag, 0, 1)
    bflag = set_bit(bflag, 0, 1)

    # Set first and second in pair, in case not set
    aflag = set_bit(aflag, 6, 1)
    bflag = set_bit(bflag, 7, 1)

    # Turn off any mate reverse flags, these should be reset
    aflag = set_bit(aflag, 5, 0)
    bflag = set_bit(bflag, 5, 0)

    # If either read is unmapped
    if aflag & 4:
        bflag = set_bit(bflag, 3, 1)  # Position 3, change to 1
    if bflag & 4:
        aflag = set_bit(aflag, 3, 1)

    # If either read on reverse strand
    if aflag & 16:
        bflag = set_bit(bflag, 5, 1)
    if bflag & 16:
        aflag = set_bit(aflag, 5, 1)

    # Set unmapped
    arname = a[1]
    apos = a[2]
    if apos == "0":  # -1 means unmapped
        aflag = set_bit(aflag, 2, 1)
        bflag = set_bit(bflag, 8, 1)

    brname = b[1]
    bpos = b[2]
    if b[2] == "0":
        bflag = set_bit(bflag, 2, 1)
        aflag = set_bit(aflag, 8, 1)

    # Set RNEXT and PNEXT
    a[5] = brname
    a[6] = bpos

    b[5] = arname
    b[6] = apos

    cdef int p1, p2
    if not (apos == "-1" or bpos == "-1"):

        if arname == brname:

            p1, p2 = int(apos), int(bpos)

            # Set proper-pair flag
            # echo(aflag, bflag, p1, p2)
            # echo((p1 < p2, (not aflag & 16), (bflag & 16)), (p2 <= p1, (not bflag & 16), (aflag & 16)))
            if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

                if abs(p1 - p2) < max_d:
                    # Check for FR or RF orientation
                    if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                        aflag = set_bit(aflag, 1, 1)
                        bflag = set_bit(bflag, 1, 1)

                        # If proper pair, sometimes the mate-reverse-strand flag is set
                        # this subsequently means the sequence should be reverse complemented
                        if aflag & 16 and not bflag & 32:
                            # Mate-reverse strand not set
                            bflag = set_bit(bflag, 5, 1)

                        if not aflag & 16 and bflag & 32:
                            # Mate-reverse should'nt be set
                            bflag = set_bit(bflag, 5, 0)
                            reverse_A = True

                        if bflag & 16 and not aflag & 32:
                            # Mate-reverse strand not set
                            aflag = set_bit(aflag, 5, 1)

                        if not bflag & 16 and aflag & 32:
                            # Mate-revsere should'nt be set
                            aflag = set_bit(aflag, 5, 0)
                            reverse_B = True

    a[0] = aflag
    b[0] = bflag

    return reverse_A, reverse_B, a, b


cdef set_supp_flags(sup, pri, bint ori_primary_reversed, bint primary_will_be_reversed):

    # Set paired
    cdef int supflag = sup[0]
    cdef int priflag = pri[0]

    # Set paired and supplementary flag
    if not supflag & 1:
        supflag = set_bit(supflag, 0, 1)
    if not supflag & 2048:
        supflag = set_bit(supflag, 11, 1)

    # If primary is on reverse strand, set the mate reverse strand tag
    if priflag & 16 and not supflag & 32:
        supflag = set_bit(supflag, 5, 1)
    # If primary is on forward srand, turn off mate rev strand
    if not priflag & 16 and supflag & 32:
        supflag = set_bit(supflag, 5, 0)

    # Turn off not-primary-alignment
    if supflag & 256:
        supflag = set_bit(supflag, 8, 0)

    cdef bint rev_sup = False

    if ori_primary_reversed:
        if not supflag & 16:  # Read on forward strand
            rev_sup = True

    elif supflag & 16:  # Read on reverse strand
        if not ori_primary_reversed:
            rev_sup = True

    sup[0] = supflag
    sup[5] = pri[1]
    sup[6] = pri[2]

    return rev_sup


cdef void set_supp_flag_single(sup, pri):
    supflag = sup[0]
    priflag = pri[0]

    if not supflag & 2048:
        supflag = set_bit(supflag, 11, 1)

    # Turn off not-primary-alignment
    if supflag & 256:
        supflag = set_bit(supflag, 8, 0)

    sup[0] = supflag


cdef add_sequence_back(item, reverse_me, template):
    # item is the alignment
    cdef int flag = item[0]
    cigar = item[4]
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    cdef int i, l
    cdef str opp
    cdef int string_length = 0
    cdef int hard_clip_length = 0
    for i in range(0, len(c), 2):
        l = int(c[i])
        opp = c[i + 1]
        if opp != "D":
            if opp == "H":
                hard_clip_length += l
            else:
                string_length += l
    cdef int cigar_length = string_length + hard_clip_length
    if flag & 64:  # Read1
        # sometimes a read seq can be replaced by "*" if the read is a duplicate of its mate, so use mate!
        if flag & 4 and not template.read1_seq:
            seq = template.read2_seq
            q = template.read2_q
        else:
            seq = template.read1_seq
            q = template.read1_q
    elif flag & 128:
        if flag & 4 and not template.read2_seq:
            seq = template.read1_seq
            q = template.read1_q
        else:
            seq = template.read2_seq
            q = template.read2_q
    else:
        seq = template.read1_seq  # Unpaired
        q = template.read1_q

    cdef int end = len(seq)
    if not seq:
        return item, False

    if cigar_length == len(seq) or cigar == "*":
        item[4] = item[4].replace("H", "S")
        item[8] = seq
        if q:
            item[9] = q
        return item, True

    return item, False


cdef list replace_sa_tags(alns):
    """
    because original sa tags will contain all original alignments, re-generate the SA tag content with only the 
    alignments remaining following dodi filtering
    :param alns: list containing all alignment info for all alignments
    :return: same info in same format, with relevant SAs removed from SA tags
    """
    # if there are no alignments, return
    if len(alns) == 0:
        return alns

    # if any of the alignments are supplementary, then re-generate the SA tags to ensure relevant info is removed
    if any([i[0] == "sup" for i in alns]):
        sa_tags = {}  # Read1: tag, might be multiple split alignments
        alns2 = []
        for i, j, k in alns:  # i = aln type (e.g. sup/pri), j = aln info as a list, k = Bool, whether to reverse (???)
            # logging.debug(f"{i}, {k}")
            # Remove any SA tags in alignment, might be wrong
            j = [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and not item.startswith("SA"))]

            # extract relevant information to re-create SA tags
            flag = j[0]
            chrom = j[1]
            pos = j[2]
            mapq = j[3]
            cigar = j[4]
            nm = 0

            # get the nm tag info
            for tg in j[10:]:
                if tg.startswith("NM"):
                    nm = tg[5:]
                    break

            # set strand information
            strand = "-" if flag & 16 else "+"

            # create the sa information string
            # in sam files, all SA tags end in ';', therefore it can be added at this stage
            sa = f"{chrom},{pos},{strand},{cigar},{mapq},{nm};"

            # determine whether paired end, and whether first in pair
            first_in_pair = flag & 64

            # debug statements for dev
            logging.debug(first_in_pair)
            logging.debug(sa)

            # use first in pair flag info (0 or 1) to correctly populate SA tags for each read in a pair
            # for single end read data, there will only be one key (all items enter the same element)
            if first_in_pair in sa_tags:
                sa_tags[first_in_pair].append(sa)
            else:
                sa_tags[first_in_pair] = [sa]

            # append the alignment list object, plus the SA info for the alignment, to the list of alns for the read
            alns2.append([i, j, k, sa, first_in_pair])

        # Now add SA tag info back in, excluding self (SA tag order will not necessarily be preserved from input file)
        out = []
        for i, j, k, sa, fip in alns2:
            # use first in pair info to access the relevant dict item
            if fip in sa_tags:
                # create SA tag string from list of SAs (excluding self), without altering the reference list
                sa_string = ""
                for each in sa_tags[fip]:
                    # if not self, then append
                    if not each == sa:
                        sa_string += each

            # if any non-self alignments, add SA tag to alignment
            if sa_string:
                # remove any newline characters from the ends of existing alignment file lines, e.g. XA tags
                if j[-1][-1] == "\n":
                    j[-1] = j[-1].rstrip('\n')
                # append the SA string to the list of alignment file fields
                j.append(f"SA:Z:{sa_string}\n")
            # append the edited alignment list to the list of all alignments for the read, to be returned
            out.append((i, j, k))
        return out
    else:
        # if there aren't any supplementary alignments for the read, then just remove the SA tag
        # (unlikely to exist anyway), and return
        return [(i, [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and not item.startswith("SA"))], ii) for i, j, ii in alns]


cdef list replace_mc_tags(alns):

    # Replace MC mate cigar tag if set
    if len(alns) <= 1:
        return alns

    cdef int i, last

    a = alns[0][1]

    if alns[1][0] != "sup":
        b = alns[1][1]
    else:
        return alns

    read1_cigar = a[4] if a[0] & 64 else b[4]
    read2_cigar = b[4] if b[0] & 128 else a[4]

    for count, (ps, a, _) in enumerate(alns):

        last = len(a) - 1
        for i in range(10, len(a)):
            if a[i].startswith("MC"):
                if i == last:
                    if a[0] & 64:
                        a[i] = f"MC:Z:{read2_cigar}\n"
                    else:
                        a[i] = f"MC:Z:{read1_cigar}\n"
                else:
                    if a[0] & 64:
                        a[i] = f"MC:Z:{read2_cigar}"
                    else:
                        a[i] = f"MC:Z:{read1_cigar}"
                break
    return alns


# cdef joint_mapq(float mapQ1, float mapQ2):
#     # Convert individual mapQs to probabilities
#     cdef double p_incorrect1 = 10 ** (-mapQ1 / 10)
#     cdef double p_incorrect2 = 10 ** (-mapQ2 / 10)
#
#     # Joint probability of both reads being correct
#     cdef double p_joint_correct = 1 - (p_incorrect1 * p_incorrect2)
#
#     # Convert back to PHRED scale
#     if p_joint_correct >= 0.999999:
#         return 60
#     return int(-10 * math.log10(1 - p_joint_correct))


cdef modify_split_read_mapq(template, primary1, primary2, out, max_d):
    if len(out) == 0:
        return

    # pq1 = float(primary1[3])
    # pq2 = float(primary2[3])

    for _, sup, _ in out:
        flag_sup = sup[0]

        if flag_sup & 64:  # supplementary is first-in-pair
            if abs(int(primary2[2]) - int(sup[2])) < max_d:
                # flag_pri = primary1[0]
                # only apply to forward-reverse?
                # pp = False
                # if primary1.pos < sup.pos:
                #     if not flag_pri & 16 and flag_sup & 16:
                #         # proper pair
                #         pp = True
                #
                # else:
                #     if not flag_sup & 16 and flag_pri & 16:
                #         # proper pair
                #         pp = True
                # if pp:

                # prob that both a wrong

                # q2 = int(sup[3])
                # joint = str(joint_mapq(pq1, q2))
                # primary2[3] = joint
                # sup[3] = joint

                # Set minimum to avg method
                # avg_mapq = int((int(primary2[3]) + int(sup[3])) / 2)
                # primary2[3] = str(max(int(primary2[3]), avg_mapq))
                # sup[3] = str(max(int(sup[3]), avg_mapq))

                # Take the max method
                max_mapq = str(max(int(primary2[3]), int(sup[3])))
                primary2[3] = max_mapq
                sup[3] = max_mapq

        else:
            if abs(int(primary1[2]) - int(sup[2])) < max_d:
                # q2 = int(sup[3])
                # joint = str(joint_mapq(pq1, q2))
                # primary1[3] = joint
                # sup[3] = joint

                # avg_mapq = int((int(primary1[3]) + int(sup[3])) / 2)
                # primary1[3] = str(max(int(primary1[3]), avg_mapq))
                # sup[3] = str(max(int(sup[3]), avg_mapq))

                max_mapq = str(max(int(primary1[3]), int(sup[3])))
                primary1[3] = max_mapq
                sup[3] = max_mapq


cdef list fixsam(Template template, Params params):

    cdef int j, row_idx

    max_d = params.default_max_d
    paired = False if template.read2_length is 0 else True
    out = []
    primary1 = None
    primary2 = None
    rev_A = False
    rev_B = False

    inputdata = template.inputdata
    tabledata = template.data_ori

    primary1_idx = template.primary1
    primary2_idx = template.primary2

    # for row_idx in template.rows:
    for row_idx in template.path_result.path:

        l = inputdata[row_idx]

        l[0] = int(l[0])  # Convert flag to int

        strand = "-1" if l[0] & 16 else "1"
        rid = str(2 if l[0] & 128 else 1)
        xs = int(tabledata[row_idx, 4])  # the biased alignment score

        if row_idx == primary1_idx:
            # set flag to pri
            flag = set_bit(l[0], 11, 0)
            l[0] = set_bit(flag, 8, 0)
            primary1 = l

        elif row_idx == primary2_idx:
            flag = set_bit(l[0], 11, 0)
            l[0] = set_bit(flag, 8, 0)
            primary2 = l
        else:
            out.append(['sup', l, False])  # Supplementary, False to decide if rev comp

    if (primary1 is None or primary2 is None) and params.paired_end:
        if primary1 is None:
            primary1 = template.inputdata[0]
            primary1[0] = int(primary1[0])
        if primary2 is None:
            primary2 = template.inputdata[template.first_read2_index]  # unmapped
            primary2[0] = int(primary2[0])


    if paired and params.paired_end:

        # rev_A/B are set to true/false indicating if the primary aligns should eb reverse complemented
        rev_A, rev_B, primary1, primary2 = set_mate_flag(primary1, primary2, max_d, template.read1_reverse, template.read2_reverse)

        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if out[i][1][0] & 64:  # First in pair  Note primary2 and primary1 order
                revsup = set_supp_flags(out[i][1], primary2, template.read1_reverse, rev_A)
            else:
                revsup = set_supp_flags(out[i][1], primary1, template.read2_reverse, rev_B)
            if revsup:
                out[i][2] = True

        # increase mapq of supplementary to match primary if possible
        if params.modify_mapq:
            modify_split_read_mapq(template, primary1, primary2, out, params.default_max_d)

    if primary1 is None:
        return []

    if params.paired_end:
        out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out
        out = set_tlen(out)

    else:
        pri_is_reverse = bool(primary1[0] & 16)
        for i in range(len(out)):
            if bool(out[i][1][0] & 16) != pri_is_reverse:
                out[i][2] = True
            set_supp_flag_single(out[i][1], primary1)
        out = [('pri', primary1, rev_A)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary

    for a_type, aln, reverse_me in out:
        if aln:  # None here means no alignment for primary2
            # Do for all supplementary
            if aln[8] == "*" or "H" in aln[4]:  # Sequence might be "*", needs adding back in
                aln, success = add_sequence_back(aln, reverse_me, template)
                if reverse_me and success:
                    aln[8] = io_funcs.reverse_complement(str(aln[8]), len(aln[8]))
                    aln[9] = aln[9][::-1]

            # Turn off not primary here
            aln[0] = set_bit(aln[0], 8, 0)

    out = replace_sa_tags(out)
    out = replace_mc_tags(out)

    # Set discordant flag on supplementary, convert flags back to string
    for j in range(len(out)):

        if out[j][1][-1][-1] != "\n":
            out[j][1][-1] += "\n"
        # Set for discordant
        if out[j][0] == "sup":
            rec = out[j][1]
            flag = rec[0]
            # Outside max dist, same chrom, or same strand
            if abs(int(rec[7])) >= max_d or rec[1] != rec[5] or bool(flag & 16) == bool(flag & 32):
                flag = set_bit(flag, 1, 0)
                out[j][1][0] = flag

        if out[j][1][2] == '0':
            flag = int(out[j][1][0])
            flag = set_bit(flag, 2, 1)  # set unmapped
            out[j][1][0] = flag
            out[j][1][1] = '*'

        out[j][1][0] = str(out[j][1][0])

    return [i[1] for i in out if i[1] != 0]

