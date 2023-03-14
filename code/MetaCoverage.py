from pybedtools import BedTool
from collections import defaultdict
from interval import interval
import numpy as np
import math
import itertools
import pandas as pd

# Methods for computing Metacoverage for all Genes

class MetaCoverage():


################### ANALYTICS #################################################################

    @staticmethod
    def get_meta_exon_bedtool(exon_bed_path):
    # Method generates Bedtool for all 'metaexons' built from all exons in exon bed.
    # Metaexons are constructed from by intersecting non-overlapping exon sequences
    # exons need to be formatted as GeneName_ExonID

        # Generate Gene Exon Dict, i.e. for each gene a dictionary with all exons for that gene
        exon_bed = BedTool(exon_bed_path)
        gene_exon_dict = defaultdict(list)

        # Sort BED Tools Features for each gene
        for b_item in exon_bed:
            gene = b_item.name.split('_')[0]
            gene_exon_dict[gene].append(b_item)

        # Compute consensus intervals for exons, i.e. for each gene the intervals that correspond to exons
        bed_tool_list = []

        # For all exons of each gene:
        #
        for gene, intervals in itertools.islice(gene_exon_dict.items(), None):

            # Define Interval for Gene
            merge_int_old = interval()
            strand = intervals[0].strand
            chrom = intervals[0].chrom

            # Loop over all intervals and generate consensus interval for gene
            # union old (merge) interval with new interval and define set as old,
            # repeat for all intervals of gene

            for i in intervals:
                int_q = interval([i.start, i.end])
                merge_int_new = merge_int_old | int_q
                merge_int_old = merge_int_new

            # Define new 'metaexons' and write to bed_tool_list
            for i_ind, i in enumerate(merge_int_old):
                i_start, i_end = i
                bed_line = (chrom, int(i_start), int(i_end), gene + '_meta_cod_' + str(i_ind), '', strand)
                bed_tool_list.append(bed_line)

        # Generate BedTool for 'metaexons' and
        gene_meta_exon_bedtool = BedTool(bed_tool_list)
        return(gene_meta_exon_bedtool)

    @staticmethod
    def compute_meta_exon_coverage_vectors(query_bam_file,
                                           gene_meta_exon_bedtool):
    # Method computes coverage vectors for each 'MetaExon' in gene_meta_exon_bedtool.
    # Method computation is done using pybedtools.

        #Call bedtools wrapper: set reverse strandedness (S) and individual positions (d)
        meta_exon_coverage = gene_meta_exon_bedtool.coverage(query_bam_file, S=True, d=True, split=True)

        gene_coverage_vectors = defaultdict(list)
        # Define vars for query meta exon and vector for storing counts
        cur_meta_ex = ''
        cur_gene = ''
        cur_strand = ''

        meta_ex_array = np.zeros(0)

        for i in itertools.islice(meta_exon_coverage, None):

            meta_ex = i.name

            # If new meta exon is reached in coverage iterator, add old to defaultdict(list for current gene)
            if meta_ex != cur_meta_ex:
                if not len(meta_ex_array) == 0:
                    # if strand == +, add no.array at end of list, else insert first and reverse order
                    if cur_strand == '+':
                        gene_coverage_vectors[cur_gene].append(meta_ex_array)
                    else:
                        gene_coverage_vectors[cur_gene].insert(0, meta_ex_array[::-1])
                # Generate new meta_ex_array with coords from i
                meta_len = int(i.end - i.start)
                meta_ex_array = np.zeros(meta_len)
                cur_meta_ex = meta_ex
                cur_gene = meta_ex.split('_')[0]
                cur_strand = i.strand

            ex_pos = int(i.fields[-2]) - 1
            ex_count = int(i.fields[-1])
            meta_ex_array[ex_pos] = ex_count
        else:
            if cur_strand == '+':
                gene_coverage_vectors[cur_gene].append(meta_ex_array)
            else:
                gene_coverage_vectors[cur_gene].insert(0, meta_ex_array[::-1])
        return gene_coverage_vectors

    @staticmethod
    def generate_scale_metacoverage_dict(gene_coverage_vectors,
                                         min_cov_per_single_pos = 5,
                                         min_coverage_vec_lenght = 100,
                                         trim5only = False):
    # Method merges coverage vectors for each gene, filters and scales to standard vector for each gene

        gene_scale_cov_vectors = dict()

        for gene, meta_ex_cov_mats in gene_coverage_vectors.items():

            # Transform to one coverage vector
            merge_arr = np.concatenate(meta_ex_cov_mats)

            # Define Max Filter Remove Genes with less than min_cov_per_single_pos counts
            try:
                merge_arr_max = np.max(merge_arr)
            except ValueError:
                continue

            if merge_arr_max < min_cov_per_single_pos:
                continue

            # Trim vector positions from both end of merged vector if they have less than 0.1x merge_arr_max

            for ind in range(len(merge_arr)):
                ind_r = len(merge_arr) - ind - 1
                if merge_arr[ind_r] >= 0.1*merge_arr_max:
                    break

            # Trim Metagene Vector from 5' only, if trim5only FLAG is set
            if trim5only:
                for ind in range(len(merge_arr)):
                    if merge_arr[ind] >= 0.1*merge_arr_max:
                        break

            merge_arr = merge_arr[ind:ind_r]

            # Filter Left Positions with less than 100 nt
            if len(merge_arr) < min_coverage_vec_lenght:
                continue

            scale_array = np.zeros(1000)

            # Define scaling factor
            scale_factor = len(merge_arr) / 1000.0

            # Project merged coverage vector to scale_array
            start_ind = 0

            for ind in range(1000):
                end_ind = math.floor((ind+1)*scale_factor)

                # Compute mean for selected position in coverage vector
                mean_cov_scale_int = np.mean(merge_arr[start_ind:end_ind])
                # Set NaN postions (if empty vector is selected) to 0
                if math.isnan(mean_cov_scale_int):
                    mean_cov_scale_int = 0
                scale_array[ind] = mean_cov_scale_int
                start_ind = end_ind

            gene_scale_cov_vectors[gene] = scale_array
        return gene_scale_cov_vectors


#------------------------------------------------------------------------------------------------------#####
#####--------------------------------------------------------------------------------------------------#####


# Add Method for processing MetaCoverage Vectors
    @staticmethod
    def get_clean_merge_metacoverage_dict(gene_clean_coverage_vectors,
                                            min_cov_per_single_pos = 5,
                                            min_coverage_vec_lenght = 100,
                                            trim5pos = True):

        # Method merges coverage vectors for each gene, filters and scales to standard vector for each gene

        gene_clean_metacov_vec = dict()

        for gene, meta_ex_cov_mats in gene_clean_coverage_vectors.items():

            # Clean individual metaExons and remove non-coverage regions
            # First MetaExon: Clean Only 3', Last Exon: clean only 5'

            metaExLClean = []

            numEx = len(meta_ex_cov_mats)

            for nEx, metaEx in enumerate(meta_ex_cov_mats):
                exMax = np.max(metaEx)

                # Check if first exon
                if nEx == 0:
                    # Trim 3' of exon
                    for ind in range(len(metaEx)):
                        ind_r = len(metaEx) - ind - 1
                        if metaEx[ind_r] >= 0.1 * exMax:
                            break
                    metaEx = metaEx[:ind_r]
                    metaExLClean.append(metaEx)
                # Check if last exon
                elif nEx == (numEx - 1):
                    # Trim 5' of exon
                    for ind in range(len(metaEx)):
                        if metaEx[ind] >= 0.1 * exMax:
                            break
                    metaEx = metaEx[ind:]
                    metaExLClean.append(metaEx)
                else:
                    metaEx = metaEx[metaEx >= 0.1 * exMax]
                    metaExLClean.append(metaEx)

            # Transform to one coverage vector
            merge_arr = np.concatenate(metaExLClean)

            # Define Max Filter Remove Genes with less than min_cov_per_single_pos counts
            try:
                merge_arr_max = np.max(merge_arr)
            except ValueError:
                continue

            if merge_arr_max < min_cov_per_single_pos:
                continue

            # Trim vector positions from both end of merged vector if they have less than 0.1x merge_arr_max
            for ind in range(len(merge_arr)):
                ind_r = len(merge_arr) - ind - 1
                if merge_arr[ind_r] >= 0.1 * merge_arr_max:
                    break

            # Trim 5' end of meta exon if FLAG trim5pos is True
            if trim5pos:
                for ind in range(len(merge_arr)):
                    if merge_arr[ind] >= 0.1 * merge_arr_max:
                        break
            else:
                ind = 0

            merge_arr = merge_arr[ind:ind_r]

            # Filter Left Positions with less than 100 nt
            if len(merge_arr) < min_coverage_vec_lenght:
                continue
            gene_clean_metacov_vec[gene] = merge_arr

        return gene_clean_metacov_vec



    # Add Method for computing scaled MetaCoverage Vectors
    @staticmethod
    def get_scaled_metacoverage_vectors(gene_clean_merged_coverage_vectors):

        gene_scale_metacov_dict = dict()

        for g, mergeMetaCov in gene_clean_merged_coverage_vectors.items():

            scale_array = np.zeros(1000)

            # Define scale factor
            scaleFac = len(mergeMetaCov) / 1000

            # Project merged coverage vector to scale_array
            start_ind = 0

            for ind in range(1000):
                end_ind = math.floor((ind + 1) * scaleFac)

                # Compute mean for selected position in coverage vector
                mean_cov_scale_int = np.mean(mergeMetaCov[start_ind:end_ind])
                # Set NaN postions (if empty vector is selected) to 0
                if math.isnan(mean_cov_scale_int):
                    mean_cov_scale_int = 0
                scale_array[ind] = mean_cov_scale_int
                start_ind = end_ind

            gene_scale_metacov_dict[g] = scale_array

        return gene_scale_metacov_dict


    @staticmethod
    def generate_scale_clean_metacoverage_dict(gene_clean_coverage_vectors,
                                                min_cov_per_single_pos = 5,
                                                min_coverage_vec_lenght = 100):
    # Method merges coverage vectors for each gene, filters and scales to standard vector for each gene

        gene_scale_cov_vectors = dict()

        for gene, meta_ex_cov_mats in gene_coverage_vectors.items():

            # Transform to one coverage vector
            merge_arr = np.concat(meta_ex_cov_mats)

            # Trim Ends of Vector such that

            # Define Max Filter Remove Genes with less than min_cov_per_single_pos counts
            try:
                merge_arr_max = np.max(merge_arr)
            except ValueError:
                continue

            if merge_arr_max < min_cov_per_single_pos:
                continue

            # Trim vector positions from both end of merged vector if they have less than 0.1x merge_arr_max
            for ind in range(len(merge_arr)):
                ind_r = len(merge_arr) - ind - 1
                if merge_arr[ind_r] >= 0.1*merge_arr_max:
                    break

            for ind in range(len(merge_arr)):
                if merge_arr[ind] >= 0.1*merge_arr_max:
                    break

            merge_arr = merge_arr[ind:ind_r]

            # Filter Left Positions with less than 100 nt
            if len(merge_arr) < min_coverage_vec_lenght:
                continue

            scale_array = np.zeros(1000)

            # Define scaling factor
            scale_factor = len(merge_arr) / 1000.0

            # Project merged coverage vector to scale_array
            start_ind = 0

            for ind in range(1000):
                end_ind = math.floor((ind+1)*scale_factor)

                # Compute mean for selected position in coverage vector
                mean_cov_scale_int = np.mean(merge_arr[start_ind:end_ind])
                # Set NaN postions (if empty vector is selected) to 0
                if math.isnan(mean_cov_scale_int):
                    mean_cov_scale_int = 0
                scale_array[ind] = mean_cov_scale_int
                start_ind = end_ind

            gene_scale_cov_vectors[gene] = scale_array
        return gene_scale_cov_vectors




############ WRITER ##############################################

    @staticmethod
    def write_metacoverage_matrix(gene_scale_cov_vectors,
                                  metacov_mat_path):
    # Output MetaCoverage Matrix at metacov_mat_path as .txt file

        gene_meta_cov_file = open(metacov_mat_path, 'w')
        for gene, meta_ex_cov_mats in gene_scale_cov_vectors.items():

            gene_cov = ["{:.1f}".format(cov) for cov in meta_ex_cov_mats]
            cov_line = gene + '\t' + '\t'.join(gene_cov) + '\n'
            gene_meta_cov_file.write(cov_line)
        gene_meta_cov_file.close()

############## WRAPPER ##################################################

    @staticmethod
    def get_metacoverage_matrix(exon_bed_path,
                                input_bam_path,
                                min_cov_per_single_pos = 5,
                                min_coverage_vec_lenght = 100,
                                matrix_path = None):
    # Method parses exon bed file and convert to BedTool.
    # These are converted to 'metaexons'
    # Then computes coverage vs. input BAM file and generates coverage vectors for each metaexon, which are
    # scaled and written to matrix_path if specified, else written

        gene_meta_exon_bedtool = MetaCoverage.get_meta_exon_bedtool(exon_bed_path=exon_bed_path)

        meta_exon_coverage_vectors = MetaCoverage.compute_meta_exon_coverage_vectors(query_bam_file=input_bam_path,
                                                                                     gene_meta_exon_bedtool=gene_meta_exon_bedtool)

        scaled_gene_metacoverge_dict = MetaCoverage.generate_scale_metacoverage_dict(gene_coverage_vectors=meta_exon_coverage_vectors,
                                                                                     min_cov_per_single_pos=min_cov_per_single_pos,
                                                                                     min_coverage_vec_lenght=min_coverage_vec_lenght)

        if not matrix_path is None:
            MetaCoverage.write_metacoverage_matrix(gene_scale_cov_vectors=scaled_gene_metacoverge_dict,
                                                   metacov_mat_path=matrix_path)
        else:
            return scaled_gene_metacoverge_dict