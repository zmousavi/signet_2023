import numpy as np
import csv
import yaml
from collections import defaultdict
import warnings
import pandas as pd
import os
import pdb
import utils_data
import utils_ppi
import gene_selection
import argparse
import re

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_args():
    # ,epilog='small call: see run.sh', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(description='Select genes')
    parser.add_argument('-c', '--config', help='config file',
                        required=False, type=str, default='./run_geneselection.yaml')
    parser.add_argument('-p', '--phenotype_list', help='phenotype list',
                        required=False, nargs='*', default=['QT', 'HR', 'QRS', 'PR', 'JT'])
    parser.add_argument('-n', '--project_dirname', help='project name',
                        required=False, type=str, default='TopMed')
    parser.add_argument('-k', '--k_closestgenes', help='number of closest genes',
                        required=False, type=int, default=1)
    parser.add_argument('-d', '--flankdistance', help='distance range from which to pick genes',
                        required=False, type=int, default=250000)
    parser.add_argument('-a', '--alpha', help='pathway score significance',
                        required=False, type=float, default=1e-10)
    parser.add_argument('-s0', '--seed_start', help='seed0',
                        required=False, type=int, default=0)
    parser.add_argument('-s1', '--seed_end', help='seed1',
                        required=False, type=int, default=1)
    parser.add_argument('--nodegcor', help='degree corrected null',
                        action='store_true')  # default is False
    parser.add_argument(
        '--plot', help='plot network', action='store_true')  # default is False
    parser.add_argument('--nooptim', help='degree corrected null',
                        action='store_true')  # default is False
    parser.add_argument('--notTF', help='transcription factor',
                        action='store_false')  # default is True
    parser.add_argument('--unannotated', help='include unannotated genes',
                        action='store_true')  # default is False
    parser.add_argument('--pseudogene', help='include pseudogene',
                        action='store_true')  # default is False
    parser.add_argument('--no_antisense', help='exclude antisense genes',
                        action='store_true')  # default is False
    parser.add_argument('--notPPI', help='not_run_ppi', action='store_false')
    parser.add_argument('-t', '--prob_threshold', help='threshold selecting gene',
                        required=False, type=float, default=0.5)

    args = parser.parse_args()
    return(args)


class Config():
    # generate run configuration object from user-provided config file in yaml format
    def __init__(self, path):
        f_h = open(path, 'r')
        config = yaml.safe_load(f_h)
        f_h.close()

        self.__dict__ = config
        logger.info('config %s', str(self.__dict__))
        return


args = get_args()


logger.info('args %s', str(args))
phenotype_list = args.phenotype_list
project_dirname = args.project_dirname
k = args.k_closestgenes
FLANKDIST = args.flankdistance
alpha = args.alpha
degree_correction = not args.nodegcor
optim_run = not args.nooptim
seed_start = args.seed_start
seed_end = args.seed_end
TF_run = args.notTF
ppi_run = args.notPPI
unannotated = args.unannotated
pseudogene = args.pseudogene
antisense = not args.no_antisense
network_plot = args.plot
prob_threshold = args.prob_threshold

config = Config(args.config)

logger.info('parameters: phenotypes %s k %d alpha %f',
            str(phenotype_list), k, alpha)

path = config.path
data_dir = config.data_dir
results_dir = config.results_dir
gene_filename = config.gene_filename
ppi_filename = config.ppi_filename
protein_gene_filename = config.protein_gene_filename
trrust_filename = config.trrust_filename
gene_file = os.path.join(data_dir, gene_filename)


protein_gene_file = os.path.join(data_dir, protein_gene_filename)
ppi_file = os.path.join(data_dir, ppi_filename)
trrust_file = os.path.join(data_dir, trrust_filename)

snploc_filename = config.snploc_filename
gwas_filename = config.gwas_filename
feature_initscore_filename = config.feature_initscore_filename

project_datadir = os.path.join(data_dir, project_dirname)

results_dir = os.path.join(results_dir, project_dirname, 'dist_' + str(
    FLANKDIST) + '_' + '_'.join(phenotype_list))


# results_dir = os.path.join(
#     results_dir, project_dirname, 'doNOTdeleteJESUS/DONT/' + 'dist_' + str(FLANKDIST) + '_' + '_'.join(phenotype_list))


if not os.path.exists(results_dir):
    os.makedirs(results_dir)

project_datadir = os.path.join(data_dir, project_dirname)
snploc_file = os.path.join(project_datadir, snploc_filename)
gwas_file = os.path.join(project_datadir, gwas_filename)
feature_initscore_file = os.path.join(
    project_datadir, feature_initscore_filename)

seedless_base_filename = 'ppi' + str(ppi_run) + '_TF' + str(TF_run) + '_degcor' + str(degree_correction) + '_optim' + str(
    optim_run) + '_unannotated' + str(unannotated) + '_pseudogene' + str(pseudogene) + '_antisense' + str(antisense)


def get_gene_snp_distance(missing_genes, gene_bp_dict, snp_bp_dict):
       # calculate distance to all genes on same chromosome
    gene_snp_distance = defaultdict(list)
    gene_closestsnp_distance = defaultdict()

    not_in_ensembl = []
    for gene in missing_genes:
        if gene not in gene_bp_dict:
            not_in_ensembl.append(gene)
            print(gene)
            continue

        gene_chromosome = gene_bp_dict[gene]['chromosome']
        gene_list = []

        for variant in snp_bp_dict.keys():
            if not snp_bp_dict[variant]['chromosome'] == gene_chromosome:
                continue
            else:
                gene_tss = int(gene_bp_dict[gene]['gene_tss'])
                #gene_end = int(gene_bp_dict[gene]['gene_end'])

                variant_bp = int(snp_bp_dict[variant]['rsid_start'])
                distance = abs(gene_tss - variant_bp)

            gene_list.append((variant, distance))

        # find snp with min_distance
        gene_snp_distance[gene] = sorted(gene_list, key=lambda x: abs(x[1]))
        gene_closestsnp_distance[gene] = min(
            gene_list, key=lambda x: abs(x[1]))

    return (gene_snp_distance, gene_closestsnp_distance)


# protein_gene_file = os.path.join(data_dir, protein_gene_filename)
# ppi_file = os.path.join(data_dir, ppi_filename)
# trrust_file = os.path.join(data_dir, trrust_filename)

# snploc_filename = config.snploc_filename
# gwas_filename = config.gwas_filename
# feature_initscore_filename = config.feature_initscore_filename

# results_dir = os.path.join(results_dir, project_dirname,
#                            'dist_' + str(FLANKDIST) + '_' + '_'.join(phenotype_list))
# if not os.path.exists(results_dir):
#     os.makedirs(results_dir)

# project_datadir = os.path.join(data_dir, project_dirname)


# snploc_file = os.path.join(project_datadir, snploc_filename)
# gwas_file = os.path.join(project_datadir, gwas_filename)
# feature_initscore_file = os.path.join(
#     project_datadir, feature_initscore_filename)



####################################################################################################################
####################################################################################################################
# Preparing files to run eMAGMA
####################################################################################################################
####################################################################################################################
# combine GWS SNP information with ensembl GRCh37 information, in the form eMAGMA requires
if 0:
    ret = utils_data.read_gwas_file(
        project_datadir, gwas_filename, phenotype_list, delimiter='\t')

    for snp in ret:
        try:  # some of the GWAS variants do not match with the ensembl database
            ret[snp]['CHR'] = snp_bp_dict[snp]['chromosome']
            ret[snp]['BP'] = snp_bp_dict[snp]['rsid_start']
        except:
            continue

    emagma_df = pd.DataFrame(columns=['SNP', 'CHR', 'BP', 'P', 'Neff'])
    for snp in ret:
        emagma_df = emagma_df.append({'SNP': snp, 'CHR': ret[snp]['CHR'],
                                      'BP': ret[snp]['BP'], 'P': ret[snp]['P'], 'Neff': ret[snp]['Neff']}, ignore_index=True)

    emagma_df_filname = os.path.join(
        results_dir, project_dirname + '_emagma_GRCh37.txt')

    emagma_df.to_csv(emagma_df_filname, sep='\t', index=False)


# eMAGMA commands:
#./magma --bfile g1000_eur --pval BreastCancer_emagma_GRCh37.txt ncol=Neff --gene-settings adap-permp=10000 --out emagma_BreastCancer --gene-annot emagma_annot_2/Breast_Mammary_Tissue.genes.annot awk '{if ($9<=3.8431e-5) print $0}' emagma_BreastCancer.genes.out > emagma_BreastCancer_signif_genes.txt


#./magma --bfile g1000_eur --pval TopMed_emagma_GRCh37.txt ncol=Neff --gene-settings adap-permp=10000 --out emagma_HeartLeftVentricle --gene-annot Heart/Heart_Left_Ventricle.genes.annot
# awk '{if ($9<=3.8431e-5) print $0}' emagma_HeartLeftVentricle.genes.out > emagma_HeartLeftVentricle_signif_genes.txt


#./magma --bfile g1000_eur --pval TopMed_emagma_GRCh37.txt ncol=Neff --gene-settings adap-permp=10000 --out emagma_HeartAtrial --gene-annot Heart/Heart_Atrial_Appendage.genes.annot
# awk '{if ($9<=3.8431e-5) print $0}' emagma_HeartAtrial.genes.out > emagma_HeartAtrial_signif_genes.txt

if 0:
    emagma_genefiles = ['/Users/z/Desktop/eMAGMA/emagma_HeartAtrial_signif_genes.txt',
                        '/Users/z/Desktop/eMAGMA/emagma_HeartLeftVentricle_signif_genes.txt']

    emagma_genefiles = [
        '/Users/z/Desktop/eMAGMA/emagma_BreastCancer_signif_genes.txt']
    utils_data.get_emagmagenes(
        emagma_genefiles, results_dir, project_dirname, delimiter=' ')
# once you get the gene ID then go to ensembl to convert to "names" save as results_dir, project_dirname + '_eMAGMA_genenames.txt'
# change gene name title to 'gene_name' also add a column with 'eMAGMA' written in it next to each gene


####################################################################################################################
####################################################################################################################


seedless_base_filename = 'ppi' + str(ppi_run) + '_TF' + str(TF_run) + '_degcor' + str(degree_correction) + '_optim' + str(
    optim_run) + '_unannotated' + str(unannotated) + '_pseudogene' + str(pseudogene) + '_antisense' + str(antisense)

snp_bp_dict = defaultdict(dict)
for phenotype in phenotype_list:
    SNPfile = os.path.join(project_datadir, phenotype + snploc_filename)
    ph_snp_bp = utils_data.read_snp_file(SNPfile)
    snp_bp_dict.update(ph_snp_bp)

# external genes:
feature_geneset = utils_data.get_feature_geneset(
    project_datadir, phenotype_list)
external_genes = set().union(*list(feature_geneset.values()))

gene_bp_dict = utils_data.read_gene_file(
    gene_file, unannotated, pseudogene, antisense, project_datadir, external_genes)


def write_genestats_network_vs_eMAGMA(results_dir):

    genestats_filename = os.path.join(
        results_dir, seedless_base_filename + '_finalprob_genestats.txt')

    genestats_df = pd.read_csv(
        genestats_filename, sep='\t', dtype=str)

    # new_df = allgenes_df.merge(selected_genes_df, how='outer').convert_dtypes()
    # new_df = selected_genes_df
    # new_df_filname = os.path.join(
    #     results_dir, 'gene_snp_locus_probselected.txt')
    # new_df.to_csv(new_df_filname, sep='\t', index=False)

    new_df = genestats_df

    eMAGMA_genenames_filname = os.path.join(
        project_datadir, project_dirname + '_eMAGMA_genenames.txt')
    eMAGMA_gene_df = pd.read_csv(eMAGMA_genenames_filname, sep='\t', dtype=str)
    eMAGMA_genes = set(eMAGMA_gene_df['gene_name'])

    gene_bp_dict = utils_data.read_gene_file(
        gene_file, unannotated, pseudogene, antisense, project_datadir, external_genes)

    not_in_ensembl = []
    for gene in eMAGMA_genes:
        if gene not in gene_bp_dict:
            not_in_ensembl.append(gene)
        # manual SEARCHING -
        # For TopMed:
        #'MIX23' is equivalent to 'CCDC58' in ensembl
        #'DYNLT5' is equivalent to VSTM2L' in ensembl
        # For BreastCancer:
        #['RUSC1-AS1', 'PANO1'] first one was replaced by RUSC1 and

    eMAGMA_gene_df = eMAGMA_gene_df[~eMAGMA_gene_df['gene_name'].isin(
        not_in_ensembl)]
    if project_dirname == 'TopMed':
        eMAGMA_gene_df = eMAGMA_gene_df.append(
            {'gene_name': 'CCDC58'}, ignore_index=True)
        eMAGMA_gene_df = eMAGMA_gene_df.append(
            {'gene_name': 'VSTM2L'}, ignore_index=True)

    if project_dirname == 'BreastCancer':
        eMAGMA_gene_df = eMAGMA_gene_df.append(
            {'gene_name': 'RUSC1'}, ignore_index=True)

    ###

    final_df = new_df.merge(
        eMAGMA_gene_df, how='outer').convert_dtypes()

    # now add distance to closest snp
    eMAMGMA_only_genes = set(
        eMAGMA_gene_df['gene_name']) - set(new_df['gene_name'])

    snp_to_genes_filename = os.path.join(results_dir, 'snp_to_genes.txt')
    snp_df = pd.read_csv(snp_to_genes_filename, sep='\t')

    eMAGMAgene_snp_distance, eMAGMAgene_closestsnp_distance = get_gene_snp_distance(
        eMAMGMA_only_genes, gene_bp_dict, snp_bp_dict)

    for gene in eMAMGMA_only_genes:

        gene_tss = int(gene_bp_dict[gene]['gene_tss'])
        gene_end = int(gene_bp_dict[gene]['gene_end'])

        closest_snp = eMAGMAgene_closestsnp_distance[gene][0]
        distance_closest_snp_str = str(
            float(eMAGMAgene_closestsnp_distance[gene][1]) / 1000.0)

        # also fill in locus columns (locus is closest_snp)
        final_df.loc[final_df['gene_name']
                     == gene, 'locus_name'] = closest_snp
        final_df.loc[final_df['gene_name']
                     == gene, 'gene_locusdist_kbp'] = distance_closest_snp_str

        final_df.loc[final_df['gene_name']
                     == gene, 'locus_chr'] = str(snp_df.loc[snp_df['snp_rsid'] == closest_snp, 'snp_chr'].values[0])

        final_df.loc[final_df['gene_name'] == gene,
                     'locus_start'] = str(snp_df.loc[snp_df['snp_rsid'] == closest_snp, 'snp_loc'].values[0])
        final_df.loc[final_df['gene_name'] == gene,
                     'locus_end'] = str(snp_df.loc[snp_df['snp_rsid'] == closest_snp, 'snp_loc'].values[0])

        final_df.loc[final_df['gene_name'] == gene,
                     'gene_closestsnp'] = closest_snp

        final_df.loc[final_df['gene_name'] == gene, 'gene_tss'] = str(gene_tss)
        final_df.loc[final_df['gene_name'] == gene, 'gene_end'] = str(gene_end)

    final_df['gene_eMAGMA'] = final_df['gene_eMAGMA'].fillna('.')

    final_df['gene_interesting'] = '.'
    final_df.loc[(final_df['gene_nonzero'] == 'nonzero') | (final_df['gene_mindist'] == 'mindist') | (
        final_df['gene_special'] != '.') | (final_df['gene_eMAGMA'] == 'eMAGMA'), 'gene_interesting'] = 'interesting'

    final_df = final_df[['locus_chr', 'locus_start', 'locus_end', 'locus_name', 'locus_ngene', 'locus_width', 'locus_pheno', 'gene_name', 'gene_prob_selected', 'gene_eMAGMA', 'gene_interesting', 'gene_special', 'gene_pattern',  'gene_nonzero', 'gene_mindist',  'gene_tss', 'gene_end', 'gene_locusdist_kbp', 'gene_closestsnp', 'gene_type',
                         'ppi_deg', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names']]

    final_df.drop_duplicates(inplace=True)

    final_df = final_df.astype({'locus_chr': 'int', 'locus_start': 'int',
                                'locus_end': 'int', 'gene_tss': 'int', 'gene_end': 'int'})
    final_df = final_df.sort_values(
        by=['locus_chr', 'locus_start', 'locus_end', 'gene_tss', 'gene_end'])

    final_eMAGMA_filname = os.path.join(
        results_dir, 'genestats_network_vs_eMAGMA.txt')
    final_df.to_csv(final_eMAGMA_filname, sep='\t', index=False)


def write_summary_genestats_network_vs_eMAGMA(prob_threshold=0.5):
    compare_eMAGMA_filname = os.path.join(
        results_dir, 'genestats_network_vs_eMAGMA.txt')
    compare_eMAGMA_df = pd.read_csv(
        compare_eMAGMA_filname, sep='\t', dtype=str)

    summary_compare_eMAGMA_df = pd.DataFrame()

    for idx in range(len(compare_eMAGMA_df)):
        row = compare_eMAGMA_df.iloc[[idx]]

        if (row['gene_name'].values[0] == row['gene_mindist'].values[0]) or (float(row['gene_prob_selected']) > prob_threshold) or (row['gene_eMAGMA'].values[0] == 'eMAGMA') or not (row['gene_special'].isnull().values[0]):
            summary_compare_eMAGMA_df = pd.concat(
                [summary_compare_eMAGMA_df, row], axis=0, ignore_index=True)

    eMAMGMA_only_genes = list(
        summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['locus_name'].isnull()]['gene_name'].values)

    # NOW CHECK FOR PPI, TF, SPECIAL:

    uniprot_to_hgnc = utils_ppi.protein_gene(protein_gene_file)
    uniprot_directed_dict = utils_ppi.read_ppi(ppi_file)
    protein_directed_dict = utils_ppi.convert_uniprot_to_hgnc(
        uniprot_directed_dict, uniprot_to_hgnc)
    protein_proteins = utils_ppi.make_undirected(protein_directed_dict)

    UniverseTFsource_target, UniverseTFsource_target_sign = utils_ppi.Trrust_dict(
        trrust_file)

    eMAGMA_ppi = [x for x in eMAMGMA_only_genes if x in protein_proteins]

    eMAGMA_TF = [
        x for x in eMAMGMA_only_genes if x in UniverseTFsource_target]

    eMAGMA_feature = defaultdict(set)

    for feature in ['omim', 'exome', 'cosmic']:
        if os.path.exists(os.path.join(project_datadir, feature + '_genes.txt')):
            eMAGMA_feature[feature] = set(eMAMGMA_only_genes).intersection(
                utils_data.get_geneset(project_datadir, feature + '_genes.txt'))

    coloc_set = set()
    for phenotype in phenotype_list:
        coloc_filename = phenotype + '_coloc_genes.txt'
        if os.path.exists(os.path.join(project_datadir, coloc_filename)):
            coloc_set = coloc_set.union(
                utils_data.get_geneset(project_datadir, coloc_filename))

    eMAGMA_feature['coloc'] = set(eMAMGMA_only_genes).intersection(coloc_set)

    for gene in eMAMGMA_only_genes:

        if gene in eMAGMA_TF:
            summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['gene_name']
                                          == gene, 'TF_indeg'] = 'Yes'
            summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['gene_name']
                                          == gene, 'TF_outdeg'] = 'Yes'

        if gene in eMAGMA_ppi:
            summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['gene_name']
                                          == gene, 'ppi_deg'] = 'Yes'

        special_str = ''
        for feature in eMAGMA_feature:
            if gene in eMAGMA_feature[feature]:
                '+'.join([special_str, feature])
        summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['gene_name']
                                      == gene, 'gene_special'] = special_str

    summary_compare_eMAGMA_df = summary_compare_eMAGMA_df.astype(
        {'locus_chr': 'int', 'locus_start': 'int', 'locus_end': 'int', 'gene_tss': 'int', 'gene_end': 'int'})
    summary_compare_eMAGMA_df = summary_compare_eMAGMA_df.sort_values(
        by=['locus_chr', 'locus_start', 'locus_end', 'gene_tss', 'gene_end'])

    summary_compare_eMAGMA_filname = os.path.join(
        results_dir, 'summary_genestats_network_vs_eMAGMA.txt')
    summary_compare_eMAGMA_df.to_csv(
        summary_compare_eMAGMA_filname, sep='\t', index=False)


# if 0:
#     emagma_file = os.path.join(
#         results_dir, seedless_base_filename + '_emagma_vs_Z_pathways.txt')
#     emagma_df = pd.read_csv(emagma_file, sep='\t')

#     emagma_df['EmagmaGenes']

#     algo_genes = set((new_df.loc[new_df['gene_prob_selected'] >= 0.5]['gene_name']))
#     len(algo_genes.intersection(eMAGMA_genes))
#     len(set(summary_compare_eMAGMA_df.loc[summary_compare_eMAGMA_df['gene_name'].isin(
#         eMAGMA_genes)]['locus_name']))


# if 0:
#     network_file = os.path.join(
#         results_dir, seedless_base_filename + '_counts_updated.txt')

#     selected_genes = set(
#         selected_genes_df[selected_genes_df['gene_prob_selected'] > prob_threshold]['gene_name'])

#     pdb.set_trace()

#     eMAMGMA_only_genes = set(eMAGMA_genes) - set(selected_genes)
#     selected_only_genes = set(selected_genes) - set(eMAGMA_genes)

#     both = set(selected_genes).intersection(set(eMAGMA_genes))
