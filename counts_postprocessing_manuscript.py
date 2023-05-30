import argparse
import yaml
import logging
import pandas as pd
from collections import defaultdict
import numpy as np
from scipy.stats import sem
import os as os
import glob as glob
import csv
import pdb
import pickle
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import utils_ppi
import utils_data
import utils_plot
import random


logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('gene_selection_simple')
logger.setLevel(logging.INFO)


def get_args():
    # ,epilog='small call: see run.sh', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(description='Select genes')
    parser.add_argument('-c', '--config', help='config file',
                        required=False, type=str, default='./signet_config.yaml')
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


protein_gene_file = os.path.join(path, data_dir, protein_gene_filename)
ppi_file = os.path.join(path, data_dir, ppi_filename)
trrust_file = os.path.join(path, data_dir, trrust_filename)

# project specific files
snploc_filename = config.snploc_filename
gwas_filename = config.gwas_filename
feature_initscore_filename = config.feature_initscore_filename


project_dirname = project_dirname
phenostr = '+'.join(phenotype_list)
flankstr = str(int(int(FLANKDIST) / 1000)) + 'kb'

results_dir = os.path.join(results_dir, '_'.join(
    [project_dirname, phenostr, flankstr]))

project_datadir = os.path.join(data_dir, project_dirname)
snploc_file = os.path.join(project_datadir, snploc_filename)
gwas_file = os.path.join(project_datadir, gwas_filename)
feature_initscore_file = os.path.join(
    project_datadir, feature_initscore_filename)


seedless_base_filename = 'SIGNET'


def write_closestgene_to_snp(results_dir):
    snp_to_genes_filename = os.path.join(results_dir, 'snp_to_genes.txt')
    df = pd.read_csv(snp_to_genes_filename, sep='\t', dtype=str)

    gene_mindistsnp = defaultdict(set)
    for idx, row in df.iterrows():
        gene = row['snp_mindistgene']

        if gene in gene_mindistsnp:
            gene_mindistsnp[gene].add(row['snp_name'])
        else:
            gene_mindistsnp[gene] = {row['snp_name']}

    for key, value in gene_mindistsnp.items():
        gene_mindistsnp[key] = '+'.join(sorted(value))

    gene_mindistsnp_df = pd.DataFrame([gene_mindistsnp]).T
    gene_mindistsnp_df.reset_index(inplace=True)
    gene_mindistsnp_df = gene_mindistsnp_df.rename(
        columns={'index': 'gene_name', 0: 'snp_mindistgene'})

    gene_mindistsnp_filename = os.path.join(
        results_dir, 'gene_snpmindist_gene.txt')
    gene_mindistsnp_df.to_csv(
        gene_mindistsnp_filename, index=False, sep='\t')


def write_finalnetworkscores(results_dir, seed_start=0, seed_end=100):
    df = pd.DataFrame(columns=['seed', 'finalnetworkscore'])

    for seed in range(seed_start, seed_end):
        base_filename = f'seed{seed:02d}'
        file = os.path.join(results_dir, base_filename +
                            '_networkscore_log.txt')
        seed_df = pd.read_csv(file, delimiter='\t')
        seed = seed_df['seed'].values[-1]
        # [-1] to only consider final value
        seed_finalnetworkscore = seed_df['network_score'].values[-1]
        df.loc[len(df.index)] = [seed, seed_finalnetworkscore]

    max_scores = np.max(df['finalnetworkscore'])
    min_scsores = np.min(df['finalnetworkscore'])
    mean_scsores = np.mean(df['finalnetworkscore'])
    std_scsores = np.std(df['finalnetworkscore'])

    finalscores_file = os.path.join(results_dir, 'final_networkscores.txt')
    df.to_csv(finalscores_file, index=False, sep='\t')

    maxnetworkscore_seed = int(df.loc[df['finalnetworkscore']
                                      == max_scores, 'seed'].values[0])
    return maxnetworkscore_seed


def write_aggregated_statsfiles(seedless_base_filename, results_dir, seed_start=0, seed_end=100):

    maxnetworkscore_seed = write_finalnetworkscores(
        results_dir, seed_start=seed_start, seed_end=seed_end)

    print(' seed %s has maxnetwork score' % maxnetworkscore_seed)
    pdb.set_trace()

    base_filename = f'seed{maxnetworkscore_seed:02d}'

    genestats_file = os.path.join(
        results_dir, base_filename + '_extended_finalpass_scores.txt')

    genestats_df = pd.read_csv(genestats_file, sep='\t', dtype=str)

    genestats_df = genestats_df[['locus_name', 'gene_name', 'gene_pattern', 'gene_mindist',
                                 'gene_locusdist_kbp', 'distance_score', 'gene_special', 'gene_nonzero', 'special_score', 'ppi_score', 'exp_score', 'ppi_deg', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'TF_score', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names']]

    counts_filename = os.path.join(
        results_dir, seedless_base_filename + '_counts_allruns_prob.txt')

    counts_df = pd.read_csv(counts_filename, sep='\t', dtype=str)
    counts_df = counts_df[['locus_name', 'gene_name', 'gene_prob_selected']]
    final_df = pd.merge(counts_df, genestats_df, on=[
                        'locus_name', 'gene_name'], validate='one_to_one')

    gene_closest_filename = os.path.join(results_dir, 'gene_closest.txt')
    gene_df = pd.read_csv(gene_closest_filename, sep='\t', dtype=str)
    gene_df = gene_df[['gene_name', 'gene_closestsnp']]
    final_df = final_df.merge(gene_df, how='left', on='gene_name')

    locus_to_genes_filename = os.path.join(results_dir, 'locus_to_genes.txt')
    locus_df = pd.read_csv(locus_to_genes_filename, sep='\t', dtype=str)

    final_df = final_df.merge(locus_df, how='left', on=[
                              'locus_name', 'gene_name', 'gene_locusdist_kbp'])

    # add a column for locus_mindist_gene

    final_df['gene_signet'] = '.'
    final_df['gene_signet+'] = '.'
    for locus in final_df['locus_name'].unique():

        locus_df = final_df[final_df['locus_name'] == locus].copy()

        mindist_genestr = locus_df.loc[locus_df['gene_mindist']
                                       == 'mindist']['gene_name'].values[0]

        locus_nspecialgene = sum(locus_df['gene_special'] != '.')

        locus_ninteractiongene = sum(locus_df['ppi_deg'].astype(
            int) + locus_df['TF_indeg'].astype(int) + locus_df['TF_outdeg'].astype(int) != 0)

        max_prob = np.max(locus_df['gene_prob_selected'])
        highest_genes = locus_df.loc[locus_df['gene_prob_selected']
                                     == max_prob]['gene_name'].values

        selected_gene = random.choice(sorted(highest_genes))

        final_df.loc[final_df['locus_name'] == locus,
                     'locus_mindistgene'] = mindist_genestr
        final_df.loc[final_df['locus_name'] == locus,
                     'locus_nspecialgene'] = locus_nspecialgene
        final_df.loc[final_df['locus_name'] == locus,
                     'locus_ninteractiongene'] = locus_ninteractiongene

        final_df.loc[(final_df['locus_name'] == locus) & (
            final_df['gene_name'] == selected_gene), 'gene_signet'] = 'signet'
        final_df.loc[(final_df['locus_name'] == locus) & ((final_df['gene_name'] == selected_gene) | (
            final_df['gene_special'] != '.')), 'gene_signet+'] = 'signet+'
        #(row['gene_prob_selected'] == max_prob) or (row['gene_special'] != '.')

    final_df['gene_interesting'] = '.'
    final_df.loc[(final_df['gene_nonzero'] == 'nonzero') | (final_df['gene_mindist'] == 'mindist') | (
        final_df['gene_special'] != '.'), 'gene_interesting'] = 'interesting'

    final_df = final_df.astype({'locus_chr': 'int', 'locus_start': 'int',
                                'locus_end': 'int', 'gene_tss': 'int', 'gene_end': 'int'})
    final_df = final_df.sort_values(
        by=['locus_chr', 'locus_start', 'locus_end', 'gene_tss', 'gene_end'])

    summary_df = final_df

    summary_df['locus_hisig_type'] = '.'
    summary_df['gene_hisig_type'] = '.'
    locus_df = pd.DataFrame()
    for locus in summary_df['locus_name'].unique():
        if locus.startswith('loc'):
            continue
        locus_df = summary_df[summary_df['locus_name'] == locus].copy()
        locus_special_str = "".join(locus_df['gene_special'])
        if 'omim' in locus_special_str:
            summary_df.loc[summary_df['locus_name']
                           == locus, 'locus_hisig_type'] = 'omim'
        elif 'exome' in locus_special_str:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'exome'
        elif 'coloc' in locus_special_str:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'coloc'
        else:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'poorsig'

        locus_geneset = set(locus_df['gene_name'].values)
        locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                      != '.']['gene_name'].values[0]

        for gene in locus_geneset:
            gene_special_str = "".join(
                locus_df.loc[locus_df['gene_name'] == gene]['gene_special'])

            if 'omim' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'omim'
            elif 'exome' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'exome'
            elif 'coloc' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'coloc'
            elif gene == locus_mindist_gene:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'mindist'
            else:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'poorsig'

    summary_df = summary_df[['locus_chr', 'locus_start', 'locus_end', 'locus_name', 'locus_ngene', 'locus_width', 'locus_pheno', 'locus_mindistgene', 'locus_nspecialgene', 'locus_ninteractiongene', 'locus_hisig_type', 'gene_tss', 'gene_end', 'gene_locusdist_kbp', 'gene_closestsnp', 'gene_name', 'gene_type', 'gene_prob_selected',
                             'gene_signet', 'gene_mindist', 'gene_special',  'gene_signet+', 'distance_score', 'special_score', 'ppi_score', 'TF_score', 'exp_score', 'gene_hisig_type', 'ppi_deg', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names']]

    summary_df_filename = os.path.join(
        results_dir, 'summary_genestats_manuscript.txt')
    summary_df.to_csv(summary_df_filename, index=False, sep='\t')


def write_finalstats_selected(results_dir, prob_threshold):

    filename = os.path.join(results_dir,  'summary_genestats_manuscript.txt')
    final_prob_genestats_df = pd.read_csv(filename, sep='\t')
    selected_df = final_prob_genestats_df.loc[final_prob_genestats_df['gene_prob_selected'] >= prob_threshold]
    selected_df = selected_df[~selected_df['locus_name'].str.startswith('loc')]

    selected_dist = selected_df['gene_locusdist_kbp'].values
    min_selecteddist = min(selected_dist)
    max_selecteddist = max(selected_dist)
    median_selecteddist = np.median(selected_dist)

    locus_width = selected_df['locus_width'].values
    min_locus_width = min(locus_width)
    max_locus_width = max(locus_width)
    median_locus_width = np.median(locus_width)

    loci_stats = {'min_selecteddist': min_selecteddist,
                  'max_selecteddist': max_selecteddist, 'median_selecteddist': median_selecteddist, 'min_locus_width': min_locus_width, 'max_locus_width': max_locus_width, 'median_locus_width': median_locus_width}

    df = pd.DataFrame([loci_stats])
    filename = os.path.join(results_dir, 'selected_dist_stats.csv')

    df.to_csv(filename, index=False, header=True)


def kegg_pathway_comparison(file1path, file1name, file2path, file2name, results_dir, seedless_base_filename):

    gene_pathwaysdf_file1 = pd.read_csv(file1path, sep='\t')
    gene_pathwaysdf_file2 = pd.read_csv(file2path, sep='\t')

    pathways_df = gene_pathwaysdf_file1.merge(
        right=gene_pathwaysdf_file2, on="Term", validate='one_to_one', suffixes=('_' + file1name, '_' + file2name))

    pathways_df[('-log10(Adjsuted_pval)_' + file1name)] = - \
        np.log10(pathways_df[('Adjusted P-value_' + file1name)])
    pathways_df[('-log10(Adjsuted_pval)_' + file2name)] = - \
        np.log10(pathways_df[('Adjusted P-value_' + file2name)])

    pathways_df['pval_diff'] = pathways_df['-log10(Adjsuted_pval)_' + file1name] - \
        pathways_df['-log10(Adjsuted_pval)_' + file2name]

    pathways_df[file1name + '_genes'] = [sorted(x.split(';'))
                                         for x in pathways_df['Genes_' + file1name].values]
    pathways_df[file2name + '_genes'] = [sorted(x.split(';'))
                                         for x in pathways_df['Genes_' + file2name].values]

    pathways_df['CommonGenes'] = sorted([set(pathways_df[file1name + '_genes'][x]).intersection(
        set(pathways_df[file2name + '_genes'][x])) for x in range(len(pathways_df))])

    pathways_df[file1name + '_only'] = [sorted(set(pathways_df[file1name + '_genes'][x]) - set(
        pathways_df[file2name + '_genes'][x])) for x in range(len(pathways_df))]

    pathways_df[file2name + '_only'] = [sorted(set(pathways_df[file2name + '_genes'][x]) -
                                               set(pathways_df[file1name + '_genes'][x])) for x in range(len(pathways_df))]

    file = os.path.join(results_dir, (seedless_base_filename + '_' +
                                      file1name + '_vs_' + file2name + '_pathways.txt'))
    pathways_df.to_csv(file, index=False, sep='\t')


def write_feature_weight_dict(results_dir):
    feature_weights = defaultdict(list)
    # if project_dirname == 'TopMed':
    feature_list = ['distfac_active',
                    'distfac_inactive', 'omim', 'exome', 'coloc']

    feature_finalvalues_df = pd.DataFrame()

    for feature in feature_list:
        feature_weights[feature] = []
        for file in glob.glob(results_dir + '/seed*' + feature + '_log.txt'):
            df = pd.read_csv(file, delimiter='\t')
            # ADDED the [-1] to only consider final value
            feature_weights[feature].append(df[feature].values[-1])

        feature_finalvalues_df[feature] = feature_weights[feature]

    file = os.path.join(
        results_dir, 'features_weights.txt')
    feature_finalvalues_df.to_csv(file, index=False, sep='\t')

    for feature in feature_list:
        print('feature: ', feature, 'mean: ', np.mean(
            feature_weights[feature]), 'std: ', np.std(feature_weights[feature]))

    feature_summary_df = pd.DataFrame(columns=['feature', 'mean', 'std'])

    for feature in feature_list:
        df = {'feature': feature, 'mean': np.mean(
            feature_weights[feature]), 'std': np.std(feature_weights[feature])}
        feature_summary_df = feature_summary_df.append(df, ignore_index=True)

    file = os.path.join(results_dir, 'features_weights_summary.txt')
    feature_summary_df.to_csv(file, index=False, sep='\t')


if 1:
    write_finalnetworkscores(results_dir, seed_end=1)
    write_aggregated_statsfiles(
        seedless_base_filename, results_dir, seed_end=1)


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)

    summary_filename = 'summary_genestats_manuscript.txt'
    summary_file = os.path.join(results_dir, summary_filename)

    summary_df = pd.read_csv(summary_file, sep='\t')
    selected_geneset = set(
        summary_df[summary_df['gene_signet'] == 'signet']['gene_name'].values)
    for anchor_gene in selected_geneset:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_filename, network_ppi, gene_special, TFsource_target)

if 0:
    write_feature_weight_dict(results_dir)


if 0:
    file1path = os.path.join(results_dir, 'KEGG_2021_signet.txt')
    file2path = os.path.join(results_dir,  'KEGG_2021_mindist.txt')
    file3path = os.path.join(results_dir,  'KEGG_2021_signetplus.txt')
    file1name = 'selected'
    file2name = 'mindist'
    file3name = 'signetplus'
    kegg_pathway_comparison(file1path, file1name, file2path,
                            file2name, results_dir, seedless_base_filename)
    kegg_pathway_comparison(file1path, file1name, file3path,
                            file3name, results_dir, seedless_base_filename)
