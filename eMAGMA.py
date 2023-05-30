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
import glob

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
project_dirname = args.project_dirname

config = Config(args.config)
path = config.path
snploc_filename = config.snploc_filename


####################################################################################################################
####################################################################################################################
# Preparing files to run eMAGMA
####################################################################################################################
####################################################################################################################
# combine GWS SNP information with ensembl GRCh37 information, in the form eMAGMA requires


project_dirname = 'AUTISM'

if project_dirname == 'AUTISM':
    phenotype_list = ['ASD']
    general_datadir = os.path.join(path, 'GeneNetwork/AUTISM/data_general')
    project_datadir = os.path.join(path, 'GeneNetwork/AUTISM/data_ASD')
    results_dir = os.path.join(path, 'GeneNetwork/AUTISM/results')
    gwas_filename = '_GWAS_rows.txt'  # 'gwascatalog_asd_signif.txt'

if 0:
    ret = utils_data.read_gwas_file(
        project_datadir, gwas_filename, phenotype_list, delimiter='\t')

    snp_bp_dict = defaultdict(dict)
    for phenotype in phenotype_list:
        SNPfile = os.path.join(
            project_datadir, phenotype + snploc_filename)
        ph_snp_bp = utils_data.read_snp_file(SNPfile)
        snp_bp_dict.update(ph_snp_bp)

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
        results_dir, project_dirname + '_emagma_GRCh38.txt')

    emagma_df.to_csv(emagma_df_filname, sep='\t', index=False)

#########
# run emagma shell script. Create emagma_TISSUE_signif.genes.txt
#########
emagma_path = os.path.join(path, 'magma_v1.10_mac')
file_list = glob.glob(emagma_path + '/*_signif_genes.txt')
utils_data.get_emagmagenes(file_list, results_dir,
                           project_dirname, delimiter=' ')

# once you get the gene ID (file: '_emagma_NCBIgenes.txt') then go to ensembl to convert to "names" save as results_dir, project_dirname + 'emagma_genes.txt'
# change gene name column header to 'gene_name' also add a column with 'eMAGMA' written in it next to each gene
