import pandas as pd
import numpy as np
import csv
import os
import pdb
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
# import counts_postprocessing
import utils_plot
import pickle
from scipy.stats import fisher_exact

results_dir = '/Users/z/Dropbox/Bader_Mousavi/SigNet_zm/results/TopMed_QT+HR+QRS+PR+JT_250kb'
data_dir = '/Users/z/Dropbox/Bader_Mousavi/SigNet_zm/data/TopMed'

summary_file = os.path.join(results_dir, 'summary_genestats_manuscript.txt')
summary_df = pd.read_csv(summary_file, sep='\t')


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)

    summary_file = 'summary_genestats_manuscript.txt'

    for anchor_gene in ['NKX2-5', 'PMP22']:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_file, network_ppi, gene_special, TFsource_target)


# ##########
# # New Table 4
# ##########
# table = np.array([[25, 4], [4, 18]])
# res = fisher_exact(table, alternative='two-sided')
# print('omim_exome:', res[1])

# table = np.array([[25, 1], [1, 36]])
# res = fisher_exact(table, alternative='two-sided')
# print('omim_coloc:', res[1])

# table = np.array([[18, 0], [0, 36]])
# res = fisher_exact(table, alternative='two-sided')
# print('exome_coloc:', res[1])


# pdb.set_trace()


##########
# TABLE 2
##########
if 0:
    omim_filename = 'omim_genes.txt'
    omim_file = os.path.join(data_dir, omim_filename)
    omim_df = pd.read_csv(omim_file, sep='\t')
    print(omim_df.groupby(['Syndrome']).count())
    allomim_geneset = set(omim_df['gene_name'])
    # remove "singleton" omim loci:
    singletonomim_geneset = set(
        summary_df[summary_df['locus_name'].str.startswith('loc')]['gene_name'])
    gwasomim_geneset = allomim_geneset - singletonomim_geneset


##########
# TABLE 1
##########
if 0:
    phenotype_list = ['QT', 'HR', 'QRS', 'PR', 'JT']
    snpset = set()
    for phenotype in phenotype_list:
        pheno_rsid_filename = phenotype + '_rsid_GRCh38p13.txt'
        pheno_rsid_file = os.path.join(data_dir, pheno_rsid_filename)
        GWAS_df = pd.read_csv(pheno_rsid_file, sep='\t').drop_duplicates(
            subset=['Variant name'])
        pheno_snpset = set(GWAS_df['Variant name'])
        snpset = snpset.union(pheno_snpset)
        print('number of unique %s variants: %d' % (phenotype, len(GWAS_df)))
    print('number of unique variants: %d' % (len(snpset)))

    print('coloc gene breakdown')
    coloc_geneset = set()
    for phenotype in phenotype_list:
        pheno_coloc_filename = phenotype + '_coloc_genes.txt'
        pheno_coloc_file = os.path.join(data_dir, pheno_coloc_filename)
        if os.path.exists(pheno_coloc_file):
            pheno_colocgeneset = set(pd.read_csv(
                pheno_coloc_file, sep='\t')['gene_name'])
            coloc_geneset = coloc_geneset.union(pheno_colocgeneset)

            print('number of %s coloc genes: %d' %
                  (phenotype, len(pheno_colocgeneset)))

    print('number of unique coloc genes: %d' % (len(coloc_geneset)))

    # Exome genes
    exome_filename = 'exome_genes.txt'
    exome_file = os.path.join(data_dir, exome_filename)
    exome_df = pd.read_csv(exome_file, sep='\t')
    exome_geneset = set(exome_df['gene_name'])


##########
# TABLE 3
##########
# Look at loci_aggregatedstats.txt from :
# gene_selection.write_loci_aggregatedstats(results_dir, locus_stats)
# look at selected_dist_stats.csv from:
# counts_postprocessing.write_finalstats_selected(results_dir, prob_threshold=0.5)


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)

    summary_filename = 'summary_genestats_manuscript.txt'

    for anchor_gene in all_selected_geneset:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_filename, network_ppi, gene_special, TFsource_target)


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)

    summary_file = 'summary_genestats_manuscript.txt'

    for anchor_gene in selected_geneset:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_file, network_ppi, gene_special, TFsource_target)

##########
# New Table 6
##########
omim_loci = set()
exome_loci = set()
coloc_loci = set()
poorsig_loci = set()

omim_loci_omim_selected = set()
omim_loci_exome_selected = set()
omim_loci_coloc_selected = set()
omim_loci_poorsig_selected = set()
omim_loci_mindist_selected = set()

exome_loci_omim_selected = set()
exome_loci_exome_selected = set()
exome_loci_coloc_selected = set()
exome_loci_poorsig_selected = set()
exome_loci_mindist_selected = set()

coloc_loci_omim_selected = set()
coloc_loci_exome_selected = set()
coloc_loci_coloc_selected = set()
coloc_loci_poorsig_selected = set()
coloc_loci_mindist_selected = set()

poorsig_loci_omim_selected = set()
poorsig_loci_exome_selected = set()
poorsig_loci_coloc_selected = set()
poorsig_loci_poorsig_selected = set()
poorsig_loci_mindist_selected = set()

locus_df = pd.DataFrame()
for locus in summary_df['locus_name'].unique():
    if locus.startswith('loc'):
        continue

    locus_df = summary_df[summary_df['locus_name'] == locus].copy()
    locus_special_str = "".join(locus_df['gene_special'])
    locus_geneset = set(locus_df['gene_name'].values)

    selected_gene = locus_df[locus_df['gene_signet']
                             == 'signet']['gene_name'].values[0]

    selected_gene_special_str = "".join(
        locus_df.loc[locus_df['gene_name'] == selected_gene]['gene_special'])
    locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                  != '.']['gene_name'].values[0]

    if 'omim' in locus_special_str:
        omim_loci.add(locus)

        # now keep count of the type of selected gene:

        if 'omim' in selected_gene_special_str:
            omim_loci_omim_selected.add(locus)
        elif 'exome' in selected_gene_special_str:  # should be 0
            omim_loci_exome_selected.add(locus)
        elif 'coloc' in selected_gene_special_str:
            omim_loci_coloc_selected.add(locus)
        elif selected_gene == locus_mindist_gene:
            omim_loci_mindist_selected.add(locus)
        else:
            omim_loci_poorsig_selected.add(locus)

        continue

    if 'exome' in locus_special_str:
        exome_loci.add(locus)
        if 'omim' in selected_gene_special_str:
            exome_loci_omim_selected.add(locus)
        elif 'exome' in selected_gene_special_str:
            exome_loci_exome_selected.add(locus)
        elif 'coloc' in selected_gene_special_str:
            exome_loci_coloc_selected.add(locus)
        elif selected_gene == locus_mindist_gene:
            exome_loci_mindist_selected.add(locus)
        else:
            exome_loci_poorsig_selected.add(locus)

        continue
    if 'coloc' in locus_special_str:
        coloc_loci.add(locus)
        if 'omim' in selected_gene_special_str:
            coloc_loci_omim_selected.add(locus)
        elif 'exome' in selected_gene_special_str:
            coloc_loci_exome_selected.add(locus)
        elif 'coloc' in selected_gene_special_str:
            coloc_loci_coloc_selected.add(locus)
        elif selected_gene == locus_mindist_gene:
            coloc_loci_mindist_selected.add(locus)
        else:
            coloc_loci_poorsig_selected.add(locus)

        continue

    else:
        poorsig_loci.add(locus)
        if 'omim' in selected_gene_special_str:
            poorsig_loci_omim_selected.add(locus)
        elif 'exome' in selected_gene_special_str:
            poorsig_loci_exome_selected.add(locus)
        elif 'coloc' in selected_gene_special_str:
            poorsig_loci_coloc_selected.add(locus)
        elif selected_gene == locus_mindist_gene:
            poorsig_loci_mindist_selected.add(locus)
        else:
            poorsig_loci_poorsig_selected.add(locus)

        continue


print('omim_loci:  %d' % len(omim_loci))
print(len(omim_loci_omim_selected))
print(len(omim_loci_exome_selected))
print(len(omim_loci_coloc_selected))
print(len(omim_loci_mindist_selected))
print(len(omim_loci_poorsig_selected))


print('exome_loci: %d' % len(exome_loci))
print(len(exome_loci_omim_selected))
print(len(exome_loci_exome_selected))
print(len(exome_loci_coloc_selected))
print(len(exome_loci_mindist_selected))
print(len(exome_loci_poorsig_selected))


print('coloc_loci: %d' % len(coloc_loci))
print(len(coloc_loci_omim_selected))
print(len(coloc_loci_exome_selected))
print(len(coloc_loci_coloc_selected))
print(len(coloc_loci_mindist_selected))
print(len(coloc_loci_poorsig_selected))


print('poorsig_loci: %d' % len(poorsig_loci))
print(len(poorsig_loci_omim_selected))
print(len(poorsig_loci_exome_selected))
print(len(poorsig_loci_coloc_selected))
print(len(poorsig_loci_mindist_selected))
print(len(poorsig_loci_poorsig_selected))


omim_df = summary_df.loc[summary_df['locus_name'].isin(omim_loci)]
omim_file = os.path.join(results_dir, 'summary_genestats_omimloci.txt')
omim_df.to_csv(omim_file, index=False, sep='\t')

exome_df = summary_df.loc[summary_df['locus_name'].isin(exome_loci)]
exome_file = os.path.join(results_dir, 'summary_genestats_exomeloci.txt')
exome_df.to_csv(exome_file, index=False, sep='\t')

coloc_df = summary_df.loc[summary_df['locus_name'].isin(coloc_loci)]
coloc_file = os.path.join(results_dir, 'summary_genestats_colocloci.txt')
coloc_df.to_csv(coloc_file, index=False, sep='\t')

poorsig_df = summary_df.loc[summary_df['locus_name'].isin(poorsig_loci)]
poorsig_file = os.path.join(results_dir, 'summary_genestats_poorsigloci.txt')
poorsig_df.to_csv(poorsig_file, index=False, sep='\t')


pdb.set_trace()


##########
# New Table 7
##########
hiomim_geneset = set()
hiomim_selected = set()
hiomim_notselected = set()
hiomim_gene_selected_omim = set()
hiomim_gene_selected_exome = set()
hiomim_gene_selected_coloc = set()
hiomim_gene_selected_mindist = set()
hiomim_gene_selected_poorsig = set()

hiexome_geneset = set()
hiexome_selected = set()
hiexome_notselected = set()
hiexome_gene_selected_omim = set()
hiexome_gene_selected_exome = set()
hiexome_gene_selected_coloc = set()
hiexome_gene_selected_mindist = set()
hiexome_gene_selected_poorsig = set()

hicoloc_geneset = set()
hicoloc_selected = set()
hicoloc_notselected = set()
hicoloc_gene_selected_omim = set()
hicoloc_gene_selected_exome = set()
hicoloc_gene_selected_coloc = set()
hicoloc_gene_selected_mindist = set()
hicoloc_gene_selected_poorsig = set()

himindist_geneset = set()
himindist_selected = set()
himindist_notselected = set()
himindist_gene_selected_omim = set()
himindist_gene_selected_exome = set()
himindist_gene_selected_coloc = set()
himindist_gene_selected_mindist = set()
himindist_gene_selected_poorsig = set()

hipoorsig_geneset = set()
hipoorsig_selected = set()
hipoorsig_notselected = set()
hipoorsig_gene_selected_omim = set()
hipoorsig_gene_selected_exome = set()
hipoorsig_gene_selected_coloc = set()
hipoorsig_gene_selected_mindist = set()
hipoorsig_gene_selected_poorsig = set()


locus_df = pd.DataFrame()
for locus in summary_df['locus_name'].unique():
    if locus.startswith('loc'):
        continue
    locus_df = summary_df[summary_df['locus_name'] == locus].copy()

    locus_geneset = set(locus_df['gene_name'].values)
    selected_gene = locus_df[locus_df['gene_signet']
                             == 'signet']['gene_name'].values[0]
    selected_gene_special_str = "".join(
        locus_df.loc[locus_df['gene_name'] == selected_gene]['gene_special'])

    for gene in locus_geneset:
        gene_special_str = "".join(
            locus_df.loc[locus_df['gene_name'] == gene]['gene_special'])
        locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                      != '.']['gene_name'].values[0]

        if 'omim' in gene_special_str:
            hiomim_geneset.add(gene)
            if gene == selected_gene:
                hiomim_selected.add(gene)
            else:
                hiomim_notselected.add(gene)
                if 'omim' in selected_gene_special_str:
                    hiomim_gene_selected_omim.add(gene)
                elif 'exome' in selected_gene_special_str:
                    hiomim_gene_selected_exome.add(gene)
                elif 'coloc' in selected_gene_special_str:
                    hiomim_gene_selected_coloc.add(gene)
                elif selected_gene == locus_mindist_gene:
                    hiomim_gene_selected_mindist.add(gene)
                else:
                    hiomim_gene_selected_poorsig.add(gene)

        elif 'exome' in gene_special_str:
            hiexome_geneset.add(gene)
            if gene == selected_gene:
                hiexome_selected.add(gene)
            else:
                hiexome_notselected.add(gene)
                if 'omim' in selected_gene_special_str:
                    hiexome_gene_selected_omim.add(gene)
                elif 'exome' in selected_gene_special_str:
                    hiexome_gene_selected_exome.add(gene)
                elif 'coloc' in selected_gene_special_str:
                    hiexome_gene_selected_coloc.add(gene)
                elif selected_gene == locus_mindist_gene:
                    hiexome_gene_selected_mindist.add(gene)
                else:
                    hiexome_gene_selected_poorsig.add(gene)

        elif 'coloc' in gene_special_str:
            hicoloc_geneset.add(gene)
            if gene == selected_gene:
                hicoloc_selected.add(gene)
            else:
                hicoloc_notselected.add(gene)
                if 'omim' in selected_gene_special_str:
                    hicoloc_gene_selected_omim.add(gene)
                elif 'exome' in selected_gene_special_str:
                    hicoloc_gene_selected_exome.add(gene)
                elif 'coloc' in selected_gene_special_str:
                    hicoloc_gene_selected_coloc.add(gene)
                elif selected_gene == locus_mindist_gene:
                    hicoloc_gene_selected_mindist.add(gene)
                else:
                    hicoloc_gene_selected_poorsig.add(gene)

        elif gene == locus_mindist_gene:  # selected_gene == locus_mindist_gene:
            himindist_geneset.add(gene)
            if gene == selected_gene:
                himindist_selected.add(gene)
            else:
                himindist_notselected.add(gene)
                if 'omim' in selected_gene_special_str:
                    himindist_gene_selected_omim.add(gene)
                elif 'exome' in selected_gene_special_str:
                    himindist_gene_selected_exome.add(gene)
                elif 'coloc' in selected_gene_special_str:
                    himindist_gene_selected_coloc.add(gene)
                elif selected_gene == locus_mindist_gene:
                    himindist_gene_selected_mindist.add(gene)
                else:
                    himindist_gene_selected_poorsig.add(gene)

        else:
            hipoorsig_geneset.add(gene)
            if gene == selected_gene:
                hipoorsig_selected.add(gene)
            else:
                hipoorsig_notselected.add(gene)
                if 'omim' in selected_gene_special_str:
                    hipoorsig_gene_selected_omim.add(gene)
                elif 'exome' in selected_gene_special_str:
                    hipoorsig_gene_selected_exome.add(gene)
                elif 'coloc' in selected_gene_special_str:
                    hipoorsig_gene_selected_coloc.add(gene)
                elif selected_gene == locus_mindist_gene:
                    hipoorsig_gene_selected_mindist.add(gene)
                else:
                    hipoorsig_gene_selected_poorsig.add(gene)


print(len(hiomim_geneset))
print(len(hiomim_selected))
# print(len(hiomim_notselected))
print(len(hiomim_gene_selected_omim))
print(len(hiomim_gene_selected_exome))
print(len(hiomim_gene_selected_coloc))
print(len(hiomim_gene_selected_mindist))
print(len(hiomim_gene_selected_poorsig))

print('exome')
print(len(hiexome_geneset))
print(len(hiexome_selected))
# print(len(hiexome_notselected))
print(len(hiexome_gene_selected_omim))
print(len(hiexome_gene_selected_exome))
print(len(hiexome_gene_selected_coloc))
print(len(hiexome_gene_selected_mindist))
print(len(hiexome_gene_selected_poorsig))

print('coloc')
print(len(hicoloc_geneset))
print(len(hicoloc_selected))
# print(len(hicoloc_notselected))
print(len(hicoloc_gene_selected_omim))
print(len(hicoloc_gene_selected_exome))
print(len(hicoloc_gene_selected_coloc))
print(len(hicoloc_gene_selected_mindist))
print(len(hicoloc_gene_selected_poorsig))

print('mindist')
print(len(himindist_geneset))
print(len(himindist_selected))
# print(len(himindist_notselected))
print(len(himindist_gene_selected_omim))
print(len(himindist_gene_selected_exome))
print(len(himindist_gene_selected_coloc))
print(len(himindist_gene_selected_mindist))
print(len(himindist_gene_selected_poorsig))

print('poorsig')
print(len(hipoorsig_geneset))
print(len(hipoorsig_selected))
# print(len(hipoorsig_notselected))
print(len(hipoorsig_gene_selected_omim))
print(len(hipoorsig_gene_selected_exome))
print(len(hipoorsig_gene_selected_coloc))
print(len(hipoorsig_gene_selected_mindist))
print(len(hipoorsig_gene_selected_poorsig))


pdb.set_trace()


###########################################
# JUST REMOVE THE SINGLETON OMIM LOCI HERE:
###########################################
#selected_geneset_allomim_incl = set(summary_df[summary_df['gene_selected']=='selected']['gene_name'].values)
#selected_geneset = selected_geneset_allomim_incl  - singletonomim_geneset

summary_df = summary_df[~summary_df['locus_name'].str.startswith('loc')]
selected_geneset = set(
    summary_df[summary_df['gene_signet'] == 'signet']['gene_name'].values)

selected_df = summary_df[summary_df['gene_signet'] == 'signet']

special_geneset = set(
    summary_df[summary_df['gene_special'] != '.']['gene_name'].values)


# number of loci with one 'special' genes, where special means omim, coloc, exome
locus_onespecial_df = summary_df.loc[summary_df['locus_nspecialgene'] == 1]
print('number of loci with one special gene: ', len(
    np.unique(locus_onespecial_df['locus_name'])))
print('identical for the %d loci with a single gene with functional' %
      len(np.unique(locus_onespecial_df['locus_name'])))

locus_onespecial_selectedgene_is_special = locus_onespecial_df.loc[(
    locus_onespecial_df['gene_name'].isin(selected_geneset)) & (locus_onespecial_df['gene_special'] != '.')]
print('which is also the selected gene for %d of these loci' %
      len(locus_onespecial_selectedgene_is_special))


# number of loci with more than one 'special' genes, where special means omim, coloc, exome
locus_multispecial_df = summary_df.loc[summary_df['locus_nspecialgene'] > 1]
print('number of loci with more than one special gene: ',
      len(np.unique(locus_multispecial_df['locus_name'])))
print('The results differ necessarily at the %d loci' %
      len(np.unique(locus_multispecial_df['locus_name'])))


selected_and_mindist = selected_df.loc[selected_df['gene_mindist'] == 'mindist']
print('Signet and Mindist agree at %d loci' % len(selected_and_mindist))
selected_and_notmindist_df = selected_df.loc[selected_df['gene_mindist'] == '.']

print(selected_and_notmindist_df[['gene_name', 'gene_special']].groupby(
    ['gene_special']).count())

locus_some_special_df = summary_df.loc[summary_df['locus_nspecialgene'] != 0]
locus_some_special_selected_df = locus_some_special_df.loc[locus_some_special_df['gene_name'].isin(
    selected_geneset)]
print('Of the loci %d have functional evidence' %
      len(locus_some_special_selected_df))


nselected_sig_and_min = len(set(locus_some_special_selected_df['gene_name']).intersection(
    set(locus_some_special_selected_df['locus_mindistgene'])))
print('Signet and Mindist agree at %d' % nselected_sig_and_min)


locus_some_special_selected_notmin_special = (set(locus_some_special_selected_df['gene_name']) - set(
    locus_some_special_selected_df['locus_mindistgene'])).intersection(special_geneset)
print('and Signet selects a more distant gene that has functional evidence at %d loci',
      len(locus_some_special_selected_notmin_special))

locus_no_special_df = summary_df.loc[summary_df['locus_nspecialgene'] == 0]
locus_no_special_selected_df = locus_no_special_df .loc[locus_no_special_df['gene_name'].isin(
    selected_geneset)]
nlocus_no_special_selected_min_special = len(set(locus_no_special_selected_df['gene_name']).intersection(
    set(locus_no_special_selected_df['locus_mindistgene'])))
print('Among the information-poor loci, Signet and Mindist agreed at %d loci' %
      nlocus_no_special_selected_min_special)

pdb.set_trace()


omim_geneset - mindist_geneset
print('Of the %d Mendelian genes in the loci, %d were not selected by Mindist' %
      (len(omim_geneset), len(omim_geneset - mindist_geneset)))


print('The results differ necessarily at the %d loci with multiple genes with functional evidence, with Signet selecting one of these genes at %d of the loci.' %
      (len(np.unique(locus_multispecial_df['locus_name'])), nlocus_multispecial_selectedspecial))


locus_somespecial_df = summary_df.loc[summary_df['locus_nspecialgene'] >= 1]

print('Of the n_loci loci, %d have functional evidence.' %
      (len(np.unique(locus_somespecial_df['locus_name']))))


selectedgenes_locus_somespecial = set(
    locus_somespecial_df.loc[locus_somespecial_df['gene_prob_selected'] >= 0.5]['gene_name'])
mindistgenes_locus_somespecial = set(locus_somespecial_df['locus_mindistgene'])


n_withspecial_agree = len(
    selectedgenes_locus_somespecial.intersection(mindistgenes_locus_somespecial))
print('Of the loci with functional evidence, Signet and Mindist agree at %d' %
      (n_withspecial_agree))


print('The number of information-poor loci, lacking strong functional evidence, is %d' %
      len(np.unique(locus_nospecial_df['locus_name'])))


selectedgenes_locus_nospecial = set(
    locus_nospecial_df.loc[locus_nospecial_df['gene_prob_selected'] >= 0.5]['gene_name'])
mindistgenes_locus_nospecial = set(locus_nospecial_df['locus_mindistgene'])
n_nospecial_agree = len(
    selectedgenes_locus_nospecial.intersection(mindistgenes_locus_nospecial))

print('Among the information-poor loci, Signet and Mindist agreed at %d loci' %
      n_nospecial_agree)


project_datadir = '/Users/z/Dropbox/Bader_Mousavi/GeneNetwork/data/TopMed'
# locus and genes
file1path = os.path.join(project_datadir, 'KEGG_2021_nospecial_selected.txt')
file2path = os.path.join(project_datadir,  'KEGG_2021_nospecial_mindist.txt')

file1name = 'nospecial_selected'
file2name = 'nospecial_mindist'

counts_postprocessing.kegg_pathway_comparison(file1path, file1name, file2path,
                                              file2name, results_dir, '')


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)
    # level 2 interactions:
    summary_file = os.path.join(
        results_dir, 'summary_genestats_manuscript.txt')
    selected_file = os.path.join(results_dir, 'summary_genestats_selected.txt')

    for anchor_gene in ['GNG3', 'FNDC4']:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_file, selected_file, network_ppi, gene_special, TFsource_target)

    # utils_plot.plot_InteractionNetwork(results_dir, summary_file, network_ppi,
    #                                    gene_special,  TFsource_target, prob_threshold_min=0.5, prob_threshold_max=1)


selected_df = pd.DataFrame(columns=summary_df.columns)
notselected_df = pd.DataFrame(columns=summary_df.columns)

for locus in summary_df['locus_name'].unique():

    locus_df = summary_df[summary_df['locus_name'] == locus].copy()
    max_prob = np.max(locus_df['gene_prob_selected'])

    for idx, row in locus_df.iterrows():

        if row['gene_prob_selected'] == max_prob:
            selected_df = selected_df.append(
                locus_df.loc[[idx]], ignore_index=True)

        else:
            notselected_df = notselected_df.append(
                locus_df.loc[[idx]], ignore_index=True)

selected_file = os.path.join(results_dir, 'summary_genestats_selected.txt')
selected_df.to_csv(selected_file, sep='\t')

notselected_file = os.path.join(
    results_dir, 'summary_genestats_notselected.txt')
notselected_df.to_csv(notselected_file, sep='\t')


selected_geneset = set(selected_df['gene_name'])
notselected_geneset = set(notselected_df['gene_name'])


pdb.set_trace()


##########
# Venn Diagrams
# New Table 4
#########


A_set = omim_geneset
A_label = 'omim'
B_set = coloc_geneset
B_label = 'coloc'
C_set = exome_geneset
C_label = 'exome'
utils_plot.plot_venn(A_set, B_set, C_set, A_label, B_label,
                     C_label,  'venn_specialgenes.pdf', results_dir)


pdb.set_trace()


##########
# others
##########


# feature stats
counts_postprocessing.write_feature_weight_dict(results_dir)
prob1 = len(summary_df.loc[(summary_df['gene_prob_selected'] == 1)])
print('final active genes were identical across all 100 runs for %d of the ntotal loci' % prob1)


# Special gene breakdown
print('special gene breakdown')
print(summary_df[['gene_name', 'gene_special']
                 ].groupby(['gene_special']).count())

ngene_special = sum(summary_df['gene_special'] != '.')
print('total number of unique special genes: %d' % (ngene_special))

# special_geneset = omim_geneset.union(exome_geneset).union(coloc_geneset) #same

special_geneset = set(
    summary_df.loc[summary_df['gene_special'] != '.']['gene_name'])
# check CREBBP, ID2, TBX5
interaction_geneset = set(summary_df.loc[(summary_df['ppi_deg'] != 0) | (
    summary_df['TF_indeg'] != 0) | (summary_df['TF_outdeg'] != 0)]['gene_name'])
gene_universe = set(summary_df['gene_name'])
mindist_geneset = set(summary_df['locus_mindistgene'])


# summary_df.loc[summary_df['gene_type']]
# FROM BEFORE FEB25
# selected_df = summary_df.loc[(summary_df['gene_prob_selected'] >= 0.5)]
# accepted_types = set({'omim', 'coloc', 'protein_coding'})
# unaccepted_df = selected_df.loc[~selected_df['gene_type'].isin(accepted_types)]

# unaccepted_selected_filename = os.path.join(
#     results_dir, 'unaccepted_selected.csv')
# unaccepted_df.to_csv(unaccepted_selected_filename, index=False, header=True)
# np.unique(selected_df['gene_type'])
# # EEEND


if 0:
    with open(results_dir + '/TFsource_target.pickle', 'rb') as handle:
        TFsource_target = pickle.load(handle)

    with open(results_dir + '/network_ppi.pickle', 'rb') as handle:
        network_ppi = pickle.load(handle)

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)
    # level 2 interactions:
    pdb.set_trace()
    for anchor_gene in selected_geneset:
        utils_plot.plot_InteractionNetwork_oup(
            anchor_gene, results_dir, summary_file, selected_file, network_ppi, gene_special, TFsource_target)


pdb.set_trace()


# of selected genes how many are special
# set(selected_df.loc[selected_df['gene_special'] != '.']['gene_name'])
special_selected_geneset = special_geneset.intersection(selected_geneset)
mindist_selected_geneset = mindist_geneset.intersection(selected_geneset)
interaction_selected_geneset = interaction_geneset.intersection(
    selected_geneset)


A_set = interaction_selected_geneset
A_label = 'PPI, GRI'
B_set = mindist_selected_geneset
B_label = 'mindist'
C_set = special_selected_geneset
C_label = 'functional'
utils_plot.plot_venn(A_set, B_set, C_set, A_label, B_label,
                     C_label, 'venn_selected_genes.pdf', results_dir)


omim_selected_geneset = omim_geneset.intersection(selected_geneset)
coloc_selected_geneset = coloc_geneset.intersection(selected_geneset)
exome_selected_geneset = exome_geneset.intersection(selected_geneset)

A_set = omim_selected_geneset
A_label = 'omim'
B_set = coloc_selected_geneset
B_label = 'coloc'
C_set = exome_selected_geneset
C_label = 'exome'

utils_plot.plot_venn(A_set, B_set, C_set, A_label, B_label,
                     C_label, 'venn_selected_specialgenes.pdf', results_dir)


print('aloooooooooooo')
pdb.set_trace()

#
omim_notselected = omim_geneset - selected_geneset

print('total number of special genes that were selected: %d' %
      (len(special_selected_geneset)))
print('of the %d genes with any functional evidence, %d were the selected gene' %
      (ngene_special, len(special_selected_geneset)))

pdb.set_trace()

special_notselected_df = notselected_df.loc[notselected_df['gene_special'] != '.']
special_notselected_geneset = set(special_notselected_df['gene_name'])
print('total number of special genes that were not selected: %d' %
      (len(special_notselected_geneset)))


pdb.set_trace()
somespecial_notselected_locuslist = np.unique(
    special_notselected_df['locus_name'])
somespecial_notselected_df = summary_df.loc[summary_df['locus_name'].isin(
    somespecial_notselected_locuslist)]
somespecial_notselected_selectedgenes = somespecial_notselected_df.loc[somespecial_notselected_df['gene_name'].isin(
    selected_geneset)]
print('Of the %d genes with functional evidence that were not selected, %d were in loci where the selected gene did have functional evidence.' %
      (len(special_notselected_geneset), len(somespecial_notselected_selectedgenes)))


# number of loci that lack 'special' genes, where special means omim, coloc, exome
locus_nospecial_df = summary_df.loc[summary_df['locus_nspecialgene'] == 0]
nospecial_file = os.path.join(results_dir, 'nospecial_genestats.txt')
locus_nospecial_df.to_csv(nospecial_file, sep='\t', index=False, header=True)
print('number of loci with no special gene: ', len(
    np.unique(locus_nospecial_df['locus_name'])))
print('and for the %d information-poor loci with no genes with functional evidence' %
      len(np.unique(locus_nospecial_df['locus_name'])))


locus_multispecial_selectedgene_is_special = locus_multispecial_df.loc[(
    locus_multispecial_df['gene_name'].isin(selected_geneset)) & (locus_multispecial_df['gene_special'] != '.')]
print('Signet selecting one of these genes at %d of the loci' %
      len(locus_multispecial_selectedgene_is_special))


#############
#############
#############
#############
locus_nospecial_summary_df = locus_nospecial_df.loc[locus_nospecial_df['gene_name'].isin(
    selected_geneset)][['locus_name', 'locus_mindistgene', 'gene_name', 'gene_prob_selected']]
nospecial_summary_file = os.path.join(
    results_dir, 'nospecial_summary_genestats.txt')
locus_nospecial_summary_df.to_csv(
    nospecial_summary_file, sep='\t', index=False, header=True)


pdb.set_trace()

# number of loci with one 'special' gene, with selected 'not' the special gene
locus_onespecial_selectednotspecial = locus_onespecial_df.loc[(
    locus_onespecial_df['gene_prob_selected'] >= 0.5) & locus_onespecial_df['gene_special'] == '.']
# len(np.unique(locus_onespecial_df.loc[(locus_onespecial_df['gene_prob_selected'] >= 0.5) & locus_onespecial_df['gene_special'] != '.']['locus_name']))


# number of loci with multi 'special' genes, with selected being a special gene
nlocus_multispecial_selectedspecial = len(np.unique(locus_multispecial_df.loc[(
    locus_multispecial_df['gene_prob_selected'] >= 0.5) & locus_multispecial_df['gene_special'] != '.']['locus_name']))
