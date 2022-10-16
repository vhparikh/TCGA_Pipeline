from urllib.parse import urlparse

import mlflow
import pandas as pd
from scipy.stats import ttest_ind

#read in subtype coefficients
luma_genes = pd.read_csv('Coefficients/LUMA_Coefficients.csv')['GENE'].head(200)
her2_genes = pd.read_csv('Coefficients/HER2_Coefficients.csv')['GENE'].head(200)
lumb_genes = pd.read_csv('Coefficients/LUMB_Coefficients.csv')['GENE'].head(200)
basal_genes = pd.read_csv('Coefficients/BASAL_Coefficients.csv')['GENE'].head(200)

#get x matrix for each subtype
x_matrix = pd.read_csv('../../../BRCA/BRCA_Dropped.csv', usecols=[*luma_genes, *her2_genes, *lumb_genes, *basal_genes, 'SUBTYPE'])

subtypes = ['LumA', 'Her2', 'LumB', 'Basal']
genes = list(x_matrix.columns)[:-1]

#splitting data by subtypes
for subtype in subtypes:
    subtype_df = x_matrix.loc[x_matrix['SUBTYPE'] == subtype]
    comparison_df = x_matrix.loc[x_matrix['SUBTYPE'] != subtype]
    subtype_impact_factor = pd.DataFrame(columns=['GENE', 'FACTOR', 'PVALUE'])
    for gene in genes:
        total_base = 0.0
        total_subtype = 0.0
        base_avg = 0.0
        subtype_avg = 0.0
        factor = 0.0

        gene_p_value = ttest_ind(comparison_df[gene], subtype_df[gene])
        # gets base subtype gene expression sum
        for sample in list(comparison_df[gene]):
            total_base = total_base + float(sample)

        # gets subtype gene expression sum
        for sample in list(subtype_df[gene]):
            total_subtype = total_subtype + float(sample)

        # compares base subtype and subtype gene expressions
        if total_base != float(0) or total_subtype != float(0):
            base_avg = total_base / len(comparison_df.index)
            subtype_avg = total_subtype / len(subtype_df.index)
            factor = base_avg / subtype_avg
            subtype_impact_factor.loc[len(subtype_impact_factor)] = [gene, factor, gene_p_value[1]]

    subtype_impact_factor = subtype_impact_factor.sort_values(by=['PVALUE'], ascending=True)
    top_20_subtype_genes = subtype_impact_factor.head(20)['GENE']
    filepath_genes = '../../../BRCA/BRCA_Subtypes/' + str(subtype) + '_top_20.csv'
    filepath_full = '../../../BRCA/BRCA_Subtypes/' + str(subtype) + '.csv'
    subtype_impact_factor.to_csv(filepath_full, index=False)
    top_20_subtype_genes.to_csv(filepath_genes, index=False)