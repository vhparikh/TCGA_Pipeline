from urllib.parse import urlparse

import mlflow
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

#get x_matrix and top_genes


#read in subtype coefficients
IDHwt_c_genes = pd.read_csv('Coefficients/IDHwt_C_Coefficients.csv')['GENE'].head(200)
IDHwt_genes = pd.read_csv('Coefficients/IDHwt_Coefficients.csv')['GENE'].head(200)
IDHwt_nc_genes = pd.read_csv('Coefficients/IDHwt_NC_Coefficients.csv')['GENE'].head(200)

#get x matrix for each subtype
x_matrix = pd.read_csv('../../LGG/LGG.csv', usecols=[*IDHwt_c_genes, *IDHwt_genes, *IDHwt_nc_genes, 'SUBTYPE'])
x_matrix.to_csv('LGG_Diff_Genes.csv')
print('done')
#print(x_matrix.head())

subtypes = ['IDHmut-non-codel', 'IDHwt', 'IDHmut-codel']

#df_dict = {}
genes = list(x_matrix.columns)[:-1]

with mlflow.start_run(run_name='Differential Expression'):

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
        filepath_genes = '../../LGG/LGG_Subtypes/' + str(subtype) + '_top_20.csv'
        filepath_full = '../../LGG/LGG_Subtypes/' + str(subtype) + '.csv'
        subtype_impact_factor.to_csv(filepath_full, index=False)
        top_20_subtype_genes.to_csv(filepath_genes, index=False)
        mlflow.log_artifact(filepath_full)
        mlflow.log_artifact(filepath_genes)
        tracking_url_type_store = urlparse(mlflow.get_tracking_uri()).scheme