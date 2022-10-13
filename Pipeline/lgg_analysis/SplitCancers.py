import pandas as pd
import numpy as np

x_matrix = pd.read_csv("../../Xy_matrices/PIK3CA/X_matrix.tsv", delimiter='\t')
y_matrix = pd.read_csv("../../Xy_matrices/PIK3CA/y_matrix.tsv", delimiter='\t', usecols=['SAMPLE_BARCODE', 'DISEASE', 'PIK3CA_snv', 'SUBTYPE'])

mutations = y_matrix['PIK3CA_snv'].to_list()
x_matrix['PIK3CA_snv'] = mutations

diseases = y_matrix['DISEASE']
x_matrix['DISEASE'] = diseases

subtypes = y_matrix['SUBTYPE']
x_matrix['SUBTYPE'] = subtypes

disease_names = ['GBM', 'OV', 'LUAD', 'LUSC', 'PRAD', 'UCEC', 'BLCA', 'ESCA', 'LIHC', 'SARC', 'LGG', 'COAD', 'STAD', 'SKCM', 'KIRC', 'CESC', 'HNSC', 'READ', 'LGG', 'UCS']

for disease in disease_names:
    df = x_matrix.loc[x_matrix['DISEASE'] == disease]
    filepath = '../../Cancers/' + disease + '.csv'
    df.to_csv(filepath, index=False)
