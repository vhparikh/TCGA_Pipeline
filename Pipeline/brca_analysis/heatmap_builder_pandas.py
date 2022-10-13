from urllib.parse import urlparse

import mlflow
import matplotlib

import pandas as pd
from matplotlib.patches import Patch
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns

print('Start')
data_short = []
data_full = []
genes = []

first = True

luma = pd.read_csv('../../BRCA/BRCA_Subtypes/LumA_top_20.csv')['GENE']
lumb = pd.read_csv('../../BRCA/BRCA_Subtypes/LumB_top_20.csv')['GENE']
her2 = pd.read_csv('../../BRCA/BRCA_Subtypes/Her2_top_20.csv')['GENE']
basal = pd.read_csv('../../BRCA/BRCA_Subtypes/Basal_top_20.csv')['GENE']

allgenes = pd.read_csv('../../BRCA/BRCA_Subtypes/LumA.csv')['GENE']

all_subtypes = pd.concat([basal, her2, luma, lumb])
unique_genes = all_subtypes.unique()
unique_genes = list(unique_genes)
dataframe_short = pd.read_csv('../../BRCA/BRCA_Dropped.csv', usecols=[*unique_genes, 'SUBTYPE'])
dataframe_short = dataframe_short.reindex(dataframe_short['SUBTYPE'].sort_values(ascending=False).index)
subtype = dataframe_short['SUBTYPE']
print(subtype.unique())

dataframe_full = pd.read_csv('../../BRCA/BRCA_Dropped.csv', usecols = [*allgenes, 'SUBTYPE'])
dataframe_full = dataframe_full.reindex(dataframe_full['SUBTYPE'].sort_values(ascending=False).index)


# reformatting data short
for index, row in dataframe_short.iterrows():
    if first:
        sample_names = row[1:]
        first = False

    else:
        genes.append(row[0])
        data_short.append(row[1:])

for index, row in dataframe_full.iterrows():
    if first:
        sample_names = row[1:]
        first = False

    else:
        genes.append(row[0])
        data_full.append(row[1:])

# creating dataframe
data_short = pd.DataFrame(dataframe_short)
data_full = pd.DataFrame(dataframe_full)

# getting gene columns
columns_short = data_short.columns
numerical_cols_short = list(columns_short)[:-1]

columns_full = data_full.columns
numerical_cols_full = list(columns_full)[:-1]

lut = dict(zip(subtype.unique(), ['MediumAquamarine', 'Salmon', 'Green', 'DeepSkyBlue']))
row_colors = subtype.map(lut)

# creating heatmap
with mlflow.start_run(run_name='Differential Expression Heatmaps png'):
    print('creating heatmap')
    sns.set_context("paper", font_scale=1.3)
    sns_plot = sns.clustermap(data_short[numerical_cols_short], xticklabels=True, yticklabels=False, row_colors=row_colors, row_cluster=False,
                              figsize=(15, 10))
    path_short = Path('heatmap_20.png')
    sns_plot.savefig("heatmap_20.png")

    print('creating legend')
    #creating legend
    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(handles, lut, title='Subtypes',
               bbox_to_anchor=(0, 0), bbox_transform=plt.gcf().transFigure, loc='lower left', handleheight=2, handlelength=7,
               fontsize=7, title_fontsize='xx-large')

    mlflow.log_artifact(path_short)

    print('creating heatmap')
    sns.set_context("paper", font_scale=1.3)
    sns_plot = sns.clustermap(data_full[numerical_cols_full], xticklabels=False, yticklabels=False,
                              row_colors=row_colors, row_cluster=False,
                              figsize=(15, 10))
    path_full = Path('heatmap_full.png')
    sns_plot.savefig("heatmap_full.png")

    print('creating legend')
    # creating legend
    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(handles, lut, title='Subtypes',
               bbox_to_anchor=(0, 0), bbox_transform=plt.gcf().transFigure, loc='lower left', handleheight=2,
               handlelength=7,
               fontsize=7, title_fontsize='xx-large')

    mlflow.log_artifact(path_full)

    tracking_url_type_store = urlparse(mlflow.get_tracking_uri()).scheme