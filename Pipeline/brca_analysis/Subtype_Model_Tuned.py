import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import Normalizer, OrdinalEncoder
from sklearn.metrics import plot_confusion_matrix, accuracy_score, r2_score, precision_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import classification_report
import mlflow
from urllib.parse import urlparse
from pathlib import Path

# reading data
print('Running tool---')
# get the data
genes = list(pd.read_csv('../../../BRCA/top_500_genes_brca.csv')['GENES']) #top 500 from ETC
brca = pd.read_csv('../../../BRCA/BRCA_Dropped.csv', usecols=[*genes, 'SUBTYPE'])
subtype_brca = pd.DataFrame(columns=['SUBTYPE'])
subtypes = brca['SUBTYPE']
print(subtypes.unique())
subtype_brca['SUBTYPE'] = subtypes

brca.drop(['SUBTYPE'], axis=1, inplace=True)

print('Transforming data...')
# Encoding y matrix
enc = OrdinalEncoder()
print(subtype_brca['SUBTYPE'].unique())
y_encoded = enc.fit_transform(subtype_brca)
y_encoded = y_encoded.ravel()
print(np.unique(y_encoded))

print('Y data transformed')

# transforming x matrix
transformer_1 = Normalizer()
transformer_2 = Normalizer()
# print update
print('X data transformed')

x_train, x_test, y_train, y_test = train_test_split(brca, y_encoded, test_size=0.1, random_state=10)
y_test_df = pd.DataFrame

x_train = transformer_1.fit_transform(x_train)
x_test = transformer_2.fit_transform(x_test)


# print update
print('Training model...')


# training model
model = LogisticRegression(C=3792.690190732246, max_iter=10000, penalty='l2', solver='sag')
model = model.fit(x_train, y_train)
y_hat = model.predict(x_test)

accuracy = accuracy_score(y_test, y_hat)
r2 = r2_score(y_test, y_hat)
precision = precision_score(y_test, y_hat, average='micro')

labels = ['Luminal A', 'HER2', 'Luminal B', 'Triple-Negative']
print(classification_report(y_test, y_hat, target_names=labels))

# Confusion Matrix Builder
cm = confusion_matrix(y_test, y_hat)
#print(cm)
matrix = plot_confusion_matrix(model, x_test, y_test, cmap=plt.cm.Reds, normalize=None)
matrix.ax_.set_title('Confusion Matrix', color='black')
plt.xlabel('Predicted Subtype', color='black')
plt.ylabel('Actual Subtype', color='black')
plt.gcf().axes[0].tick_params(colors='black')
plt.gcf().axes[1].tick_params(colors='black')
plt.gcf().set_size_inches(10, 6)
plt.rcParams['font.size'] = '30'
#plt.show()
path = Path('Visuals/Tuned_Confusion_Matrix.png')
plt.savefig('Visuals/Tuned_Confusion_Matrix.png')

values = model.coef_
coefficients_LUMA = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[0]})
coefficients_LUMA = coefficients_LUMA.reindex(
    coefficients_LUMA['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_LUMA = Path('Coefficients/LUMA_Coefficients.csv')
coefficients_LUMA.to_csv(path_LUMA, index=False)

coefficients_HER2 = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[1]})
coefficients_HER2 = coefficients_HER2.reindex(
    coefficients_HER2['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_HER2 = Path('Coefficients/HER2_Coefficients.csv')
coefficients_HER2.to_csv(path_HER2, index=False)

coefficients_LUMB = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[2]})
coefficients_LUMB = coefficients_LUMB.reindex(
    coefficients_LUMB['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_LUMB = Path('Coefficients/LUMB_Coefficients.csv')
coefficients_LUMB.to_csv(path_LUMB, index=False)

coefficients_BASAL = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[3]})
coefficients_BASAL = coefficients_BASAL.reindex(
    coefficients_BASAL['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_BASAL = Path('Coefficients/BASAL_Coefficients.csv')
coefficients_BASAL.to_csv(path_BASAL, index=False)