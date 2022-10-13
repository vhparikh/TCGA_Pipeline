import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import Normalizer, OrdinalEncoder
from sklearn.metrics import plot_confusion_matrix, accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import classification_report

# reading data
print('Running tool---')
# get the data
genes = list(pd.read_csv('../../BRCA/top_500_genes_brca.csv')['GENES']) #top 500 from ETC
brca = pd.read_csv('../../BRCA/BRCA_Dropped.csv', usecols=[*genes, 'SUBTYPE'])
subtype_brca = pd.DataFrame(columns=['SUBTYPE'])
subtypes = brca['SUBTYPE']
print(subtypes.unique())
subtype_brca['SUBTYPE'] = subtypes

brca.drop(['SUBTYPE'], axis=1, inplace=True)
#'SAMPLE_BARCODE', 'DISEASE', 'PIK3CA_snv', 'log10_mut',
#^^ if not using top 500 genes csv need to drop

print('Transforming data...')
# Encoding y matrix
enc = OrdinalEncoder()
print(subtype_brca['SUBTYPE'].unique())
y_encoded = enc.fit_transform(subtype_brca)
y_encoded = y_encoded.ravel()
print(np.unique(y_encoded))

print('Y data transformed')

# transforming x matrix
transformer = Normalizer()
x_encoded = transformer.fit_transform(brca)

# print update
print('X data transformed')

x_train, x_test, y_train, y_test = train_test_split(x_encoded, y_encoded, test_size=0.2, random_state=10)
y_test_df = pd.DataFrame

# print update
print('Training model...')

# training model
model = LogisticRegression(C=3792.690190732246, max_iter=5000, penalty='l2', solver='sag')
model.fit(x_train, y_train)
print('Model trained')

y_hat = model.predict(x_test)

print(accuracy_score(y_test, y_hat))
labels = ['Luminal A', 'HER2', 'Luminal B', 'Triple-Negative']
print(classification_report(y_test, y_hat, target_names=labels))

#Getting Model Coefficients
values = model.coef_
print(values)

#Placing Coefficients into df and sorting by abs
coefficients_LUMA = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[0]})
coefficients_LUMA = coefficients_LUMA.reindex(coefficients_LUMA['COEFFICIENT'].abs().sort_values(ascending=False).index)
coefficients_LUMA.to_csv('../../LGG/LUMA_Coefficients.csv', index=False)

coefficients_HER2 = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[1]})
coefficients_HER2 = coefficients_HER2.reindex(coefficients_HER2['COEFFICIENT'].abs().sort_values(ascending=False).index)
coefficients_HER2.to_csv( '../../LGG/HER2_Coefficients.csv', index=False)

coefficients_LUMB = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[2]})
coefficients_LUMB = coefficients_LUMB.reindex(coefficients_LUMB['COEFFICIENT'].abs().sort_values(ascending=False).index)
coefficients_LUMB.to_csv( '../../LGG/LUMB_Coefficients.csv', index=False)

coefficients_BASAL = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[3]})
coefficients_BASAL = coefficients_BASAL.reindex(coefficients_BASAL['COEFFICIENT'].abs().sort_values(ascending=False).index)
coefficients_BASAL.to_csv('../../LGG/BASAL_Coefficients.csv', index=False)

#Confusion Matrix Builder
cm = confusion_matrix(y_test, y_hat)
print(cm)
matrix = plot_confusion_matrix(model, x_test, y_test, cmap=plt.cm.Reds, normalize=None)
matrix.ax_.set_title('Confusion Matrix', color='black')
plt.xlabel('Predicted Subtype', color='black')
plt.ylabel('Actual Subtype', color='black')
plt.gcf().axes[0].tick_params(colors='black')
plt.gcf().axes[1].tick_params(colors='black')
plt.gcf().set_size_inches(10, 6)
plt.rcParams['font.size'] = '30'
plt.show()


