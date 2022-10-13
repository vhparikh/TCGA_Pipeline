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
from urllib.parse import urlparse


# reading data
print('Running tool---')
# get the data
genes = list(pd.read_csv('../../LGG/top_500_genes_lgg.csv')['GENES']) #top 500 from ETC
lgg = pd.read_csv('../../LGG/LGG.csv', usecols=[*genes, 'SUBTYPE'])
subtype_lgg = pd.DataFrame(columns=['SUBTYPE'])
subtypes = lgg['SUBTYPE']
print(subtypes.unique())
subtype_lgg['SUBTYPE'] = subtypes

lgg.drop(['SUBTYPE'], axis=1, inplace=True)
#'SAMPLE_BARCODE', 'DISEASE', 'PIK3CA_snv', 'log10_mut',
#^^ if not using top 500 genes csv need to drop

print('Transforming data...')
# Encoding y matrix
enc = OrdinalEncoder()
print(subtype_lgg['SUBTYPE'].unique())
y_encoded = enc.fit_transform(subtype_lgg)
y_encoded = y_encoded.ravel()
print(np.unique(y_encoded))

print('Y data transformed')

# transforming x matrix
transformer_1 = Normalizer()
transformer_2 = Normalizer()
#x_encoded = transformer.fit_transform(lgg)

# print update
print('X data transformed')

x_train, x_test, y_train, y_test = train_test_split(lgg, y_encoded, test_size=0.2, random_state=10)
y_test_df = pd.DataFrame

x_train = transformer_1.fit_transform(x_train)
x_test = transformer_2.fit_transform(x_test)
# print update
print('Training model...')

#solver = 'sag'
#penalty = 'l2'
#max_iter = 5000
#C = 3792.690190732246
d
# training model
model = LogisticRegression()

param_grid = [
    {'penalty': ['l1', 'l2', 'elasticnet'],
     'C': np.logspace(-4, 4, 20),
     'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'],
     'max_iter': [5000, 7500, 10000]
    }
]

clf = GridSearchCV(model, param_grid=param_grid, cv=5, verbose=2, n_jobs=-1)
best_clf = clf.fit(x_train, y_train)
print('score: ', clf.best_score_)
print('parameters: ', clf.best_params_)
model.fit(x_train, y_train)
print('Model trained')