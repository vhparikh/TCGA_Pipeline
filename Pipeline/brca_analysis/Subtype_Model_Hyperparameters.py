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

x_train, x_test, y_train, y_test = train_test_split(brca, y_encoded, test_size=0.2, random_state=10)
y_test_df = pd.DataFrame

x_train = transformer_1.fit_transform(x_train)
x_test = transformer_2.fit_transform(x_test)

# print update
print('Training model...')

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