from pathlib import Path

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

# reading data
print('Running tool---')
# get the data
brca = pd.read_csv('../../BRCA/BRCA_Dropped.csv')
subtype_brca = pd.DataFrame(columns=['SUBTYPE'])
subtypes = brca['SUBTYPE']
print(subtypes.unique())
subtype_brca['SUBTYPE'] = subtypes

brca.drop(['SUBTYPE', 'SAMPLE_BARCODE', 'DISEASE', 'PIK3CA_snv', 'log10_mut'], axis=1, inplace=True)

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

#with mlflow.start_run(run_name='Basic Model'):
#training model
model = LogisticRegression(solver='sag')
model.fit(x_train, y_train)
print('Model trained')

y_hat = model.predict(x_test)

accuracy = accuracy_score(y_test, y_hat)
r2 = r2_score(y_test, y_hat)
precision = precision_score(y_test, y_hat, average='micro')

labels = ['Luminal A', 'HER2', 'Luminal B', 'Triple-Negative']
print(classification_report(y_test, y_hat, target_names=labels))

#Confusion Matrix Builder
cm = [[34, 0, 2, 0],
      [4, 0, 8, 0],
      [0, 0, 110, 0],
      [0, 0, 26, 4]]
print(cm)
matrix = plot_confusion_matrix(model, x_test, y_test, cmap=plt.cm.Reds, normalize=None)
matrix.ax_.set_title('Confusion Matrix', color='black')
plt.xlabel('Predicted Subtype', color='black')
plt.ylabel('Actual Subtype', color='black')
plt.gcf().axes[0].tick_params(colors='black')
plt.gcf().axes[1].tick_params(colors='black')
plt.gcf().set_size_inches(10, 6)
plt.rcParams['font.size'] = '30'
path = Path('Visuals/Confusion_Matrix.svg')
plt.savefig('Visuals/Confusion_Matrix.svg')