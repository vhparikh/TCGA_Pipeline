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
import mlflow.sklearn
from urllib.parse import urlparse
from pathlib import Path


# reading data
print('Running tool---')
# get the data
genes = list(pd.read_csv('../../LGG/top_500_genes_lgg.csv')['GENES']) #top 500 from ETC
lgg = pd.read_csv('../../LGG/LGG.csv', usecols=[*genes, 'SUBTYPE'])
lgg.to_csv('LGG_500.csv')
print('done')
subtype_lgg = pd.DataFrame(columns=['SUBTYPE'])
subtypes = lgg['SUBTYPE']
print(subtypes.unique())
subtype_lgg['SUBTYPE'] = subtypes

lgg.drop(['SUBTYPE'], axis=1, inplace=True)

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

# print update
print('X data transformed')

x_train, x_test, y_train, y_test = train_test_split(lgg, y_encoded, test_size=0.2, random_state=10)
y_test_df = pd.DataFrame

x_train = transformer_1.fit_transform(x_train)
x_test = transformer_2.fit_transform(x_test)

# print update
print('Training model...')

# training model
model = LogisticRegression(C=78.47599703514607, max_iter=5000, penalty='l2', solver='lbfgs')
model = model.fit(x_train, y_train)
y_hat = model.predict(x_test)

accuracy = accuracy_score(y_test, y_hat)
r2 = r2_score(y_test, y_hat)
precision = precision_score(y_test, y_hat, average='micro')

labels = ['IDHmut-non-codel' 'IDHwt' 'IDHmut-codel']
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
path = Path('Visuals/Tuned_Confusion_Matrix.svg')
plt.savefig('Visuals/Tuned_Confusion_Matrix.svg')
#mlflow.log_artifact(path)

values = model.coef_
coefficients_idhwt_nc = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[0]})
coefficients_idhwt_nc = coefficients_idhwt_nc.reindex(
    coefficients_idhwt_nc['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_idhwt_nc = Path('Coefficients/IDHwt_NC_Coefficients.csv')
coefficients_idhwt_nc.to_csv(path_idhwt_nc, index=False)
#mlflow.log_artifact(path_idhwt_nc, 'Coefficients')

coefficients_idhwt = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[1]})
coefficients_idhwt = coefficients_idhwt.reindex(
    coefficients_idhwt['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_idhwt = Path('Coefficients/IDHwt_Coefficients.csv')
coefficients_idhwt.to_csv(path_idhwt, index=False)
#mlflow.log_artifact(path_idhwt, 'Coefficients')

coefficients_idhwt_c = pd.DataFrame({'GENE': genes, 'COEFFICIENT': values[2]})
coefficients_idhwt_c = coefficients_idhwt_c.reindex(
    coefficients_idhwt_c['COEFFICIENT'].abs().sort_values(ascending=False).index)
path_idhwt_c = Path('Coefficients/IDHwt_C_Coefficients.csv')
coefficients_idhwt_c.to_csv(path_idhwt_c, index=False)
#mlflow.log_artifact(path_idhwt_c, 'Coefficients')

#tracking_url_type_store = urlparse(mlflow.get_tracking_uri()).scheme