import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.preprocessing import OrdinalEncoder, Normalizer
from sklearn.model_selection import train_test_split

# reading data
print('Running tool---')

# get the data
brca = pd.read_csv('../../BRCA/LGG.csv')
subtype_brca = pd.DataFrame(columns=['SUBTYPE'])
subtypes = brca['SUBTYPE']
print(subtypes.unique())
subtype_brca['SUBTYPE'] = subtypes
brca.drop(['SAMPLE_BARCODE', 'log10_mut', 'PIK3CA_snv', 'DISEASE', 'SUBTYPE'], axis=1, inplace=True)

genes = brca.columns
print(genes)
print('Transforming data...')
# Encoding y matrix
enc = OrdinalEncoder()
y_encoded = enc.fit_transform(subtype_brca)
y_encoded = y_encoded.ravel()

print('Y data transformed')

# transforming x matrix
transformer = Normalizer()
x_encoded = transformer.fit_transform(brca)

# print update
print('X data transformed')

x_train, x_test, y_train, y_test = train_test_split(x_encoded, y_encoded, test_size=0.1, random_state=10)

etc = ExtraTreesClassifier(n_estimators=50)
etc.fit(x_encoded, y_encoded)
feature_importances = etc.feature_importances_
feature_importances_df = pd.DataFrame({'GENES': genes, 'IMPORTANCE': feature_importances})
feature_importances_df = feature_importances_df.sort_values(by='IMPORTANCE', ascending=False)
top_500_genes = feature_importances_df.head(500)['GENES']
top_500_genes.to_csv('../../LGG/top_500_genes_brca.csv', index=False)