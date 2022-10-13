import pandas as pd

# reading data
print('Running tool---')
# get the data
brca = pd.read_csv('../../Cancers/BRCA.csv')
genes = brca.columns
subtypes = brca['SUBTYPE']
print(len(subtypes))
print('LumA', len(brca[brca['SUBTYPE'] == 'LumA']))
print('LumB', len(brca[brca['SUBTYPE'] == 'LumB']))
print('Her2', len(brca[brca['SUBTYPE'] == 'Her2']))
print('Basal', len(brca[brca['SUBTYPE'] == 'Basal']))
print(subtypes.unique())
#print(len(genes))
count = 0

#print(len(brca[brca['SUBTYPE']=='LumA']))
#print(brca.shape)
#brca_dropped = gbm[gbm.SUBTYPE != 'Normal']
#print(brca_dropped.shape)
#print(len(brca_dropped[brca_dropped['SUBTYPE']=='LumA']))

#for row in iter(brca['SUBTYPE']):
    #if row == 'Normal':
        #print(count)
        #brca = brca.drop(brca.index[count])
    #count += 1
#brca_dropped.to_csv('../../LGG/BRCA_Dropped.csv', index=False)