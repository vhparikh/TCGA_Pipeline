import pandas as pd

# reading data
print('Running tool---')
# get the data
lgg = pd.read_csv('../../Cancers/LGG.csv')
genes = lgg.columns
subtypes = lgg['SUBTYPE']
print(len(subtypes))
print(len(lgg[lgg['SUBTYPE'] == 'IDHmut-codel']))
print(len(lgg[lgg['SUBTYPE'] == 'IDHmut-non-codel']))
print(len(lgg[lgg['SUBTYPE'] == 'IDHwt']))
print(subtypes.unique())
#print(len(genes))
count = 0

#print(brca.shape)
#print(len(brca[brca['SUBTYPE']=='LumA']))
#brca_dropped = gbm[gbm.SUBTYPE != 'Normal']
#print(brca_dropped.shape)
#print(len(brca_dropped[brca_dropped['SUBTYPE']=='LumA']))

#for row in iter(brca['SUBTYPE']):
    #if row == 'Normal':
        #print(count)
        #brca = brca.drop(brca.index[count])
    #count += 1
#brca_dropped.to_csv('../../LGG/BRCA_Dropped.csv', index=False)