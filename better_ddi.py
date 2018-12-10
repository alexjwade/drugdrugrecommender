import pandas

df = pandas.read_csv('smol_twosides.tsv', sep='\t', header=0)
df.drop(['drug1', 'drug2', 'event_umls_id', 'event_name', 'proportional_reporting_ratio', 'pvalue', 'drug1_prr', 'drug2_prr', 'observed', 'expected'], axis=1)
df = df[df.confidence < 5]
print(df)

print df.confidence