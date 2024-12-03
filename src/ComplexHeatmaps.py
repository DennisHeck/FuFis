import os,sys
import matplotlib.pylab as plt
import pickle
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
plt.rcParams['font.family']='sans serif'
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype']=42
# sys.path.append(os.path.expanduser("~/Projects/Github/PyComplexHeatmap/"))
import PyComplexHeatmap
from PyComplexHeatmap import *

data = pd.DataFrame([['TCGA-05-4384-01', 'KRAS', 1, 2, 3],
                    ['TCGA-05-4384-01', 	'HRAS', 3, 2, 1],
                    ['TCGA-05-4384-01', 	"KRAS", 1, 1, 1]], columns=['Sample', "Gene", 'A', 'B', 'C'])
row_vc=data.groupby('Gene').apply(lambda x:x.loc[:,['A', 'B', 'C']].sum())

data=pd.read_csv("../data/tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt",sep='\t',index_col=0)
data=data.iloc[:,:-1]
data=data.stack().reset_index()
data.columns=['SampleID','Genes','Variants']
data.Variants.replace({'  ':np.nan},inplace=True)

unique_variants=[]
for var in data.Variants.dropna().unique():
    for v1 in var.split(';'):
        v1=v1.strip()
        if v1=='':
            continue
        if v1 not in unique_variants:
            unique_variants.append(v1)
print(unique_variants)
for var in unique_variants:
    data[var]=data.Variants.fillna('').apply(lambda x:1 if var in x else 0)

cols=['AMP','HOMDEL','MUT']
colors=["red","blue","#008000"]

# calculate genes (row) mutation frequencies.
row_vc=data.groupby('Genes').apply(lambda x:x.loc[:,cols].sum())
# calculate samples (cols) mutation frequencies.
col_vc=data.groupby('SampleID').apply(lambda x:x.loc[:,cols].sum())

#Samples with variants at KRAS
kras_samples=data.loc[(data.Genes=='KRAS') & (data.loc[:,cols].sum(axis=1)>0)].SampleID.unique().tolist()
df_col_split=pd.DataFrame(index=data.SampleID.unique(),data={'KRAS':['No KRAS Var']*data.SampleID.nunique()})
df_col_split.loc[kras_samples,'KRAS']='KRAS Var'

top_annotation=HeatmapAnnotation(axis=1,
                                KRAS=anno_simple(df_col_split.KRAS,add_text=True,height=6),
                                Col=anno_barplot(col_vc,colors=colors,legend=False,height=10,linewidth=0.1),
                                verbose=0, label_side='left', label_kws={'horizontalalignment': 'right','visible':False})
right_annotation = HeatmapAnnotation(axis=0,orientation='right',
                                Row=anno_barplot(row_vc,colors=colors,legend=False,height=10,linewidth=0.1),
                                verbose=0, label_side='top', label_kws={'horizontalalignment': 'left','rotation':45,'visible':False})

plt.figure(figsize=(12,8))
op=oncoPrintPlotter(data=data,y='Genes',x='SampleID',
                    values=cols,colors=colors,subplot_gap=3,label='Alteration',
                    top_annotation=top_annotation,right_annotation=right_annotation,
                    col_split=df_col_split.KRAS,col_split_order=['KRAS Var','No KRAS Var'],col_split_gap=3,
                    legend_hpad=0,show_rownames=True,show_colnames=False) #xticklabels_kws={'labelsize':3}
plt.savefig("/Users/dennis/Dev/CADlinc/Plots/TEST.pdf",bbox_inches='tight')
plt.show()

