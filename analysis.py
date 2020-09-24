import pandas as pd
import seaborn as sns

def length_filter(df,length,cnv_type=False):
    variants={'<1kb':[0,1e3],'1-50kb':[1e3,5e4],'50-500kb':[5e4,5e5],'>500kb':[5e5,1e9]}
    
    if cnv_type:
        return df[(df['length']>variants[length][0])&(df['length']<=variants[length][1])&(df['type']==cnv_type)]
    else:    
        return df[(df['length']>variants[length][0])&(df['length']<=variants[length][1])]

def targets_filter(df,length,cnv_type=False):
    variants={'1':[0,1],'2-3':[1,3],'4-10':[3,10],'>10':[10,200]}
    
    if cnv_type:
        return df[(df['targets']>variants[length][0])&(df['targets']<=variants[length][1])&(df['type']==cnv_type)]
    else:    
        return df[(df['targets']>variants[length][0])&(df['targets']<=variants[length][1])]
    
def targets_filter2(df,length):
    variants={'1':[0,1],'2-3':[1,3],'4-7':[3,7],'>7':[7,200]}
    d=len(df[(df>variants[length][0])&(df<=variants[length][1])])
    return float(d)/len(df)
    
    
def add_alg(ress, alg):
    df=pd.read_csv('exon_int/'+alg,sep='\t',header=None, 
                      names=['chr','start','end','name','0','strand','gene',
                             'c','s','e','type_'+alg,'length_'+alg,'targets'])

    df= df.drop_duplicates(subset=['name'])
    df = df.reset_index(drop=True)
    df[alg]=df['type_'+alg].map(lambda x: 0 if x=='.'else 1)
    ress=pd.merge(ress,df[['name',alg]],on="name",how='left')
    return ress
    
    
def plot_v3(data):
    g = sns.FacetGrid(data,col="Size",sharex=True,size=4)

    g.map(sns.barplot,
          "Algorithm", "% of  predicted CNV",
          palette=sns.cubehelix_palette(17, start=.5, rot=-.75))
        
    axes = np.array(g.axes.flat)
    return plt.gcf(), axes

def set_style():
  
    sns.set(font='serif')
    sns.set_style("whitegrid", {
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "serif"]
    })
    
#####    
    
alg_list=['canoes','clamms','cnmops','cnvkit','codex','conifer','contra','deanncnv',
            'excavator2','exomecopy','exomedepth','exondel',
            'fishingcnv','hmzdelfinder','patterncnv','xhmm']

excavator2=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/excavator2.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
exomedepth=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/exomedepth.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
xhmm=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/xhmm.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
conifer=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/conifer.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
contra=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/contra.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
canoes=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/canoes.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
cnmops=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/cnmops.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
clamms=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/clamms.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
cnvkit=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/cnvkit.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
codex=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/codex.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
deanncnv=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/deanncnv.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
exomecopy=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/exomecopy.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
exondel=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/exondel.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
fishingcnv=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/fishingcnv.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
hmzdelfinder=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/hmzdelfinder.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
patterncnv=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/patterncnv.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])
refcnv=pd.read_csv('/home/gordeeva/./comparasion_study/calling_tools/alg_res/refcnv.sort',sep='\t',header=None, 
                      names=['chr','start','end','type','length','targets'])


ALGORITHMS=[canoes,clamms,cnmops, cnvkit, codex,conifer,contra,deanncnv,
            excavator2,exomecopy,exomedepth,exondel,
            fishingcnv,hmzdelfinder,patterncnv,xhmm]
#####

results=pd.DataFrame(index=['CANOES', 'CLAMMS', 'cn.MOPS', 'CNVkit', 'CODEX', 'CoNIFER',
       'CONTRA', 'DeAnnCNV', 'EXCAVATOR2', 'exomeCopy', 'ExomeDepth',
       'ExonDel', 'FishingCNV', 'HMZDelFinder', 'PatternCNV','XHMM'])
results['CNVs']=[len(i) for i in ALGORITHMS]
results['Loss']=[len(i[i['type']=='Loss']) for i in ALGORITHMS]
results['Gain']=[len(i[i['type']=='Gain']) for i in ALGORITHMS]
results['Targets']=[i['targets'].sum() for i in ALGORITHMS]
results['Length']=[i['length'].sum() /1000. for i in ALGORITHMS]
results['Mean_length,kb']=[i['length'].mean()/1000 for i in ALGORITHMS]
results['Median_length,kb']=[i['length'].median()/1000 for i in ALGORITHMS]
results['Min_length,kb']=[i['length'].min()/1000. for i in ALGORITHMS]
results['Max_length,kb']=[i['length'].max()/1000. for i in ALGORITHMS]
results['Mean_targets']=[i['targets'].mean() for i in ALGORITHMS]
results['Median_targets']=[i['targets'].median() for i in ALGORITHMS]
results['Min_targets']=[i['targets'].min() for i in ALGORITHMS]
results['Max_targets']=[i['targets'].max() for i in ALGORITHMS]

for l in ['<1kb','1-50kb','50-500kb','>500kb']:
    results[l]=[len(length_filter(i,l)) for i in ALGORITHMS]
    
    for t in ['Loss', 'Gain']:
        results[l+'_'+t]=[len(length_filter(i,l,t)) for i in ALGORITHMS] 
        
for l in ['1','2-3','4-10','>10']:
    results[l]=[len(targets_filter(i,l)) for i in ALGORITHMS]

results[['CNVs','Loss','Gain',
         'Targets','Length','Mean_targets',
         'Mean_length,kb','Min_length,kb','Max_length,kb']].to_csv('/home/gordeeva/./comparasion_study/results/CNV_description.csv',sep="\t", index=False)

#Plot size of predicted CNV calls

col_order=["PatternCNV", "CONTRA", "CLAMMS",'HMZDelFinder','FishingCNV','ExomeDepth','CODEX','CoNIFER',
                  'cn.MOPS','ExonDel','XHMM','DeAnnCNV','EXCAVATOR2','CANOES','CNVkit','exomeCopy']

res=pd.DataFrame(columns=['Algorithm','Size','% of  predicted CNV'])
for i in col_order:
    
    for l in ['<1kb','1-50kb','50-500kb','>500kb']:       
        res=res.append({'Algorithm':i,
                        'Size':l,
                        '% of  predicted CNV':results.loc[i,l]/float(results.loc[i,'CNVs'])},
                       ignore_index=True)


fig,axes=plot_v3(res)
for i, ax in enumerate(axes):
    for tick in ax.get_xticklabels():
            tick.set_rotation(90)
            
    ax.set_xlabel("")
    labels = ['<1kb','1-50kb','50-500kb','>500kb']
    ax.set_title(labels[i],fontweight='bold',size=12)
    set_style()
    sns.set_context("paper",font_scale=1.3)
    
fig.savefig('/home/gordeeva/./comparasion_study/result/length_distributiond.png', format='png', dpi=500)

res=pd.DataFrame(columns=['Algorithm','Size','% of  predicted CNV'])
for i in col_order:
    
    for l in ['1','2-3','4-10','>10']:
        res=res.append({'Algorithm':i,
                        'Size':l,
                        '% of  predicted CNV':results.loc[i,l]/float(results.loc[i,'CNVs'])},
                       ignore_index=True)
        
fig,axes= plot_v3(res)
for i, ax in enumerate(axes):
    for tick in ax.get_xticklabels():
            tick.set_rotation(90)
            
    ax.set_xlabel("")
    labels = ['1','2-3','4-10','>10']
    ax.set_title(labels[i],fontweight='bold',size=12)
    
set_style()
sns.set_context("paper",font_scale=1.3)
fig.savefig('/home/gordeeva/./comparasion_study/result/target_distributiond.png', format='png', dpi=500)

#########
#Comparasion at exon level
exons=pd.read_csv("/home/gordeeva/./human_genome/gencode_exons",sep='\t',header=None,
                  names=['chr','start','end','name','0','strand','gene'])

for i in alg_list:
    exons=add_alg(exons,str(i))
exons['sum']=target[alg_list].sum(axis=1)

df=pd.read_csv('/home/gordeeva/./comparasion_study/standard.bed',sep='\t',header=None,
              names=['c','s','e','name','gene','st'])
exons=pd.merge(exons,df[['name','st']], on='name',how='left')

df=pd.read_csv('/home/gordeeva/./comparasion_study/exon_int/target',sep='\t',header=None,
              names=['chr','start','end','name','0','strand','gene','c','s','e','target'])
df= df.drop_duplicates(subset=['name']).reset_index(drop=True)
df['target']=df['target'].map(lambda x: 0 if x=='.' else 1)

exons=pd.merge(exons,df[['name','target']],on="name",how='left')
exons=exons[(exons['chr']!='X')&(exons['chr']!='Y')&(exons['target']>0)]


inter=exons[alg_list+['sum']]
alg_comb=[]

for i,row in inter[alg_list].iterrows():
    alg_comb.append('_'.join(set(row.index[row.values>0].values)))
inter['alg']=alg_comb

l=pd.DataFrame(inter['alg'].value_counts()[1:])
l['count']=l.index.map(lambda x: len(x.split('_')))
l.sort_values([ 'count','alg'], ascending=[True, False])
l.to_csv('/home/gordeeva/./comparasion_study/result/call_intersect.csv', sep="\t")


stat=pd.DataFrame(np.zeros((len(alg_list),len(alg_list))))

for i in range (len(alg_list)):
    stat.iloc[i,i]=1
    for k in range (i+1, len(alg_list)):
       
        c=inter[(inter['alg'].str.contains(alg_list[i])) & (inter['alg'].str.contains(alg_list[k]))].shape[0]
        stat.iloc[i,k]=float(c)/(inter[(inter['alg'].str.contains(alg_list[i]))].shape[0]+inter[(inter['alg'].str.contains(alg_list[k]))].shape[0] -c)
        
    
stat.columns=alg_list
stat.index=alg_list
(stat.T).to_csv('/home/gordeeva/./comparasion_study/result/paired_intersect.csv', sep="\t")


resu=pd.DataFrame(col_order, columns=['alg'])

for l in ['1','2-3','4-7','>7']:
    resu[l]=[targets_filter(inter[inter[i]==1]['sum'],l) for i in col_order]

resu['0']=resu['1']
resu['1-2']=resu['1']+resu['2-3']
resu['3-6']=resu['1-2']+resu['4-7']
resu['>6']=resu['3-6']+resu['>7']


ax=sns.barplot(x="alg", y=">6", data=resu, 
            color=sns.cubehelix_palette(50, start=.7, rot=-.75, dark=.2)[45],label=">6")
ax=sns.barplot(x="alg", y="3-6", data=resu, 
            color=sns.cubehelix_palette(50, start=.7, rot=-.75, dark=.2)[30],label="3-6")
ax=sns.barplot(x="alg", y="1-2", data=resu, 
            color=sns.cubehelix_palette(50, start=.7, rot=-.75, dark=.2)[15],label="1-2")
ax=sns.barplot(x="alg", y="0", data=resu,
            color=sns.cubehelix_palette(50, start=.7, rot=-.75, dark=.2)[0],label="0")
ax.set_ylabel("")
ax.set_xlabel("")
ax.set_title('CNV-exon confirmed by other tools')
ax.legend(bbox_to_anchor=(1.05, 1.05))
handles, labels = ax.get_legend_handles_labels()
for tick in ax.get_xticklabels():
            tick.set_rotation(85)
        
ax.savefig('/home/gordeeva/./comparasion_study/result/confirmed_cnv.png', format='png', dpi=500)


####Precision-recall
from sklearn import metrics
exons2=exons[alg_list+['st']][(exons.st.notnull())]
exons22.columns=['CANOES', 'CLAMMS', 'cn.MOPS', 'CNVkit', 'CODEX', 'CoNIFER',
       'CONTRA', 'DeAnnCNV', 'EXCAVATOR2', 'exomeCopy', 'ExomeDepth',
       'ExonDel', 'FishingCNV', 'HMZDelFinder', 'PatternCNV','XHMM', 'st'] 
results=pd.DataFrame(columns=['algorithm','precision','recall'])
for i in  col_order: 
    
 
    results=results.append({'algorithm':i,
                            'precision': metrics.precision_score(exons2['st'],exons2[i]),
                            'recall': metrics.recall_score(exons2['st'],exons2[i])} ,ignore_index=True)
results.to_csv('/home/gordeeva/./comparasion_study/result/precision-recall.csv', sep="\t", index=False)


pal=[[0.5488483550919352, 0.7555792405897566, 0.5823941160739373],
     [0.3256286372570367, 0.5824294714811111, 0.551260440725878],
     [0.24073795863422207, 0.34003177164489373, 0.4786370052053605],
     [0.17250549177124488, 0.11951843162770595, 0.24320155229883056]]*4
mark=["v","v","v","v","o","o","o","o","s","s","s","s","^","^","^","^"] 

g=sns.lmplot(x="recall", y="precision", hue="algorithm", 
                       data=results,fit_reg=False,palette=pal,
                      size=8,scatter_kws={"s": 200},legend=True, markers=mark)

set_style()
sns.set_context("paper",font_scale=1.5)

g.savefig('/home/gordeeva/./comparasion_study/result/precision-recall.png', format='png', dpi=500)

