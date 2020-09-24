def add_standard(df,x):


    table_col=['chr','start','end','exon','0','strand','gene','chr1','start1','end1',str(x)]
    rt = pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets/exon_level',str(x)),
                    sep="\t",header=None, names=table_col)

    rt=rt.drop_duplicates(subset='exon', keep="last")
    df=pd.merge(df,rt[['exon',str(x)]],on = 'exon')

    return df

sets=['conrad2006','pilot','metasv','pinto2007','wang2007','conrad2010','hapmap',
    'mccarroll2006','pacbio','redon2006','cooper2008','phase3',
    'lumpy','mccarroll2008','phase1','svclassify']

df=pd.read_csv("/home/gordeeva/./human_genome/gencode_exons.bed",
               sep="\t",header=None,
               names=['chr','start','end','exon','0','strand','gene'])

for st in sets:
    df=add_standard(df,st)
    
df=df.replace('.',np.nan)
df[st]=df[sets].astype(float)
df.to_csv(os.path.join(DATA_FOLDER1,'exon-CNVset'),
                        sep="\t",header=False,index=False)  

population=pd.read_csv(os.path.join(DATA_FOLDER1,'population_cnv'),sep="\t", header=None,
                      names=['chr1','start','end','exon','0','strand','gene','chr2','start2','end2','European']) 

population=population.drop_duplicates(subset='exon')
population=population[population['European']!='.']
x=population['European'].astype(float)
x[x>1000]=1000

import scipy
a,b,sd,res=scipy.stats.beta.fit(x)
from scipy import special
from scipy.stats import beta



df['pro']=(df[st]==1).sum(axis=1)   
df['con']=(df[st]==0).sum(axis=1) 
df['low_pacbio']= df['pacbio'].map(lambda x: 1 if x == 0.1 else 0)
df['low_metasv']= df['metasv'].map(lambda x: 1 if x == 0.1 else 0)
df['NA']=df[st].isnull().sum(axis=1)


df['exon_rate']=map(lambda x,y,w,z: np.nan if x+y+z+w<1 else special.betaincinv(x+a+w*358/6546.+z*2256/24306,y+b,5/21.), 
                    df['pro'],df['con'],df['low_pacbio'],df['low_metasv'])

threshold=df['exon_rate'].max()*0.5

df['our_standard']= map(lambda x: 1 if (x>=threshold) else (0  if (x<threshold) else '-'),df['exon_rate'])
df[['chr','start','end','exon','gene','our_standard']].to_csv('/home/gordeeva/./comparasion_study/standard.bed',
                                                                                  sep="\t", index=False, header=None)
