#Convert GVF to BED

#grep "NA12878" file.gvf | grep -v "#" - | awk '{  if ($3=="copy_number_loss") print $1"\t"$4-1"\t"$5"\t1"; else if #($3=="copy_number_gain") print $1"\t"$4-1"\t"$5"\t2";}' - |sort -k 1,1 -k2,2n - | bedtools merge -i - -nms > file

##############

import pandas as pd
import os

###Hapmap  
df= pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','hm3_cnv_submission.txt'), sep ='\t')

df=df[['chr','start','end','NA12878']].dropna()
df.NA12878=df.NA12878.map(lambda x: 0 if x == 2 else 1) 
df.NA12878=df.NA12878.astype(int)

df.to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','hapmap_hg18'), sep ='\t',header=False,index=False)



###McCarroll
genotypes=pd.read_excel(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','41588_2008_BFng238_MOESM25_ESM.xls'),sheet_name='CEU')
location=pd.read_excel(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','41588_2008_BFng238_MOESM24_ESM.xls'),sheet_name='cnp1319_summary.txt',skiprows=6)

df=pd.merge(location[['CNP_id','chr','start_hg18','end_hg18']],
            genotypes[['CNP_id','NA12878']],on='CNP_id',how='outer')
df=df.dropna()
df.NA12878=df.NA12878.map(lambda x: 0 if x == 2 else 1)

df.to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','mccarroll_hg18'), sep ='\t',header=False,index=False)

##liftover hapmap and mccarol to hg19 (https://genome.ucsc.edu/cgi-bin/hgLiftOver)



####Conrad
genotypes=pd.read_excel(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','41586_2010_BFnature08516_MOESM10_ESM.xls'),sheet_name='CEU')
location=pd.read_excel(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','41586_2010_BFnature08516_MOESM10_ESM.xls'),sheet_name='Genotype Map')

df=pd.merge(location[['CNV','chr','start','end']], genotypes[['CNV','NA12878']],on='CNV',how='outer')
df=df.dropna()
df.NA12878=df.NA12878.map(lambda x: 0 if x == 2 else 1).astype(int)
df.chr[df.chr ==23] = 'X'

df.to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','conrad-table.bed'), sep ='\t',header=False,index=False)


###Pilot
df=pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets', 'MasterValidation.Pilot2.all.leftmost.061510.txt'),sep='\t')

df=df[['CHR','START_CI_BKPT','END_CI_BKPT','SV_TYPE','SAMPLES','VALIDATION_STATUS']]
df=df[(df.SV_TYPE !='INV') & (df.SV_TYPE != 'INS') 
    & (df.VALIDATION_STATUS !='unvalidated')] 

df=df[df.SAMPLES.str.contains(r'NA12878')].sort_values(by=['CHR','START_CI_BKPT','END_CI_BKPT']) 
df.NA12878=df['VALIDATION_STATUS'].map(lambda x: 1 if x == 'validated' else 0).astype(int)

df[['CHR','START_CI_BKPT','END_CI_BKPT','NA12878']].to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets',
                                                                        'pilot_hg19'), sep='\t',header=False,index=False)


###1KG - Phase3
df=pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets', 
                           'ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf'),sep='\t')
df=df[['NA12878','#CHROM','POS','REF','ALT','INFO']]
df['POS']=r['POS']-1 
df['TYPE']=map( lambda i: i.split('SVTYPE=')[1].split(';')[0], r.INFO)


df2=df[(df.TYPE == 'DEL_ALU' )|(df.TYPE == 'DEL_LINE1' )|(df.TYPE == 'DEL_SVA')|(df.TYPE == 'DEL_HERV')] 
df2['END']= map( lambda i: i.split(';END=')[1].split(';')[0], df2.INFO)
#ALT: A,T,C,G,CN0
df2['SV']= map( lambda i: 0 if i == '0|0' else 1 , df2['NA12878'])

df1=df[(df.TYPE == 'CNV' )|(df.TYPE == 'DEL' )|(df.TYPE == 'DUP')]
df1['END']= map( lambda i: i.split(';END=')[1].split(';')[0], df1.INFO) 

# alternative allele - nucleotide
rAltNuc=df1[(df1.ALT == 'T')|(df1.ALT == 'A')|(df1.ALT == 'G')|(df1.ALT == 'C')]
rAltNuc['SV']= map( lambda i: 0 if i == '0|0' else 1 , rAltNuc['NA12878']) #all deletions

#alternative alllele - copy number
rAltCN=df1[-((df1.ALT == 'T')|(df1.ALT == 'A')|(df1.ALT == 'G')|(df1.ALT == 'C'))] 
rAltCN2=rAltCN[rAltCN['NA12878']=='0|0'] # 2 референсные копии
rAltCN2['SV']= map( lambda i: 0 , rAltCN2['NA12878'])

#CHROM  POS ID  REF ALT INFO
#1   25688926    DUP_gs_CNV_1_25688926_25700415  G   <CN0>,<CN2> .   PASS    SVTYPE=CNV;END=25700415;CS=DUP_gs
#0|0:  copy number 2
#0|1: copy number 1 (1 + 0)
#0|2: copy number 3 (1+ 2)
#1|2: copy number 2 (0 + 2)
#2|2: copy number 4 (2 + 2)

rAltCNN=rAltCN[rAltCN['NA12878']!='0|0']
rAltCNN['AL1']= map( lambda GT,ALT :int((ALT.split(',')[int(GT.split('|')[0])-1]).split('CN')[1].split('>')[0]) if GT.split('|')[0] != '0' else 1,
                   rAltCNN['NA12878'],rAltCNN.ALT)

rAltCNN['AL2']= map( lambda GT,ALT :int((ALT.split(',')[int(GT.split('|')[1])-1]).split('CN')[1].split('>')[0]) if GT.split('|')[1] != '0' else 1,
                   rAltCNN['NA12878'],rAltCNN.ALT)

rAltCNN['SV']=map(lambda x,y: 0 if x+y==2 else 1, rAltCNN.AL1,rAltCNN.AL2)

col=['#CHROM','POS','END','SV']
GS= pd.concat([df2[col], rAltNuc[col],
               rAltCN2[col], rAltCNN[col]],ignore_index=True).sort_values(by=['#CHROM','POS','END']) 

GS.to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets', '1KG'),sep='\t',index=False, header=False)



## conrad and pilot, phase3
#awk ' { if ($4 == "1") print $0}' $set | sort -k1,1 -k2,2n - | bedtools merge -i -  | awk '{print $0"\t1"}' - > $set_cnv
#awk ' { if ($4 == "0") print $0}' $set | sort -k1,1 -k2,2n - | bedtools merge -i -  | awk '{print $0"\t0"}' - > $set_noncnv
#bedtools multiinter -i $set_cnv  $set_noncnv  -header -names 1  0 | awk '{ if ($4 == 1) print $1"\t"$2"\t"$3"\t"$5; else print  $1"\t"$2"\t"$3"\t"1}'>  $set_hg18
#liftover  conrad and pilot to hg19 (https://genome.ucsc.edu/cgi-bin/hgLiftOver)

##lumpy
awk '{print $1"\t"$2"\t"$6"\t1"}' /home/gordeeva/./comparasion_study/cnv_sets/3717462611446476_add4.bedpe > lumpy

###pacbio
df=pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','NA12878.sorted.vcf'),sep="\t")

df['END']= map(lambda i: df.iloc[i,7].split(';END=')[1].split(';')[0], range(len(df)))
df['POS']=df['POS']-1
df=df.replace(['PASS', 'lt3'], [1,0.1]) #Quality PASSdetected by at least 3 algoritnms

df[['#CHROM','POS','END','FILTER']].to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','pacbio'),
                                          sep="\t",header=False,index=False)

### metasv                                          
df=pd.read_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','NA12878_svs.vcf'),sep="\t")

df['END']= map( lambda i: df.iloc[i,7].split('END=')[1].split(';')[0], range(len(df)))
df['TYPE']=map( lambda i: df.iloc[i,7].split('SVTYPE=')[1].split(';')[0], range(len(df)))
df=df.replace(['PASS', 'LowQual'], [1,0.1])#PASS - two or more algorithms in MetaSV confirm
df['POS']=df['POS']-1

df[['#CHROM','POS','END','FILTER']].to_csv(os.path.join('/home/gordeeva/./comparasion_study/cnv_sets','metasv'),
                                           sep="\t",header=False,index=False)