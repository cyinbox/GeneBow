# -*- coding: utf-8 -*-
'''
# Program detecting antibiotic resistance genes
#
# Changchuan Yin, Ph.D.
# Dept. of Mathematics, Statistics and Computer Science
# University of Illinois at Chicago
# Chicago, IL 60607
# USA
#
# Email cyin1@uic.edu
# Last update 06/16/2017
#
# Citation
# Yin, C. (2017) GeneBow: a pipeline for detecting antibiotic resistance genes.
#
'''
# Program to construct ARG sequences 
# Inputs:ar_genes.tab contains ARG protein accessId,resistance_profile.tab is for resistance profile
# output: ar_genes_cds.fasta: ARG genes and their DNA sequences
import pandas as pd
import GenePod as gene

fastaFileName='./ARGData/ar_genes_cds.fasta'
handle = open(fastaFileName, 'w')

# resistance_profile.tab is downloaded from https://ardb.cbcb.umd.edu/
df = pd.read_csv('./ARGData/resistance_profile.tab',sep='\t',header=None)
# df's first index is column, second index is row number of tab data file
#print(df[0][0])
#print(df[0][1])
numRow=df[0].count()

antibioticProfiles=dict() #dictionary of all antibiotic profile in ARDB

for row in range(1,numRow):
  profileId=df[0][row]
  profile=df[1][row]

  if profileId in antibioticProfiles:
      antibioticProfiles[profileId].append(profile)
  else:
      antibioticProfiles[profileId]=[profile]

#print('Profile dictionary:',antibioticProfiles)

#file ar_genes.tab: downloaded from https://ardb.cbcb.umd.edu/
df = pd.read_csv('./ARGData/ar_genes.tab',sep='\t',header=None)
numRow=df[0].count()

print(df[0][2])
print(df[0].count())#row 
numRow=df[0].count()

proteinIds=[]
profiles=[]

#for row in range(10641,numRow):
IdProfiles=dict()
for row in range(0,numRow):
  proteinId=df[0][row]
  profile=df[1][row]
  proteinIds.append(proteinId)
  IdProfiles[proteinId]=profile #Store id and profile relationship

''' For dev test
print('Test2:',proteinIds)
print('IdProfile',IdProfiles)
proteinId='AAC44793'
profileNameByProteinId=IdProfiles[proteinId]
print('Test2',profileNameByProteinId)
antibioticProfileByProfileId=antibioticProfiles[profileNameByProteinId]
print('Test3',antibioticProfileByProfileId)
print('proteinIds:',proteinIds)
'''
#proteinIds = ['ABS70977','CAD38267','BAB71949XX','CAG44462']#,'AAS64386','ABS70977','CAD38267','NP_334681','YP_905232','ABR04609','ABL89492','ABG06386'];

cnt=0
for proteinAccessId in proteinIds:
 try:
  print('\n'+'starting to process protein:',proteinAccessId)
 
  [genBankId,cdsDirection,startPos,endPos,desc] = gene.getGenBankByProteinId(proteinAccessId)
  profileNameByProteinId=IdProfiles[proteinAccessId]
  antibioticProfileByProfileId=antibioticProfiles[profileNameByProteinId]
  #print('antibioticProfile:',antibioticProfileByProfileId)
  
  print('protein access Id:',proteinAccessId)
  print('gene access Id:',genBankId)
  print('cds direction:',cdsDirection)
  print('cds starts:',startPos)
  print('cds ends:',endPos)
  print('descs:',desc)
  print('antibiotic resistance:',antibioticProfileByProfileId[0])
 
  cds=gene.getGenBankCDS(genBankId,cdsDirection,startPos,endPos) 
  header=proteinAccessId+'|'+genBankId+'|'+str(cdsDirection)+'|'+str(startPos)+'|'+str(endPos)+'|'+desc+'|'+antibioticProfileByProfileId[0]
  print('Header:',header)
  print('seq:',cds)
  gene.writeFastaFileOneRecord(header,cds,handle)
  cnt=cnt+1
 
 except Exception as err:
    print('Error:',err)
    continue  

handle.close()
print('Number of ARGs:',cnt)
print('All protein Ids have been processed.') 


 
