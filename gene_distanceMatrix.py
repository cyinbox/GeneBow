# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 11:55:34 2018

@author: Changchuan Yin
         University of Illinois at Chicago
         Email cyinbox@gmail.com | cyin1@uic.edu
"""
import gene_ncbi as ncbi

#------------------------------------------------------------------------------
# Compute the pair-wise distances for phylogenetic tree analysis
#------------------------------------------------------------------------------
def removeNonATCG(seq):
 seq = str(seq)
 for ch in seq:
   if ch not in ['A','T','C','G']:
      seq = seq.replace(ch,'')
 
 return seq

#------------------------------------------------------------------------------
# 1. The list of GenBank　Ids and gene or genome names.
# The gene/genome name should be a single word, and unique
genes = {'V00662':' Human',
            'D38116':'Common_chimpanzee',
            'D38113':'Pigmy_chimpanzee',
            'D38114':'Gorilla',
            'NC_002764':'Macaca_thibetana',
            'D38115':'Bornean_Orangutan',
            'NC_002083':'Sumatran_Orangutan',
            'X99256':'Gibbon',
            'AY863426':'Vervet_Monkey',
            'AB245427':'Deer',
            'NC_020679':'Antelope',
            'AP003424':'Giraffe',
            'EU177871':'Plains_bison',
            'NC_010302':'Manatee',
            'Y18001':'Baboon',
            'X79547':'Horse',
            'X97337':'Donkey',  
            'X97336':'Indian_rhinoceros',
            'Y07726':'White_rhinoceros',
            'AJ428947':'Tapir', # %close to rhino
            'X63726':'Harbor_seal',
            'X72004':'Gray_seal',
            'KJ868120':'Gray_Kangaroo',
            'U20753':'Cat',
            'U96639':'Dog',
            'X61145':'Fin_whale',
            'X72204':'Blue_whale',
            'EU557091':'Dolphin',
            'GU475464':'Sea_lion',
            'AJ428576':'Walrus',
            'NC_016470':'Cougar',
            'V00654':'Cow',
            'AF010406':'Sheep',
            'AF533441':'Goat',
            'AY488491':'Buffalo',
            'AJ002189':'Pig',
            'AJ010957':'Hippopotamus',
            'EF212882':'Giant_panda',
            'DQ402478':'Black_bear',
            'AF303111':'Polar_bear',
            'AF303110':'Brown_bear',
            'EU442884':'Wolf',
            'EF551002':'Leonard',
            'EF551003':'Tiger',
            'KF907306':'Lion',
            'AM181037':'Fox',
            'NC_008093':'Coyote',
            'NC_015108':'Beaver',
            'NC_018367':'Marmot ',
            'NC_016008':'Pangolin', 
            'AF061340':'Neotropical_fruit_bat',
            'AJ224821':'African_elephant',
            'Y18475':'Aardvark',
            'Y11832':'Armadillo',
            'NC_009126':'Raccoon',
            'AJ001588':'Rabbit',
            'AJ537415':'Pika', #%close to Rabbit
            'AJ222767':'Guinea_pig',
            'AJ001562':'Fat_dormouse',
            'X14848':'Rat',
            'V00711':'Mouse',
            'AJ238588':'Squirrel',
            'X88898':'Hedgehog', 
            'NC_019626':'Shrew_Gymnure',
            'Z29573':'Opossum',
            'X83427':'Platypus',
            'Y10524':'Wallaroo',
            'MH142215': 'Anteater',
            'AY960979':'Sloth',
            'NC_002391':'Mole'}

#----------------------------------------------------------------------------
# 2. Make sure the gene/genome names are unique
geneIds = []
geneNames = []

for geneId, geneName in genes.items():
  if geneName not in geneNames:
   geneIds.append(geneId)
   geneNames.append(geneName)
  else:
   print('Not unique name:',geneName)

# 3. Retrieve Genbank sequence data and write sequences to a fasta file
headers = []
sequences = []
m = len(geneIds)

from pathlib import Path

fileName ='genes.fasta'
geneFile = Path(fileName)

if geneFile.is_file():
    print('Gene fasta file exists')
else:
 for i in range(0,m):
  geneId = geneIds[i]
  geneName = geneNames[i]
  [header,seq] = ncbi.getGenBank(geneId) 
  headers.append(geneName)
  seq = removeNonATCG(seq) # Remove non-ATCG nucleotides in a sequence
  sequences.append(seq)
  
  ncbi.writeFastas(headers,sequences,fileName)
  print('Completed in getting genomes')

#-----------------------------------------------------------------------------
# 4. Compute pair-wise distances of genes/genomes
file = open('distances.txt', 'w')

# MATLAB program phylogeneticTreeByDistances.m uses this distance file to construct phylogenetic tree
# https://github.com/cyinbox/NGS/blob/master/phylogeneticTreeByDistances.m

from Bio import SeqIO
import gene_CGR as cgr

records = list(SeqIO.parse("genes.fasta", "fasta"))
m = len(records)

for i in range(0,m):
   geneNameA = records[i].id
   seqA = records[i].seq
   
   for j in range(i+1,m):
    geneNameB = records[j].id
    seqB = records[j].seq
    
    dist = cgr.getDist_CGR(seqA,seqB) # or any distance function
    strDist = geneNameA+','+geneNameB+','+str(dist)+'\n'
    print(strDist)
    file.write(strDist)
    
print('Distance matrix computation completed') 
file.close()

