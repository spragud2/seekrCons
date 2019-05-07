import pandas as pd
import pickle
import numpy as np
from itertools import product
import pyBigWig as bw
from seekr.fasta_reader import Reader
from collections import defaultdict

import gtfparse

import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing
import argparse

def dSeekrCoords(df,w,s):
    hitDict = {}
    for row in df.iterrows():
        loc,data = row[0],row[1]
        seekrGenomicCoords = []
        for idx,queryHitLocs in enumerate(data):
            if isinstance(queryHitLocs,list):
                continuousCoords =  [np.arange(start=i*s,stop=i*s+w,step=1) for i in queryHitLocs]
                mergedCoords = np.concatenate(continuousCoords)
                flatCoords = np.unique(mergedCoords)
                seekrGenomicCoords.append(flatCoords)
        hitDict[loc] = np.unique(np.concatenate(seekrGenomicCoords))
    return hitDict

def getCons(sCoords,gtf):
    info,coords = sCoords
    print(info,coords)
    infoSplit = info.split('|')
    coordsArray = np.array(coords,dtype=int)
    gene,transcript,strand,start,end = infoSplit[0][1:],infoSplit[1],int(infoSplit[2]),int(infoSplit[3]),int(infoSplit[4])
    if strand == 1:
        coordsArray+=start
    elif strand == -1:
        coordsArray = end - coordsArray
    if not mouseGTF[mouseGTF['gene_id']==gene].empty:
        chrom = mouseGTF[mouseGTF['gene_id']==gene]['seqname'].unique()[0]
        exons = mouseGTF[(mouseGTF['gene_id']==gene)&(mouseGTF['feature']=='exon')]
        if not exons.empty:
            exonCoords = list(zip(exons['start'],exons['end']))
            exonCoords = [np.arange(istart,iend) for istart,iend in exonCoords]
            exonCoords = np.unique(np.concatenate(exonCoords))
        else:
            exonCoords = np.array([])
        allTranscriptCoords = np.arange(start,end)
        setAllTranscriptCoords,setExonCoords,setSeekrCoords = set(allTranscriptCoords),set(exonCoords),set(coordsArray)
        exonsInTranscript = setAllTranscriptCoords.intersection(exonCoords)
        intronsNoSeekrHits = setAllTranscriptCoords - exonsInTranscript - setSeekrCoords
        seekrNoExons = setSeekrCoords - exonsInTranscript
        SEEKR,INTRONS,EXONS = sorted(list(seekrNoExons)),sorted(list(intronsNoSeekrHits)),sorted(list(exonsInTranscript))

        if (SEEKR) and (len(SEEKR)>1):
            seekrPhyloP = np.array([phyloP.values(chrom,start=i,end=i+1) for i in SEEKR])
        if (INTRONS) and (len(INTRONS)>1):
            intronsPhyloP = np.array([phyloP.values(chrom,start=i,end=i+1) for i in INTRONS])
        if (EXONS) and (len(EXONS)>1):
            exonsPhyloP = np.array([phyloP.values(chrom,start=i,end=i+1) for i in EXONS])

        return gene, (seekrPhyloP,intronsPhyloP,exonsPhyloP)

parser = argparse.ArgumentParser()
parser.add_argument("--phyloP",type=str,help='path to phyloP bigwig')
parser.add_argument('--GTF', type=str,help='path to gtf file')
parser.add_argument('--hits', type=str,help='path to pickle file of seekr hits')
parser.add_argument('-n', type=int, help='Number of processors,default = number cpus avail',
                    default=multiprocessing.cpu_count()-1)
parser.add_argument('-w', type=int, help='Window for tile size', default=1000)
parser.add_argument(
    '-s', type=int, help='How many bp to slide tiles', default=100)
args = parser.parse_args()


phyloP = bw.open(args.phyloP)

mouseGTF = gtfparse.read_gtf(args.GTF)
mouseGTF['gene_id'] = [i[:-2] for i in mouseGTF['gene_id']]
hits = pickle.load(open(args.hits,'rb'))

df = pd.DataFrame.from_dict(hits,orient='index')
query_dict = defaultdict(list)
for row in df.iterrows():
    query_dict = defaultdict(list)
    name = row[0]
    data = row[1:][0][0]
    for i in data:
        query_dict[i[1]].append(i[0])
    hits[name] = query_dict.copy()
    query_dict.clear()

parsedDf = pd.DataFrame.from_dict(hits,orient='index')

queryFasta = Reader('./queries.fa')
queryHeaders = queryFasta.get_headers()

queryMap = dict(zip(list(range(13)),queryHeaders))
parsedDf.rename(columns=queryMap,inplace=True)
queryCons = {}

for query in parsedDf:
    seekrGenomicCoords = dSeekrCoords(parsedDf,args.w,args.s)
    with multiprocessing.Pool(args.n) as pool:
        ha = pool.starmap(getCons, product(
            *[[list(seekrGenomicCoords.items())],[mouseGTF]))
        merge = dict(ha)
    queryCons[query] = merge
    print(queryCons)
    1/0


df = pd.DataFrame([seekrAgg,intronsAgg,exonsAgg])
df = df.T
df.columns=['seekr','introns','exons']


sns.set_context('talk')
plt.figure(figsize=(4,7))
sns.boxplot(data=df,showfliers=False)
plt.axhline(y=0,linestyle='--',color='black')
plt.tight_layout()
plt.savefig('./proteincoding_seekrhits_phyloP_new.pdf',bbox_inches='tight')
