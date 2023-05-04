import sys
import igraph as ig
from igraph import Graph, mean
from igraph import VertexSeq
from igraph import EdgeSeq
from igraph import summary
from igraph import plot
from igraph import GraphBase
from igraph import VertexClustering
from igraph import clustering

import time
import re
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch


def sanitize(string):
	string=re.sub(r'[^\w]','',string)	
	string=re.sub(r'[\d_]','',string)		
	string=string.lower()
	return string

def listtostring(list, sep):
	string=""
	if len(list)>0:
		string+=list[0]#.decode('utf-8')
	for item in list[1:]:
		string+=sep+item#.decode('utf-8')
	return string

def readFasta(file, keepFile="NA"):
	import re
	sequences=[]
	kp=set()
	with open(file,'r') as input:
		if keepFile!="NA":
			with open(keepFile,'r') as keep:
				for line in keep:
					kp.add(line[:-1])
		#print kp
		for line in input:
			if line[0]!=">":
				if len(kp)>0:
					if sequences[-1][0] in kp:
						#print "keep"
						sequences[-1][1]+=sanitize(line)
					else:
						#print sequences[-1], "not on keep list"
						continue
				else:
					sequences[-1][1]+=sanitize(line)
			else:
				if len(sequences)>0 and sequences[-1][1]=="":
					del sequences[-1]
				sequences.append( [re.sub('[^0-9]', '',line[1:]),""] )
	print(len(sequences), sequences[0])
	return sequences

def levenshteinDistance(s1,s2):
		
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return distances[-1]

def levenshteinDistanceThresh(s1,s2,k):
        from Bio.Align import substitution_matrices as matlist
        matrix = matlist.load('BLOSUM62')
        if len(s1) > len(s2):
                s1,s2 = s2,s1
        distances = range(len(s1) + 1)
        for index2,char2 in enumerate(s2):
                newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if abs(index1-index2) > k:
                newDistances.append(sys.maxsize)
            elif char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],distances[index1+1],newDistances[-1])))	
        distances = newDistances

        if min(distances)>k:
            return sys.maxsize		
        return distances[-1]

def makeGraphEfficiently(species, repeatmatrix, filename, sparse=False):
	
	G = Graph(directed=True)
	G.add_vertices(len(species))
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))
	for row in	 range(len(species)):
		if sparse==True:
			import scipy
			m=min(list(scipy.sparse.find(repeatmatrix.tocsc()[:,row])[2]) + list(scipy.sparse.find(repeatmatrix.tocsc()[row,:])[2]))
			edgeTo = getSparseMinED(repeatmatrix, row, m)
		else:
			m=min(l for l in repeatmatrix[row] if l > 0)
			edgeTo = getMinED(repeatmatrix, row, m)
		for e in edgeTo:
			G.add_edge(e[0], e[1], weight=m)
		
	G.simplify(combine_edges=min)
	summary(G)
	G.write(filename+"DWRNgraphEff.gml","gml")
	return G

def makeGraph(species, RepeatMatrix, filename):

	G_new = Graph(directed=False)
	G_new.add_vertices(len(species))
	
	G = Graph()
	G.add_vertices(len(species))
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))

	
	its=len(species)*len(species)
	it=0
	for i in range(len(species)):
		minED = min(l for l in RepeatMatrix[i] if l > 0)
		
		for j in range(len(species)):
			if (RepeatMatrix[i][j] == minED):
				G_new.add_edge(i,j,weight=RepeatMatrix[i][j])
			if i>j and RepeatMatrix[i][j]<2000000:
				G.add_edge(i,j,weight=RepeatMatrix[i][j])
			it+=1
		print(round((it*1.0)/its,3))
	summary(G_new)
	summary(G)
	
	
	G_new.write(filename+"DWN-base.gml","gml")
	G.write(filename+"Thresh.gml","gml")
 
	return G

def getMinED(repeatmatrix, row, minVal):
	
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
		if repeatmatrix[row][col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row][col] < minVal and repeatmatrix[row][col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row][col] 
	
	return ret
	
def getSparseMinED(repeatmatrix, row, minVal):
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
		if repeatmatrix[row,col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row,col] < minVal and repeatmatrix[row,col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row,col] 
	
	return ret

def createMinMatrix(species, useBoundedLD=False):
    from Bio.Align import substitution_matrices as matlist
    from Bio import pairwise2
    matrix = matlist.load('BLOSUM62')
    skippedCells=0	
    repeatmatrix=[]
    its=(len(species)*len(species))/2
    it=len(species)
    for row in range(len(species)):
        repeatmatrix.append([])
        repeatmatrix[row].extend([sys.maxsize]*len(species))
    for row in range(len(species)):
        if row>0:
            print("Percent complete:", round((it*1.0)/its,3))
            maxED=[]
            minED=[]
            if row>1:
                rowMin=min(repeatmatrix[row][0:row-1])
            else:
                rowMin=sys.maxsize
			
            for col in range(row+1,len(species)):
                minED.append(abs(repeatmatrix[0][col]-repeatmatrix[0][row]))
                maxED.append(repeatmatrix[0][col]+repeatmatrix[0][row])
			
            if len(maxED)>0:
                lowestMax = min(maxED)
            else:
                lowestMax = sys.maxsize
			
            for col in range(row+1,len(species)):
                it+=1
                colMin = min(repeatmatrix[col])
                if (minED[col-(row+1)] > lowestMax or minED[col-(row+1)] > rowMin) and (minED[col-(row+1)] > colMin): 
                    repeatmatrix[row][col]=sys.maxsize
                    repeatmatrix[col][row]=sys.maxsize
                    skippedCells+=1
                else:
                    if useBoundedLD == True:
                        repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0], max(colMin,rowMin))
                        if repeatmatrix[row][col] > max(colMin,rowMin):
                            repeatmatrix[row][col]=sys.maxsize
                    else:
                        repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                    repeatmatrix[col][row]=repeatmatrix[row][col]
                    if repeatmatrix[row][col] < rowMin:
                        rowMin=repeatmatrix[row][col]

        else:
            for col in range(len(species)):
                repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                repeatmatrix[col][row]=repeatmatrix[row][col]
            repeatmatrix[0][0] = sys.maxsize

    return repeatmatrix

def createRepeatMatrix(species, thresh=2):
    from Bio.Align import substitution_matrices as matlist
    from Bio import pairwise2
    matrix = matlist.load('BLOSUM62')

    repeatmatrix=[]
    for row in range(len(species)):
        repeatmatrix.append([ listtostring(labelss[row],";") ])
        repeatmatrix[row].extend([0]*len(species))
    its=len(species)*len(species)/2
    i=0
    for row in range(len(species)):
        for col in range(row,len(species)):
            repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0],thresh)
            repeatmatrix[col][row]=repeatmatrix[row][col]
            
            i+=1
        print("Percent complete:", round((i*1.0)/its,3))

    return repeatmatrix


def MatrixtoCSV(mat,file):
	with open(file, 'wb') as csv:
		csv.write(u'\ufeff'.encode('utf-8'))
		for row in mat:
			for col in row:
				csv.write(str(col)+",")
			csv.write("\n")


def makeNetwork(filename,thresh,bthresh,type="", species=None,BL=""):
    graph=""
    if type=="PB":
        print("Prune + bound*****************")
        print(filename)
        start = time.time()
        C = createMinMatrix(species, False)
        C = np.asarray(C)
        end = time.time()
        print("Pruning and Bounding:", (end-start))
        graph = makeGraphEfficiently(species, C, filename+"minNet")
        summary(graph)
    
    return graph

style = {}
style["edge_curved"] = False
style["margin"] = 50
style["edge_arrow_size"]=0.1
style["vertex.label.cex"]=0.4
style["vertex.label.family"]="Helvetica"
style["vetrex.label.color"]='black'
style["edge_width"]=0.5
style["edge_color"]="black"
style["vertex_size"]=9
style["vertex.label.font"]=1.5

USA_unique_T10 = pd.read_csv('SARS2+1_data.csv',header=0)

uni_seq_T10 = USA_unique_T10.iloc[:,1]
uni_seq_T10 = uni_seq_T10.values.tolist()
uni_labels_T10 = USA_unique_T10.iloc[:,0]
uni_labels_T10 = uni_labels_T10.values.tolist()
uni_acc_T10 = USA_unique_T10.iloc[:,2]
uni_acc_T10 = uni_acc_T10.values.tolist()
uni_date_T10 = USA_unique_T10.iloc[:,3]
uni_date_T10 = uni_date_T10.values.tolist()
uni_count_T10 = USA_unique_T10.iloc[:,4]
uni_count_T10 = uni_count_T10.values.tolist()

time_seq = USA_unique_T10.iloc[:,[4]]
time_seq = time_seq.values.tolist()
df = pd.DataFrame()
df.insert(0, 'label', (range(len(uni_seq_T10))))
df = df.applymap(str)
x_labels = df.iloc[:,[0]]
labelss=x_labels.values.tolist()
rfile="USATime10"

#create DiWANN network
diwann_usa_time10=Graph()				
diwann_usa_time10=makeNetwork(rfile,thresh="",bthresh="",type="PB",species=time_seq)
ig.plot(diwann_usa_time10,'diwann_USA_Time10.png',**style)

#clustering
diwann_usa_undirected=diwann_usa_time10.as_undirected()
communities= diwann_usa_undirected.community_multilevel()
clusters_cov=communities
ig.plot(clusters_cov, 'community_diwann_sars1+2_len.png',mark_groups=True, **style, vertex_label = USA_unique_T10.iloc[:,3], vertex_label_size = 9)
clust_list_sars_combo=list(clusters_cov)

#label propagation
communities= diwann_usa_undirected.community_label_propagation()
clusters_cov=communities
ig.plot(clusters_cov, 'community_diwann_sars1+2_lp_len.png',mark_groups=True, **style, vertex_label = USA_unique_T10.iloc[:,3], vertex_label_size = 9)


#community detection in diwann using betweenness
communities = diwann_usa_undirected.community_edge_betweenness(clusters = 38, directed = False)
clusters_cov=communities.as_clustering()
plot(clusters_cov, 'community_diwann_sars1+2_3b_len.png',mark_groups=True, **style, vertex_label = USA_unique_T10.iloc[:,3], vertex_label_size = 9)
clust_list_cov_b=list(clusters_cov)



