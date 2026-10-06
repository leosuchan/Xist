import numpy as np
import igraph as ig
import pandas as pd
from timeit import default_timer as timer
import math
import csv
import subprocess
import os.path
import scipy.sparse as sp
import sys
import os
from PIL import Image
import leidenalg
import kahip
import pymetis
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import SpectralClustering
os.environ["R_LIBS_USER"] = ""
os.environ["R_LIBS_SITE"] = ""
from rpy2 import robjects
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
#from rpy2.rlike.container import NamedList, TaggedList
from rpy2.robjects.vectors import ListVector
from collections.abc import Mapping

#### General helper methods:

def reduce_to_lcc(df, preserve=False, reduce=True):
    if len(df.columns) > 2:
        weightlst = df.iloc[:,2]
    else:
        weightlst = [1] * len(df)
    df_graph = ig.Graph(edges=df.iloc[:,0:2].values.tolist(), edge_attrs={'weight': weightlst})
    if reduce:
        lcc_set = max(df_graph.connected_components(mode='weak'), key=len)
    else:
        lcc_set = df_graph.vs
    dfg = df_graph.subgraph(lcc_set).simplify(combine_edges={'weight': max})
    if preserve:
        return pd.DataFrame([list(sorted(e.tuple)) + [e['weight']] for e in dfg.es]), [preserve[i.index] for i in lcc_set]
    else:
        return pd.DataFrame([list(sorted(e.tuple)) + [e['weight']] for e in dfg.es])

def transform_to_assignment(S):
    S_all = [0] * sum(len(i) for i in S)
    for i in range(1, len(S)):
        for j in S[i]:
            S_all[j] = i
    return S_all

def classification_rate(partition, T):
    if len(partition)==0:
        return math.nan
    elif type(partition[0]) is list:
        cr_best = -1
        for idx in range(len(partition)):
            S = [1 if i == idx else 0 for i in transform_to_assignment(partition)]
            cr_current = max(sum(S[i]==T[i] for i in range(len(S))), sum(S[i]==1-T[i] for i in range(len(S))))
            if cr_current > cr_best:
                cr_best = cr_current
        #return cr_current / len(T)
        return cr_best / len(T)
    else:
        return max(sum(partition[i]==T[i] for i in range(len(partition))), sum(partition[i]==1-T[i] for i in range(len(partition)))) / len(T)

def ncut_value(g, partition):
    if len(partition) == len(g.vs) and (len(partition) <= 2 or partition[2] in [0,1]):
        S = [i for i in range(len(partition)) if partition[i]]
    else:
        S = partition
    S_C = list(set(range(len(g.vs))) - set(S))
    mc_value = sum(g.es.select(_between = (S, S_C))['weight'])
    vol_S = 2 * sum(g.es.select(_within = S)['weight']) + mc_value
    vol_S_C = 2 * sum(g.es.select(_within = S_C)['weight']) + mc_value
    if vol_S and vol_S_C:
        return mc_value / (vol_S * vol_S_C)
    else:
        return math.inf

def generate_tiff_edges(imgpath, m, t=1, sigm=math.nan, weights_by_sample=True):
    xedgevec, yedgevec, weightvec = [], [], []
    img = Image.open(imgpath)
    img_arr = np.array(img)[8:,:504]
    img_arr = img_arr / np.sum(img_arr)
    ndim_img = np.shape(img_arr)
    xm = int(ndim_img[0]/m)
    ym = int(ndim_img[1]/m)
    sample = np.array([[np.sum(img_arr[i*xm:(i+1)*xm,j*ym:(j+1)*ym]) for j in range(m)] for i in range(m)])
    #if sample[-1,-1]==0:
    #    print(sample)
    ndim = np.shape(sample)
    
    for x1 in range(ndim[0]):
        for x2 in range(ndim[1]):
            if sample[x1][x2] > 0:
                for y1 in range(max(0, x1 - math.floor(t)), min(ndim[0], x1 + math.floor(t) + 1)):
                    for y2 in range(x2, min(ndim[1], x2 + math.floor(math.sqrt(t ** 2 - (x1 - y1) ** 2)) + 1)):
                        if (x2 != y2 or x1 < y1) and sample[y1][y2] > 0:
                            xedgevec.append(x2 + x1 * ndim[1])
                            yedgevec.append(y2 + y1 * ndim[1])
                            if math.isnan(sigm) and weights_by_sample:
                                current_weight = sample[x1][x2] * sample[y1][y2]
                            elif math.isnan(sigm):
                                current_weight = 1
                            elif weights_by_sample:
                                current_weight = sample[x1][x2] * sample[y1][y2] * math.exp(
                                    -math.sqrt((x1 - y1) ** 2 + (x2 - y2) ** 2) / sigm)
                            else:
                                current_weight = math.exp(-math.sqrt((x1 - y1) ** 2 + (x2 - y2) ** 2) / sigm)
                            weightvec.append(current_weight)
    edge_df, _ = reduce_to_lcc(pd.DataFrame({"id_1": xedgevec, "id_2": yedgevec, "weight": weightvec}), preserve=[0] * len(xedgevec), reduce=False)
    return edge_df, sample


###### FUNCTIONS FOR UNWEIGHTED GRAPHS

#### Algorithm-specific helper methods:

def adjacency_for_metis(edgearray, start_at=0):
    adj = [ [] for i in range(len(np.unique(edgearray.iloc[:,0:2]))) ]
    if len(edgearray.columns) == 2:
        edgearray.columns = ['id_1', 'id_2']
        for idx, row in edgearray.iterrows():
                adj[int(row['id_1'])] += [int(row['id_2']) + start_at]
                adj[int(row['id_2'])] += [int(row['id_1']) + start_at]
    elif len(edgearray.columns) == 3:
        edgearray.columns = ['id_1', 'id_2', 'weight']
        for idx, row in edgearray.iterrows():
            adj[int(row['id_1'])] += [int(row['id_2']) + start_at, row['weight']]
            adj[int(row['id_2'])] += [int(row['id_1']) + start_at, row['weight']]
    else:
        return []
    return adj

def adjacency_for_kahip_unweighted(edgearray):
    adj = adjacency_for_metis(edgearray.iloc[:,0:2])
    xadj = np.cumsum([len(a) for a in adj])
    adjncy = [x for a in adj for x in a]
    return [list(np.concatenate(([0], xadj))), adjncy, edgearray.iloc[:,2]]

def adjacency_for_scoreplus(df, num_nodes=None):
    if num_nodes== None:
        num_nodes = df[['id_1', 'id_2']].to_numpy().max() + 1  # number of nodes
    adj = np.zeros((num_nodes, num_nodes))                    # initialize with zeros
    for _, row in df.iterrows():
        i, j = int(row['id_1']), int(row['id_2'])
        adj[i, j] = 1
        adj[j, i] = 1 # alternative: w, but we only use score on unweighted matrizes  
    return adj


#### Unweighted algorithms:

def xist_unweighted(df):
    adj = adjacency_for_metis(df)
    st = timer()
    df_ncut = []
    df_ncut_value = math.inf
    df_graph = ig.Graph(edges=df.values.tolist())
    deg_lst = df_graph.degree()
    locmax_01 = [1] * len(deg_lst)
    for i in range(len(deg_lst)):
        for j in adj[i]:
            if deg_lst[i] > deg_lst[j]:
                locmax_01[j] = 0
            elif deg_lst[i] < deg_lst[j]:
                locmax_01[i] = 0
    locmax = [i for i in range(len(locmax_01)) if locmax_01[i]]
    if len(locmax) == 1:
        return [math.inf, []]
    tau_lst = [0] * len(locmax)
    for i in range(1, len(locmax)):
        stmc = df_graph.mincut(source=locmax[i], target=locmax[tau_lst[i]])
        if i in stmc.partition[0]:
            part_i = stmc.partition[0]
        else:
            part_i = stmc.partition[1]
        vol_adj_cut_S = sum(deg_lst[j] for j in part_i)
        vol_adj_cut_S_C = sum(deg_lst) - vol_adj_cut_S
        stmc_ncut_value = stmc.value / (vol_adj_cut_S * vol_adj_cut_S_C)
        if stmc_ncut_value < df_ncut_value:
            df_ncut = stmc.partition
            df_ncut_value = stmc_ncut_value
        for j in range(i + 1, len(locmax)):
            if locmax[j] in part_i and tau_lst[j] == tau_lst[i]:
                tau_lst[j] = i
    et = timer()
    return [df_ncut_value, et - st, df_ncut]

def leidenoracle_unweighted(g):
    best_ncut_value, best_partition = math.inf, []
    for eps in range(0.01,1,0.01):
        lpart = leidenalg.find_partition(dfg, leidenalg.CPMVertexPartition, resolution_parameter=eps)
        lpart_ncut_value = ncut_value2(lpart[0])
        if lpart_ncut_value < best_ncut_value:
            best_ncut_value, best_partition = lpart_ncut_value, lpart[0]
    return(best_ncut_value, best_partition)

def ncut_kahip_best_unweighted(df, prnt=False):
    minpart = [math.inf, math.inf]
    xadj_m, adjncy_m = adjacency_for_kahip(df)
    st = timer()
    for i in range(1,101):
        kahip_part = ncut_kahip_unweighted(xadj_m, adjncy_m, mode=5, imbalance=i/100)
        if kahip_part[0] < minpart[0]:
            minpart = kahip_part
            if prnt:
                print("Iteration", i, "with imbalance", i/100, "and NCut", kahip_part[0])
    et = timer()
    if print:
        print("Best KaHIP NCut for musae_squirrel:", minpart, "in total time of", et-st, "seconds")
    return [minpart, et-st]

def ncut_kahip_unweighted(xadj, adjncy, mode=5, imbalance=0.03, seed=0):
    vwgt = [1] * (len(xadj) - 1)
    adjcwgt = [1] * len(adjncy)
    st = timer()
    adj_mcut_value, adj_cut = kahip.kaffpa(vwgt, xadj, adjcwgt, adjncy, 2, imbalance, 0, seed, mode)
    et = timer()
    vol_adj_cut_S = sum((xadj[i+1] - xadj[i]) * adj_cut[i] for i in range(len(adj_cut)))
    vol_adj_cut_S_C = sum((xadj[i+1] - xadj[i]) * (1-adj_cut[i]) for i in range(len(adj_cut)))
    if vol_adj_cut_S > 0 and vol_adj_cut_S_C > 0:
        return [adj_mcut_value / (vol_adj_cut_S * vol_adj_cut_S_C), et - st]
    else:
        return [math.inf, et - st, adj_cut]
# Mode for KaHIP: FAST=0, ECO=1, STRONG=2, FASTSOCIAL=3, ECOSOCIAL=4, STRONGSOCIAL=5

def ncut_metis_unweighted(edgearray):
    adj = adjacency_for_metis(edgearray)
    st = timer()
    adj_mcut_value, adj_cut = pymetis.part_graph(2, adjacency=adj)
    et = timer()
    vol_adj_cut_S = sum(len(adj[i]) * adj_cut[i] for i in range(len(adj_cut)))
    vol_adj_cut_S_C = sum(len(adj[i]) * (1 - adj_cut[i]) for i in range(len(adj_cut)))
    return [adj_mcut_value / (vol_adj_cut_S * vol_adj_cut_S_C), et - st]

def ncut_chaco_unweighted(df, name, numit=1, print_chaco_output=True):
    chaco_best_ncut_value = math.inf
    adj_m = [[x+1 for x in a] for a in adjacency_for_metis(df)]
    with open("/home/path/Chaco/Chaco-2.2/exec/{0}_chaco.graph".format(name), "w") as f:
        f.write("%s\t%s\n" %(max(map(max,adj_m)), int(sum([len(a) for a in adj_m])/2)))
        wr = csv.writer(f, delimiter="\t")
        wr.writerows(adj_m)
    if(print_chaco_output):
        call(["/home/path/Chaco/Chaco-2.2/exec/bash_chaco.sh", "{0}_chaco.graph".format(name), str(numit)])
    else:
        check_output(["/home/path/Chaco/Chaco-2.2/exec/bash_chaco.sh", "{0}_chaco.graph".format(name), str(numit)])
    df_graph = ig.Graph(edges=df.values.tolist()).simplify()
    for k in range(numit):
        chaco_p1 = open("/home/path/Chaco/Chaco-2.2/exec/ChacoOutput/cout{0}.txt".format(k), "r").read().split('\n')
        chaco_p = [int(i) for i in chaco_p1[:-1]]
        chaco_ncut_p = ncut_value(df_graph, chaco_p)
        if chaco_ncut_p < chaco_best_ncut_value:
            chaco_imbalance = k/100
            chaco_best_ncut_value = chaco_ncut_p
            chaco_best_ncut = chaco_p
    return [chaco_best_ncut_value, chaco_imbalance, chaco_best_ncut]



# ScorePlus is designed for unweighted graphs, so we just use it here

# Load your R script
robjects.r['source']('SCOREplus.R')

# Get the R function
SCOREplus = robjects.globalenv['SCOREplus']

def ncut_scoreplus(df):
    df.columns = ['id_1', 'id_2', 'weight']
    adj= adjacency_for_scoreplus(df)
    #st= timer()
    with localconverter(default_converter + numpy2ri.converter):
        st= timer()
        result = SCOREplus(adj, 2) #k=2
        et= timer()
    #et= timer()
    # Handle both NamedList and TaggedList returns
    #if isinstance(result, (NamedList, TaggedList)):
    #    py_result = dict(zip(result.names(), list(result)))
        # Convert numeric components to NumPy arrays
    #    labels = np.array(py_result['labels'])
    #else:
    #    print("Unexpected result type:", type(result))
    #if isinstance(result, ListVector):
    #    py_result = dict(zip(result.names, list(result)))
    #    labels = np.array(py_result['labels'])
    #else:
    #    print("Unexpected result type:", type(result))
    #    raise TypeError(type(result))
    if isinstance(result, Mapping):
        labels = np.array(result['labels'])
    else:
        py_result = dict(zip(result.names, list(result)))
        labels = np.array(py_result['labels'])
    adj_cut= labels-1 # get a list of 0 and 1s instead of 1 and 2s. 
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (adj_cut[int(e['id_1'])] and not adj_cut[int(e['id_2'])])
                         or (adj_cut[int(e['id_2'])] and not adj_cut[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S
    return adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)),et-st, adj_cut



###### FUNCTIONS FOR WEIGHTED GRAPHS

#### Algorithm helper methods:

def adjacency_for_kahip(edgearray, toint=False):
    edgearray.columns = ['id_1', 'id_2', 'weight']
    adj = [ [] for i in range(len(np.unique(edgearray.iloc[:,0:2]))) ]
    eadj = [ [] for i in range(len(np.unique(edgearray.iloc[:,0:2]))) ]
    for idx, row in edgearray.iterrows():
        adj[int(row['id_1'])] += [int(row['id_2'])]
        adj[int(row['id_2'])] += [int(row['id_1'])]
        eadj[int(row['id_1'])] += [row['weight']]
        eadj[int(row['id_2'])] += [row['weight']]
    xadj = np.cumsum([len(a) for a in adj])
    adjncy = [x for a in adj for x in a]
    eweights = [x for a in eadj for x in a]
    if toint and not all(isinstance(x, int) for x in eweights):
        eweights_normalizer = 2147483647 / sum(eweights)
        eweights = [round(w * eweights_normalizer) for a in eadj for w in a]
    return list(np.concatenate(([0], xadj))), adjncy, eweights


#### Weighted algorithms:

def xvst(df):
    st = timer()
    df_ncut = []
    df_ncut_value = math.inf
    df_graph = ig.Graph(edges=df.iloc[:,0:2].values.tolist(), edge_attrs={'weight': df.iloc[:,2]})
    deg_lst = df_graph.strength(weights='weight')
    for i in range(len(deg_lst)-1):
        for j in range(i+1, len(deg_lst)):
            stmc = df_graph.mincut(source=i, target=j, capacity='weight')
            if i in stmc.partition[0]:
                part_i = stmc.partition[0]
            else:
                part_i = stmc.partition[1]
            vol_part_i = sum([deg_lst[i] for i in part_i])
            stmc_ncut_value = stmc.value / (vol_part_i * (sum(deg_lst) - vol_part_i))
            if stmc_ncut_value < df_ncut_value:
                df_ncut = stmc.partition
                df_ncut_value = stmc_ncut_value
    et = timer()
    return [df_ncut_value, et - st, df_ncut]

def xist(df):
    st = timer()
    df_ncut = []
    df_ncut_value = math.inf
    df_graph = ig.Graph(edges=df.iloc[:,0:2].values.tolist(), edge_attrs={'weight': df.iloc[:,2]})
    deg_lst = df_graph.strength(weights='weight')
    locmax_01 = [True] * len(deg_lst)
    for e in df_graph.es:
        if deg_lst[e.source] < deg_lst[e.target]:
            locmax_01[e.source] = False
        elif deg_lst[e.source] > deg_lst[e.target]:
            locmax_01[e.target] = False
    locmax = [i for i in range(len(locmax_01)) if locmax_01[i]]
    if len(locmax) == 1:
        et = timer()
        return [math.inf, et - st, []]
    tau_lst = [0] * len(locmax)
    for i in range(1,len(locmax)):
        stmc = df_graph.mincut(source=locmax[i], target=locmax[tau_lst[i]], capacity='weight')
        if locmax[i] in stmc.partition[0]:
        #if i in stmc.partition[0]:
            part_i = stmc.partition[0]
        else:
            part_i = stmc.partition[1]
        vol_part_i = sum([deg_lst[i] for i in part_i])
        vol_part_i_C = sum(deg_lst) - vol_part_i
        if vol_part_i > 0 and vol_part_i_C > 0:
            stmc_ncut_value = stmc.value / (vol_part_i * vol_part_i_C)
            if stmc_ncut_value < df_ncut_value:
                df_ncut = stmc.partition
                df_ncut_value = stmc_ncut_value
        for j in range(i+1, len(locmax)):
            if locmax[j] in part_i and tau_lst[j]==tau_lst[i]:
                tau_lst[j] = i
    et = timer()
    return [df_ncut_value, et - st, df_ncut]




def leidenoracle(df, exponential_resolution_scaling=False, classif_truth=None):
    best_ncut_value, best_part_n, best_eps_n, best_time_n = math.inf, [], math.nan, math.nan
    best_classif, best_part_c, best_eps_c, best_time_c = -math.inf, [], math.nan, math.nan
    df_graph = ig.Graph(edges=df.iloc[:,0:2].values.tolist(), edge_attrs={'weight': df.iloc[:,2]})
    deg_lst = df_graph.strength(weights='weight')
    if exponential_resolution_scaling:
        eps_lst = [(1.1**x - 1)/100000 for x in range(1,101)]
    else:
        eps_lst = [5*x/100 * 10**(-5) for x in range(1,101)] # OR: eps_lst = [x/100 for x in range(1,101)]
    st_total = timer()
    for eps in eps_lst:
        st = timer()
        lpart = leidenalg.find_partition(df_graph, leidenalg.CPMVertexPartition, resolution_parameter=eps)
        et = timer()
        if len(lpart) > 1 and len(lpart) < 10: #len(df_graph.vs):
            S_01 = [0] * len(df_graph.vs)
            for i in range(1,len(lpart)):
                for j in lpart[i]:
                    S_01[j] = i
            vol_lpart = [sum(deg_lst[i] for i in S) for S in lpart]
            total_deg = sum(deg_lst)
            lpart_ncut_value = sum(sum(e['weight'] for e in df_graph.es if (S_01[e.source] == i and S_01[e.target] != i)
                                       or (S_01[e.source] != i and S_01[e.target] == i)) / (vol_lpart[i] * (total_deg - vol_lpart[i])) for i in range(len(lpart))) / 2
            if classif_truth is not None:
                lpart_classif = classification_rate(lpart, classif_truth)
            if lpart_ncut_value < best_ncut_value:
                best_ncut_value, best_part_n, best_eps_n, best_time_n = lpart_ncut_value, lpart, eps, et - st
            if classif_truth is not None and lpart_classif > best_classif:
                best_classif, best_part_c, best_eps_c, best_time_c = lpart_classif, lpart, eps, et - st
    et_total = timer()
    if classif_truth is not None:
        return [[best_ncut_value, best_part_n, best_eps_n, best_time_n], [best_classif, best_part_c, best_eps_c, best_time_c], et_total - st_total]
    else:
        return [[best_ncut_value, best_part_n, best_eps_n, best_time_n], et_total - st_total]

def ncut_kahip(df, mode=5, imbalance=0.03, seed=8472):
    xadj, adjncy, adjcwgt = adjacency_for_kahip(df, toint=True)
    vwgt = [1] * (len(xadj) - 1)
    st = timer()
    _, adj_cut = kahip.kaffpa(vwgt, xadj, adjcwgt, adjncy, 2, imbalance, 0, seed, mode)
    et = timer()
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (adj_cut[int(e['id_1'])] and not adj_cut[int(e['id_2'])])
                         or (adj_cut[int(e['id_2'])] and not adj_cut[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S
    if adj_inc_S + adj_out_S > 0 and adj_inc_S_C + adj_out_S_C > 0:
        return [adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)), et - st, adj_cut]
    else:
        return [math.inf, et - st, adj_cut]
# Mode for KaHIP: FAST=0, ECO=1, STRONG=2, FASTSOCIAL=3, ECOSOCIAL=4, STRONGSOCIAL=5

def ncut_metis(df):
    df.columns = ['id_1', 'id_2', 'weight']
    xadj, adjncy, adjcwgt = adjacency_for_kahip(df, toint=True)
    st = timer()
    _, adj_cut = pymetis.part_graph(2, xadj=xadj, adjncy=adjncy, eweights=adjcwgt)
    et = timer()
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (adj_cut[int(e['id_1'])] and not adj_cut[int(e['id_2'])])
                         or (adj_cut[int(e['id_2'])] and not adj_cut[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S
    return [adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)), et - st, adj_cut]

def ncut_chaco(df, name, numit=1, print_chaco_output=True):
    st, et = 0, math.inf
    chaco_best_ncut_value = math.inf
    chaco_imbalance = math.nan
    chaco_best_ncut = []
    adj_m = adjacency_for_metis(df, start_at=1)
    with open("/home/path/Chaco/Chaco-2.2/exec/{0}_chaco.graph".format(name), "w") as f:
        f.write("%s\t%s\t001\n" %(max(map(max,adj_m)), int(sum([len(a) for a in adj_m])/4)))
        wr = csv.writer(f, delimiter="\t")
        wr.writerows(adj_m)
    if(print_chaco_output):
        subprocess.call(["/home/path/Chaco/Chaco-2.2/exec/bash_chaco.sh", "{0}_chaco.graph".format(name), str(numit)])
        exit_status = 0
    else:
        try:
            st = timer()
            subprocess.check_output(["/home/path/Chaco/Chaco-2.2/exec/bash_chaco.sh", "{0}_chaco.graph".format(name), str(numit)], timeout=30)
            et = timer()
            exit_status = 0
        except subprocess.CalledProcessError as e:
            print("Error:", e.returncode, "\n")
            exit_status = 1
        except subprocess.TimeoutExpired as e:
            print("Timeout after", e.timeout, "seconds\n")
            exit_status = -1
    if not exit_status:
        df_graph = ig.Graph(edges=df.iloc[:,0:2].values.tolist(), edge_attrs = {'weight': df.iloc[:,2]}).simplify(combine_edges={'weight': max})
        for k in range(numit):
            if os.path.isfile("/home/path/Chaco/Chaco-2.2/exec/ChacoOutput/cout{0}.txt".format(k)):
                chaco_p1 = open("/home/path/Chaco/Chaco-2.2/exec/ChacoOutput/cout{0}.txt".format(k), "r").read().split('\n')
                chaco_p = [int(i) for i in chaco_p1[:-1]]
                chaco_ncut_p = ncut_value(df_graph, chaco_p)
                if chaco_ncut_p < chaco_best_ncut_value:
                    chaco_imbalance = k/100
                    chaco_best_ncut_value = chaco_ncut_p
                    chaco_best_ncut = chaco_p
    return [chaco_best_ncut_value, chaco_imbalance, chaco_best_ncut, et - st]

def spectral_clustering(df, delta=0.2):
    df.columns = ['id_1', 'id_2', 'weight']
    maxind = int(np.max(df.iloc[:,-2]))
    df_simmat = np.zeros((maxind+1, maxind+1))
    for _, e in df.iterrows():
        df_simmat[int(e['id_1']),int(e['id_2'])] = df_simmat[int(e['id_2']),int(e['id_1'])] = np.exp(-e['weight']**2 / (2*delta**2))
        # This is the suggested method of computing the similarity matrix from the sklearn documentation.
    st = timer()
    clust = SpectralClustering(n_clusters=2, affinity='precomputed', random_state=8472).fit(df_simmat).labels_
    et = timer()
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (clust[int(e['id_1'])] and not clust[int(e['id_2'])])
                         or (clust[int(e['id_2'])] and not clust[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if clust[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if clust[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S

    return adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)), et - st, clust

def spectral_clustering_sparse(df, delta=0.2):
    df.columns = ['id_1', 'id_2', 'weight']
    maxind = int(np.max(df.iloc[:,-2]))
    df_simmat = np.zeros((maxind+1, maxind+1))
    for _, e in df.iterrows():
        df_simmat[int(e['id_1']),int(e['id_2'])] = df_simmat[int(e['id_2']),int(e['id_1'])] = np.exp(-e['weight']**2 / (2*delta**2))
        # This is the suggested method of computing the similarity matrix from the sklearn documentation.
    st = timer()
    clust = SpectralClustering(n_clusters=2, affinity='precomputed', random_state=8472).fit(df_simmat).labels_
    et = timer()
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (clust[int(e['id_1'])] and not clust[int(e['id_2'])])
                         or (clust[int(e['id_2'])] and not clust[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if clust[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if clust[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S

    return adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)), et - st, clust



def ncut_xcut(df, name, print_xcut_output=True, weighted= True):
    df.columns = ['id_1', 'id_2', 'weight']
    df_ori= df.copy()
    if weighted and not all(isinstance(x, int) for x in df.iloc[:, 2]): 
        weights_normalizer = 2147483647 /(df.iloc[:, 2].sum()+1)
        df.iloc[:, 2] = (df.iloc[:, 2] * weights_normalizer).astype(int)+1
    #df = df[df.iloc[:, 2] != 0]
    df.iloc[:, 0] = df.iloc[:, 0] + 1
    df.iloc[:, 1] = df.iloc[:, 1] + 1
    df.to_csv(f"/home/path/code/helpdata/{name}.edges",
              sep="\t", header=False, index=False)
    if(print_xcut_output):
        st= timer()
        subprocess.call(["/home/path/xcut/build/app/xcut", "/home/path/code/helpdata/{0}.edges".format(name), "2"], cwd="/home/path/code/helpdata")
        et= timer()
    else:
        st= timer()
        subprocess.check_output(["/home/path/xcut/build/app/xcut", "/home/path/code/helpdata/{0}.edges".format(name), "2"], cwd="/home/path/code/helpdata")
        et = timer()
    adj_cut = np.loadtxt(f"/home/path/code/helpdata/{name}.part2",
                delimiter="\t")
    df.iloc[:, 0] = df.iloc[:, 0] - 1
    df.iloc[:, 1] = df.iloc[:, 1] - 1
    adj_mcut_value = sum(e['weight'] for _, e in df_ori.iterrows() if (adj_cut[int(e['id_1'])] and not adj_cut[int(e['id_2'])])
                         or (adj_cut[int(e['id_2'])] and not adj_cut[int(e['id_1'])]))
    adj_inc_S = sum(df_ori['weight'][i] for i in range(len(df)) if adj_cut[int(df_ori['id_1'][i])])
    adj_out_S = sum(df_ori['weight'][i] for i in range(len(df)) if adj_cut[int(df_ori['id_2'][i])])
    adj_inc_S_C = sum(df_ori['weight']) - adj_inc_S
    adj_out_S_C = sum(df_ori['weight']) - adj_out_S
    return adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C)),et-st, adj_cut


# wrapper for C++ implemention of xist and xvst

import xist_dinic_faster
import xvst_dinic_faster



def _prepare_edges_for_cpp(df):
    """Return list of (u,v,w) with integer weights suitable for C++ capacity fields.
    Scales small/float weights to reasonable integer range to avoid DP algorithm hangs.
    """
    # Extract columns efficiently
    ucol = df.iloc[:, 0].values
    vcol = df.iloc[:, 1].values
    wcol = df.iloc[:, 2].values.astype(float)

    # Check if weights are already small integers
    min_w = wcol.min()
    max_w = wcol.max()
    
    #if min_w >= 1 and max_w < 100000:
    if np.all(wcol % 1 == 0):
        # Already reasonable, use as-is
        int_weights = wcol.astype(np.int64)
    else:
        eweights_normalizer = 2147483647 / wcol.sum()

        int_weights = np.round(wcol * eweights_normalizer).astype(np.int64)
        # alternative add int_weights = np.maximum(int_weights, 1)

    edges = [(int(u), int(v), int(w)) for u, v, w in zip(ucol, vcol, int_weights)]
    return edges
    

def xist_dinic_faster_wrapper(df):
    """Call C++ xist_dinic_faster and compute ncut in Python with timing."""
    edges = _prepare_edges_for_cpp(df)
    df.columns = ['id_1', 'id_2', 'weight'] # this is new
    edges_objects = [xist_dinic_faster.EdgeInput(u, v, w) for u, v, w in edges]
    n = int(max(max(u, v) for u, v in [[e[0], e[1]] for e in edges]) + 1) if edges else 0
    st= timer()
    cpp_result = xist_dinic_faster.xist(n, edges_objects)
    et= timer()
    #ncut_value = cpp_result[0]
    elapsed= et-st  
    partition = cpp_result[1]
    if not partition:
        return [np.inf, elapsed, partition]
    adj_cut = np.zeros(n, dtype=bool)
    adj_cut[partition[0]] = True
    adj_mcut_value = sum(e['weight'] for _, e in df.iterrows() if (adj_cut[int(e['id_1'])] and not adj_cut[int(e['id_2'])])
                         or (adj_cut[int(e['id_2'])] and not adj_cut[int(e['id_1'])]))
    adj_inc_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_1'][i])])
    adj_out_S = sum(df['weight'][i] for i in range(len(df)) if adj_cut[int(df['id_2'][i])])
    adj_inc_S_C = sum(df['weight']) - adj_inc_S
    adj_out_S_C = sum(df['weight']) - adj_out_S
    ncut_value = adj_mcut_value / ((adj_inc_S + adj_out_S) * (adj_inc_S_C + adj_out_S_C))
    return [ncut_value, elapsed, partition]

def xvst_dinic_faster_wrapper(df):
    """Call the C++ xvst implementation and return [ncut_value, elapsed, partition]."""
    edges = _prepare_edges_for_cpp(df)
    edges_objects = [xvst_dinic_faster.EdgeInput(u, v, w) for u, v, w in edges]
    n = int(max(max(u, v) for u, v, _ in edges) + 1) if edges else 0

    st = timer()
    cpp_result = xvst_dinic_faster.xvst(n, edges_objects)
    et = timer()

    ncut_value = float(cpp_result[0])
    elapsed = et - st
    partition = cpp_result[1]

    return [ncut_value, elapsed, partition]


