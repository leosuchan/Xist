from xist import *
from tqdm import tqdm
import warnings
from concurrent.futures import ProcessPoolExecutor
import networkx as nx


#### Applying Xist to a dataset:

artists = reduce_to_lcc(pd.read_csv('Datasets/artist_edges.csv'))
artists_xist = xist(artists)
artists_xist = xist_dinic_faster_wrapper(artists)
print("Xist (NCut) computed a value of", artists_xist[0], "for the artists dataset in", artists_xist[1],"seconds.")


#### NCut value and classification rate comparison between Xist and other algorithms on a Gaussian mixture both for weighted and unweighted graphs:


def generate_mixed_gaussian_edges(BigN, delta, k, r, sigma, weighted= True):
    Nsize = np.random.binomial(n=BigN, p=1/2)
    smpl = np.concatenate((np.random.multivariate_normal([0,0], [[1,0],[0,1]], Nsize), np.random.multivariate_normal([delta,delta], [[1,0],[0,1]], BigN-Nsize)))
    kNN_model = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(smpl)
    dstncs, indices = kNN_model.kneighbors(smpl)
    if weighted:
        kNN_edges = [[i,indices[i][j],math.exp(-math.sqrt(sum((smpl[i] - smpl[indices[i][j]])**2))/sigma)] for j in range(1,k+1) for i in range(BigN)]
        exp_edges = [[i,j,math.exp(-math.sqrt(sum((smpl[i] - smpl[j])**2))/sigma)] for i in range(BigN-1) for j in range(i+1,BigN) if sum((smpl[i] - smpl[j])**2) <= r]
    else:
        kNN_edges = [[i,indices[i][j], 1] for j in range(1,k+1) for i in range(BigN)]
        exp_edges = [[i,j, 1] for i in range(BigN-1) for j in range(i+1,BigN) if sum((smpl[i] - smpl[j])**2) <= r]
    mixed_gaussian_edges, mg_assign = reduce_to_lcc(pd.DataFrame(kNN_edges + exp_edges), preserve=[0] * Nsize + [1] * (BigN - Nsize), reduce=False)
    return mixed_gaussian_edges, mg_assign, smpl


# Runs a single iteration for one delta and returns (ncut values, rate values) for all algorithms
def _single_iteration_mixed_gaussian_unweighted(BigN, delta, k, r, sigma, seed, name='mixed_gaussian'):
    np.random.seed(int(seed))
    mge_df, mge_truth, _ = generate_mixed_gaussian_edges(
        BigN, delta, k, r, sigma, weighted=False
    )

    mge_xist = xist_dinic_faster_wrapper(mge_df)
    #mge_xist = xist(mge_df)
    rate_xist = classification_rate(mge_xist[2], mge_truth)
    ncut_xist = mge_xist[0]

    mge_scoreplus = ncut_scoreplus(mge_df)
    rate_scoreplus = classification_rate(mge_scoreplus[2], mge_truth)
    ncut_scoreplus_val = mge_scoreplus[0]

    return (
        [ncut_xist, ncut_scoreplus_val],
        [rate_xist, rate_scoreplus]
    )


def mixed_gaussian_alg_test_unweighted(N_iter, BigN, delta_lst, k, r, sigma, num_cpus=25):
    gauss_rate_trafo, gauss_ncut_trafo = [], []
    gauss_rate_trafo_mean, gauss_ncut_trafo_mean = [], []
    algnames = ["Xist", "Scoreplus"]

    for delta in tqdm(delta_lst):
        # Run N_iter independent computations in parallel
        rng = np.random.default_rng()
        seeds = rng.integers(low=0, high=2**32 - 1, size=N_iter, dtype=np.uint32)
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            results = list(
                executor.map(
                    _single_iteration_mixed_gaussian_unweighted,
                    [BigN] * N_iter,
                    [delta] * N_iter,
                    [k] * N_iter,
                    [r] * N_iter,
                    [sigma] * N_iter,
                    seeds.tolist(),
                    [f"mixed_gaussian_{i}" for i in range(N_iter)],
                )
            )

        # Separate ncut and rate lists
        all_ncut = [res[0] for res in results]   # list of length 3
        all_rate = [res[1] for res in results]   # list of length 3


        all_ncut = np.array(all_ncut)   # shape (N_iter, 3)
        all_rate = np.array(all_rate)   # shape (N_iter, 3)

        for i, alg in enumerate(algnames):
            gauss_ncut_trafo.append([delta, np.nanmedian(all_ncut[:, i]), alg])
            gauss_rate_trafo.append([delta, np.nanmedian(all_rate[:, i]), alg])
            gauss_ncut_trafo_mean.append([delta, np.nanmean(all_ncut[:, i]), alg])
            gauss_rate_trafo_mean.append([delta, np.nanmean(all_rate[:, i]), alg])            

    df_rate = pd.DataFrame(gauss_rate_trafo, columns=["x", "values", "alg"])
    df_ncut = pd.DataFrame(gauss_ncut_trafo, columns=["x", "values", "alg"])
    df_rate_mean = pd.DataFrame(gauss_rate_trafo_mean, columns=["x", "values", "alg"])
    df_ncut_mean = pd.DataFrame(gauss_ncut_trafo_mean, columns=["x", "values", "alg"])

    return df_rate, df_ncut, df_rate_mean, df_ncut_mean



def ncut_truth(mge_df, mge_truth):
    G = ig.Graph(edges=mge_df.iloc[:, 0:2].values.tolist(), edge_attrs={'weight': mge_df.iloc[:, 2]})
    deg = G.strength(weights='weight')
    S = [i for i, c in enumerate(mge_truth) if c == 0]
    T = [i for i, c in enumerate(mge_truth) if c == 1]
    cut_value = 0.0
    for e in G.es:
        u, v = e.tuple
        if mge_truth[u] != mge_truth[v]:
            cut_value += e['weight']
    vol_S = sum(deg[i] for i in S)
    vol_T = sum(deg[i] for i in T)
    true_ncut = cut_value / (vol_S * vol_T)
    return true_ncut


def _single_iteration_mixed_gaussian_weighted(BigN, delta, k, r, sigma, seed, name='mixed_gaussian'):
    np.random.seed(int(seed))
    mge_df, mge_truth, _ = generate_mixed_gaussian_edges(BigN, delta, k, r, sigma, weighted=True)

    mge_xist = xist_dinic_faster_wrapper(mge_df)
    #mge_xist = xist(mge_df)
    rate_xist = classification_rate(mge_xist[2], mge_truth)
    ncut_xist = mge_xist[0]

    mge_leiden = leidenoracle(mge_df, exponential_resolution_scaling=True, classif_truth=mge_truth)
    rate_leiden= mge_leiden[1][0]
    ncut_leiden= mge_leiden[0][0]

    mge_kahip = ncut_kahip(mge_df)
    rate_kahip= classification_rate(mge_kahip[2], mge_truth)
    ncut_kahip_val= mge_kahip[0]

    mge_metis = ncut_metis(mge_df)
    rate_metis= classification_rate(mge_metis[2], mge_truth)
    ncut_metis_val= mge_metis[0]

    mge_xcut= ncut_xcut(mge_df, name, print_xcut_output= False)
    rate_xcut= classification_rate(mge_xcut[2], mge_truth)
    ncut_xcut_val= mge_xcut[0]

    rate_truth= classification_rate(mge_truth, mge_truth)
    ncut_truth_val= ncut_truth(mge_df, mge_truth)


    return (
        [ncut_xist, ncut_leiden, ncut_kahip_val, ncut_metis_val, ncut_xcut_val, ncut_truth_val], 
        [rate_xist, rate_leiden, rate_kahip, rate_metis, rate_xcut, rate_truth]  
    )


# spectral and chaco clustering without parallizing
def mixed_gaussian_alg_test_spectral_chaco(N_iter, BigN, delta_lst, k, r, sigma, seeds_by_delta):
    gauss_rate_trafo, gauss_ncut_trafo = [], []
    gauss_rate_trafo_mean, gauss_ncut_trafo_mean = [], []
    algnames = ["SpecClust", "Chaco"] 
    #for delta in tqdm(delta_lst):
    for delta_idx, delta in enumerate(tqdm(delta_lst)):  
        seeds= seeds_by_delta[delta_idx]
        iter_rate = [[] for name in algnames]
        iter_ncut = [[] for name in algnames]
        #for i in range(N_iter):
        for seed in seeds:
            np.random.seed(int(seed))
            mge_df, mge_truth, _ = generate_mixed_gaussian_edges(BigN, delta, k, r, sigma, weighted=True)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mge_spec = spectral_clustering(mge_df)
            iter_rate[0].append(classification_rate(mge_spec[2], mge_truth))
            iter_ncut[0].append(mge_spec[0])


            mge_chaco = ncut_chaco(mge_df, "mixed_gaussian", print_chaco_output=False)
            iter_rate[1].append(classification_rate(mge_chaco[2], mge_truth))
            iter_ncut[1].append(mge_chaco[0])
        
        gauss_ncut_trafo.extend(([delta, np.nanmedian(np.array(iter_ncut[i])[np.isfinite(iter_ncut[i])]), algnames[i]] for i in range(len(algnames))))
        gauss_rate_trafo.extend(([delta, np.nanmedian(np.array(iter_rate[i])[np.isfinite(iter_rate[i])]), algnames[i]] for i in range(len(algnames))))
        gauss_ncut_trafo_mean.extend(([delta, np.nanmean(np.array(iter_ncut[i])[np.isfinite(iter_ncut[i])]), algnames[i]] for i in range(len(algnames))))
        gauss_rate_trafo_mean.extend(([delta, np.nanmean(np.array(iter_rate[i])[np.isfinite(iter_rate[i])]), algnames[i]] for i in range(len(algnames))))

    return pd.DataFrame(gauss_rate_trafo, columns=["x", "values", "alg"]), pd.DataFrame(gauss_ncut_trafo, columns=["x", "values", "alg"]), pd.DataFrame(gauss_rate_trafo_mean, columns=["x", "values", "alg"]), pd.DataFrame(gauss_ncut_trafo_mean, columns=["x", "values", "alg"])


def mixed_gaussian_alg_test_weighted(N_iter, BigN, delta_lst, k, r, sigma, num_cpus=25):
    gauss_rate_trafo, gauss_ncut_trafo = [], []
    gauss_rate_trafo_mean, gauss_ncut_trafo_mean = [], []
    algnames = ["Xist", "Leiden", "KaHIP", "Metis", "Xcut", "true_results"]
    rng = np.random.default_rng(20000)
    seeds_by_delta = [
    rng.integers(
        0,
        2**32 - 1,
        size=N_iter,
        dtype=np.uint32
    )
    for _ in delta_lst
    ]
    for delta_idx, delta in enumerate(tqdm(delta_lst)):
    #for delta in tqdm(delta_lst):
        # Run N_iter independent computations in parallel
        #rng = np.random.default_rng()
        #seeds = rng.integers(low=0, high=2**32 - 1, size=N_iter, dtype=np.uint32)
        seeds= seeds_by_delta[delta_idx]
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            results = list(
                executor.map(
                    _single_iteration_mixed_gaussian_weighted,
                    [BigN] * N_iter,
                    [delta] * N_iter,
                    [k] * N_iter,
                    [r] * N_iter,
                    [sigma] * N_iter,
                    seeds.tolist(),
                    [f"mixed_gaussian_{i}" for i in range(N_iter)],
                )
            )

        # Separate ncut and rate lists
        all_ncut = [res[0] for res in results]   # list of length 6
        all_rate = [res[1] for res in results]   # list of length 6


        all_ncut = np.array(all_ncut)   # shape (N_iter, 6)
        all_rate = np.array(all_rate)   # shape (N_iter, 6)

        for i, alg in enumerate(algnames):
            gauss_ncut_trafo.append([delta, np.nanmedian(all_ncut[:, i]), alg])
            gauss_rate_trafo.append([delta, np.nanmedian(all_rate[:, i]), alg])
            gauss_ncut_trafo_mean.append([delta, np.nanmean(all_ncut[:, i]), alg])
            gauss_rate_trafo_mean.append([delta, np.nanmean(all_rate[:, i]), alg])



    # add spec_clust
    df_rate_spec, df_ncut_spec, df_rate_spec_mean, df_ncut_spec_mean= mixed_gaussian_alg_test_spectral_chaco(N_iter, BigN, delta_lst, k, r, sigma, seeds_by_delta)

    df_rate_all_but_spec = pd.DataFrame(gauss_rate_trafo, columns=["x", "values", "alg"])
    df_ncut_all_but_spec = pd.DataFrame(gauss_ncut_trafo, columns=["x", "values", "alg"])
    df_rate_all_but_spec_mean = pd.DataFrame(gauss_rate_trafo_mean, columns=["x", "values", "alg"])
    df_ncut_all_but_spec_mean = pd.DataFrame(gauss_ncut_trafo_mean, columns=["x", "values", "alg"])

#    return df_rate_all_but_spec, df_ncut_all_but_spec, df_rate_all_but_spec_mean, df_ncut_all_but_spec_mean

    return pd.concat([df_rate_spec, df_rate_all_but_spec], ignore_index=True), pd.concat([df_ncut_spec, df_ncut_all_but_spec], ignore_index=True), pd.concat([df_rate_spec_mean, df_rate_all_but_spec_mean], ignore_index=True), pd.concat([df_ncut_spec_mean, df_ncut_all_but_spec_mean], ignore_index=True)



#num_cpus = 40
gauss_rate_trafo, gauss_ncut_trafo, gauss_rate_trafo_mean, gauss_ncut_trafo_mean  = mixed_gaussian_alg_test_unweighted(100, 100, [x/10 for x in range(20,51)], 5, 0.2, 0.2, num_cpus=40)
#print("Rate:\n", gauss_rate_trafo, "\nNCut:\n", gauss_ncut_trafo)
print("Rate mean:\n", gauss_rate_trafo_mean, "\nNCut mean:\n", gauss_ncut_trafo_mean)
#gauss_rate_trafo.to_csv("gauss_rate_trafo_100_unweighted.csv", index=False)
#gauss_ncut_trafo.to_csv("gauss_ncut_trafo_100_unweighted.csv", index=False)
gauss_rate_trafo_mean.to_csv("gauss_rate_trafo_100_unweighted_mean.csv", index=False)
gauss_ncut_trafo_mean.to_csv("gauss_ncut_trafo_100_unweighted_mean.csv", index=False)
  
gauss_rate_trafo, gauss_ncut_trafo, gauss_rate_trafo_mean, gauss_ncut_trafo_mean = mixed_gaussian_alg_test_weighted(100, 100, [x/10 for x in range(20,51)], 5, 0.2, 0.2, num_cpus=40)
#print("Rate:\n", gauss_rate_trafo, "\nNCut:\n", gauss_ncut_trafo)
print("Rate mean:\n", gauss_rate_trafo_mean, "\nNCut mean :\n", gauss_ncut_trafo_mean)
#gauss_rate_trafo.to_csv("gauss_rate_trafo_100_weighted.csv", index=False)
#gauss_ncut_trafo.to_csv("gauss_ncut_trafo_100_weighted.csv", index=False)
gauss_rate_trafo_mean.to_csv("gauss_rate_trafo_100_weighted_mean.csv", index=False)
gauss_ncut_trafo_mean.to_csv("gauss_ncut_trafo_100_weighted_mean.csv", index=False)




#### Comparing the runtime of Xist to that of other algorithms on discretized NIH3T3 cell images (again weighted and unweighted)

# Runs all algorithms for a single (m, image index) and returns: list of [m, time, alg] entries
def _single_image_runtime_unweighted(m, imgidx, t, sigm, weights_by_sample):
    mge_df, _ = generate_tiff_edges(
        f"NIH3T3_Data/s_C001Tubulin{imgidx}.tif",
        m=m, t=t, sigm=sigm, weights_by_sample=weights_by_sample
    )

    out = []

    mge_xist = xist_dinic_faster_wrapper(mge_df)
    #mge_xist = xist(mge_df)
    out.append([m, mge_xist[1], "Xist"])

    mge_scoreplus = ncut_scoreplus(mge_df)
    out.append([m, mge_scoreplus[1], "Scoreplus"])

    mge_xvst = xvst_dinic_faster_wrapper(mge_df)
    #mge_xvst = xvst(mge_df)
    out.append([m, mge_xvst[1], "Xvst"])


    return out


def algs_runtime_comparison_parallel_unweighted(imgindices, mlist=[8,9,12,14,18,21,24,28], t=1, sigm=math.nan, weights_by_sample=True, num_cpus=25):
    results_all = []
    for m in tqdm(mlist):
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            # Each result is a list of 4 entries
            out = list(
                executor.map(
                    _single_image_runtime_unweighted,
                    [m] * len(imgindices),
                    imgindices,
                    [t] * len(imgindices),
                    [sigm] * len(imgindices),
                    [weights_by_sample] * len(imgindices),
                )
            )

        # Flatten
        for entry in out:
            results_all.extend(entry)


    return pd.DataFrame(results_all, columns=["m", "time", "alg"])





def _single_image_runtime_weighted(m, imgidx, t, sigm, weights_by_sample):

    mge_df, _ = generate_tiff_edges(
        f"NIH3T3_Data/s_C001Tubulin{imgidx}.tif",
        m=m, t=t, sigm=sigm, weights_by_sample=weights_by_sample
    )

    out = []

    mge_xist = xist_dinic_faster_wrapper(mge_df)
    #mge_xist = xist(mge_df)
    out.append([m, mge_xist[1], "Xist"])

    mge_leiden = leidenoracle(mge_df, exponential_resolution_scaling=True)
    out.append([m, mge_leiden[0][3], "Leiden"])

    mge_kahip = ncut_kahip(mge_df)
    out.append([m, mge_kahip[1], "KaHIP"])

    mge_metis = ncut_metis(mge_df)
    out.append([m, mge_metis[1], "METIS"])

    name_xcut= "nih3t3_" + str(imgidx) + "_" + str(m)
    mge_xcut = ncut_xcut(mge_df, name_xcut, print_xcut_output=False)
    out.append([m, mge_xcut[1], "Xcut"])

    mge_xvst = xvst_dinic_faster_wrapper(mge_df)
    #mge_xvst = xvst(mge_df)
    out.append([m, mge_xvst[1], "Xvst"])

    return out


def algs_runtime_comparison_parallel_weighted(imgindices, mlist=[8,9,12,14,18,21,24,28], t=1, sigm=math.nan, weights_by_sample=True, num_cpus=1):
    results_all = []
    for m in tqdm(mlist):
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            out = list(
                executor.map(
                    _single_image_runtime_weighted,
                    [m] * len(imgindices),
                    imgindices,
                    [t] * len(imgindices),
                    [sigm] * len(imgindices),
                    [weights_by_sample] * len(imgindices),
                )
            )

        # Flatten
        for entry in out:
            results_all.extend(entry)
            
    # spectral_clustering and chaco without parallizing: 
    runtimes = []
    for m in mlist:
        for i in tqdm(range(len(imgindices))):
            imgidx= imgindices[i]
            mge_df, _ = generate_tiff_edges(f"NIH3T3_Data/s_C001Tubulin{imgidx}.tif", m=m, t=t, sigm=sigm, weights_by_sample=weights_by_sample)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mge_spec = spectral_clustering(mge_df)
            runtimes.append([m, mge_spec[1], "SpecClust"])
            mge_chaco = ncut_chaco(mge_df, "nih3t3", print_chaco_output=False)
            runtimes.append([m, mge_chaco[3], "Chaco"])
    df_spec= pd.DataFrame(runtimes, columns=["m", "time", "alg"])
    df_all_but_spec= pd.DataFrame(results_all, columns=["m", "time", "alg"])

    return pd.concat([df_spec, df_all_but_spec], ignore_index=True)


nih3t3_runtimes_unweighted = algs_runtime_comparison_parallel_unweighted( [i for i in range(1, 22) if i != 17], weights_by_sample=False, num_cpus=20)   
print(nih3t3_runtimes_unweighted)
nih3t3_runtimes_unweighted.to_csv("nih3t3_timecomp_unweighted.csv", index=False)

nih3t3_runtimes_weighted = algs_runtime_comparison_parallel_weighted( [i for i in range(1, 22) if i != 17], weights_by_sample=True, num_cpus=2)
print(nih3t3_runtimes_weighted)
nih3t3_runtimes_weighted.to_csv("nih3t3_timecomp_weighted.csv", index=False)



#### Applying the algorithms to large network datasets:

# functions for modifying the datasets:
def remove_top_k_percent_hubs(G, k=0.002):
    degrees = dict(G.degree())
    n_remove = int(len(degrees) * k)
    #print(n_remove)

    # sort nodes by degree descending
    top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:n_remove]

    G_new = G.copy()
    G_new.remove_nodes_from(top_nodes)

    return G_new



def graph_to_edge_dataframe(G, weight_attr='weight'):
    rows = []

    for u, v, data in G.edges(data=True):
        weight = data.get(weight_attr, 1)  
        rows.append((u, v, weight))

    df = pd.DataFrame(rows, columns=['node1', 'node2', 'weight'])
    return df

def make_dataset_star(dataset,k=0.002):
    G = nx.from_pandas_edgelist(dataset,source=0,target=1,edge_attr=2,create_using=nx.Graph())
    G_no_hubs = remove_top_k_percent_hubs(G, k=k)
    no_hubs= graph_to_edge_dataframe(G_no_hubs)
    no_hubs_lcc= reduce_to_lcc(no_hubs)
    return no_hubs_lcc



musae_squirrel = reduce_to_lcc(pd.read_csv('Datasets/musae_squirrel_edges.csv'))
ca_hepph = reduce_to_lcc(pd.read_csv('Datasets/CA-HepPh.txt', sep='\t', header=3))
musae_facebook = reduce_to_lcc(pd.read_csv('Datasets/musae_facebook_edges.csv'))
email_enron = reduce_to_lcc(pd.read_csv('Datasets/Email-Enron.txt', sep='\t', header=3))
artists = reduce_to_lcc(pd.read_csv('Datasets/artist_edges.csv'))
twitch_gamers = reduce_to_lcc(pd.read_csv('Datasets/large_twitch_edges.csv'))
twitch_star = make_dataset_star(twitch_gamers, k=0.002)

datasets = [musae_squirrel, ca_hepph, musae_facebook, email_enron, artists, twitch_gamers, twitch_star]
dataset_names = ["musae_squirrel", "ca_hepph", "musae_facebook", "email_enron", "artists", "twitch_gamers", "twitch_star"]



for i in range(len(datasets)): 
    print("Xist NCut for", dataset_names[i], "is given by", xist_dinic_faster_wrapper(datasets[i])[0:2])
    print("Leiden Oracle NCut for", dataset_names[i], "is given by", leidenoracle(datasets[i])[0:2])
    print("KaHIP NCut for", dataset_names[i], "is given by", ncut_kahip(datasets[i])[0:2])
    print("METIS NCut for", dataset_names[i], "is given by", ncut_metis(datasets[i])[0:2])
    df_chaco = ncut_chaco(datasets[i], dataset_names[i], print_chaco_output=False)
    print("Best Chaco NCut for", dataset_names[i], "is given by", df_chaco[0], "with imbalance", df_chaco[1], "in total time", df_chaco[3])
    print("Xcut NCut for", dataset_names[i], "is given by", ncut_xcut(datasets[i], dataset_names[i])[0:2])
    if dataset_names[i] != "twitch_gamers" and dataset_names[i] != "twitch_star":
        print("Scoreplus NCut for", dataset_names[i], "is given by", ncut_scoreplus(datasets[i])[0:2])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dataset_spec = spectral_clustering(datasets[i])
        print("SpecClust NCut for", dataset_names[i], "is given by", dataset_spec[0:2])
    print("-----")


