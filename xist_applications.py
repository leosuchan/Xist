from xist import *
from tqdm import tqdm
import warnings

#### Applying Xist to a dataset:

artists = reduce_to_lcc(pd.read_csv('Datasets/artist_edges.csv'))
artists_xist = xist(artists)
print("Xist (NCut) computed a value of", artists_xist[0], "for the artists dataset in", artists_xist[1],"seconds.")


#### NCut value and classification rate comparison between Xist and other algorithms on a Gaussian mixture:

def generate_mixed_gaussian_edges(BigN, delta, k, r, sigma):
    Nsize = np.random.binomial(n=BigN, p=1/2)
    smpl = np.concatenate((np.random.multivariate_normal([0,0], [[1,0],[0,1]], Nsize), np.random.multivariate_normal([delta,delta], [[1,0],[0,1]], BigN-Nsize)))
    kNN_model = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(smpl)
    dstncs, indices = kNN_model.kneighbors(smpl)
    kNN_edges = [[i,indices[i][j],math.exp(-math.sqrt(sum((smpl[i] - smpl[indices[i][j]])**2))/sigma)] for j in range(1,k+1) for i in range(BigN)]
    exp_edges = [[i,j,math.exp(-math.sqrt(sum((smpl[i] - smpl[j])**2))/sigma)] for i in range(BigN-1) for j in range(i+1,BigN) if sum((smpl[i] - smpl[j])**2) <= r]
    mixed_gaussian_edges, mg_assign = reduce_to_lcc(pd.DataFrame(kNN_edges + exp_edges), preserve=[0] * Nsize + [1] * (BigN - Nsize), reduce=False)
    return mixed_gaussian_edges, mg_assign, smpl

def mixed_gaussian_alg_test(N_iter, BigN, delta_lst, k, r, sigma):
    gauss_rate_trafo, gauss_ncut_trafo = [], []
    algnames = ["Xist", "Leiden (oracle)", "KaHIP", "METIS", "Chaco", "SpecClust"]
    for delta in tqdm(delta_lst):
        iter_rate = [[] for name in algnames]
        iter_ncut = [[] for name in algnames]
        for i in range(N_iter):
            mge_df, mge_truth, _ = generate_mixed_gaussian_edges(BigN, delta, k, r, sigma)
            mge_xist = xist(mge_df)
            iter_rate[0].append(classification_rate(mge_xist[2], mge_truth))
            iter_ncut[0].append(mge_xist[0])

            mge_leiden = leidenoracle(mge_df, exponential_resolution_scaling=True, classif_truth=mge_truth)
            iter_rate[1].append(mge_leiden[1][0])
            iter_ncut[1].append(mge_leiden[0][0])

            mge_kahip = ncut_kahip(mge_df)
            iter_rate[2].append(classification_rate(mge_kahip[2], mge_truth))
            iter_ncut[2].append(mge_kahip[0])

            mge_metis = ncut_metis(mge_df)
            iter_rate[3].append(classification_rate(mge_metis[2], mge_truth))
            iter_ncut[3].append(mge_metis[0])

            mge_chaco = ncut_chaco(mge_df, "mixed_gaussian", print_chaco_output=False)
            iter_rate[4].append(classification_rate(mge_chaco[2], mge_truth))
            iter_ncut[4].append(mge_chaco[0])

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mge_spec = spectral_clustering(mge_df)
            iter_rate[5].append(classification_rate(mge_spec[2], mge_truth))
            iter_ncut[5].append(mge_spec[0])

        gauss_ncut_trafo.extend(([delta, np.nanmedian(iter_ncut[i]), algnames[i]] for i in range(len(algnames))))
        gauss_rate_trafo.extend(([delta, np.nanmedian(iter_rate[i]), algnames[i]] for i in range(len(algnames))))

    return pd.DataFrame(gauss_rate_trafo, columns=["x", "values", "alg"]), pd.DataFrame(gauss_ncut_trafo, columns=["x", "values", "alg"])
# This method might result in a segmentation fault. The reason for that is Chaco. If the problem persist, disable the Chaco computation.
# Also, don't mind the timeouts and errors; these are also only produced by Chaco and are caught and automatically ignored.

gauss_rate_trafo, gauss_ncut_trafo = mixed_gaussian_alg_test(100, 100, [x/10 for x in range(10,51)], 5, 0.2, 0.2)
print("Rate:\n", gauss_rate_trafo, "\nNCut:\n", gauss_ncut_trafo)
gauss_rate_trafo.to_csv("gauss_rate_trafo_100.csv", index=False)
gauss_ncut_trafo.to_csv("gauss_ncut_trafo_100.csv", index=False)


#### Comparing the runtime of Xist to that of other algorithms on discretized NIH3T3 cell images

def algs_runtime_comparison(imgindices, mlist=[8,9,12,14,18,21,24,28], t=1, sigm=math.nan, weights_by_sample=True):
    runtimes = []
    for m in mlist:
        for i in tqdm(range(len(imgindices))):
            mge_df, _ = generate_tiff_edges("NIH3T3_Data/s_C001Tubulin{0}.tif".format(imgindices[i]), m=m, t=t, sigm=sigm, weights_by_sample=weights_by_sample)

            mge_xist = xist(mge_df)
            runtimes.append([m, mge_xist[1], "Xist"])

            mge_leiden = leidenoracle(mge_df, exponential_resolution_scaling=True)
            runtimes.append([m, mge_leiden[0][3], "Leiden"])

            mge_kahip = ncut_kahip(mge_df)
            runtimes.append([m, mge_kahip[1], "KaHIP"])

            mge_metis = ncut_metis(mge_df)
            runtimes.append([m, mge_metis[1], "METIS"])

            mge_chaco = ncut_chaco(mge_df, "nih3t3", print_chaco_output=False)
            runtimes.append([m, mge_chaco[3], "Chaco"])

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mge_spec = spectral_clustering(mge_df)
            runtimes.append([m, mge_spec[1], "SpecClust"])

            mge_xvst = xvst(mge_df)
            runtimes.append([m, mge_xvst[1], "Xvst"])

    return pd.DataFrame(runtimes, columns=["m", "time", "alg"])

nih3t3_runtimes = algs_runtime_comparison([i for i in range(1,22) if i!=17])
print(nih3t3_runtimes)
nih3t3_runtimes.to_csv("nih3t3_timecomp.csv", index=False)


#### Applying the algorithms to large network datasets:

musae_squirrel = reduce_to_lcc(pd.read_csv('Datasets/musae_squirrel_edges.csv'))
ca_hepph = reduce_to_lcc(pd.read_csv('Datasets/CA-HepPh.txt', sep='\t', header=3))
musae_facebook = reduce_to_lcc(pd.read_csv('Datasets/musae_facebook_edges.csv'))
email_enron = reduce_to_lcc(pd.read_csv('Datasets/Email-Enron.txt', sep='\t', header=3))
artists = reduce_to_lcc(pd.read_csv('Datasets/artist_edges.csv'))
twitch_gamers = reduce_to_lcc(pd.read_csv('Datasets/large_twitch_edges.csv'))

datasets = [musae_squirrel, ca_hepph, musae_facebook, email_enron, artists, twitch_gamers]
dataset_names = ["musae_squirrel", "ca_hepph", "musae_facebook", "email_enron", "artists", "twitch_gamers"]

for i in range(len(datasets)):
    print("Xist NCut for", dataset_names[i], "is given by", xist(datasets[i])[0:2])
    print("Leiden Oracle NCut for", dataset_names[i], "is given by", leidenoracle(datasets[i])[0:2])
    print("KaHIP NCut for", dataset_names[i], "is given by", ncut_kahip(datasets[i])[0:2])
    print("METIS NCut for", dataset_names[i], "is given by", ncut_metis(datasets[i])[0:2])
    df_chaco = ncut_chaco(datasets[i], dataset_names[i], print_chaco_output=False)
    print("Best Chaco NCut for", dataset_names[i], "is given by", df_chaco[0], "with imbalance", df_chaco[1])
    print("-----")