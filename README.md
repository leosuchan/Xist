# Content

This repository contains five related sets of files:

1. Python, R and C++ code supplementary to the paper "Xist: Scalable graph cut clustering with statistical guarantees" (Li, Munk, Suchan, Kratz; 2023), see [here on arXiv](https://arxiv.org/abs/2308.09613):
   - `xist.py`
   - `xist_application.py`
   - `bash_chaco.sh`
   - `figures.ipynb`

2. R code supplementary to the same paper:
   - `xist.R`
   - `xist_applications.R`

3. C++ code supplementary to the same paper:
    - `xist_dinic_faster.cpp`
    - `xvst.cpp`

4. R code supplementary to the paper "Distributional limits of graph cuts on discretized grids" (Suchan, Li, Munk; 2024), to appear on arXiv:
   - `graph_cut_limits.R`
   - `graph_cut_limit_applications.R`

5. The NIH 3T3 dataset:
   - The folder `NIH3T3_Data` containing 21 mouse embryo stem cell images. This data belongs to Ulrike Rölleke and Sarah Köster (University of Göttingen).

In particular, the Xist algorithm is implemented in `xist.py` for Python as well as in `xist.R` for R. It is further implemented in C++ in `xist_dinic_faster.cpp' with a Python wrapper that can be found in `xist.py`


# Usage

## 1. Installation of the Python code supplement to "A scalable algorithm to approximate graph cuts"

1. Install KaHIP for Python (https://github.com/KaHIP/KaHIP - follow the installation instructions in their README under section "Using KaHIP in Python", we used pip install kahip
2. Install the Chaco algorithm (https://www3.cs.stonybrook.edu/~algorith/implement/chaco/implement.shtml)
3. Download `xist.py`, `xist_application.py`, `bash_chaco.sh`, `xist_dinic_faster.cpp` and `xvst.cpp` from this repository. Put them in a directory of your choosing and put `bash_chaco.sh` into the `exec` folder of your Chaco installation. Compile the C++ files.
4. Edit the path fragments `/home/kratz10/project_xist/Chaco/Chaco-2.2` inside the functions `ncut_chaco_unweighted` and `ncut_chaco` in `xist.py` to point towards your Chaco installation directory.
5. If your Chaco installation directory is not `~/Chaco-2.2/`, change this expression in lines 5 and 6 of `bash_chaco.sh` so that it points towards your Chaco installation directory.
6. Install xcut ( https://gitlab.com/vietaa/xcut) and follow the instructions to in their README. We used the preset "gcc- release"
7. Edit the path fragment `/home/kratz10/project_xist/xcut/build/`  in "ncut_xcut" to point towards your xcut installations and the `/home/kratz10/project_xist/code/helpdata` to point towards your helpdata folder.
8. Download Scoreplus (https://cran.r-project.org/src/contrib/Archive/ScorePlus/). Put the file SCOREplus.R also into the same directory as `xist.py`. Install the necessary R packages (igraph, Rspectral)
9. Install the following Python packages via pip:
   - `numpy`
   - `igraph`
   - `pandas`
   - `PIL`
   - `leidenalg`
   - `pymetis`
   - `tqdm`
   - `skicit-learn`
   - `rpy2`
   - `concurrent.futures`
   - `networkx`
   - `collections.abc`
   - `timeit`
   - `math`
   - `csv`
   - `subprocess`
   - `os.path` and `os`
   - `scipy,sparse`
   - `sys`
10. Run `xist.py`.
11. (Optional) If you desire to work with the NIH 3T3 Dataset, download the `NIH3T3_Data` folder from this repository and place it into your Python working directory.
12. (Optional) If you desire to work with the large datasets used in the paper, download them from [the SNAP database](https://snap.stanford.edu/data/). Create a `Datasets` folder inside the `deploy` folder of your KaHIP installation, and put the following files into it:
   - `musae_squirrel_edges.csv` (from https://snap.stanford.edu/data/wikipedia-article-networks.html)
   - `CA-HepPh.txt` (from https://snap.stanford.edu/data/cit-HepPh.html)
   - `musae_facebook_edges.csv` (from https://snap.stanford.edu/data/facebook-large-page-page-network.html)
   - `Email-Enron.txt` (from https://snap.stanford.edu/data/email-Enron.html)
   - `artist_edges.csv` (from https://snap.stanford.edu/data/gemsec-Facebook.html)
   - `large_twitch_edges.csv` (from https://snap.stanford.edu/data/twitch_gamers.html)
13. Done! You are now ready to run any part of `xist_application.py` and should therefore be able to reproduce the results from "Xist: Scalable graph cut clustering with statistical guarantees" (Li, Munk, Suchan, Kratz; 2023).


## 2. Usage of the R code supplement to "A scalable algorithm to approximate graph cuts"

1. Download `xist.R` and `xist_applications.R` from this repository.
2. Run `xist.R`.
3. (Optional) If you desire to work with the NIH 3T3 Dataset, download the `NIH3T3_Data` folder from this repository and place it into your Python working directory. Edit `xist_applications.R` to replace the two occurences of `setwd("/path/to/NIH3T3/data")` by the appropriate path to the NIH 3T3 dataset folder.
4. (Optional) Similarly, if you desire to work with the large datasets used in the paper, download them from [the SNAP database](https://snap.stanford.edu/data/), putting the files listed above in section 1.9. into a folder. Then replace `setwd("/path/to/SNAP/data")` in `xist_applications.R` with the path to your newly created folder.
5. Done! You should now be able to run any part of `graph_cut_limit_applications.R`.

All the time comparisons and algorithm comparisons in "A scalable algorithm to approximate graph cuts" (Suchan, Li, Munk; 2023) have been done using Python and are not present in the R supplement. This is because some SOTA algorithms used in the paper, namely Chaco, KaHIP, and METIS, do not have an R implementation.

Notice that the entirety of `xist.R` is heavily commented to aid the user. Please read through the comments if the use of some functions is not immediately obvious.


## 3. Usage of the R code supplementary to "Distributional limits of graph cuts on discretized grids"

1. Download `graph_cut_limits.R` and `graph_cut_limit_applications.R` from this repository.
2. Run `graph_cut_limits.R`.
3. Done! You are now ready to run any part of `graph_cut_limit_applications.R` and should therefore be able to reproduce the results from "Distributional limits of graph cuts on discretized grids" (Suchan, Li, Munk; 2024).


## 4. Citing the NIH 3T3 Dataset

The NIH 3T3 Dataset can be found in the folder `NIH3T3_Data` in this repository. It should be cited using `NIH3T3_Data/CITATION.cff`
