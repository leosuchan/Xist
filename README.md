# Content

This repository contains four related sets of files:

1. Python code supplementary to the paper "A scalable algorithm to approximate graph cuts" (Suchan, Li, Munk; 2023), see [here on arXiv](https://arxiv.org/abs/2308.09613):
 - `xist.py`
 - `xist_applications.py`
 - `bash_chaco.sh`

2. R code supplementary to the same paper:
 - `xist.R` (Contains the main functions, but not all application examples from the paper.)

3. R code supplementary to the paper "Distributional limits of graph cuts on discretized grids" (Suchan, Li, Munk; 2024), to appear on arXiv:
 - `graph_cut_limits.R`
 - `graph_cut_limit_applications.R`

4. The NIH 3T3 dataset:
 - The folder `NIH3T3_Data` containing 21 mouse embryo stem cell images. This data belongs to Ulrike Rölleke and Sarah Köster (University of Göttingen).

In particular, the Xist algorithm is implemented in `xist.py` for Python as well as in `xist.R` for R.


# Usage

## 1. Installation of the Python code supplement to "A scalable algorithm to approximate graph cuts" for Linux-based systems

1. Install KaHIP for Python (https://github.com/KaHIP/KaHIP - follow the installation instructions in their README under section "Using KaHIP in Python")
2. Install the Chaco algorithm (https://www3.cs.stonybrook.edu/~algorith/implement/chaco/implement.shtml)
3. Download `xist.py`, `xist_applications.py` and `bash_chaco.sh` from this repository. Put both Python files into the `deploy` folder of your KaHIP installation, and put `bash_chaco.sh` into the `exec` folder of your Chaco installation.
4. Edit the path fragments `/home/lsuchan/Chaco-2.2` inside the functions `ncut_chaco_unweighted` and `ncut_chaco` in `xist.py` to point towards your Chaco installation directory.
5. If your Chaco installation directory is not `~/Chaco-2.2/`, change this expression in lines 5 and 6 of `bash_chaco.sh` so that it points towards your Chaco installation directory.
6. Install the following Python packages via pip:
 - `numpy`
 - `igraph`
 - `pandas`
 - `PIL`
 - `leidenalg`
 - `pymetis`
 - `tqdm`
7. Run `xist.py`
8. (Optional) If you desire to work with the NIH 3T3 Dataset, download the `NIH3T3_Data` folder from this repository and place it into your Python working directory.
9. (Optional) If you desire to work with the large datasets used in the paper, download them from [the SNAP database](https://snap.stanford.edu/data/). Create a `Datasets` folder inside the `deploy` folder of your KaHIP installation, and put the following files into it:
 - `musae_squirrel_edges.csv` (from https://snap.stanford.edu/data/wikipedia-article-networks.html)
 - `CA-HepPh.txt` (from https://snap.stanford.edu/data/cit-HepPh.html)
 - `musae_facebook_edges.csv` (from https://snap.stanford.edu/data/facebook-large-page-page-network.html)
 - `Email-Enron.txt` (from https://snap.stanford.edu/data/email-Enron.html)
 - `artist_edges.csv` (from https://snap.stanford.edu/data/gemsec-Facebook.html)
 - `large_twitch_edges.csv` (from https://snap.stanford.edu/data/twitch_gamers.html)
10. Done! You are now ready to run any part of `xist_applications.py` and should therefore be able to reproduce the results from "A scalable algorithm to approximate graph cuts" (Suchan, Li, Munk; 2023).


## 2. Usage of the R code supplement to "A scalable algorithm to approximate graph cuts"

1. Download `xist.R` from this repository.
2. (Optional) If you desire to work with the NIH 3T3 Dataset, download the `NIH3T3_Data` folder from this repository and place it into your Python working directory. Edit `xist.R` to replace `setwd("/path/to/NIH3T3/data")` (around line 750 and again around line 820) by the appropriate path to the NIH 3T3 dataset folder.
3. (Optional) Similarly, if you desire to work with the large datasets used in the paper, download them from [the SNAP database](https://snap.stanford.edu/data/), putting the files listed above in section 1.9. into a folder. Then replace `setwd("/path/to/SNAP/data")` in `xist.R` (around line 790) with the path to your newly created folder.
2. Run the first part of the R file, i.e. everything until `### PART II: APPLICATION EXAMPLES ###` (appears roughly around line 1350).
3. Done! You should now be able to run the application examples found below the aforementioned line.

The entirety of `xist.R` is heavily commented to aid the user. Please read through the comments if the use of some functions is not immediately obvious.


## 3. Usage of the R code supplementary to "Distributional limits of graph cuts on discretized grids"

1. Download `graph_cut_limits.R` and `graph_cut_limit_applications.R` from this repository.
2. Run `graph_cut_limits.R`.
3. Done! You are now ready to run any part of `graph_cut_limit_applications.R` and should therefore be able to reproduce the results from "Distributional limits of graph cuts on discretized grids" (Suchan, Li, Munk; 2024).


## 4. Citing the NIH 3T3 Dataset

The NIH 3T3 Dataset can be found in the folder `NIH3T3_Data` in this repository. It should be cited using `NIH3T3_Data/CITATION.cff`
