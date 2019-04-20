# Simultaneous Dimensionality and Complexity Model Selection for Spectral Graph Clustering

Congyuan Yang, Carey E. Priebe, Youngser Park, David J. Marchette

https://arxiv.org/abs/1904.02926

## Abstract

Our problem of interest is to cluster vertices of a graph by identifying its underlying community structure. Among various vertex clustering approaches, spectral clustering is one of the most popular methods, because it is easy to implement while often outperforming traditional clustering algorithms. However, there are two inherent model selection problems in spectral clustering, namely estimating the embedding dimension and number of clusters. This paper attempts to address the issue by establishing a novel model selection framework specifically for vertex clustering on graphs under a stochastic block model. The first contribution is a probabilistic model which approximates the distribution of the extended spectral embedding of a graph. The model is constructed based on a theoretical result of asymptotic normality of the informative part of the embedding, and on a simulation result of limiting behavior of the redundant part of the embedding. The second contribution is a simultaneous model selection framework. In contrast with the traditional approaches, our model selection procedure estimates embedding dimension and number of clusters simultaneously. Based on our proposed distributional model, a theorem on the consistency of the estimates of model parameters is stated and proven. The theorem provides a statistical support for the validity of our method. Heuristic algorithms via the simultaneous model selection framework for vertex clustering are proposed, with good performance shown in the experiment on synthetic data and on the real application of connectome analysis.

## Code and Results

All the experiments are written in `R`, and here is the instruction for getting the results in the paper. 

* To run the synthetic data simulation, please follow these steps:

	1. (optional) Run `main_simulation.R` 600 times with different parameter pair `(i, mc)`, where `i` ranges from 1 to 6 and `mc` ranges from 1 to 100. `batch_simulation.sh` is a script for running these jobs in our cluster. Each run produces a `.RData` file as a partial result (600 files total). The precomputed `.RData` files are stored in Google Drive and loaded inside `loadData.R`.

	2. Run `plot_simulation.R` to get Figure 5, 6, 7 in the paper. The results are shown in [here](http://www.cis.jhu.edu/~parky/dhatKhat/plot_simulation.html)

* To run the connectome data experiment, please follow these steps:

	1. (optional) Run `main_DS01216.R` 114 times with different parameter `fileIndex`, which  ranges from 1 to 114. `batch_realdata.sh` is a script for running these jobs in our cluster. Each run will produces a `.RData` file as a partial result (114 files total). The precomputed `.RData` files are also stored in Google Drive and loaded inside `loadData.R`.

	2. Run `plot_DS01216.R` to get Figure 8, 9 and the results of Table 1 in the paper.  The results are shown in [here](http://www.cis.jhu.edu/~parky/dhatKhat/plot_DS01216.html).
