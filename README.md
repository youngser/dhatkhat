# Simultaneous Dimensionality and Complexity Model Selection for Spectral Graph Clustering

Congyuan Yang, Carey E. Priebe, Youngser Park, David J. Marchette

https://arxiv.org/abs/1904.02926

## Abstract

Our problem of interest is to cluster vertices of a graph by identifying underlying community structure. Among various vertex clustering approaches, spectral clustering is one of the most popular methods because it is easy to implement while often outperforming more traditional clustering algorithms. However, there are two inherent model selection problems in spectral clustering, namely estimating both the embedding dimension and number of clusters. This paper attempts to address the issue by establishing a novel model selection framework specifically for vertex clustering on graphs under a stochastic block model. The first contribution is a probabilistic model which approximates the distribution of the extended spectral embedding of a graph. The model is constructed based on a theoretical result of asymptotic normality of the informative part of the embedding, and on a simulation result providing a conjecture for the limiting behavior of the redundant part of the embedding. The second contribution is a simultaneous model selection framework. In contrast with the traditional approaches, our model selection procedure estimates embedding dimension and number of clusters simultaneously. Based on our conjectured distributional model, a theorem on the consistency of the estimates of model parameters is presented, providing support for the validity of our method. Algorithms for our simultaneous model selection for vertex clustering are proposed, demonstrating superior performance in simulation experiments. We illustrate our method via application to a collection of brain graphs.

## Code and Results

All the experiments are written in `R`, and here is the instruction for getting the results in the paper. 

### Demo

To run a simple working live demo, please run these lines
```
if (!require(RCurl)) install.packages("RCurl")
library(RCurl)

script <- getURL("https://raw.githubusercontent.com/youngser/dhatkhat/master/R/demo.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
```
and the result should look like [this](http://www.cis.jhu.edu/~parky/dhatKhat/demo.html).  
(This may take around 10 minutes on a typical laptop.)

This is exactly the same setting as Figures 8 and 10 in the paper except that the graph is smaller, e.g., $n = 100$ as opposed to $n = 500$ in the paper.

### Simulations

To run the synthetic data simulation in the paper, please follow these steps:

1. (optional) Run `main_simulation.R` 600 times with different parameter pair `(i, mc)`, where `i` ranges from 1 to 6 and `mc` ranges from 1 to 100. `batch_simulation.sh` is a script for running these jobs in our cluster. Each run produces a `.RData` file as a partial result (600 files total). The precomputed `.RData` files are loaded inside `loadData.R`.

2. Run `simulation1.R` to get Figures 1, 2, 3, 5, 6, 7, 9 in the paper. The results are shown in [here](http://www.cis.jhu.edu/~parky/dhatKhat/simulation1.html).

3. Run `simulation2.R` to get Figures 4, 8, 10, 11 in the paper. The results are shown in [here](http://www.cis.jhu.edu/~parky/dhatKhat/simulation2.html).

### Real-data example

To run the connectome data experiment in the paper, please follow these steps:

1. (optional) Run `main_DS01216.R` 114 times with different parameter `fileIndex`, which  ranges from 1 to 114. `batch_realdata.sh` is a script for running these jobs in our cluster. Each run will produces a `.RData` file as a partial result (114 files total). The precomputed `.RData` files are also stored in Google Drive and loaded inside `loadData.R`.

2. Run `plot_DS01216.R` to get Figures 12, 13 and the results of Table 1 in the paper.  The results are shown in [here](http://www.cis.jhu.edu/~parky/dhatKhat/plot_DS01216.html).
