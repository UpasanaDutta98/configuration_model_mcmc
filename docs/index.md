# User Guide

## Introduction

**configuration_model_mcmc** is a tool for sampling networks from the Configuration model, given a network and a graph space. This code package builds upon the Double Edge Swap MCMC Graph Sampler by Fosdick et al. [[1]](https://epubs.siam.org/doi/pdf/10.1137/16M1087175). It detects convergence in the Double Edge Swap MCMC and samples networks from the MCMC's stationary distribution, so that the samples are uniform random draws from the Configuration model.

[1] Bailey K. Fosdick, Daniel B. Larremore, Joel Nishimura, Johan Ugander (2018) Configuring Random Graph Models with Fixed Degree Sequences. SIAM Review 60(2):315–355.

## Why use configuration_model_mcmc?

The [configuration_model](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.generators.degree_seq.configuration_model.html) function of the [networkx](https://pypi.org/project/networkx/) package in Python works best only for the loopy multigraph space because the graph returned is a pseudograph, i.e., the graph is allowed to have self-loops and parallel edges(multi-edges). Often times practitioners remove the self-loops and collapse the multi-edges in the network returned by the function to get a *simple network*, but this modification changes the degree sequence of the network. It also introduces a bias in the network generating process because the high-degree nodes are more likely to have self-loops and multi-edges attached to them than are the low-degree nodes. Therefore, the network generated is a biased sample. The **configuration_model_mcmc** package lets you sample an unbiased sample from the Configuration model on eight different graph spaces parameterized by self-loops/no self-loops, multi-edges/no multi-edges and stub-labeled/vertex-labeled.

## Installing

`pip install configuration_model_mcmc`

This package supports Python 3.7 and 3.8, and requires the packages Numpy>=1.17.1, Networkx>=2.4, Scipy>=1.4.1, Numba==0.49.1. These dependencies are automatically installed while installing the package.
s


## Notes

The package will not work for weighted networks, directed networks, hypergraphs, or simplical complexes. Please refer to the [website](https://UpasanaDutta98.github.io/configuration_model_mcmc/) for complete documentation.

## Feedback and bugs

If you find a bug or you have any feedback, please contact upasana.dutta@colorado.edu.

## License

GNU General Public License v3+