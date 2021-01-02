# configuration_model_mcmc

## What is it?

**configuration_model_mcmc** is a tool for sampling networks from the Configuration model, given a network and a graph space. This code package builds upon the Double Edge Swap MCMC Graph Sampler by Fosdick et al. [1]. It detects convergence in the Double Edge Swap MCMC and samples networks from the MCMC's stationary distribution, such that the samples are uniform random draws from the Configuration model.

[1]Bailey K. Fosdick, Daniel B. Larremore, Joel Nishimura, Johan Ugander (2018) Configuring Random Graph Models with Fixed Degree Sequences. SIAM Review 60(2):315–355.(https://epubs.siam.org/doi/pdf/10.1137/16M1087175)

### Why use configuration_model_mcmc?

The [configuration_model](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.generators.degree_seq.configuration_model.html) function of the [networkx](https://pypi.org/project/networkx/) package in Python works best only for the loopy multigraph space, where a graph is allowed to have self-loops and parallel edges (multi-edges), since the graph returned by this function is a pseudograph(can have multi-edges and self loops). Often times practitioners remove the self-loops and collapse the multi-edges in the network returned by the function to get a *simple network*, but this modification changes the degree sequence of the network. It also introduces a bias in the network generating process because the high-degree nodes are more likely to have self-loops and multi-edges attached to them than are the low-degree nodes. Therefore, the network generated is a biased sample. The **configuration_model_mcmc** package lets you sample an unbiased sample from the Configuration model on eight different graph spaces parameterized by self-loops/no self-loops, multi-edges/no multi-edges and stub-labeled/vertex-labeled.

## Installing

`pip install configuration_model_mcmc`

This package supports Python>=3.7.x, and requires the packages Numpy>=1.17.1, Networkx>=2.4, Scipy>=1.4.1, Numba>=0.49.1. These dependencies are automatically installed while installing the package.

## An example

1.
```python
import configuration_model_mcmc as CM
import networkx as nx

# read the Karate Club network
G = nx.karate_club_graph()

# obtain a new graph (G_2) from the Configuration model
allow_loops = False
allow_multi = False
is_vertex_labeled = False
G_2 = CM.get_graphs_from_ConfigModel(G, allow_loops, allow_multi, is_vertex_labeled)
```

In the above example, `G_2` is a network with the sample degree sequence as the Karate Club network. Details on how to choose the graph space can be found in [1]. `G_2` is sampled from the stub-labeled simple graph space.

2.
```python
import configuration_model_mcmc as CM
import networkx as nx

# read the network
G = nx.read_edgelist("JazzMucisians.txt")
allow_loops = False
allow_multi = False
is_vertex_labeled = True
list_of_graphs = CM.get_graphs_from_ConfigModel(G, allow_loops, allow_multi, 
is_vertex_labeled, count=5)

# Output:
# The network does not satisfy the density criterion for automatic selection of sampling gap.
# Running the Sampling Gap Algorithm. This will take a while.....
```

The above code reads the edgelist of a network of Jazz musicians and then samples 5 graphs from the vertex-labeled simple graph space. The network does not satisfy the contraints for the use of an automated sampling gap, so you have to wait till the Sampling Gap Algorithm is run. Once the run is over, the variable `list_of_graphs` contains all the 5 graphs from the Configuration model.

## Notes

The package will not work for weighted networks, directed networks, hypergraphs, or simplical complexes. To know more on how to use the package, please refer to the documentation.

## Feedback and bugs

If you find a bug or you have any feedback, please contact upasana.dutta@colorado.edu.

## License

GNU General Public License v3+