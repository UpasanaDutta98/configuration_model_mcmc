# configuration_model_mcmc

## What is it?

**configuration_model_mcmc** is a tool for sampling networks from the Configuration model, given a network and a graph space. This code package builds upon the Double Edge Swap MCMC Graph Sampler by Fosdick et al. [1]. It detects convergence in the Double Edge Swap MCMC and samples networks from the MCMC's stationary distribution, so that the samples are uniform random draws from the Configuration model.

[[1]](https://epubs.siam.org/doi/pdf/10.1137/16M1087175) Bailey K. Fosdick, Daniel B. Larremore, Joel Nishimura, Johan Ugander (2018) Configuring Random Graph Models with Fixed Degree Sequences. SIAM Review 60(2):315–355.

### Why use configuration_model_mcmc?

The [configuration_model](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.generators.degree_seq.configuration_model.html) function of the [networkx](https://pypi.org/project/networkx/) package in Python works best only for the loopy multigraph space because the graph returned is a pseudograph, i.e., the graph is allowed to have self-loops and parallel edges(multi-edges). Often times practitioners remove the self-loops and collapse the multi-edges in the network returned by the function to get a *simple network*, but this modification changes the degree sequence of the network. It also introduces a bias in the network generating process because the high-degree nodes are more likely to have self-loops and multi-edges attached to them than are the low-degree nodes. Therefore, the network generated is a biased sample. The **configuration_model_mcmc** package lets you sample an unbiased sample from the Configuration model on eight different graph spaces parameterized by self-loops/no self-loops, multi-edges/no multi-edges and stub-labeled/vertex-labeled.

## Installing

`pip install configuration_model_mcmc`

This package supports Python>=3.7.x, and requires the packages Numpy>=1.17.1, Networkx>=2.4, Scipy>=1.4.1, Numba>=0.49.1. These dependencies are automatically installed while installing the package.

## Examples

1. A simple example.
```python
import configuration_model_mcmc as CM
import networkx as nx

# read the Karate Club network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = False
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# obtain a new graph (G_2) from the Configuration model
G2 = mcmc_object.get_graphs()
```

In the above example, `G_2` is sampled from the stub-labeled simple graph space. `G_2` has the same degree sequence as the Karate Club network. Details on how to choose the graph space can be found in [1].

2. Multiple graphs can also be sampled from the Configuration model with/without a loop
```python
# read the Karate Club network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = False
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# Obtain 10 graphs from the Configuration model over 2 iterations.
for i in range(2):
    list_of_graphs = mcmc_object.get_graphs(count=5)
    for each_graph in list_of_graphs:
        print(round(nx.degree_pearson_correlation_coefficient(each_graph),4), end = " ")
    print()
    
# Output:
# -0.2735 -0.2984 -0.2735 -0.3496 -0.3126
# -0.2646 -0.3723 -0.2864 -0.3131 -0.2646
```

In the above code, by specifying ```count=5``` in the function parameter, we sample 10 networks by running the loop twice. The default value of ```count``` is 1. The degree assortativity values of the samples networks are printed as output.

3. If the network does not satisfy the conditions under which our Sampling Gap heuristic can be applied, the Sampling Gap algorithm is run.
```python
# read the network
G = nx.read_edgelist("JazzMusicians.txt")

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# Obtain 5 graphs from the Configuration model
list_of_graphs = mcmc_object.get_graphs(count=5)

# Output:
# The network does not satisfy the density criterion for automatic selection of sampling gap.
# Running the Sampling Gap Algorithm. This will take a while.....
```

The above code reads the edgelist of a network of Jazz musicians and samples 5 graphs from the vertex-labeled simple graph space. The network does not satisfy the contraints for the use of an automated sampling gap, so the Sampling Gap Algorithm is run. Once the run is over, the variable `list_of_graphs` contains the 5 graphs sampled from the Configuration model.

4. The warning messages in the output can be muted by specifying ```warning_outputs = False``` in the function parameter.
```python
# read the network
G = nx.read_edgelist("JazzMusicians.txt")

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# Obtain 5 graphs from the Configuration model
list_of_graphs = mcmc_object.get_graphs(count=5, warning_outputs = False)
```

5. One can run the Sampling Gap Algorithm as follows:
```python
# read the network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops=False, allow_multi=False, is_vertex_labeled=False)

# run the Sampling Gap Algorithm
sampling_gap = mcmc_object.run_sampling_gap_algorithm()
print("Sampling gap obtained = ", sampling_gap)

# Output:
# Sampling gap obtained = 31
```

Again, the warning messages may be muted using 
```python
sampling_gap = mcmc_object.run_sampling_gap_algorithm(warning_outputs = False)
```

6. One can also specify a bespoke sampling gap one wants to run the convergence detection test with, using the ```sampling_gap``` function parameter.
```python
# read the network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops=False, allow_multi=False, is_vertex_labeled=False)

# Specify the sampling gap 
gap = 100 # any user-defined value
list_of_graphs = mcmc_object.get_graphs(count=5, sampling_gap = gap)
for each_graph in list_of_graphs:
    print(round(nx.degree_pearson_correlation_coefficient(each_graph),4), end = " ")

# Output:
# -0.2614 -0.4342 -0.3118 -0.2864 -0.2939
```


## Notes

The package will not work for weighted networks, directed networks, hypergraphs, or simplical complexes.

## Feedback and bugs

If you find a bug or you have any feedback, please contact upasana.dutta@colorado.edu.

## License

GNU General Public License v3+