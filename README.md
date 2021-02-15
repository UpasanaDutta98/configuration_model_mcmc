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

### A simple example.

Here is the most basic example.

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
G2 = mcmc_object.get_graph()
```

In the above example, `G_2` is sampled from the stub-labeled simple graph space. `G_2` has the same degree sequence as the Karate Club network. Details of how to choose the graph space can be found in [[1]](https://epubs.siam.org/doi/pdf/10.1137/16M1087175).

By default, the chosen graph space is the simple vertex-labeled graph space. In the example below, the graph `G_2` is obtained from the simple vertex-labeled graph space.

```python
# read the Karate Club network
G = nx.karate_club_graph()

# initialise the MCMC
mcmc_object = CM.MCMC(G)

# obtain a new graph (G_2) from the Configuration model
G2 = mcmc_object.get_graph()
```


### Sample multiple graphs

Multiple graphs can also be sampled from the Configuration model with/without a loop

```python
# read the Karate Club network
G = nx.karate_club_graph()

# initialise the MCMC in simple stub-labeled space
allow_loops = False
allow_multi = False
is_vertex_labeled = False
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# Obtain 10 graphs from the Configuration model over 2 iterations.
for i in range(2):
    list_of_graphs = mcmc_object.get_graph(count=5)
    for each_graph in list_of_graphs:
        print(round(nx.degree_pearson_correlation_coefficient(each_graph),4), end = " ")
    print()
```
Output:
```
-0.2735 -0.2984 -0.2735 -0.3496 -0.3126
-0.2646 -0.3723 -0.2864 -0.3131 -0.2646
```

In the above code, by specifying ```count=5``` in the function parameter, we sample 10 networks by running the loop twice. The default value of ```count``` is 1. The degree assortativity values of the samples networks are printed as output.

### Sampling Gap heuristics

If the network does not satisfy the conditions under which our Sampling Gap heuristic can automatically choose a sampling gap, the Sampling Gap algorithm is run, which might take a while.

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
```

Output:
```
The network does not satisfy the density criterion for automatic selection of sampling gap.
Running the Sampling Gap Algorithm. This will take a while.....
----- Running initial burn-in -----
76%|████████████████████████████         | 7568/10000 [00:33<00:10, 229.00it/s]
----- Initial burn-in complete -----
fs = 0.601, α = 0.05, time = 0.01 mins.
fs = 0.44, α = 0.05, time = 0.51 mins.
....
....
....
fs = 0.052, α = 0.05, time = 11.23 mins.
```

The above code reads the edgelist of a network of Jazz musicians and samples 5 graphs from the vertex-labeled simple graph space. The network does not satisfy the contraints necessary for the automatic selection of the sampling gap, so the Sampling Gap Algorithm is run. A progress bar is displayed during the burn-in period of the MCMC walk. The printed lines thereafter show the values of the variables fs and α as the run proceeds. When fs becomes less than α, the run is completed. The variable `list_of_graphs` contains the 5 graphs sampled from the Configuration model.

The warning messages and the print statements in the output can be muted by specifying ```warnings = False``` while creating the MCMC object. The deafult value is ```warnings = True```.

```python
# read the network
G = nx.read_edgelist("JazzMusicians.txt")

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
warnings = False
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled, warnings)

# Obtain 5 graphs from the Configuration model
list_of_graphs = mcmc_object.get_graphs(count=5)
```

### Running the Sampling Gap Algorithm

If you want to call the Sampling Gap Algorithm to obatin a bespoke sampling gap for your graph, you may do so as follows:

```python
# read the network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# run the Sampling Gap Algorithm
sampling_gap = mcmc_object.run_sampling_gap_algorithm()
print("Sampling gap obtained = ", sampling_gap)
```

Output:
~~~
Sampling gap obtained = 31
~~~

Again, warning messages will be muted if specified so while creating the MCMC object. The warning messages of the sampling gap algorithm in particular can be muted using the following code, even when it was not muted while creating the MCMC object.
```python
sampling_gap = mcmc_object.run_sampling_gap_algorithm(warnings = False)
```

The default significance level of the the autocorrelation hypothesis tests = 5% and the default number of parallel MCMC chains run for the Sampling Gap Algorithm = 100. However, you can change them by specifying in the function parameter. For example, we can set the significance level as 10% and run 50 parallel MCMC chains as follows:

```python
sampling_gap = mcmc_object.run_sampling_gap_algorithm(alpha_c = 0.1, numchains = 50)
```


### Customise sampling gap

You can also specify a custom sampling gap that you want to run the convergence detection test with, using the ```sampling_gap``` function parameter.

```python
# read the network
G = nx.karate_club_graph()

# initialise the MCMC
allow_loops = False
allow_multi = False
is_vertex_labeled = True
mcmc_object = CM.MCMC(G, allow_loops, allow_multi, is_vertex_labeled)

# Specify the sampling gap 
gap = 100 # any user-defined value
list_of_graphs = mcmc_object.get_graph(count=5, sampling_gap = gap)
for each_graph in list_of_graphs:
    print(round(nx.degree_pearson_correlation_coefficient(each_graph),4), end = " ")
```
Output:
```
-0.2614 -0.4342 -0.3118 -0.2864 -0.2939
```


## Notes

The package will not work for weighted networks, directed networks, hypergraphs, or simplical complexes.

## Feedback and bugs

If you find a bug or you have any feedback, please contact upasana.dutta@colorado.edu.

## License

GNU General Public License v3+