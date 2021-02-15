# Examples

## A simple example.

Here's the most basic example.
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

## Sample multiple graphs

Multiple graphs can also be sampled from the Configuration model with/without a loop

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

## Sampling Gap heuristic conditions
If the network does not satisfy the conditions under which our Sampling Gap heuristic can be applied, the Sampling Gap algorithm is run.
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

The warning messages in the output can be muted by specifying ```warning_outputs = False``` in the function parameter.
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

## Calling the Sampling Gap Algorithm
One can run the Sampling Gap Algorithm as follows:
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

## Customise sampling gap
One can also specify a bespoke sampling gap one wants to run the convergence detection test with, using the ```sampling_gap``` function parameter.
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
