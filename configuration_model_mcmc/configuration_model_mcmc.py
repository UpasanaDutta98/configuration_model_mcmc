"""
Created on Thu Jan 25 2021

@author: Upasana Dutta
"""

import networkx as nx
from scipy import stats
import copy
from pathlib import Path
import sys
import numpy as np

sys.path.append(str(Path(__file__).parent))
import dbl_edge_mcmc as mcmc

def graphs_after_McmcConvergence(G, spacing, allow_loops, allow_multi, is_vertex_labeled, count, has_converged):
    '''
    Input --------------------
    G (networkx_class) : Graph with node labels from integer 0 to n-1. It can be a simple graph or a multi-graph, but not weighted or directed. The graph returned will have the same degree sequence as this graph.
    spacing : Sampling gap for the network.
    allow_loops (boolean): True only if loops allowed in the sampling procedure
    allow_multi (boolean): True only if multiedges are allowed in the sampling procedure
    is_v_labeled (boolean): True only if choice of graph labelling is vertex-labeled, inconsequential when graph is simple
    count : Number of graphs to sample from the Configuration model
    has_converged: Indicates if this MCMC is already in its converged region.

    Returns -------------------
    A graph or a list of graphs with the given degree sequence, chosen uniformly at random from the Configuration model.

    '''
    graph_of_graphs = mcmc.MCMC_class(G, allow_loops, allow_multi, is_vertex_labeled) # initialise the MCMC

    if has_converged == False: # We test for convergence only if the MCMC is not in its converged region already.
        G_degree = list(nx.degree(G))
        m = G.number_of_edges()
        n = G.number_of_nodes()      

        S1 = 2*m
        S2 = 0
        S3 = 0
        for i in range(n):
            S2 += (G_degree[i][1])**2
            S3 += (G_degree[i][1])**3

        denominator = S1*S3 - (S2**2) # We calculate the denominator of r only once.

        SL = 0
        for e in G.edges():
            SL += 2*G_degree[e[0]][1]*G_degree[e[1]][1]
        numerator = S1*SL - (S2**2)
        r = float(numerator)/denominator # The starting r value is also calculated only once.
        initial_deg_assort = r 

        datapoints_1 = []
        datapoints_2 = []
        graphs_from_configModel = []

        datapoints_1.append(initial_deg_assort)

        window_size = 29*spacing # Needed to populate each list with 30 points.

        converged_KSTest = 0
        swaps_done = 0
        for i in range(window_size): 
            G2, swap_list = graph_of_graphs.step_and_get_graph()
            delta_r = 0
            if len(swap_list) != 0:
                numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]            
                delta_r = float(2*numerator*2*m)/denominator

            calculated_new_degAssort = initial_deg_assort + delta_r

            if (i+1)%spacing == 0: # To add 29 more points to datapoints_1
                datapoints_1.append(calculated_new_degAssort)


            initial_deg_assort = calculated_new_degAssort
            swaps_done+=1
        while converged_KSTest == 0:
            datapoints_2 = []
            for i in range(window_size+1):        
                G2, swap_list = graph_of_graphs.step_and_get_graph()

                delta_r = 0
                if len(swap_list) != 0:
                    numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]            
                    delta_r = float(2*numerator*2*m)/denominator

                calculated_new_degAssort = initial_deg_assort + delta_r

                if i%spacing == 0: # To add 30 points to datapoints_2
                    datapoints_2.append(calculated_new_degAssort)


                initial_deg_assort = calculated_new_degAssort
                swaps_done+=1

            # Now doing the KS Test here
            KSstatistic, p_value = stats.ks_2samp(datapoints_1,datapoints_2)
            if converged_KSTest == 0 and p_value >= 0.05 :
                converged_KSTest = 1

            datapoints_1 = copy.deepcopy(datapoints_2)

    # At this point, the convergence has been detected.
    if has_converged == False:
        newgraph = copy.deepcopy(G2)
        graphs_from_configModel.append(newgraph)
        count = count - 1
        has_converged = True
    else:
        G2 = G
        graphs_from_configModel = []

    for i in range(count):
        for swaps in range(spacing):
            G2, swap_list = graph_of_graphs.step_and_get_graph()
        newgraph = copy.deepcopy(G2)
        graphs_from_configModel.append(newgraph)

    return graphs_from_configModel, has_converged   

    
def convert_graph(G, allow_loops, allow_multi, is_vertex_labeled):
    '''
    Input --------------------------
    G (networkx_class) : A graph, it can be simple graph or a multi-graph, but not weighted or directed.
    allow_loops (boolean): True only if loops allowed in the sampling procedure
    allow_multi (boolean): True only if multiedges are allowed in the sampling procedure
    is_v_labeled (boolean): True only if choice of graph labelling is vertex-labeled, inconsequential when graph is simple
    
    Returns ------------------------
    converted_G (networkx_class) : A graph having node labels from 0 to n-1, where n is number of nodes in G. It has the same number of nodes and the same edges as present in graph G.
    hash_map (dictionary) : key = node label in G and value = the corresponding node number between 0 and n-1.
    reverse_hash_map (dictionary) : key = node number between 0 and n-1 and value = the corresponding node label in G.
    
    Notes -------------------
    Prints error and returns false if G is a weighted or a directed network.
    
    '''
    
    if nx.is_weighted(G):
        raise ValueError("Cannot apply double edge swap to weighted networks. Exiting.")
    if G.is_directed():
        raise ValueError("Cannot apply double edge swap to directed networks. Exiting.")

    if allow_loops and allow_multi:
        converted_G = nx.MultiGraph()
    elif allow_loops and not allow_multi:
        converted_G = nx.Graph()
    elif not allow_multi and not allow_loops:
        converted_G = nx.Graph()
    else:
        converted_G = nx.MultiGraph()
        
        
    hash_map = {} # A dictionary where key = node label in the original network G and value = the corresponding node number between 0 and n-1
    reverse_hash_map = {} # A dictionary where key = node number between 0 and n-1 and value = the corresponding node label in the original network G    
    i = 0
    for node in G.nodes():
        hash_map[node] = i
        reverse_hash_map[i] = node
        converted_G.add_node(i)
        i+=1

    # Add edges to the new graph
    for eachedge in list(G.edges()):
        converted_G.add_edge(hash_map[eachedge[0]], hash_map[eachedge[1]])

    return converted_G, hash_map, reverse_hash_map

def revert_back(G, reverse_hash_map, allow_loops, allow_multi, is_vertex_labeled):
    '''
    Input -----------------------
    G (networkx_class) : A graph, it can be simple graph or a multi-graph, but not weighted or directed.
    Node labels of G must be from integer 0 to n-1.
    reverse_hash_map (dictionary) : key = node number between 0 and n-1 and value = the corresponding node label in the graph to be returned.
    allow_loops (boolean): True only if loops allowed in the sampling procedure
    allow_multi (boolean): True only if multiedges are allowed in the sampling procedure
    is_v_labeled (boolean): True only if choice of graph labelling is vertex-labeled, inconsequential when graph is simple
    
    Returns ------------------------
    reverted_G (networkx_class) : A graph with node labels as in the reverse_hash_map.
    
    Notes -------------------
    Prints error and returns false if G is a weighted or a directed network.
    
    '''
    if nx.is_weighted(G):
        raise ValueError("Cannot apply double edge swap to weighted networks. Exiting.")
    if G.is_directed():
        raise ValueError("Cannot apply double edge swap to directed networks. Exiting.")

    if allow_loops and allow_multi:
        reverted_G = nx.MultiGraph()
    elif allow_loops and not allow_multi:
        reverted_G = nx.Graph()
    elif not allow_multi and not allow_loops:
        reverted_G = nx.Graph()
    else:
        reverted_G = nx.MultiGraph()
        
        
    n = G.number_of_nodes()
    for i in range(n):
        reverted_G.add_node(reverse_hash_map[i])
        
    for eachedge in list(G.edges()):
        reverted_G.add_edge(reverse_hash_map[eachedge[0]], reverse_hash_map[eachedge[1]])
    
    return reverted_G

def check_density_criterion(G, allow_loops):
    m = G.number_of_edges()
    n = G.number_of_nodes()
    if allow_loops == True:
        rho = m/(n*n/2)
    else:
        rho = m/((n*(n-1))/2)
    density_factor = 2*rho - rho*rho
    if density_factor <= 0.25:
        return 1
    else:
        return 0
    
def check_maxDegree_criterion(G):
    m = G.number_of_edges()
    degrees = list(nx.degree(G))
    degreeList = []
    for eachtuple in degrees:
        degreeList.append(eachtuple[1])
    sorted_degrees = sorted(degreeList, reverse = True)
    if sorted_degrees[0]*sorted_degrees[0] < (2*m/3): 
        return 1
    else:
        return 0
    
def is_disconnected(deg_seq):
    '''
    This function will return True if the degree_seq correponds to a disconnected graph space of 
    loopy-graphs.
    '''
    n = 0
    deg_list = [0]*(max(deg_seq)+1)
    n_unique = 0
    for d in deg_seq:
        if d>0:
            if deg_list[d]==0:
                n_unique+=1 
            n+= 1
            deg_list[d]+=1

    return is_dis_main(deg_list,n_unique,n,deg_seq)
    
def is_dis_main(deg_list,n_unique,n,deg_seq):
    
    deg_list_ori = list(deg_list)
    
    if n==0:
        return False
    
    while deg_list[-1]==0:
        deg_list.pop()


    if n<=2 or len(deg_list)<3:
        return False

    if deg_list[2] == n:
        return True
    
    if len(deg_list)>= n and deg_list[n-1] == n:
        return True

        
    for i in range(0,n+2):
        if deg_list[i]>0:
            min_deg=i
            break

    while n>0:
        
        if n_unique<=2:
            a = min_deg
            b = len(deg_list)-1
            n_a = deg_list[a]
            n_b = deg_list[b]


            if a>=3 and n_a>=3 and a==b-2 and a==n_a+n_b-1 :
                return True

            if a>=3 and n_a>=3 and b-2 == n_a+n_b-1 and a-2==n_b :
                return True
                         
    
        deg_list[min_deg] += -1
        deg_list_ori[min_deg] = 0
        if deg_list[min_deg]==0:
            n_unique += -1
        n += -1
        if min_deg> n:
            return False
    
        s_deg = min_deg 
        prev_s = 0
        i = -1
        while deg_list[i]<= s_deg and s_deg>0:
    
            cur_s = deg_list[i]
            deg_list[i]= prev_s
            if cur_s>0 and prev_s==0:
                n_unique += -1
            if cur_s==0 and prev_s>0:
                n_unique += 1
            s_deg += -cur_s
            prev_s = cur_s
            i += -1
            
         
        if deg_list[i] == 0 and prev_s>0:
            n_unique += 1
            
        deg_list[i] += -s_deg
        deg_list[i] += prev_s
        
        if s_deg>0 and deg_list[i-1]==0:
            n_unique += 1
        deg_list[i-1] += s_deg
     
     
        n += -deg_list[0]
        if deg_list[0]>0:
            n_unique += -1
        deg_list[0] = 0
     
        
        # find new min_deg
        for i in range(0,n+2):              
            if deg_list[i]>0:
                min_deg = i
                break
            
            
        while deg_list[-1]==0:
            deg_list.pop()
            
        if n==0 or len(deg_list)<=3:
            return False

    return False
    
def autocorrelation_function(series, ax=None, **kwds):
    n = len(series)
    data = np.asarray(series)
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        corr = ((data[: n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
        return corr

    x = np.arange(n) + 1
    y = [r(loc) for loc in x]

    z95 = 1.959963984540054
    z99 = 2.5758293035489004
    confidence_limit = z95 / np.sqrt(n)
    return y, confidence_limit

            
class MCMC:
    '''
    Input --------------------
    G (networkx_class) : Graph
    allow_loops (boolean): True only if loops allowed in the sampling procedure
    allow_multi (boolean): True only if multiedges are allowed in the sampling procedure
    is_v_labeled (boolean): True only if choice of graph labelling is vertex-labeled, inconsequential when graph is simple
    count : Number of graphs to sample from the Configuration model
    
    Returns -------------------
    None
    
    Notes -------------------
    Prints error and returns false if G is a weighted or a directed network, or if the graph of graphs in the loopy space is disconnected for the degree sequence of G.
    
    '''
    def __init__(self,G,allow_loops,allow_multi,is_vertex_labeled = True):
        self.allow_loops = allow_loops
        self.allow_multi = allow_multi
        self.is_vertex_labeled = is_vertex_labeled
        self.has_converged = False
        if self.allow_loops == True and self.allow_multi == False: # Loopy graph space
            G_degree = list(nx.degree(G))
            degreeList = []
            for eachpair in G_degree:
                degreeList.append(eachpair[1])
            if is_disconnected(degreeList) == True:
                raise ValueError("The loopy graph space for this degree sequence is disconnected. Exiting.")
                return

        returned = convert_graph(G, allow_loops, allow_multi, is_vertex_labeled)
        self.G = returned[0]
        self.hash_map = returned[1]
        self.reverse_hash_map  = returned[2]
            
    def run_sampling_gap_algorithm(self, warning_outputs=True):
        if warning_outputs == True:
            #Add time display
            pass
            
        graph_of_graphs_list = []
        current_graphs = []
        G_degree = list(nx.degree(self.G))
        m = self.G.number_of_edges()
        n = self.G.number_of_nodes()      
            
        S1 = 2*m
        S2 = 0
        S3 = 0
        for i in range(n):
            S2 += (G_degree[i][1])**2
            S3 += (G_degree[i][1])**3

        denominator = S1*S3 - (S2**2) # We calculate the denominator of r only once.

        if self.allow_multi == False:
            base_sampling_gap = int(m/3 + 300)
        else:
            if self.is_vertex_labeled == True: 
                base_sampling_gap = int(m/2.7 + 100)
            else:
                base_sampling_gap = int(m/3 + 300)

        graph_of_graphs = mcmc.MCMC_class(self.G, self.allow_loops, self.allow_multi, self.is_vertex_labeled)
        for j in range(1000*m): # Burn-in period
            G2, swap_list = graph_of_graphs.step_and_get_graph()

        for i in range(100): # Initialising 100 MCMC walks all starting from graph G2 obatined after burn-in
            graph_of_graphs = mcmc.MCMC_class(G2, self.allow_loops, self.allow_multi, self.is_vertex_labeled)
            graph_of_graphs_list.append(graph_of_graphs)

        inspection_datapoints = [[] for i in range(100)] # Stores r values from 100 MCMC walks.
        SL = 0
        for e in G2.edges():
            SL += 2*G_degree[e[0]][1]*G_degree[e[1]][1]
        numerator = S1*SL - (S2**2)
        r = float(numerator)/denominator

        for i in range(100):    
            inspection_datapoints[i].append(r)

        max_swaps = 250*base_sampling_gap
        for i in range(100): # We run each of the 100 MCMC walks for 250*base_sampling_gap swaps.
            swaps_done = 0
            while swaps_done < max_swaps:
                G2, swap_list = graph_of_graphs_list[i].step_and_get_graph()
                delta_r = 0
                initial_deg_assort = inspection_datapoints[i][-1]
                if len(swap_list) != 0:
                    numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]
                    delta_r = float(2*numerator*2*m)/denominator

                calculated_new_degAssort = initial_deg_assort + delta_r
                inspection_datapoints[i].append(calculated_new_degAssort)
                swaps_done += 1
                initial_deg_assort = calculated_new_degAssort

        found_eta0_200 = 0
        found_eta0_250 = 0
        eta0_200 = 0
        eta0_250 = 0

        gap = base_sampling_gap
        while found_eta0_200 == 0 and found_eta0_250 == 0:
            f_s_200 = 0
            f_s_250 = 0
            for i in range(100): #for every new gap, 100 new fracs need to be calculated and then averaged
                List_of_r_200 = []
                j = 0
                for k in range(200):
                    List_of_r_200.append(inspection_datapoints[i][j])
                    j += gap
                autocorrelation_returned_200 = autocorrelation_function(List_of_r_200)
                Rh_values = autocorrelation_returned_200[0]
                intervalLimit = autocorrelation_returned_200[1]
                outsideInterval = 0
                for eachval in Rh_values:
                    if abs(eachval) > intervalLimit:
                        outsideInterval+=1
                f_s_200 += (outsideInterval/len(Rh_values))
                List_of_r_250 = []
                j = 0
                for k in range(250):
                    List_of_r_250.append(inspection_datapoints[i][j])
                    j += gap

                autocorrelation_returned_250 = autocorrelation_function(List_of_r_250)
                Rh_values = autocorrelation_returned_250[0]
                intervalLimit = autocorrelation_returned_250[1]
                outsideInterval = 0
                for eachval in Rh_values:
                    if abs(eachval) > intervalLimit:
                        outsideInterval+=1   
                f_s_250 += (outsideInterval/len(Rh_values))

            f_s_200 = f_s_200/100
            f_s_250 = f_s_250/100
            if f_s_200 <= 0.05:
                found_eta0_200 = 1
                eta0_200 = gap
            if f_s_250 <= 0.05:
                found_eta0_250 = 1
                eta0_250 = gap
            gap += 1

            if found_eta0_200 == 1 or found_eta0_250 == 1:
                break

            max_swaps = 250
            for i in range(100): # We run each of the 100 MCMC walks for additional 250 swaps
                swaps_done = 0
                while swaps_done < max_swaps:
                    G2, swap_list = graph_of_graphs_list[i].step_and_get_graph()
                    delta_r = 0
                    initial_deg_assort = inspection_datapoints[i][-1]
                    if len(swap_list) != 0:
                        numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]
                        delta_r = float(2*numerator*2*m)/denominator

                    calculated_new_degAssort = initial_deg_assort + delta_r
                    inspection_datapoints[i].append(calculated_new_degAssort)
                    swaps_done += 1
                    initial_deg_assort = calculated_new_degAssort


        while found_eta0_250 == 0:
            max_swaps = 250
            for i in range(100): # We run each of the 100 MCMC walks for additional 250 swaps
                swaps_done = 0
                while swaps_done < max_swaps:
                    G2, swap_list = graph_of_graphs_list[i].step_and_get_graph()
                    delta_r = 0
                    initial_deg_assort = inspection_datapoints[i][-1]
                    if len(swap_list) != 0:
                        numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]
                        delta_r = float(2*numerator*2*m)/denominator

                    calculated_new_degAssort = initial_deg_assort + delta_r
                    inspection_datapoints[i].append(calculated_new_degAssort)
                    swaps_done += 1
                    initial_deg_assort = calculated_new_degAssort

            f_s_250 = 0
            for i in range(100): #for every new gap, 100 new fracs need to be calculated and then averaged
                List_of_r_250 = []
                j = 0
                for k in range(250):
                    List_of_r_250.append(inspection_datapoints[i][j])
                    j += gap

                autocorrelation_returned_250 = autocorrelation_function(List_of_r_250)
                Rh_values = autocorrelation_returned_250[0]
                intervalLimit = autocorrelation_returned_250[1]
                outsideInterval = 0
                for eachval in Rh_values:
                    if abs(eachval) > intervalLimit:
                        outsideInterval+=1   
                f_s_250 += (outsideInterval/len(Rh_values))

            f_s_250 = f_s_250/100

            if f_s_250 <= 0.05:
                found_eta0_250 = 1
                eta0_250 = gap

            gap += 1

        while found_eta0_200 == 0:
            max_swaps = 200
            for i in range(100): # We run each of the 100 MCMC walks for additional 200 swaps, not 250.
                swaps_done = 0
                while swaps_done < max_swaps:
                    G2, swap_list = graph_of_graphs_list[i].step_and_get_graph()
                    delta_r = 0
                    initial_deg_assort = inspection_datapoints[i][-1]
                    if len(swap_list) != 0:
                        numerator = G_degree[swap_list[0][0]][1]*G_degree[swap_list[0][1]][1] + G_degree[swap_list[1][0]][1]*G_degree[swap_list[1][1]][1] - G_degree[swap_list[2][0]][1]*G_degree[swap_list[2][1]][1] - G_degree[swap_list[3][0]][1]*G_degree[swap_list[3][1]][1]
                        delta_r = float(2*numerator*2*m)/denominator

                    calculated_new_degAssort = initial_deg_assort + delta_r
                    inspection_datapoints[i].append(calculated_new_degAssort)
                    swaps_done += 1
                    initial_deg_assort = calculated_new_degAssort

            f_s_200 = 0
            for i in range(100): #for every new gap, 100 new fracs need to be calculated and then averaged
                List_of_r_200 = []
                j = 0
                for k in range(200):
                    List_of_r_200.append(inspection_datapoints[i][j])
                    j += gap

                autocorrelation_returned_200 = autocorrelation_function(List_of_r_200)
                Rh_values = autocorrelation_returned_200[0]
                intervalLimit = autocorrelation_returned_200[1]
                outsideInterval = 0
                for eachval in Rh_values:
                    if abs(eachval) > intervalLimit:
                        outsideInterval+=1
                f_s_200 += (outsideInterval/len(Rh_values))

            f_s_200 = f_s_200/100

            if f_s_200 <= 0.05:
                found_eta0_200 = 1
                eta0_200 = gap

            gap += 1

        return int((eta0_200+eta0_250)/2)

    def get_sampling_gap(self, warning_outputs=True):
        m = self.G.number_of_edges()
        if self.allow_multi == False:
            density_criterion_satisfied = check_density_criterion(self.G, self.allow_loops)
            if density_criterion_satisfied == 1 or m <= 1000:
                #print("Sampling gap of the network = ", int(m/3 + 300), "(m/3 + 300).")
                return int(m/3 + 300)
            else:
                if warning_outputs==True:
                    print("The network does not satisfy the density criterion for automatic selection of sampling gap.")
                    print("Running the Sampling Gap Algorithm. This will take a while.....")
                sampling_gap = self.run_sampling_gap_algorithm(warning_outputs)
                #print("Sampling gap of the network = ", sampling_gap)
                return sampling_gap
        else:
            if self.is_vertex_labeled == False: 
                #print("Sampling gap of the network = ", int(m/3 + 300), "(m/3 + 300).")
                return int(m/3 + 300)
            else:
                maxDegree_criterion_satisfied = check_maxDegree_criterion(self.G)
                if maxDegree_criterion_satisfied == 1 or m <= 1000:
                    #print("Sampling gap of the network = ", int(m/2.7 + 100), "(m/2.7 + 100).")
                    return int(m/2.7 + 100)
                else:
                    if warning_outputs==True:
                        print("The network does not satisfy the maximum degree criterion for automatic selection of sampling gap.")
                        print("Running the Sampling Gap Algorithm. This will take a while.....")
                    sampling_gap = self.run_sampling_gap_algorithm(warning_outputs)
                    #print("Sampling gap of the network = ", sampling_gap)
                    return sampling_gap

    def get_graphs(self, count=1, warning_outputs = True, sampling_gap = -999):
        # Testing the loop memory now - Done
        # Letting user to provide sampling gap - Done
        # Letting user turn of warning statements - Done
        # Letting user to run the main sampling gap algorithm - Done
        # Giving the user an estimate of time that'll be needed - Not yet done
        '''
        Input --------------------
        count (int) : Number of graphs to sample from the Configuration model. Default value is 1.
        warning_outputs (boolean) = True if user wants warning messages printed when get_sampling_gap() function runs. Set to 'False' is warnings are not desired.
        sampling_gap (int) : Number of double-edge swaps applied (accepted or rejected) between two MCMC states sampled for the convergence testing.
        If user provides a sampling_gap, it is used, otherwise the sampling gap is calculated by calling the get_sampling_gap() function
        
        Returns -------------------
        A list of graph(s) with the given degree sequence, chosen uniformly at random from the desired graph space.
        '''
        if sampling_gap == -999: # if user-defined sampling gap has not been provided
            self.spacing = self.get_sampling_gap(warning_outputs)
        else:
            self.spacing = sampling_gap
            
        list_of_networks, has_converged = graphs_after_McmcConvergence(self.G, self.spacing, self.allow_loops, self.allow_multi, self.is_vertex_labeled, count, self.has_converged)
        self.has_converged = has_converged 
        self.G = list_of_networks[-1]
        networks_from_configModel = []
        for each_net in list_of_networks:
            reverseHashed_net = revert_back(each_net, self.reverse_hash_map, self.allow_loops, self.allow_multi, self.is_vertex_labeled)
            networks_from_configModel.append(reverseHashed_net)
        if len(networks_from_configModel) == 1:
            return networks_from_configModel[0]
        return networks_from_configModel