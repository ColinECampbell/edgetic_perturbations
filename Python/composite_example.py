"""
For analyzing the example network of composite_example.txt. 

We identify the attractors of the original network, then perturb the expanded 
network and identify the new attractors.

Author: Colin Campbell
Date: November 2018
Language: Python 3.6
"""

import networkx as nx
import custom_functions as cf

f = open(r'../Network Rules/composite_example.txt','rt')
lines = f.readlines()
f.close() 


# ANALYSIS OF ORIGINAL NETWORK ================================================
# STEP 1 - FORM NETWORKS -------------------------------------------------
G_expanded = cf.form_expanded_network(lines)
G_original = cf.form_network(lines)

# STEP 2 - REDUCE THE NETWORK GIVEN OUR INITIAL CONDITIONS ---------------
status, G_original_reduced, frozen, removed  = cf.reduce_network(G_original,SCC=False)
G_expanded_reduced = cf.form_expanded_network_reduced(G_expanded, fdict=frozen, rdict=removed)

# STEP 3 - FIND THE FIXED-POINT SCCS -------------------------------------
# Require SCCs where (1) a node and a complementary node are -not- both present
# and (2) the inputs of any composite nodes are also included

SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_expanded_reduced)

# Print the SCCs?
for SCC in SCC_list:
    print('Found a SCC in the original graph:')
    print(SCC)

# Here we find no SCCs so we don't need to reduce the network and repeat on
# non-frozen components.

# STEP 4 - FIND THE OSCILLATING COMPONENTS -------------------------------  
osc_components = cf.find_oscillations(G_expanded,SCC_list)

# Print edges in SCCs?
for i,osc in enumerate(osc_components):
    print('Oscillating components {0} of {1} in the original graph:'.format(i+1,len(osc_components)))
    for j in osc:
        if ' ' not in j and '~' not in j: print(j)    



# MODIFICATION OF NETWORK =====================================================
# STEP 1 - SET UP PERTURBATION ------------------------------------------------
perturbation = [{'A':{'node':'B','value':0}}]
G_perturbed = cf.perturb_network(G_original,perturbation)
G_expanded_perturbed = cf.perturb_expanded_network(G_expanded,perturbation)

# Print edges?
print('Edges in perturbed network:')
for i,j in G_expanded_perturbed.edges():
    print('{0}\t{1}'.format(i,j))
    
# STEP 2 - REDUCE THE NETWORK -------------------------------------------------

# Now repeat analysis from above to ID oscillations and fixed points
SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_expanded_perturbed)

for SCC in SCC_list:
    print('Found a SCC in the perturbed graph:')
    print(SCC)

# Here we find a SCC, so reduce the network, then begin identifying frozen 
# nodes (function calls itself recursively if necessary. Output printed to 
# console.)
fixed_attractors, frozen_states = cf.iterate_SCC(G_perturbed,G_expanded_perturbed,SCC_list)    

# We could proceed to look for oscillating components, but we stop here since
# all nodes are frozen in the one SCC.

#
