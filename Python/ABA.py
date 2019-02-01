"""
For analysis of the ABA network from Albert et al. PLOS Bio 2017. This script 
sequentially walks through the main components of the expanded graph analysis 
described in the report. See ABA_STN.py for additional analysis on the ABA
network.

Author: Colin Campbell
Date: November 2018
Language: Python 3.6
"""

import networkx as nx
import custom_functions as cf
from itertools import product
import time
import pickle


TIME0 = time.time()

f = open(r'../Network Rules/ABA2017.txt','rt')
lines = f.readlines()
f.close() 

# ANALYSIS OF ORIGINAL NETWORK ================================================
# STEP 1 - FORM NETWORKS -------------------------------------------------
G_expanded = cf.form_expanded_network(lines)
G_original = cf.form_network(lines)

# Print edges in expanded network?
#print('Edges in the expanded network:')
#for i,j in G_expanded.edges():
#    print(i,';',j)

TIME1 = time.time()
print('STEP 1 completed in {0:.2f} minutes.'.format((TIME1-TIME0)/60.0))

# STEP 2 - REDUCE THE NETWORK GIVEN OUR INITIAL CONDITIONS ---------------
status, G_original_reduced, frozen, removed  = cf.reduce_network(G_original,SCC=False)
G_expanded_reduced = cf.form_expanded_network_reduced(G_expanded,fdict=frozen,rdict=removed)

# We refer to the expanded reduced network in ABA_STN.py; save a copy
f = open(r'ABA_expanded_reduced_network.data','wb')
pickle.dump(G_expanded_reduced,f)
f.close()

TIME2 = time.time()
print('STEP 2 completed in {0:.2f} minutes ({1:.2f} total).'.format((TIME2-TIME1)/60.0,(TIME2-TIME0)/60.0))

# STEP 3 - FIND THE FIXED-POINT SCCS -------------------------------------
# Require SCCs where (1) a node and a complementary node are -not- both present
# and (2) the inputs of any composite nodes are also included

SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_expanded_reduced)

for SCC in SCC_list:
    print('Found a SCC:')
    print(SCC)
fixed_attractors, frozen_states = cf.iterate_SCC(G_original_reduced,
                                                 G_expanded_reduced,SCC_list,
                                                 frozen=frozen, removed=removed, 
                                                 dump_to_screen = False)

TIME3 = time.time()
print('STEP 3 completed in {0:.2f} minutes ({1:.2f} total).'\
      .format((TIME3-TIME2)/60.0,(TIME3-TIME0)/60.0))

'''
We find 24 frozen_states, which correspond to all possible orderings of the 
4 stable motifs. They all yield the same state, as reported in Albert et al. 
(Table S6).
'''

# STEP 4 - FIND THE OSCILLATING COMPONENTS -------------------------------
    
# Now find valid oscillations
osc_components = cf.find_oscillations(G_expanded,SCC_list)

for i,osc in enumerate(osc_components):
    print('Oscillating components {0} of {1}:'.format(i+1,len(osc_components)))
    for j in osc:
        if ' ' not in j: print(j)   
   
TIME4 = time.time()
print('STEP 4 completed in {0:.2f} minutes ({1:.2f} total).'\
      .format((TIME4-TIME3)/60.0,(TIME4-TIME0)/60.0))
     
'''
Here we find none, as expected from Albert et al.
'''

# BEGIN PERTURBATION ANALYSIS =================================================
# STEP 5 - ANALYZE PARENT SCCS --------------------------------------
'''
We want to ID all minimal cuts that destroy the constituent stable SCCs.
Each of the 4 stable SCCs has its own parent SCC, and we would like to 
investigate destroying each of them independently or in combination (sets of 2,
 3, or all 4)

We implement this by explicitly choosing the combination below; we can make
different choices and re-run the file to investigate all possibilities.
'''
MHS_perturbations = set()

# Which parent_SCCs will we consider?
# We do this explicitly below, but first ensure that the analysis above gives 
# us what we assume
set_check= [set(x) for x in parent_SCCs]
if not all([set(x) in set_check for x in [['PLDd', 'PA', 'RBOH', 'ROS', 'NO'], 
            ['Microtubule'], ['MAPK912'], ['CPK213']]]) or len(set_check)!=4:
    raise RuntimeError("parent_SCCs does not have the expected population")

# Choose -one- of the below pairs of lines to consider in the remaining analysis

# INDIVIDUAL
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO']                          # cycle only
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS']]
#
#parent_SCC_nodes = ['Microtubule']
#SCC_list = [['Microtubule']]
#  
#parent_SCC_nodes = ['MAPK912']
#SCC_list = [['MAPK912']]
#
#parent_SCC_nodes = ['CPK213']
#SCC_list = [['CPK213']]

##PAIRS  
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','Microtubule']            # cycle, Microtubule
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['Microtubule']]
#
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','MAPK912']                # cycle, MAPK912
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['MAPK912']]
#
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','CPK213']                 # cycle, CPK213
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['CPK213']]
#
#parent_SCC_nodes = ['Microtubule','MAPK912']
#SCC_list = [['Microtubule'], ['MAPK912']]
#
#parent_SCC_nodes = ['Microtubule','CPK213']
#SCC_list = [['Microtubule'], ['CPK213']]
#
#parent_SCC_nodes = ['MAPK912','CPK213']
#SCC_list = [['MAPK912'], ['CPK213']]

##TRIPLETS
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','Microtubule','MAPK912']  # cycle, Microtubule, MAPK912  
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['Microtubule'], ['MAPK912']]
#
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','Microtubule','CPK213']   # cycle, Microtubule, CPK213
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['Microtubule'], ['CPK213']]
#
#parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','MAPK912','CPK213']       # cycle, MAPK912, CPK213
#SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['MAPK912'], ['CPK213']]

##ALL FOUR
parent_SCC_nodes = ['PLDd', 'PA', 'RBOH', 'ROS', 'NO','Microtubule','MAPK912','CPK213']
SCC_list = [['PLDd', 'PA', 'RBOH', 'ROS'], ['Microtubule'], ['MAPK912'], ['CPK213']]
# -----
# -----

# Form the parent SCC expanded graph
parent_graph = G_expanded_reduced.copy()
remove_nodes = set(parent_graph.nodes())-set(parent_SCC_nodes)
parent_graph.remove_nodes_from(remove_nodes)

# Determine which partial cycles exist entirely in the parent SCC
ppc = []
for cycle in partial_cycles:
    if set(cycle) <= set(parent_graph.nodes()): ppc.append(cycle)
    
# To ensure complete coverage, we use as source every node that is frozen in 
# the motif
source_nodes = []
for scc in SCC_list:
    source_nodes += [x for x in scc if ' ' not in x]
source_nodes = list(set(source_nodes))

# We are now ready to look for the valid circuits
path_store = []
for source in source_nodes:
    print('Considering source = {0}'.format(source))
    count = 0
    for path in cf.parent_SCC_search(parent_graph, source):
        # We have a new path, but we want to only store it if it is a unique 
        # addition. We also don't want long circuits that contain a complete 
        # valid circuit. (If we destroy the smaller we must destroy the bigger, 
        # but the reverse is not true.)
        
        # We could have [a,b,c] and [b,c,a] - i.e. different orderings, so we 
        # first check the lengths of the lists: if so then we dbl one and see 
        # if the other is a sub-component
        new = True
        
        for p in path_store:
            # Transfer these to strings to make substring searching easier       
            p_store = [str(x) for x in p]
            p_store = ''.join(p_store)
            # double so when we search we avoid 'wraparound' permutations
            p_store *= 2  
            
            p_new = [str(x) for x in path]
            p_new = ''.join(p_new)
            
            if p_new in p_store:
                new = False
                break
            
        if new:
            path_store.append(path[:])
        count += 1

    print('\tcount = {0}'.format(count))
print('All done; unique paths = {0}'.format(len(path_store)))

TIME5 = time.time()
print('STEP 5 completed in {0:.2f} minutes ({1:.2f} total).'\
      .format((TIME5-TIME4)/60.0,(TIME5-TIME0)/60.0))

# STEP 6 - FIND MHS --------------------------------------
'''
At this point we have our unique circuits in this parent SCC. Now we need
to find the edges that can be destroyed to destroy each. 
We have edges of the following types ('node' is node or antinode):
    1. node -> composite
    2. node -> node
    3. composite -> node

We want to find edges of the types (2) and (3). Moreover, we want to ID the 
minimum hitting sets of these edges: the smallest  collections of edges that 
collectively covers each of the edges in the circuits.

This is the minimum hitting set problem. We implement the MMCS algorithm
from Murakami and Uno (Discrete Applied Mathematics, 2014). 

Before applying it here, we transform our circuits from lists of sequential
nodes to sets of relevant edges.
'''

# Transform circuits from list of nodes to set of edges (only types 2 and 3)
edge_sets = []
for path in path_store:
    new_edge = set()
    for i,j in zip(range(-1,len(path)-1),range(0,len(path))):
        # Include the -1 to 0 index pair since we have a cycle but don't
        # include the start/end node in the last position
        if ' ' in path[j]:
            # Skip this because the target node is a composite (type 1)
            continue
        if ' ' in path[i]:
            # composite -> node; add each constituent of the composite
            cnodes = path[i].split()
            for cn in cnodes: new_edge.add((cn,path[j]))
        else:
            # node -> node; store this edge as a tuple of nodes
            new_edge.add((path[i],path[j]))
    edge_sets.append(new_edge)

# Also need a list of all nodes in all sets (i.e. all included edges)
edge_collection = edge_sets[0].copy()
for i in range(1,len(edge_sets)):
    edge_collection = edge_collection | edge_sets[i]

# Now we need to perform the MHS search
MHS = cf.MMCS(edge_sets,edge_collection)

# Some cleaning before printing to screen...
MHS = [list(x) for x in MHS]
print('Minimum Hitting Sets:')
for i in range(len(MHS)):
    print('# {0} of {1}:'.format(i+1,len(MHS)))
    print(MHS[i])
    print('\n')
    
TIME6 = time.time()
print('STEP 6 completed in {0:.2f} minutes ({1:.2f} total).'\
      .format((TIME6-TIME5)/60.0,(TIME6-TIME0)/60.0))


# STEP 7 - ANALYZE MHS --------------------------------------
'''
Now each of the cuts proposed above needs to be analyzed to make sure we do 
not, for instance, introduce a new attractor.

So below we directly check the new version of the expanded network for the 
creation of a new viable loop (for stable SCCs) or a new oscillation.
'''

# First we need to update the expanded graph in response to the proposed
# edge cuts

perturbation_counter = 0 
failed_combinations = 0

for mhs in MHS: 
    """
    ID all source nodes and all target nodes in the mhs
    If any source nodes are composites, we can consider -any- internal node to 
    destroy the composite.
    If we have multiple composites, we have combinations to consider. Example:
    [AB -> C ; BD -> E] as two edges in the mhs means we can consider AB, AD, 
    BB, and BD as OFF for C and E, respectively.
    
    In addition, each of those composites is set to the value that destroys 
    the node (in A -> B, A is OFF; in ~A -> B, A is ON; in AB -> C, A is OFF, 
    in ~AB ->C, A is ON ... so always the opposite state.
    
    To take care of the above we form a list of lists; each internal list 
    corresponds to a target and contains entries of the form 
    ["possible source","this target"]
    """
    master = []
    for entry in mhs:
        all_sources = entry[0].split()
        master_entry = []
        for x in all_sources:
            master_entry.append(tuple([x,entry[1],0]))
        if master_entry != []: master.append(tuple(master_entry))
    
    if master == []: continue

    for product_iter in product(*master):
        # Using the above example, we make mhs_cur [(A,C),(B,E)] on the first 
        # pass, then [(A,C),(D,E)] -- i.e. one of the valid combinations
        mhs_cur = []
        for x in range(len(product_iter)): mhs_cur.append(tuple([*product_iter[x]])) 
        mhs_cur_check = [x[:2] for x in mhs_cur]
        if len(set(mhs_cur_check)) != len(mhs_cur):
            # We have something like [(a,b,0),(a,b,1)] in mhs_cur (possible if 
            # a is in two composite nodes targeting b). Skip it
            continue
        
        # Make a copy of the expanded graph and reduced graph for the particular
        # set of changes we are considering here
        G_erc = G_expanded_reduced.copy()
        G_orc = G_original_reduced.copy()
    
        # Use set() because you could have e.g. 
        # mhs = [(ab,d),(ac,d)] -> (a,d),(a,d) 
        for edge_cut in set(mhs_cur):           
            # Identify the target and the complement of the target
            target = edge_cut[1][:]
            if target[0] == '~': rev_target = target[1:]
            else: rev_target = '~'+target[:]
            # Identify the source and the complement of the source
            
            source = edge_cut[0][:]           
            if source[0] == '~': rev_source = source[1:]
            else: rev_source = '~' + source[:]
            
            # Update the copy of the original network in response to this change
            if target[0] == '~':
                target_node = target[1:]
            else:
                target_node = target[:]
            
            if source[0] == '~':
                source_node = source[1:]
            else:
                source_node = source[:]
            
            # Modify the edges, update_rules of the original version of the network 
            
            # are we setting the node True or False? Depends on if we're 
            # setting the full node (negation or reg) 0 or 1
            if (edge_cut[2] == 0 and source[0] == '~')\
            or (edge_cut[2] == 1 and source[0] != '~'):
                # Node becomes ON: We are here setting the node negation OFF 
                # -or- the node itself ON
                rule_split = G_orc.node[source_node]['update_rules'].split()  
                new_rule = [x.replace(source_node,'True') if x.strip('() ') 
                            == source_node else x for x in rule_split]
                new_rule = ' '.join(new_rule)
                G_orc.node[target_node]['update_rules'] = new_rule
            elif(edge_cut[2] == 1 and source[0] == '~' )\
            or (edge_cut[2] == 0 and source [0] != '~'):
                # Node becomes OFF: We are here setting the node negation ON 
                # -or- the node itself OFF
                rule_split = G_orc.node[source_node]['update_rules'].split()  
                new_rule = [x.replace(source_node,'False') if x.strip('() ')
                            == source_node else x for x in rule_split]
                new_rule = ' '.join(new_rule)
                G_orc.node[target_node]['update_rules'] = new_rule
            else:
                # Should never see this
                raise RuntimeError("Invalid command in forcing node values.")

            # Also need to remove the edge
            
            G_orc.remove_edge(source_node,target_node)

            # Modify the edges in the expanded version of the network -----
            # Need to treat both the target node and the rev_target node 

            for targ in [target,rev_target]:
                preds = [x for x in G_erc.predecessors(targ)]
                for pred in preds:
                    # this predecessor can be either a node or a composite node
                    pred_list = pred.split()
                    # pred_list is a list version of each node (single entry 
                    # or 2+ if composite)
                    
                    # Below block fires if we have a node -> node interaction
                    if len(pred_list) == 1:
                        # source is a predecessor
                        if source in pred_list:
                            # source = OFF, so remove edge
                            if edge_cut[2] == 0: G_erc.remove_edge(source,targ)
                            # source = ON, so all other edges go away and targ 
                            # picks up a self-loop
                            else: 
                                cur_preds = [x for x in G_erc.predecessors(targ)]
                                G_erc.remove_edges_from(zip(cur_preds,[targ]*len(cur_preds)))
                                G_erc.add_edge(targ,targ)
                                # No need to consider other predecessors; stop
                                break
                        # rev_source (=ON) is a predecessor; this fixes the target
                        # to ON
                        elif rev_source in pred_list:
                            #rev_source = OFF, so remove ege
                            if edge_cut[2] == 1: G_erc.remove_edge(rev_source,targ)
                            #rev_source = ON, so all other edges go away and 
                            # targ picks up a self-loop
                            else:
                                cur_preds = [x for x in G_erc.predecessors(targ)]
                                G_erc.remove_edges_from(zip(cur_preds,[targ]*len(cur_preds)))
                                G_erc.add_edge(targ,targ)
                                # No need to consider other predecessors; stop
                                break
                    
                    # below block fires if we have a composite -> node interaction
                    else:
                        if source in pred_list:
                            # source (=OFF) is in a predecessor composite node;
                            # this kills the entire composite node
                            if edge_cut[2] == 0: 
                                G_erc.remove_edge(pred,targ)
                                if G_erc.out_degree(pred) == 0: G_erc.remove_node(pred)
                            # source (=ON) is in a predecessor composite node;
                            # this reduces the composite
                            else:
                                G_erc = cf.remove_from_composite(G_erc,pred,source,targ)
                        if rev_source in pred_list:
                            # rev_source (=OFF) is in a predecessor composite node;
                            # this kills the entire composite node
                            if edge_cut[2] == 1: 
                                G_erc.remove_edge(pred,targ)
                                if G_erc.out_degree(pred) == 0: G_erc.remove_node(pred)
                            # rev_source (=ON) is in a predecessor composite node;
                            # this reduces the composite
                            else:
                                G_erc = cf.remove_from_composite(G_erc,pred,rev_source,targ)
            
        
        # Once we get here, the expanded graph has been revised in response to
        # the edge cuts we are here considering. Now we need to look at the 
        # topology around this expanded graph to ID the attractors.
        
        # Essentially we are repeating steps 3 and 4 of the original analysis, above
        
        SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_erc)
        
        fixed_attractors, frozen_states = cf.iterate_SCC(G_orc,G_erc,SCC_list,
                                                         frozen={},removed={}, 
                                                         dump_to_screen=False,
                                                         output=False)
        
        # And look at oscillations
        osc_components = cf.find_oscillations(G_erc,SCC_list)
        
        # Now summarize the result of this set of edge cuts
        print('\nConsidering cuts:',product_iter)
        
        # Some of the frozen states could be duplicates; remove them      
        frozen_states_culled = []
        for fs in frozen_states:
            if fs not in frozen_states_culled: frozen_states_culled.append(fs)
        
        # Final reporting
        print('\t# oscillating components = {0}, # frozen states = {1}'\
              .format(len(osc_components),len(frozen_states_culled)))
        
        for i,fs in enumerate(frozen_states_culled):
            if 'Closure' in fs: 
                print('\tfrozen state {0} of {1} freezes {2} nodes; Closure = {3}'
                      .format(i+1,len(frozen_states_culled),len(fs),fs['Closure']))
            else: 
                print('\tfrozen state {0} of {1} freezes {2} nodes; Closure oscillates'
                      .format(i+1,len(frozen_states_culled),len(fs)))
        
TIME7 = time.time()
print('STEP 7 completed in {0:.2f} minutes ({1:.2f} total).'
      .format((TIME7-TIME6)/60.0,(TIME7-TIME0)/60.0))
