"""
For analysis of the reduced T-LGL network. 

This script sequentially walks through the main components of the expanded 
graph analysis described in the report.

Author: Colin Campbell
Date: November 2018
Language: Python 3.6
"""

import networkx as nx
import custom_functions as cf
from itertools import product

f = open(r'../Network Rules/LGL.txt','rt')
lines = f.readlines()
f.close() 


# ANALYSIS OF ORIGINAL NETWORK ================================================
# STEP 1 - FORM NETWORKS -------------------------------------------------
G_expanded = cf.form_expanded_network(lines)
G_original = cf.form_network(lines)

# print edges in expanded network?
#print('Edges in the expanded network:')
#for i,j in G_expanded.edges():
#    print(i,';',j)

# STEP 2 - REDUCE THE NETWORK GIVEN OUR INITIAL CONDITIONS ---------------
status, G_original_reduced, frozen, removed  = cf.reduce_network(G_original,SCC=False)
G_expanded_reduced = cf.form_expanded_network_reduced(G_expanded,fdict=frozen,rdict=removed)

# STEP 3 - FIND THE FIXED-POINT SCCS -------------------------------------
# Require SCCs where (1) a node and a complementary node are -not- both present
# and (2) the inputs of any composite nodes are also included

SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_expanded_reduced)

for SCC in SCC_list:
    print('Found a SCC:')
    print(SCC)
fixed_attractors, frozen_states = cf.iterate_SCC(G_original_reduced,
                                                 G_expanded_reduced, SCC_list,
                                                 frozen=frozen, removed = removed, 
                                                 dump_to_screen = True)

# Here we find three SCCs; the corresponding attractors are the known/original
# attractors of the network (healthy vs cancerous)

# STEP 4 - FIND THE OSCILLATING COMPONENTS -------------------------------
    
# Now find valid oscillations
osc_components = cf.find_oscillations(G_expanded,SCC_list)

for i,osc in enumerate(osc_components):
    print('Oscillating components {0} of {1}:'.format(i+1,len(osc_components)))
    for j in osc:
        if ' ' not in j: print(j)    

# We find none, as expected.


# BEGIN PERTURBATION ANALYSIS =================================================
# STEP 5 - ANALYZE PARENT SCCS --------------------------------------
MHS_perturbations = set()

# Form the parent SCC expanded graph
parent_graph = G_expanded_reduced.copy()
remove_nodes = set(parent_graph.nodes())-set(parent_SCCs[1])
parent_graph.remove_nodes_from(remove_nodes)

# To ensure complete coverage, we use as source every node that is frozen in a 
# constituent stable SCC.
source_nodes = []
for scc in SCC_list:
    source_nodes += [x for x in scc if ' ' not in x]
source_nodes = list(set(source_nodes))
# Don't care about Apoptosis since we are looking at the other parent SCC
source_nodes.remove('Apoptosis')               

# We are now ready to look for the valid circuits
path_store = []
for source in source_nodes:
    print('Considering source = {0}'.format(source))
    count = 0
    for path in cf.parent_SCC_search(parent_graph, source):
        # We have a new path, but we want to only store it if it is a unique 
        # addition. We also don't want long circuits that contain a complete 
        # valid circuit (if we destroy the smaller we must destroy the bigger, 
        # but the reverse is not true)
        
        # We could have [a,b,c] and [b,c,a] - i.e. different orderings
        # so we first check the lengths of the lists: if so then we dbl one and
        # see if the other is a sub-component
        new = True
        
        for p in path_store:
            # Transfer these to strings to make substring searching easier
            
            p_store = [str(x) for x in p]
            p_store = ''.join(p_store)
            p_store *= 2  # double so when we search we avoid 'wraparound' permutations
            
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
        # include the -1 to 0 index pair since we have a cycle but don't
        # include the start/end node in the last position
        if ' ' in path[j]:
            # skip this because the target node is a composite (type 1)
            continue
        # This bit new in test4 cf test3: split composites
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
print('Minimum Hitting Sets (omitting Apoptosis):')
for i in range(len(MHS)):
    print('# {0} of {1}:'.format(i+1,len(MHS)))
    print(MHS[i])
    print('\n')

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
        if 'Apoptosis' in entry[1]: continue  # We don't consider modifying Apoptosis
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
            # we have something like [(a,b,0),(a,b,1)] in mhs_cur (possible if 
            # a is in two composite nodes targeting b). Skip it
            continue
        
        # Make a copy of the expanded graph and reduced graph for the particular
        # set of changes we are considering here
        G_erc = G_expanded_reduced.copy()
        G_orc = G_original_reduced.copy()
    
        for edge_cut in set(mhs_cur):   # use set() because you could have e.g.
            # mhs = [(ab,d),(ac,d)] -> (a,d),(a,d)         
            # identify the target and the complement of the target
            target = edge_cut[1][:]
            if target[0] == '~': rev_target = target[1:]
            else: rev_target = '~'+target[:]
            # identify the source and the complement of the source
            
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
            
            # Modify the edges, update_rules of the original version of the network -----
            
            # are we setting the node True or False? Depends on if we're setting 
            # the full node (negation or reg) 0 or 1
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
                # Node becomes OFF: We are here setting the node negation ON or
                # the node itself OFF
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
                    # pred_list is a list version of each node (single entry or 
                    # 2+ if composite)
                    
                    # below block fires if we have a node -> node interaction
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
                            #rev_source = ON, so all other edges go away and targ picks up a self-loop
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
        # the edge cuts we are here considering. Now we need to look at the topology
        # around this expanded graph to ID the attractors.
        
        # Essentially we are repeating steps 3 and 4 of the original analysis, above
        
        SCC_list, parent_SCCs, partial_cycles = cf.find_fixed_SCCs_iterative(G_erc)
        # Here we store the output of iterate_SCC to see all of the attractors
        # as a list of node states (e.g. ['10?','000','??1'])
        
        # To carry through iterate_SCC below (and see how the SCCs propagate through
        # the network), we need an updated version of the non-expanded networkk
        fixed_attractors, frozen_states = cf.iterate_SCC(G_orc,G_erc,SCC_list,
                                                      frozen={},removed={}, 
                                                      dump_to_screen=False,
                                                      output=True)
        
        # And look at oscillations
        osc_components = cf.find_oscillations(G_erc,SCC_list)
        
        # For this network, a valid mhs set of edge cuts results in only 1 fixed point
        # attractor and no oscillations. If we meet these criteria, announce it
        
        if osc_components == [] and len(fixed_attractors) == 1:
            print('Found a valid MHS!')
            print('MHS =',mhs)
            print('Specific interaction changes =',mhs_cur)
            print('Results in the following fixed-point attractor:')
            print(''.join(fixed_attractors[0]))
            print('\n')
            
            stop_perturbation = False
            # Store this for comparison to explicit enumeration
            
            new_pert = set()
            for s in set(mhs_cur):
                if s[1][0] == '~': tgt = s[1][1:]
                else: tgt = s[1][:]
                if s[0][0] == '~':
                    src = s[0][1:]
                    val = (s[2]+1)%2   # we set e.g. ~A to 1, so A is set to 0
                else:
                    src = s[0][:]
                    val = s[2]
                # Just here at the end, we exclude Apoptosis as source or target
                if 'Apoptosis' in src or 'Apoptosis' in tgt: 
                    stop_perturbation = True
                    break
                
                new_pert.add(tuple([src,tgt,val]))
            if stop_perturbation:
                print('(Omitting from final count because Apoptosis is involved)')
                continue
            perturbation_counter += 1
            MHS_perturbations.add(frozenset(new_pert))
        else:
#            print('reduction not successful:',product_iter,master)
            failed_combinations += 1

print('Found {0} valid pertubations.'.format(perturbation_counter))


