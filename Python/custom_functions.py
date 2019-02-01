'''
Custom functions for carrying out the edgetic perturbation analysis. To be 
imported into the main analysis files.

Brief contents:
    form_network(): Make a NetworkX DiGraph (standard format) from a .txt of
    rules.
    
    form_expanded_network(): Make a NetworkX DiGraph (expanded format) from a 
    .txt of rules.
    
    form_expanded_network_reduced(): Reduce an expanded network when some nodes
    are frozen or removed.
    
    check_valid(): Evaluates if a given cycle on an expanded graph meets the
    criteria for a stable SCC.
    
    perturb_network(): Updates a NetworkX DiGraph (standard format) in response
    to a given edgetic perturbation.
    
    perturb_expanded_network(): Updates a NetworkX DiGraph (expanded format) in
    response to a given edgetic perturbation.
    
    find_fixed_SCCs_iterative(): Identifies stable SCCs on an expanded graph.
    
    reduce_network(): Identify frozen nodes and those that can be removed 
    (chains and sinks) from an expanded graph.
    
    find_oscillations(): Identify the oscillating components in an expanded
    graph.
    
    parent_SCC_search(): Generator function to find valid circuits on a parent
    graph.
    
    remove_from_composite(): Update an expanded graph by removing a node from
    a composite node and updated the graph accordingly.
    
    iterate_SCC(): Freezes nodes and allows the effects to propagate through an
    expanded graph.
    
    MMCS(): Implementation of the MMCS algorithm for determining the minimum 
    hitting set on a collection of sets.
    
Author: Colin Campbell
Date: November 2018
Language: Python 3.6
'''

import networkx as nx
import time
from random import randint

# FUNCTION DEFINITIONS ========================================================

def form_network(rules):
    """Create a NetworkX DiGraph (standard rep.) from a list of rules.
    
    Keyword arguments:
        rules:  Lines from a text file specifying the update rules for a 
                network.
   
    Returns:
        g:      A networkx DiGraph of the network. The update rule for each 
                node is accessed via g.node[n]['update_rules'] = <string>
    """
    
    # Remove comments and blank lines
    stream = [x.rstrip('\n') for x in rules if x != '\n' and x[0]!='#']         

    g = nx.DiGraph()

    for line in stream:
        node,rule = line.split('* = ') 
        # Here we don't care about the node rule negation 
        if node[0] == '~': continue                                             
    
        g.add_node(node)
    
        input_nodes = rule.replace(' AND ',' ')
        input_nodes = input_nodes.replace(' OR ', ' ')
        input_nodes = input_nodes.replace('NOT ', ' ')
        input_nodes = rule.replace(' and ',' ')
        input_nodes = input_nodes.replace(' or ', ' ')
        input_nodes = input_nodes.replace('not ', ' ')
        input_nodes = input_nodes.replace('(', '')
        input_nodes = input_nodes.replace(')', '')

        input_nodes = input_nodes.split() 

        for input_node in input_nodes:
            if input_node != 'False' and input_node != 'True':
                g.add_edge(input_node,node)
    
        g.node[node]['update_rules'] = rule.strip()
    
    return g
    

def form_expanded_network(rules):
    """Create a NetworkX DiGraph (expanded rep.) from a list of rules.
    
    Keyword arguments:
        rules:  Lines from a text file specifying the update rules for a 
                network.
   
    Returns:
        g:      A networkx DiGraph of the expanded network.
    """

    # Remove comments and blank lines
    stream = [x.rstrip('\n') for x in rules if x != '\n' and x[0]!='#']         

    g = nx.DiGraph()

    for line in stream:
        node,rule = line.split('* = ') 
        
        g.add_node(node)
        if rule.lower().strip() == 'true' or rule.lower().strip() == 'false': 
            # If the node state is fixed, then there is no need to assess 
            # incoming edges
            continue
        
        rule = rule.replace(' AND ',' and ')
        rule = rule.replace(' OR ',' or ')
        rule = rule.replace('NOT ','not ')

        terms = rule.split(' or ')                                         

        # Each term separated by 'or' will be a series of ANDs (with perhaps 
        # some node negations inside). It forms a composite node.        
        for term in terms:
            t = term.strip('() ')
            if '(' in t or ')' in t:
                print('bug')
            t= t.split(' and ')
            composite = ''
            for cnode in t:
                if cnode[:4] == 'not ': composite += '~'+cnode[4:]
                else: composite += cnode
                composite += ' '
            composite = composite.strip()
            # Add composite node -> node edge
            g.add_edge(composite,node)  
            # Also need the input nodes to the composite 
            # (if it truly is composite)
            if len(t) > 1:
                for inp_node in composite.split():
                    g.add_edge(inp_node,composite)
                        
    return g

def form_expanded_network_reduced(input_g, fdict = {}, rdict = {}):
    """Reduce an expanded network in response to some nodes being frozen/lost.
    
    Nodes in fdict and rdict are removed from the expanded network; remaining
    nodes and edges are updated (e.g. A->B->C becomes A->C if B is removed).
    
    Keyword arguments:
        input_g:    The expanded graph.
        fdict:      A dictionary of {node:state} pairs of frozen nodes.
        rdict:      A dictionary of {node:details} pairs of removed nodes. See
                    iterate_SCC() for information on the form of details.
                    
    Returns: 
        g:          The updated expanded network.
    """
    
    g = input_g.copy()

    remove_self_loop = set()
    preserve_self_loop = set()

    # FIRST LOOK AT NODES AND NODE NEGATIONS -----

    # Handle nodes in rdict: connect target to predecessor
    # First deal with the node; following block deals with ~node.
    # (This is somewhat repetitive...)
    for node in rdict.keys():
        if node not in g.nodes(): continue
        if rdict[node] == 'removed_as_sink':
            if g.out_degree(node) != 0: 
                raise RuntimeError('node removed as sink has nonzero out degree!')
            else: g.remove_node(node)
        else:
            if g.out_degree(node) != 1 or g.in_degree(node) != 1: 
                # The below Exception is worth including unless we are removing 
                # interactions via the MHS algorithm. If we -are-, some nodes 
                # flagged for removal b/c of a chain A-B-C structure can have 
                # nodes A and/or B artificially removed already.
#                raise RuntimeError('node removed from chain has in-degree and/or out-degree != 1!')
                
                # If we are allowing for mhs cuts, then if A and/or C are 
                # missing, we can just remove the node.
                g.remove_node(node)
            else: 
                # Remove node and add edge from dependency to target
                succ = list(g.successors(node))[0]
                g.remove_node(node)
                g.add_edge(rdict[node][0],succ)
     
    for node in rdict.keys():
        if '~'+node not in g.nodes(): continue
        if rdict[node] == 'removed_as_sink':
            if g.out_degree('~'+node) != 0: 
                raise RuntimeError('node removed as sink has nonzero out degree!')
            else: g.remove_node('~'+node)   
        else: 
            if g.out_degree('~'+node) != 1 or g.in_degree('~'+node) != 1: 
                # See comment in the block above. The exception below should 
                # only be included if not considering MHS edge cuts
#                raise RuntimeError('node removed from chain has in-degree and/or out-degree != 1!')
                g.remove_node('~'+node)
            else:
                # Remove node and add edge from dependency to target
                succ = list(g.successors('~'+node))[0]
                g.remove_node('~'+node)
                g.add_edge('~'+rdict[node][0],succ)

    # Remove frozen nodes and their negations   
    g.remove_nodes_from(fdict.keys())
    for node in fdict.keys():
        if '~'+node in g.nodes: g.remove_node('~'+node)
    
    
    # NOW UPDATE COMPOSITE NODES -----
    update_dict = {}
    nodes = list(g.nodes())
    for node in nodes:
        if ' ' not in node: continue    # composite nodes all have a space
        
        # isolated composite nodes can be removed
        if len([x for x in g.predecessors(node)]) == 0 \
        or len([x for x in g.successors(node)])==0:         
            g.remove_node(node)
            continue
        
        cnodes = node.split()
        
        deleted_node = False
        new_label = ''    
        for cnode in cnodes:
            # First check to see if a member of this composite node has its 
            # state fixed, in which case we can simplify the composite node
            if cnode[0] == '~' and cnode[1:] in fdict:
                # cnode is e.g. ~A and A is in fdict
                # We have an "AND 1", which we can remove from the label
                if (fdict[cnode[1:]] + 1)%2 == 1: 
                    continue   
                # We have an "AND 0", which kills the node
                else: 
                    g.remove_node(node)                 
                    deleted_node = True
                    break
    
            elif cnode in fdict:
                # cnode is e.g. A and A is in fdict
                # We have an "AND 1", which we can remove from the node label
                if fdict[cnode] == 1: 
                    continue          
                # We have an "AND 0", which kills the node
                else:
                    g.remove_node(node)                 
                    deleted_node = True
                    break
            
            # Now check to see if a member of this composite node has been 
            # removed with its state dependent on some other node: the 
            # composite node is updated with the upstream reference
            elif cnode in rdict:
                if rdict[cnode] == 'removed_as_sink': 
                    raise RuntimeError('Node removed as sink exists in composite node!')
                else: 
                    # A node in the composite node has been removed. Need to 
                    # follow the dependency chain until we find a node that is 
                    # not in rdict
                    neg_counter = rdict[cnode][1]
                    cur_targ = rdict[cnode][0]
                    while  cur_targ in rdict.keys():
                        cur_targ = rdict[cur_targ][0]
                        neg_counter += rdict[cur_targ][1]
                    # Now we have the node that the original composite node 
                    # depends upon (that has not also been removed).
                    # Is it frozen?
                    if cur_targ in fdict.keys():
                        # It is frozen; need to determine if we have 'AND 1' 
                        # or 'AND 0'
                        if (fdict[cur_targ] == 1 and (neg_counter +1)%2 == 1)\
                        or (fdict[cur_targ] == 0 and (neg_counter +1)%2 == 0):
                            # The target is frozen ON and the original 
                            # composite node is therefore ON. Or, the target is
                            # frozen OFF and the original composite is 
                            # therefore ON. We effectively have 'AND 1' and 
                            # don't need to append anything to the node label
                            continue
                        else:
                            # If the above don't fire, we have the opposite 
                            # cases: 'AND 0' destroys the node
                            g.remove_node(node)
                            deleted_node = True
                            break
                    else:
                        # It is not frozen; do we include 'node' or '~node'?
                        if (neg_counter +1)%2 == 1:
                            # 0, 2, 4... 'nots' in the dependency chain, so 
                            # direct relationship
                            new_label += cur_targ
                        else:
                            new_label += '~'+cur_targ
            elif cnode[0] == '~' and cnode[1:] in rdict:
                # Repeating above block for case of node negation
                if rdict[cnode[1:]] == 'removed_as_sink': 
                    raise RuntimeError('Node removed as sink exists in composite node!')
                else: 
                    # A node in the composite node has been removed. Need to 
                    # follow the dependency chain until we find a node that is 
                    # not in rdict
                    neg_counter = rdict[cnode[1:]][1]
                    cur_targ = rdict[cnode[1:]][0]
                    while  cur_targ in rdict.keys():
                        cur_targ = rdict[cur_targ][0]
                        neg_counter += rdict[cur_targ][1]
                    # Now we have the node that the original composite node 
                    # depends upon (that has not also been removed).
                    # Is it frozen?
                    if cur_targ in fdict.keys():
                        # It is frozen; need to determine if we have 'AND 1' 
                        # or 'AND 0'. cCntrary to above block we are dealing 
                        # with the original composite being a negation
                        if (fdict[cur_targ] == 0 and (neg_counter +1)%2 == 1)\
                        or (fdict[cur_targ] == 1 and (neg_counter +1)%2 == 0):
                            # The target is frozen ON and the original 
                            # composite node is therefore ON. Or, the target is 
                            # frozen OFF and the original composite is 
                            # therefore ON. We effectively have 'AND 1' and 
                            # don't need to append anything to the node label.
                            continue
                        else:
                            # If the above don't fire, we have the opposite 
                            # cases: 'AND 0' destroys the node.
                            g.remove_node(node)
                            deleted_node = True
                            break
                    else:
                        # it is not frozen; do we include 'node' or '~node'?
                        if (neg_counter +1)%2 == 1:
                            # 0, 2, 4... 'nots' in the dependency chain, so 
                            # direct relationship. But we're -starting- from 
                            # a ~ composite node, so we end with one too
                            new_label += '~' + cur_targ
                        else:
                            new_label += cur_targ
            
            # If none of the above, just add the current composite
            else:
                new_label += cnode
            new_label += ' '
        if not deleted_node:
            new_label = new_label.strip()
            if ' ' not in new_label: 
                # This composite node has become a regular node and without 
                # intervention we will add a self-loop. However, we only want 
                # to create a self-loop when it is a genuine  interaction s.t. 
                # we had e.g. AB -> A  (and not only A -> AB, which is just a 
                # logical requirement). The only exception would be if we had 
                # AB -> A and also A->A (shouldn't happen in practice)
                if not g.has_edge(node,new_label) and not g.has_edge(new_label,new_label):
                    remove_self_loop.add(new_label)
                if g.has_edge(node,new_label):
                    # But this could happen more than once for a 'new_label', 
                    # and we could come to differing conclusions about whether
                    # the self-loop shoud be preserved. Note "definitely keep" 
                    # here.
                    preserve_self_loop.add(new_label)
            else:
                # Rarely because of removed nodes we can end up with a 
                # composite node that has two instances of the same node (e.g. 
                # '~A ~A B'). To avoid that we here remove duplicate while 
                # preserving order.
                new_label_spl = new_label.split()
                new_label = []
                for x in new_label_spl:
                    if x not in new_label: new_label.append(x)
                new_label = ' '.join(new_label)
            update_dict[node] = new_label
        
    # Carry out the renaming of the composite node labels
    g = nx.relabel_nodes(g,update_dict, copy=False)

    # Carry out the elimination of self-loops noted above   
    for node in remove_self_loop:
        if g.has_edge(node,node): g.remove_edge(node,node)
    # But preserve those noted above
    for node in preserve_self_loop:
        g.add_edge(node,node)
    
    return g

def check_valid(cycle):
    """See if a particular cycle meets our criteria."""
    # First check non-composite nodes 
    c = cycle[:]
    while len(c) > 0:
        n = c.pop()
        if ' ' in n: continue   # skipping composites here
        if '~' == n[0]: 
            if n[1:] in cycle:
                # Invalid: We have e.g. 'A' and '~A'
                return 'invalid node combination'
               
        elif '~'+n in cycle:
            # Invalid: We have e.g. 'A' and '~A'
            return 'invalid node combination'
    
    # Now check composite nodes
    for i in range(len(cycle)):
        n = cycle[i]
        if ' ' in n:
            # a composite node
            parents = n.split(' ')
            for p in parents:
                if p not in cycle:
                    # Invalid: the parents of a comp. node are not in the
                    # path
                    return 'missing parents'
    
    # If none of the above return statements fire, we have a valid cycle
    return 'valid' 

def perturb_network(g_inp, pert=[]):
    """Updates a network by implementing an edgetic perturbation.
    
    Keyword arguments:
        g_inp:      The input network (A NetworkX DiGraph)
        pert:       A list of dicts in the format 
                    [{k:{'node':v,'value':<0,1>},...] where k is a node (not a 
                    composite or negation) whose incoming edges are modified 
                    such that node v is assumed to take on the indicated value.
    
    The function removes the edge from v to k and updates the rule for k in
    g.node[n]['update_rules'] = <string>
    accordingly.
    
    Returns:
        g:          The updated network.
        
    See also perturb_expanded_network(), which does the same operation on an
    expanded network.
    """
    
    g = g_inp.copy()
    
    for entry in pert:
        # Each entry only has one k:{} pair, so this loop fires 1x
        for node,d in entry.items():                                            
            # Remove the edge
            g.remove_edge(d['node'],node)
            # Update the rule for the node
            g.node[node]['update_rules'] = \
            g.node[node]['update_rules'].replace(d['node'],str(d['value']))     
      
    return g

def perturb_expanded_network(g,pert={}):
    """Updates an expanded network by implementing an edgetic perturbation.
    
    Keyword arguments:
        g_inp:      The input network (A NetworkX DiGraph)
        pert:       A list of dicts in the format 
                    [{k:{'node':v,'value':<0,1>},...] where k is a node (not a 
                    composite or negation) whose incoming edges are modified 
                    such that node v is assumed to take on the indicated value.
     
    The function updates the incoming edges for both node k and its negation.
     
    If a node k's input rules reduce to 0 as a result (meaning ~k = 1):
        We could add an edge from k to ~k, but we don't (irrelevant for
        stable/oscillating motif analysis) 
        
        A self loop is added from ~k to ~k
        (being in state k=0 is persistent)
    If a node k's input rules reduce to 1, the opposite of the above occurs
    (and likewise for ~k)
    
    Returns:
        g_copy:     The updated expanded network
        
    See also perturb_network(), which does the same operation on a standard 
    network.
    """
    
    def analyze_node(k):
        """Update node k"""
        for pred in g2.predecessors(k):                  
            if ' ' in pred:
                # if a composite node, look at each component
                src_nodes = pred.split(' ')       
                # First check if we have (... AND 1...) in this composite node,
                # meaning we can just remove the node v from the composite node
                check_against_composites = False
                if '~'+v in src_nodes and state == 0:
                    new_pred_nodes = [x for x in src_nodes if x!='~'+v]
                    new_pred = [x + ' ' for x in new_pred_nodes]
                    new_pred = ''.join(new_pred).rstrip()
                    if ' ' in new_pred and not g_copy.has_node(new_pred):
                        # If this is a new composite node, need to add its incoming edges
                        for node in new_pred_nodes: g_copy.add_edge(node,new_pred)
                    g_copy.remove_edge(pred,k)
                    # Remove the composite node if unused
                    if len(list(g_copy.successors(pred))) == 0: deletion_list.append(pred)  
                    check_against_composites = True
                    
                elif v in src_nodes and state == 1:
                    new_pred_nodes = [x for x in src_nodes if x!=v]
                    new_pred = [x + ' ' for x in new_pred_nodes]
                    new_pred = ''.join(new_pred).rstrip()
                    if ' ' in new_pred and not g_copy.has_node(new_pred):
                        # If this is a new composite node, need to add its incoming edges
                        for node in new_pred_nodes: g_copy.add_edge(node,new_pred)
                    g_copy.remove_edge(pred,k)
                    if len(list(g_copy.successors(pred))) == 0: deletion_list.append(pred)
                    check_against_composites = True
                
                # Then check if we have (...AND 0...) in this composite node, 
                # meaning the entire composite node dies
                elif ('~'+v in src_nodes and state == 1) or (v in src_nodes and state == 0):
                    g_copy.remove_edge(pred,k)
                    if len(list(g_copy.successors(pred))) == 0: deletion_list.append(pred)
    
                if check_against_composites:
                    # We've fired one of the first two blocks above, meaning 
                    # our original composite node has simplified. As a result, 
                    # a separate composite node parent of node k may have 
                    # become superfluous (and should be removed if so)
                    pred_comp_nodes = [x for x in g_copy.predecessors(k) if ' ' in x]
                    for pcn in pred_comp_nodes:
                        pcn_split = pcn.split()
                        if all([x in pcn_split for x in new_pred_nodes]):
                            # If the new composite node is a subset of this 
                            # composite node, this composite node is logically 
                            # superfluous and can be removed
                            g_copy.remove_edge(pcn,k)
                            if g_copy.out_degree(pcn) == 0: g_copy.remove_node(pcn)
                                    
                    # Finally, add the new edge
                    if new_pred == '~Apoptosis' and k == '~Apoptosis':
                        print('debug')
                    g_copy.add_edge(new_pred,k)
    
            elif pred == v or pred == '~' + v:
                # this input is just a singleton node
                # see if we have an (... OR 1...) or (...OR 0...)
                if ('~' not in pred and state == 1) or ('~' in pred and state == 0):
                    # this means we simply have a self loop
                    all_preds = [x for x in g_copy.predecessors(k)]
                    for x in all_preds: g_copy.remove_edge(x,k)
                    g_copy.add_edge(k,k)
                    break
                elif ('~' not in pred and state == 0) or ('~' in pred and state == 1):
                    # we just remove the edge ... if all are removed we have
                    # a node with no incoming edges, meaning it is never activated
                    g_copy.remove_edge(pred,k)
                
    g_copy = g.copy()
    g2 = g.copy()
    
    for entry in pert:
        # Each entry only has one k:{} pair, so this loop fires 1x
        for k,d in entry.items():                     
            v = d['node']
            state = d['value']
            
            # Treat the original node -----
            deletion_list = []
            analyze_node(k)
            # Remove unused composite nodes
            for node in set(deletion_list): g_copy.remove_node(node)  
        
            # Treat the negation of the node -----
            deletion_list = []
            analyze_node('~'+k)
            for node in set(deletion_list): g_copy.remove_node(node)        
        
        g2 = g_copy.copy()
        
    return g_copy


def find_fixed_SCCs_iterative(g):
    """Identify stable SCCs in an expanded graph.
    
    Viable SCCs obey the properties:
        1. A node and a complementary node are -not- both present
        2. The inputs of any composite nodes are also included.    
    
    The algorithm works by considering combinations of cycles on the expanded
    graph.
    
    Keyword arguments:
        g:              An expanded graph.
    
    Returns:
        valid_SCCs:     A list of lists. Each interior list contains the nodes 
                        (and composite nodes) present in a smallest viable SCC.
                        A smallest viable SCC means no subset of the nodes are
                        themselves a viable SCC.
        parent_SCCs:    A list of lists. Each interior list contains the nodes
                        (and composite nodes) present in a parent SCC. A parent
                        SCC means that its nodes are not a subset of another,
                        larger SCC.
        invalid_cycles: A list of cycles on the expanded graph that were 
                        rejected from forming a viable SCC because parents of
                        an included composite node were not present. 
    """
    def check_missing_parents(cycle):
        """
        Looks at composite nodes; outputs list of necessary input nodes  that 
        are missing.
        """
        out = []
        for i in range(len(cycle)):
            n = cycle[i]
            if ' ' in n:
                # a composite node
                parents = n.split(' ')
                for p in parents:
                    if p not in cycle:
                        out.append(p)
        return out
    
    def compare_parents(x):
        """
        x is a list of non-composites from two cycles. We here check to make 
        sure we don't have node/negation pairs.
        """
        c = x[:]
        while len(c) > 0:
            n = c.pop()
            
            if '~' == n[0]: 
                if n[1:] in x:
                    return 'invalid node combination'
                   
            elif '~'+n in x:
                return 'invalid node combination'
        return 'valid'
    
    valid_SCCs = []
    invalid_cycles = []
    parent_SCCs = []
    
    # First, look at individual cycles
    counter = 0
    TIME0 = time.time()
    for cycle in nx.simple_cycles(g):
        counter += 1
        if counter % 100000 == 0:
            TIME1 = time.time()
            print('Starting cycle {0} at {1:.2f} minutes elapsed.'
                  .format(counter,(TIME1-TIME0)/60.0))
        
        if len(cycle) >= 40: continue
        
        status = check_valid(cycle)
        if status == 'valid': 
            # Found a stable motif - make sure it isn't a duplicate or superset 
            # of another motif
            if not any([set(cycle) >= set(x) for x in valid_SCCs]):
                # must also consider reverse case to avoid duplicates
                valid_SCCs = [x for x in valid_SCCs if not set(x).issuperset(cycle)]
                valid_SCCs.append(cycle)


            # Now check to see if this SCC is a superset of an existing valid 
            # SCC. First make sure it isn't a subset of a different parent
            if not any([set(cycle) <= set(x) for x in parent_SCCs]):
                # This is a new "big one."
                # We'll track it as a parent, but remove any other parents 
                # that are a subset of this one.
                parent_SCCs = [x for x in parent_SCCs if not set(x).issubset(cycle)]
                parent_SCCs.append(cycle)
            
        elif status == 'missing parents': 
            invalid_cycles.append(cycle)

    
    # Now we consider invalid cycles due to missing parents of composite nodes.
    # At first we consider all pairs of cycles.
    
    SCCs_check = []
    for c1 in range(len(invalid_cycles)-1):
        missing_parents1 = check_missing_parents(invalid_cycles[c1])
        present_parents1 = [x for x in invalid_cycles[c1] if ' ' not in x]
        for c2 in range(c1+1,len(invalid_cycles)):
            
            # If one is the subset of the other, no chance for improvement here
            if set(invalid_cycles[c2]) <= set(invalid_cycles[c1])\
            or set(invalid_cycles[c1]) <= set(invalid_cycles[c2]): continue   
              
            missing_parents2 = check_missing_parents(invalid_cycles[c2])
            present_parents2 = [x for x in invalid_cycles[c2] if ' ' not in x]
            
            status = compare_parents(present_parents1 + present_parents2)  
            
            if status == 'valid':
                # We get here if there are no node/antinode pairs in this combination
                # Below is non-empty if we have found a new parent
                overlap1 = set(missing_parents1) & set(present_parents2) 
                overlap2 = set(missing_parents2) & set(present_parents1) 
                if len(overlap1) > 0 or len(overlap2) > 0:
                    # We've added a parent, yay -- now need to make sure we 
                    # still have a SCC
                    g_copy = g.copy()
                    excess_nodes = set(g.nodes()) - (set(invalid_cycles[c1]) | set(invalid_cycles[c2]))
                    g_copy.remove_nodes_from(excess_nodes)
                    node_set = set(g_copy.nodes())                    
            
                    if nx.is_strongly_connected(g_copy): 
                        # Eval this cycle
                        new_SCC = list(node_set)
                        status = check_valid(new_SCC)
                        if status == 'valid':
                            # Found a stable SCC - make sure it isn't a 
                            # duplicate or superset of another SCC
                            if not any([set(new_SCC) >= set(x) for x in valid_SCCs]):
                                # Must also consider reverse case to avoid duplicates
                                valid_SCCs = [x for x in valid_SCCs if not set(x).issuperset(new_SCC)]
                                valid_SCCs.append(new_SCC)

                            # Now check to see if this SCC is a superset of an 
                            # existing valid SCC. First make sure it isn't a 
                            # subset of a different parent
                            if not any([set(new_SCC) <= set(x) for x in parent_SCCs]):
                                # This is a new "big one." We'll track it as a 
                                # parent, but remove any other parent that is
                                # a subset of this one
                                parent_SCCs = [x for x in parent_SCCs if not set(x).issubset(new_SCC)]
                                parent_SCCs.append(new_SCC)  
                        elif status == 'missing parents':
                            # Still some missing parents; check it against the 
                            # base invalid cycles again
                            SCCs_check.append(new_SCC)
    
    # Now repeat with SCCs_check against the base set of invalid cycles
    # But iteratively until we don't modify SCCs_check any more
    new_SCCs_check = []
    while True:
        for c1 in SCCs_check:
            missing_parents1 = check_missing_parents(c1)
            present_parents1 = [x for x in c1 if ' ' not in x]
            for c2 in invalid_cycles:
                # If one is the subset of the other, no chance for improvement 
                # here
                if set(c1) <= set(c2) or set(c2) <= set(c1): continue   
                  
                missing_parents2 = check_missing_parents(c2)
                present_parents2 = [x for x in c2 if ' ' not in x]
                
                status = compare_parents(present_parents1 + present_parents2)  
                
                if status == 'valid':
                    # We get here if there are no node/antinode pairs in this 
                    # combination
                    # Below is non-empty if we have found a new parent
                    overlap1 = set(missing_parents1) & set(present_parents2) 
                    overlap2 = set(missing_parents2) & set(present_parents1) 
                    if len(overlap1) > 0 or len(overlap2) > 0:
                        # We've added a parent, yay -- now need to make sure we 
                        # still have a SCC
                        g_copy = g.copy()
                        excess_nodes = set(g.nodes()) - (set(c1) | set(c2))
                        g_copy.remove_nodes_from(excess_nodes)
                        node_set = set(g_copy.nodes())                    
                
                        if nx.is_strongly_connected(g_copy): 
                            # Eval this cycle
                            new_SCC = list(node_set)
                            status = check_valid(new_SCC)
                            if status == 'valid':
                                # found a new stable SCC - make sure it isn't a
                                # duplicate or superset of another SCC
                                if not any([set(new_SCC) >= set(x) for x in valid_SCCs]):
                                    # Must also consider reverse case to avoid duplicates
                                    valid_SCCs = [x for x in valid_SCCs if not set(x).issuperset(new_SCC)]
                                    valid_SCCs.append(new_SCC)

                                # Now check to see if this SCC is a superset of
                                # an existing valid SCC. First make sure it 
                                # isn't a subset of a different parent
                                if not any([set(new_SCC) <= set(x) for x in parent_SCCs]):
                                    # This is a new "big one." We'll track it 
                                    # as a parent, but remove any other parent 
                                    # that is a subset of this one
                                    parent_SCCs = [x for x in parent_SCCs if not set(x).issubset(new_SCC)]
                                    parent_SCCs.append(new_SCC)  
                            elif status == 'missing parents':
                                # Still some missing parents; check it against 
                                # the base invalid cycles again
                                new_SCCs_check.append(new_SCC)

    
        # If we've added -any- new SCCs that are not valid, we need to continue
        # seeing if other cycles can be combined with it to make it become valid.
        # But if we haven't found any, we can safely stop.
        if len(new_SCCs_check) == 0:
            break
        elif len(new_SCCs_check) > 1E6: 
            # Safety valve
            raise RuntimeError("STOP")
        else:
            SCCs_check = new_SCCs_check[:]
            new_SCCs_check = []
                
    return valid_SCCs, parent_SCCs, invalid_cycles

def reduce_network(g,SCC,remove=True):
    """Identify nodes that are frozen or can be removed.
    
    Given some source nodes that are to be fixed to the  ON or OFF state, the 
    graph is then reduced further:
        Step 1: 'freeze' the state of the source nodes, and update the rules 
        for other nodes accordingly. If they can only have 1 value, freeze them
        as well.
  
        Step 2: remove sink nodes and mediator nodes as in node B in the chain 
        A->B->C, but add the connection A->C (A replaces B in the update rule 
        for C)
    
    Keyword arguments:
        g:      A graph in the standard representation.
        SCC:    A list of nodes & node negations -or- False. If a list, 
                including a node fixes the node ON; including the negation 
                fixes the node OFF. If False, the frozen nodes are those for
                which the node property 'update rules' is True or False.
        remove: True/False; determines if we reduce linear chains and eliminate
                sink nodes.
        
    Returns:
        <'done','not done'>:    If all nodes are frozen or not.
        g_copy:                 A version of the input graph where the 
                                update_rules are updated such that frozen nodes
                                become 1 or 0 as appropriate.
        frozen:                 A dictionary with frozen nodes as keys and 
                                values of <0,1>
        removed:                A dictionary with removed nodes as keys and 
                                either (1) the nodes they depend upon as vals 
                                (if simplifed as a chain) or (2) the string 
                                'removed_as_sink' if it was removed as a sink 
                                node and depends upon all if its inc. edges.
    """
    
    def check_outputs(s):
        """
        Takes a logical string (in the appropriate format: (.. and ..) or ... )
        and considers all combinations of True/False for its inputs. 
        
        Returns 'True' if all combinations yield 'True', 'False' if all
        combinations yield 'False', or otherwise 'mixed'.
        """
        # First acquire a list of the non-operator words
        words = s.replace(' AND ',' and ')
        words = words.replace(' OR ',' or ')
        words = words.replace('NOT ','not ')
        
        words = words.replace('(','')
        words = words.replace(')','')
        words = words.replace(' and ',' ')
        words = words.replace('not ','')
        words = words.replace(' or ',' ')
                                           
        # Unify use of True/False vs. 1/0 in updates to the rules
        words = words.replace('True','1')
        words = words.replace('False','0')

        words = words.split()  

        while '0' in words: words.remove('0')
        while '1' in words: words.remove('1')
        
        # Remove duplicates since ea. node is set to 1 value below
        words = list(set(words))  
        
        # If the below fires, all input nodes are frozen
        if words == []: return str(eval(s))   

        outputs = set()
        # Now we need to iterate through all combinations of inputs
        for i in range(2**len(words)):
            state = bin(i)
            state = state[2:]
            # Prepend with 0s as needed
            state = '0'*(len(words) - len(state)) + state   
            cur_state = s[:]
            for i,c in enumerate(state):
                cur_state = cur_state.split()
                cur_state = [x.replace(words[i],c) if x.strip('() ')
                            == words[i] else x for x in cur_state]
                cur_state = ' '.join(cur_state)

            output = eval(cur_state)
            if output == False or output == 0:
                outputs.add('False')
            elif output == True or output == 1:
                outputs.add('True')
            else:
                raise RuntimeError('unexpected output in eval()')
            
            # if we get a True and a False, we know the output is 'mixed'
            if len(outputs) > 1: return 'mixed'
        
        # If we eval all outcomes and have just one output, return it in 
        # string format.
        return str(outputs.pop())
                
    
    G_copy = g.copy()
    nodes = list(G_copy.nodes())
    frozen = {}
    removed = {}
    
    # STEP 1 (initial)
    if SCC != False:
        # first set the nodes in the SCC to frozen
        for i in SCC:
            if ' ' not in i:            # omit composite nodes
                if '~' not in i:        # node is frozen ON
                   frozen[i] = 1
                else:                   # node is frozen OFF
                    frozen[i[1:]] = 0
                    
    else:
        # we instead use the source network
        for i in g.nodes():
            if g.node[i]['update_rules'] == 'True':
                frozen[i] = 1
            elif g.node[i]['update_rules'] == 'False':
                frozen[i] = 0
                
    curr_frozen = len(frozen)
    
    # Use of updated allows us to check (1) if a newly frozen node changes the 
    # rule for a node or (2) network reduction changes the rule for a node 
    # (as in A->B->C going to A->C)
    updated = set()  
    
    while True:        
        # Step 1 (subsequent)
        # now need to check if any -other- nodes have become frozen
        for i in G_copy.nodes():
            if i in frozen: continue
            for frozen_node in frozen.keys():
                # Here we need to be careful about frozen_node being a 
                # substring of some other node. Need to do below line for 
                # iterative updates (i.e. if more than 1 frozen node in the rule)
                rule_split = G_copy.node[i]['update_rules'].split()   
                if any([frozen_node == x.strip('() ') for x in rule_split]):
                    updated.add(i)
                    new_rule = [x.replace(frozen_node,str(frozen[frozen_node])) 
                                if x.strip('() ') == frozen_node else x
                                for x in rule_split]
                    new_rule = ' '.join(new_rule)
                    G_copy.node[i]['update_rules'] = new_rule
            
            if i in updated:
                # Rule has been updated so we have 1s and 0s in place of the 
                # frozen nodes. Still in form of (A AND B AND ...) OR (...) ...
                # We will check to see if all combinations of logical inputs 
                # eval to same output
                term = G_copy.node[i]['update_rules'][:]

                check_frozen = check_outputs(term)
                if check_frozen == 'mixed':
                    updated.remove(i)

                elif check_frozen == 'False' or check_frozen == '0':
                    # The rule always goes to 0: the node is frozen to False
                    frozen[i] = 0
                else:
                    # The rule always goes to 1: the node is frozen to True
                    frozen[i] = 1
             
        # STEP 2
        remove_node = False
        # Do we bother?
        if remove:
            # If yes, go through every node
            for node in nodes:
                if node in frozen: continue
                if node != 'Closure' and G_copy.out_degree(node) == 0:
                    G_copy.remove_node(node)
                    removed[node] = 'removed_as_sink'
                    remove_node = True
                elif G_copy.in_degree(node) == 1 and G_copy.out_degree(node) == 1:
                    # Make sure we don't have a loop A->B->A
                    pred = list(G_copy.predecessors(node))[0]                              
                    succ = list(G_copy.successors(node))[0]
                    if pred != succ:
                        # Delete B, add A->C, update the rule for C by changing 'B'
                        # to B's rule (either 'A' or 'NOT A')
                        G_copy.add_edge(pred,succ)
                        
                        words = G_copy.node[succ]['update_rules'].split()
                        words = [G_copy.node[node]['update_rules'] if x == node else x for x in words]
                        words = ' '.join(words)
                        G_copy.node[succ]['update_rules'] = words[:]
                        updated.add(succ)
                        if 'not ' in G_copy.node[node]['update_rules']:
                            # Record the presence of a 'not ' reversal as the 
                            # second entry in the tuple
                            removed[node] = (pred,1)
                        else:
                            removed[node] = (pred,0)
                            
                        remove_node = True
                        G_copy.remove_node(node)
        
        # See if we're done
        if remove_node == False and len(frozen) == curr_frozen:
            # No nodes removed and no change in the # of frozen nodes
            break    
        
        # otherwise update the new number of frozen nodes and repeat
        curr_frozen = len(frozen)
        
    # Now we need to see if there are other nodes that remain unfrozen
    if len(list(set(G_copy.nodes()) - set(frozen.keys()))) == 0:
        # Every node is frozen
        return 'done', G_copy, frozen, removed
    else:
        return 'not done', G_copy, frozen, removed


def find_oscillations(g,stable_SCCs):
    """Find the oscillating subgraphs in an expanded graph.
    
    Identify the SCCs with the properties required for oscillations:
        1. Contain the complementary node of every included normal node and 
           vice versa.
        2. Contain input nodes of all included composite nodes.
        3. Does not contain a stable motif consisting only of nodes and
           complementary nodes.
        
    Keyword arguments:
        g:              an expanded graph.
        stable_SCCs:    a list of SCCs.
        
    Returns:
        osc_SCCs:       a list of lists; each interior list contains the 
                        oscillating nodes in an oscillating SCC.
    """
    
    def check_osc_valid(s):
        """Check if a list of nodes forms a valid oscillating SCC."""
        invalid = False
        
        # First check all composite nodes
        for i in range(len(s)):
            n = s[i]
            if ' ' in n:
                # a composite node
                parents = n.split(' ')
                for p in parents:
                    if p not in s:
                        invalid = True
                        break
        
        # If necessary, now look at non-composite nodes
        if not invalid:
            c = s[:]
            while len(c) > 0:
                n = c.pop()
                # skip composite nodes
                if ' ' in n: continue
                
                if '~' == n[0]: 
                    # We get here if we have a complementary node
                    if n[1:] in c:
                        #  If the corresponding node exists, remove them both
                        c.remove(n[1:])
                    else:
                        # Otherwise this is not a valid SCC
                        invalid = True    
                        break
                # Below check fires if we have a regular node
                elif '~'+n in c:
                        # If the corresponding complementary node exists, 
                        # remove them both
                        c.remove('~'+n)   
                # Below only fires if we have a regular node and -not- the 
                # corresponding complementary node
                else:    
                    invalid = True    
                    break
        
        return invalid
    
    osc_SCCs = []
    
    # Viable stable SCCs have no composite nodes -- composite nodes have a space
    # in their name. So use this as a filter point
    viable_stable_SCCs = []
    for i in stable_SCCs:
        t= ''.join(i)
        if ' ' not in t: viable_stable_SCCs.append(i)
    
    # Now look at the SCCs of the expanded graph (here we want the largest
    # SCCs, so we don't need to build up cycles)
    SCCs = [x for x in nx.strongly_connected_component_subgraphs(g)]
    while len(SCCs) > 0:
        SCC = SCCs.pop()
        # First check if this SCC satisfies conditions 1 & 2
        SCC_nodes = list(SCC.nodes())
#        if len(SCC_nodes) == 1:
#            print('CHECK')
        invalid = check_osc_valid(SCC_nodes)
        # Now check condition 3
        if not invalid:
            SCC_node_set = set(SCC.nodes())
            for v in viable_stable_SCCs:
                if set(v) <= SCC_node_set:
                    # does not meet condition 3 b/c of this fixed-point SCC
                    invalid = True
                    # remove these nodes and look for smaller SCCs
                    SCC_copy = SCC.copy()
                    SCC_copy.remove_nodes_from(v)
                    for x in nx.strongly_connected_component_subgraphs(SCC_copy):
                        SCCs.append(x)
            if not invalid:
                osc_SCCs.append(SCC)

    return osc_SCCs
      

def parent_SCC_search(g,source,cutoff = 30):
    """ Generator function to identify valid circuits on a parent SCC.
    
    This function is modified from all_simple_paths() in the NetworkX module.
    Here we are looking for a circuit that does not include node/antinode pairs
    and includes all of the parent nodes for any included composite node. We
    additionally avoid walking across the same cycle 2x to avoid redundancy.
    
    Keyword arguments:
        g:          A NetworkX DiGraph of the parent SCC.
        cutoff:     The maximum path length to consider. Default is 30.
        
    Yields:
        visited:    A valid circuit on the parent graph.
    """   
    
    
    if cutoff < 1:
        return
    visited = [source]
    stack = [iter(g[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.pop()
        elif len(visited) < cutoff:
            # Make sure we haven't walked across the same cycle 2x 
            # (e.g. [C ... C ...] + new C)
            if visited.count(child) >= 2:
                # Last index of child
                end_entry_index = len(visited) - 1 - visited[::-1].index(child) 
                # The most recent child-child cycle
                potential_cycle= visited[end_entry_index:]  
                potential_cycle_length = len(potential_cycle)
                start_index = visited.index(child)
                repeat = False
                while child in visited[start_index:]:
                    cur_path = visited[start_index:
                                       start_index+potential_cycle_length]
                    if cur_path == potential_cycle:
                        # we do indeed have a repeat; stop
                        repeat = True
                        break
                    start_index += 1         
                if repeat:
                    stack.pop()
                    visited.pop()
                    continue         
            
            if child == source:
                status = check_valid(visited+[source])
                if status == 'valid':
                    # We have a valid circuit
                    yield visited
        
            else:
                # Continue this walk
                visited.append(child)
                stack.append(iter(g[child]))
                    
        else:  
            # beyond the cutoff, so stop
            stack.pop()
            visited.pop()
 
    
    
def remove_from_composite(g,composite,strip_node,target):
    """Update an expanded graph by removing a node from a composite node.
    
    The edge composite -> target is removed and replaced by a new edge from
    the composite node modified such that the strip_node is removed.
    
    If the original composite node has no outgoing edges, it is removed.
    
    If the new composite node is a subset of another composite node with an
    edge to target, then the larger composite node is also removed.
    Example:
        AB -> T
        AC -> T

        If we remove B from AB, we have
        A -> T
        AC -> T
        Which reads T* = (A) or (A and C) = A
        
    Keyword arguments:
        g:          an expanded graph 
        composite:  a composite node in g
        strip_node: a node to be removed from the composite 
        target:     a node targeted by the composite node
        
    Returns:
        g_out:      the updated expanded graph
    """
    
    g_out = g.copy()
    composite_nodes = composite.split()
    composite_nodes.remove(strip_node)
    new_composite_str = ''
    for x in composite_nodes: new_composite_str += x + ' '
    new_composite_str = new_composite_str.rstrip()
    # Remove the old edge
    g_out.remove_edge(composite,target)
    # If this predecessor composite has no other outgoing edges
    # we can remove it
    if g_out.out_degree(composite) == 0: g_out.remove_node(composite)
    
    # Look at other composite node predecessors of target
    pred_comp_nodes = [x for x in g_out.predecessors(target) if ' ' in x]
    for pcn in pred_comp_nodes:
        pcn_split = pcn.split()
        if all([x in pcn_split for x in composite_nodes]):
            # If the new composite node is a subset of this composite node,
            # this composite node is logically superfluous and can be removed
            g_out.remove_edge(pcn,target)
            if g_out.out_degree(pcn) == 0: g_out.remove_node(pcn)
    
    # Add the new edge
    g_out.add_edge(new_composite_str,target)
    
    return g_out



def iterate_SCC(g, g_expanded, SCCs, recursion = 0,frozen = {}, removed = {}, 
                dump_to_screen=True, remove = True, output=False):
    """Consider freezing nodes in a SCC and allow the effects to propagate.
    
    This function considers the effect of freezing the state of a set of 
    nodes (starting from a stable SCC) on the expanded graph. For each SCC
    provided, the function "freezes" the state of the constituent nodes and
    identifies any successor nodes that become frozen as a result.
    
    The function then calls itself recursively with the new set of updated
    nodes until every node is frozen or no additional nodes can be frozen.
    
    The output is dumped directly to console unless output == True, in which
    case the global variable fixed_attractors stores the unique outcomes as
    lists of N-length strings (for N nodes): '0' and '1' to store fixed node 
    states and '?' to indicate non-fixed node states.
    
    The approach involves recursive self-calls; we embed the main algorithm
    in an interior function to preserve the namespace independence of the
    output objects. (This is a bit clunky and could be streamlined!)
    
    Keyword arguments:
        g:                  The NetworkX graph.
        g_expanded:         The expanded representation of g.
        SCCs:               A list of strongly-connected components.
        recursion:          Tracks the number of iterative function calls.
        frozen:             A dictionary of frozen nodes: {node:state}
        remove:             True/False; do we simplify chains and remove sink
                            nodes?
        removed:            A dictionary of removed nodes: {node:details}
        dump_to_screen:     True or False; controls amount of information
                            output to console. Default is True.
        output:             True or False; if True, fixed_point attractors are
                            written to the global list fixed_attractors.
                            Default is False.
    
    Returns:
        fixed_attractors:   Text
        frozen_states:      A dictionary of frozen nodes in {node:state} form
    """
    
    # Interior function definition --------------------------------------------
    def interior(g, g_expanded, SCCs, recursion = 0,frozen = {}, removed = {}, 
                 dump_to_screen=True, remove=True, output=False):
        """Performs the main loop of iterate_SCC()"""
        # Now reduce the network for each SCC
        for SCC in SCCs:
            if dump_to_screen:
                print('recursion = {0} ========================'.format(recursion))
                print('Beginning SCC = ')
                print(SCC)
            # For each SCC we start over the process of adding to the frozen nodes
            frozen_cur = frozen.copy()  
            # If we have frozen nodes, include them with the SCC
            removed_cur = removed.copy()  # Likewise for removed nodes
            SCC_append = SCC.copy()
            for i,j in frozen_cur.items():
                if i in SCC_append: continue
                if j == 0: SCC_append.append('~'+i)
                else: SCC_append.append(i)
            status, g_updated, frozen_new, removed_new = reduce_network(g,SCC_append,remove)
            frozen_cur.update(frozen_new)
            removed_cur.update(removed_new)
                 
            # If some nodes aren't frozen, we need to pull out the unfrozen nodes
            # and repeat.
            if status != 'done':
                g_expanded_reduced = form_expanded_network_reduced(g_expanded,frozen_cur, removed_cur)   
                
                SCC_list_reduced, parent_SCCs_reduced, partial_cycles_reduced \
                = find_fixed_SCCs_iterative(g_expanded_reduced)
                
                if len(SCC_list_reduced) == 0:
                    # Some nodes remain unfrozen, but we have no more viable SCCs.
                    # Conclusion: some nodes oscillate, and we can stop.
                    
                    final_frozen_counter = 0
                    final_frozen_dict = {}
                    
                    if dump_to_screen: 
                        print('Partial Fixed Point!')
                        print('Frozen nodes, recursion = {0}'.format(recursion))
                        
                    for key,val in frozen_cur.items():                    
                        final_frozen_dict[key] = val
                        final_frozen_counter += 1
                        
                    for key,val in removed_cur.items():
                        val_node = val[0]  # just the node association
    
                        # We can have a case where frozen_cur[node1] = node2
                        # where node2 is also a key in frozen_cur, etc., eventually
                        # pointing to a node in frozen_cur or not.
                        # Count how many A* = not B instances we acquire in the chain
                        neg_count = val[1]  
                        while val_node not in frozen_cur.keys() and val_node in removed_cur.keys():
                            val_node = removed_cur[val_node][0]
                            neg_count += removed_cur[val_node][1]
                        # val now points to a node that is (1) in frozen (meaning 
                        # key is frozen too) and/or (2) not in removed_cur (meaning
                        # key is key is not frozen if condition (1) is not met)
                        if val_node in frozen_cur: 
                            # do we want the value or the opposite of the value?
                            # depends on whether we find an odd number of negatives
                            # above
                            if neg_count % 2 == 1:
                                state = (frozen_cur[val_node] + 1)%2
                            else:
                                state = frozen_cur[val_node]
                            if key not in final_frozen_dict: 
                                final_frozen_counter += 1
                            final_frozen_dict[key] = state
                            
                    if dump_to_screen: 
                        print('{0} nodes frozen in total.'.format(final_frozen_counter))                  
                    
                    frozen_states.append(final_frozen_dict)
                
                else:
                    # We have more SCCs, and a new recursion into this function
                    # is necessary
                    interior(g_updated,g_expanded_reduced,SCC_list_reduced, 
                                recursion = recursion + 1,frozen = frozen_cur, 
                                removed = removed_cur, dump_to_screen = dump_to_screen, 
                                output=output)
                    
            # If status == 'done', every node is frozen and we're at a fixed point.    
            else:      
                if dump_to_screen:
                    print('Fixed Point!')
                    print('Frozen nodes, recursion = {0}'.format(recursion))
                    for key,val in frozen_cur.items():
                        print('{0:15}:\t{1}'.format(key,val))
                    print('All nodes frozen; terminating.')
                if output:
                    out = ['?']*len(frozen_cur)
                    nodes = sorted(frozen_cur.keys())
                    for i in range(len(nodes)):
                        out[i] = str(frozen_cur[nodes[i]])
                    # store a list of node states, '0','1', or '?' if oscillating
                    # but only store unique copies
                    if out not in fixed_attractors and output:
                        fixed_attractors.append(out)
                
                frozen_states.append(frozen_cur)
    # End interior function definition ----------------------------------------  
    
    fixed_attractors, frozen_states = [],[]  
    # Just pass through all arguments
    interior(g=g, g_expanded=g_expanded, SCCs=SCCs, recursion = recursion, 
             frozen=frozen, removed=removed, dump_to_screen=dump_to_screen, 
             remove=remove, output=output)
    return fixed_attractors, frozen_states

def MMCS(edge_sets,edge_collection):
    """ The MMCS algorithm for the minimum hitting set problem.
    
    For a collection of unordered sets C, a hitting set is a set of items that
    includes at least one item from each set in C. A hitting set is minimal if 
    no proper subset is itself a hitting set. Our sets are collections of edges
    from different circuits on the expanded graph.
    
    We here implement the MMCS algorithm from Murakami and Uno (Discrete
    Applied Mathematics, 2014). The algorithm involves recursive self-calls,
    and we rely on the following global variables:
        crit:   a dictionary of edge:set() key:val pairs, across all edges
        uncov:  a list of the edge sets
        CAND:   a set of all edges
        MHS:    an empty list; will store identified minimum hitting sets
    
    Because the function involves recursive self-calls, we embed the algorithm
    in an internal loop (to avoid namespace confusion). This could definitely
    be streamlined!
    
    The keyword argument S is used internally and should not be changed.
    """
    global crit, uncov, CAND, MHS
    
    crit = {u:set() for u in edge_collection}
    uncov = edge_sets.copy()
    CAND = edge_collection.copy()
    MHS = []

    # Interior function definition --------------------------------------------
    def interior(S=set()):
        global crit, uncov, CAND, MHS
        
        # step 1
        if uncov == []:
            MHS.append(S)
            return S
        # step 2
        F = uncov[randint(0,len(uncov)-1)]
        # step 3
        C = CAND & F
        CAND = CAND - C
        # step 4
        for v in C:
            # Step 5 is to call update_crit_uncov, which we just embed here
            # Next two lines are in prep for step 7 (recovering changes done in
            # step 5)
            uncov_copy = uncov[:]
            crit_copy = crit.copy()
            # Step 5-1
            for F in edge_sets:
                if v not in F: continue
                # Step 5-2
                for u in S:
                    if F <= crit[u]:
                        crit[u] = crit[u] - F
                        break
                # Step 5-3
                if F in uncov:
                    uncov.remove(F)
                    crit[v] = crit[v] | F
            # Step 6
            overall_recursion = True
            for f in S:
                inner_recursion = False
                for f_check in edge_sets:    
                    if (S|set([v])) & f_check == set([f]):
                        # get here if crit(f,S|v) != {}
                        inner_recursion = True
                        break
                if inner_recursion == False:
                    # For this f, crit(f,S|v) = {}. Thus the overall check across
                    # all possible f in S is false and we do not call MMCS again
                    overall_recursion = False
                    break
            if overall_recursion:
                interior(S = S | set([v]))
            CAND = CAND | set([v])
            # Step 7
            uncov = uncov_copy[:]
            crit = crit_copy.copy()
    # End interior function definition ----------------------------------------
    
    interior()
    
    return MHS
    
    