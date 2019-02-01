'''
State transition network (STN) analysis of the ABA network. See ABA.py for
the expanded graph analysis.

Note: This script is designed to be run with either "Option 1" or "Option 2"
active; the other should be commented out. Option 1 generates the output file
"ABA dynamics analysis.data"; Option 2 loads and analyzes it to generate the 
bargraph in the form of Figure 9.

The version of ABA dynamics analysis.data used in the report is available on
GitHub.

Author: Colin Campbell
Date: November 2018
Language: Python 3.6
'''

import networkx as nx
import pickle
import glob
import time
import numpy as np
from numpy.random import randint, shuffle
from numpy import std
import matplotlib.pyplot as plt

 
# FUNCTION DEFINITIONS ========================================================

def form_network(rules):
    """
    Takes as input a list of rules in the format of sample_network.txt.

    Outputs a networkx DiGraph with node properties:
        'update_nodes': a list of regulating nodes
        'update_rules': a dictionary with binary strings as keys, corresponding
                        to the possible states of the update_nodes, and integers
                        as values, corresponding to the state of the node given
                        that input.

    Note that nodes are identified by their position in a sorted list of
    user-provided node names.
    """
    
    def clean_states_inf(x):
        """Clean binary representation of node input states."""
        # Strip leading 0b
        out=x[2:]                                                              
        # Append leading 0's as needed
        return '0'*(len(g.node[n]['update_nodes'])-len(out))+out               

    # Remove comments, blank lines, and rules for node negations
    stream = [x.rstrip('\n') for x in rules if x != '\n' and x[0]!='#' and x[0]!='~']  
    stream = [x.lower() for x in stream]
    # Generate a sorted list of node names
    nodes = sorted([x.split(' ',1)[0][:-1] for x in stream])                    

    g = nx.DiGraph()
    # At creation, no node is flagged for knockout or overexpression
    g.graph['knockout'] = None                                                  
    g.graph['express'] = None

    for n in range(len(stream)):
        rule = stream[n].split('* = ')[1]
        # Force decap of logical operators so as to work with eval()
        rule = rule.replace(' AND ',' and ')                                    
        rule = rule.replace(' OR ',' or ')
        rule = rule.replace('NOT ','not ')
        # For always ON or always OFF nodes
        if stream[n].find('true') >= 0 or stream[n].find('false') >= 0:         
            # We refer to nodes by their location in a sorted list of the 
            # user-provided node names
            g.add_node(n)                                                       
            g.node[n]['update_nodes'] = []
            g.node[n]['update_rules'] = {'':str(int(eval(rule.title())))}
            continue

        # Strip down to just a list of influencing nodes
        inf = rule.split(' ')                                                   
        inf = [x.lstrip('(') for x in inf]
        inf = [x.rstrip(')') for x in inf]
        
        # If we've fed in modified rules, we might have '0's or '1's that need 
        # to be ignored
        while '0' in inf: inf.remove('0')
        while '1' in inf: inf.remove('1')

        # The sort ensures that when we do text replacement (<node string> ->
        # 'True' or 'False') below in this fn, we avoid problems where node 1 
        # is a substring of node 2 (e.g. NODE1_phosph and NODE1)
        inf = sorted([x for x in inf if x not in ['','and','or','not']],key=len,reverse=True)

        # Add edges from all influencing nodes to target node
        for i in inf: g.add_edge(nodes.index(i),n)                              
        g.node[n]['update_nodes'] = list(set([nodes.index(i) for i in inf]))
        g.node[n]['update_rules'] = {}

        bool_states = map(bin,range(2**len(g.node[n]['update_nodes'])))
        bool_states = map(clean_states_inf,bool_states)
        for j in bool_states:
            rule_mod = rule[:]
            for k in range(len(j)):
                rule_words = rule_mod.split(' ')
                node_name = nodes[g.node[n]['update_nodes'][k]].lower()
                if j[k]=='0':
                    rule_mod = [x.replace(node_name,'False') if x.strip('() ') 
                                == node_name else x for x in rule_words]
                    rule_mod = ' '.join(rule_mod)
                else: 
                    rule_mod = [x.replace(node_name,'True') if x.strip('() ') 
                                == node_name else x for x in rule_words]
                    rule_mod = ' '.join(rule_mod)

            # Store outcome for every possible input
            g.node[n]['update_rules'][j] = int(eval(rule_mod))                                  

    return g,nodes


def clean_states(x):
    """Clean binary representation of node input states."""
    # Strip leading 0b
    out=x[2:]                                                                   
    # Append leading 0's as needed
    return '0'*(len(nodes)-len(out))+out                                          


def update_state_asynch(x,state):
    """Acquire string of current states of node x's input nodes."""
    input_state = ''
    for j in G.node[x]['update_nodes']: input_state += str(state[j])
    return str(G.node[x]['update_rules'][input_state])


def probe_asynch(G):
    """Probe the STN for asynchronous updates.

    Samples M (set internally) random states and, for each, updates the state P 
    (set internally) times. Stores the state of each node for the final 10 
    states. Upon completion of the M trials, determines the average state for 
    each node.
    
    We follow Albert et al. PLOS BIO 2017 in using a random order update scheme
    (all nodes are updated in random order to ID the next state; repeat). 
    
    Keyword Arguments:
        G:      The NetworkX graph.
    
    Returns:
        master: A dictionary with node # as keys and a tuple of the average 
                values and stard deviations (avg, std) as vals.
    """
    
    master = {x: [] for x in G.nodes()}
    
    N = nx.number_of_nodes(G)
    M = 1E4
    P = 40
    for i in range(int(M)):      
        # Generate a random sequence of 0,1 vals
        start_state = list(randint(0,2,N))
        start_state = ''.join([str(x) for x in start_state])
        
        trajectory = [start_state]
        # For this start state, we update P times
        cur_state = start_state[:]               # cur_state is where we are
        for update in range(P):
            # "state" is the state we're moving -to-; begin with copy of where 
            # we -are-
            state = cur_state[:]   

            # Now determine the order of updates for this round
            order = list(range(N))
            shuffle(order)
            
            # Now update each of the nodes to determine the next state                 
            for node in order:
                next_node_state = update_state_asynch(node,state)
                state=list(state)
                state[node] = next_node_state
                state=''.join(state)
            
            # Store our new state
            trajectory.append(state)
            # The current state becomes what we just moved to
            cur_state = state[:]                    
       
        # Look at the last 10 states and store the state
        for j in range(N):
            for x in range(-10,0): master[j].append(int(trajectory[x][j]))
    
    # Now we have a series of 10*M states for each of the N nodes. Compute
    # the average and standard deviation for each of the nodes
    for i in range(N):
        avg = float(sum(master[i]))/len(master[i])
        stdev = std(master[i])
        master[i] = tuple([avg,stdev])
    
    return master

# END FUNCTION DEFINITIONS ====================================================

# OPTION 1 - iterate through networks; save to file ------------
# Caution: Running this block will overwrite the default version of the
# .data file.
#dynamics_dict = {}
#TIME0 = time.time()
#for fname in glob.glob(r'../Data Files/ABA dynamics/*.txt'):
#    f = open(fname,'rt')
#    lines = f.readlines()
#    f.close()
#    
#    # call out the file name but cut the identical beginnings and .txt ending
#    file_id = fname[fname.rfind('_')+1:-4]
#    
#    # give update on timing
#    TIME1 = time.time()
#    print('Beginning file {0} after {1:.2f} minutes elapsed.'
#          .format(file_id,(TIME1-TIME0)/60.0))
#    
#    G, nodes = form_network(lines)
#    dynamics_dict[file_id]  = probe_asynch(G)
#    
#    # Save the results so far (updates as we go)
#    f = open('ABA dynamics analysis.data','wb')
#    pickle.dump(dynamics_dict,f)
#    f.close()


# OPTION 2 - analyze the output generated above ---------------
f = open(r'ABA dynamics analysis.data','rb')
dynamics_dict = pickle.load(f)
f.close()

# Need to acquire the list of sorted nodes. Use the network rules in the 
# original/standard format
f = open(r'../Network Rules/ABA2017_original.txt','rt')
lines = f.readlines()
f.close()

# Form the network
G,nodes = form_network(lines)

# Obtain a copy of the expanded reduced network (refer to nodes therein below)
f = open(r'ABA_expanded_reduced_network.data','rb')
G_expanded_reduced = pickle.load(f)
f.close()

G_er_nodes = [x.lower() for x in G_expanded_reduced.nodes()]
# Just actual node names; not composite or complement
G_er_nodes = [x for x in G_er_nodes if ' ' not in x and '~' not in x]

# Where is Closure?
Closure_id = nodes.index('closure')

# Print summary to screen and store information to make a barplot --
# First determine the order we'd like to plot (matches Table II in main text)
key_order = ['MAPK912', 'CPK213', 'microtubule', 'cycle1', 'cycle2', 'cycle3', 
             'cycle4', 'cycle5', 'all1', 'all2', 'all3', 'all4', 'all5']
closure = []
nodes_fixed = []

# Now print information to screen and store in above containers
for key in key_order:
    val = dynamics_dict[key]
    osc_count = 0
    debug_count = 0
    for node, result in val.items():
        if nodes[node] in G_er_nodes and result[0] != 1.0 and result[0] != 0.0: 
            osc_count += 1
    print('For network {0}:\n\tmean closure state is = {1:.3f}\n\t# nodes fixed ON or OFF = {2} of {3} ({4:.2f}%)'
          .format(key, val[Closure_id][0], osc_count, 
                  len(G_er_nodes),100.0*float(osc_count)/len(G_er_nodes)))
    closure.append(100.0*val[Closure_id][0])
    nodes_fixed.append(100.0*float(osc_count)/len(G_er_nodes))
    
# And create a bar plot
fig, ax = plt.subplots()

index = np.arange(len(key_order))
bar_width = 0.35

opacity = 1.0
error_config = {'ecolor': '0.3'}

rects1 = ax.bar(index, closure, bar_width,
                alpha=opacity, color='#FF6666', edgecolor='k',
                label='<Closure>')

rects2 = ax.bar(index + bar_width, nodes_fixed, bar_width,
                alpha=opacity, color='#3399FF', edgecolor='k', hatch='\\',
                label='oscillating nodes')

ax.set_xlabel('Perturbations')
ax.set_ylabel('%')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(('1','2','3','4','5','6','7','8','1-3,4','1-3,5','1-3,6','1-3,7','1-3,8'),rotation=45)
ax.legend()

fig.tight_layout()
plt.show()


# In addition to the above, we can analyze the expanded reduced network to
# generate a .txt of the rules (provided as SI)
#
#f = open(r'ABA_expanded_reduced_network.data','rb')
#G_expanded_reduced = pickle.load(f)
#f.close()
#
#nodes = sorted(G_expanded_reduced.nodes())
#nodes = [x for x in nodes if ' ' not in x and x[0] != '~']
#
#lines = []
#
#for cur_node in nodes:
#    # Look at both nodes and node negations
#    for node in (cur_node,'~'+cur_node):
#        preds = list(G_expanded_reduced.predecessors(node))
#        # Make composite nodes explicitly use 'and' instead of a space
#        preds = [x.replace(' ',' and ') for x in preds]
#        # Make negations in form '~' instead in form 'not ' for each predecessor
#        preds = [x.replace('~','not ') for x in preds]
#        # Wrap terms separated by 'or' in parentheses
#        preds = ['('+x+')' for x in preds]
#        # Now join together the predecessors as a string of 'or' statements
#        rule = ' or '.join(preds)
#        
#        # Store this rule in the appropriate format
#        complete_rule = node+'* = '+rule
#        lines.append(complete_rule+'\n')
#
## And write to file (Final version has header comment manually inserted)
#f = open('ABA2017_simplified.txt','wt')
#f.writelines(lines)
#f.close()

    
        