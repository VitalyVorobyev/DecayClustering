# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 23:12:41 2016

@author: Vitaly Vorobyev
"""
import graphviz     as gv
import node_styling as nsty

stable_particles = (11,13,22,211,321,14,12,16,130,2112,2212)

def add_nodes(graph, nodes):
    """ Adds nodes to gviz graph
    Args:
        graph: gviz graph
        nodes: tuple of nodes
    Returns:
        Modified gviz graph
    """
    [graph.node(n[0], **n[1]) for n in nodes]
    return graph

def add_edges(graph, edges):
    """ Adds edges to gviz graph
    Args:
        graph: gviz graph
        edges: tuple of edges
    Returns:
        Modified gviz graph
    """
    [graph.edge(e[0],e[1]) for e in edges]
    return graph

def is_stable_particle(idhep):
    return (abs(idhep) in stable_particles)
    
def make_node_and_edges(idx,idhep,daF,daL):
    """ Retrieves node and edges for a line of EvtGen table
    Args:
        idx: index of the EvtGen table line
        idhep: MC particle code
        daF: index of the first descendant
        daL: index of the last  descendant
    Returs:
        tuple of node and tuple of edges
    """
    sty = nsty.node_style(idhep)
    node = (str(idx),sty)
    if (abs(idhep) not in stable_particles) and (daF > 0):
        edges = tuple([(node[0],str(x)) for x in range(daF,daL+1)])
    else:
        edges = ()
    return (node, edges)

def decay_tree_gv(idhep,daF,daL):
    """ Retrieves gviz graph using the EvtGen table
    Args:
        idhep: array of MC particle codes
        daF: array of the first descendant indices
        daL: array of the last  descendant indices
    Returns:
        gviz graph
    """
    decay_tree = gv.Digraph(comment='Decay Tree')
    nodes = []
    edges = []
    for i in range(0,len(idhep)):
        if idhep[i] == 911 :
            break
        if daF[i] == -1:
            continue
        n,es = make_node_and_edges(i+1,idhep[i],daF[i],daL[i])
        nodes.append(n)
        [edges.append(e) for e in es]
    decay_tree = add_edges(add_nodes(decay_tree,nodes),edges)
    return decay_tree