import networkx as nx
from networkx.algorithms.approximation import steiner_tree, min_weighted_dominating_set
from parse import *
from utils import is_valid_network, average_pairwise_distance_fast
from heapq import heappop, heappush
import sys
import os
import queue
import csv
import datetime
import random
import math

# Considers each leaf vertex in the tree g and greedily prune the tree
# if removing the vertex will reduce the routing costs.
# Setting r to True randomizes the order of checking leaf vertices.
def contractGreedy(G, g, r=False):
    d = average_pairwise_distance_fast(g)
    q = queue.SimpleQueue()
    nodes = list(g.nodes)
    if r:
        random.shuffle(nodes)
    for n in nodes:
        if g.degree[n] == 1:
            q.put(n)
    
    while not q.empty() and len(g) > 1:
        n = q.get()
        if n in g:
            rm = g.copy()
            rm.remove_node(n)
            if nx.is_dominating_set(G, rm.nodes):
                d_rm = average_pairwise_distance_fast(rm)
                if d_rm <= d:
                    for nb in g[n]:
                        if rm.degree[nb] == 1:
                            q.put(nb)
                    g = rm
                    d = d_rm
    return g

# Prune all leaves from the tree that aren't necessary for the tree to be dominating.
def contractAll(G, g):
    q = queue.SimpleQueue()
    for n in g:
        if g.degree[n] == 1:
            q.put(n)
    
    while not q.empty() and len(g) > 1:
        n = q.get()
        if n in g:
            rm = g.copy()
            rm.remove_node(n)
            if nx.is_dominating_set(G, rm.nodes):
                for nb in g[n]:
                    if rm.degree[nb] == 1:
                        q.put(nb)
                g = rm
    return g

# When costy=True, combines the two algorithms and randomizes the contraction order for 
# a few times. Otherwise, use the first contraction algorithm.
def contract(G, g, costy=False):
    if costy:
        ts = [contractAll(G, g), contractGreedy(G, g)]
        for _ in range(9):
            ts.append(contractGreedy(G, g, True))
        return min(ts, key=average_pairwise_distance_fast)
    else:
        return contractGreedy(G, g)

# For all vertices neighboring to the tree g, add them to the tree if doing so would
# reduce the average routing costs.
def expand(G, g):
    d = average_pairwise_distance_fast(g)
    frontier = []
    for u in g:
        for v, d2 in G.adj[u].items():
            if v not in g:
                heappush(frontier, (d2['weight'], v, u))
    while len(frontier) > 0:
        w, u, prev = heappop(frontier)
        if u not in g:
            ng = g.copy()
            ng.add_node(u)
            ng.add_edge(prev, u, weight = w)
            d_ng = average_pairwise_distance_fast(ng)
            if d_ng <= d:
                g = ng
                d = d_ng
                for v, d2 in G.adj[u].items():
                    if v not in g:
                        heappush(frontier, (d2['weight'], v, u))
    return g

# Given a tree g, consider each edge in g to check if replacing the edge with another
# in the graph G would reduce the routing costs.
def localSearch(G, g):
    d = average_pairwise_distance_fast(g)
    flag = True
    while flag:
        flag = False
        t = g.copy()
        for e in g.edges(data=True):
            t.remove_edge(e[0], e[1])
            A, B = list(nx.connected_components(t))
            for u in A:
                for v, d2 in G.adj[u].items():
                    if v in B:
                        t2 = t.copy()
                        t2.add_edge(u, v, weight=d2['weight'])
                        dt = average_pairwise_distance_fast(t2)
                        if dt < d:
                            d = dt
                            e = (u, v, {'weight': d2['weight']})
                            flag = True
            t.add_edges_from([e])
        g = t
    return g

# Given a graph G and a vertex r, compute the (partial) MST starting from r with Prim's.
# Returns when the partial MST is a dominating tree. 
# Change h for different heuristics (not necessarily yieding MSTs).
def msdt(G, root, h=0):
    g = nx.Graph()
    g.add_node(root)
    heap = []
    heuristics = []
    heuristics.append(lambda w, v: w)
    heuristics.append(lambda w, v: w / (sum([nb in g for nb in G[v]]) + 1))
    for v, d in G.adj[root].items():
        heappush(heap, (heuristics[h](d['weight'], v), d['weight'], v, root))
    while len(heap) > 0 and not nx.is_dominating_set(G, g):
        _, tw, u, prev = heappop(heap)
        if u not in g:
            g.add_node(u)
            g.add_edge(prev, u, weight=tw)
            for v, d in G.adj[u].items():
                if v not in g:
                    heappush(heap, (heuristics[h](d['weight'], v), d['weight'], v, u))
    return g

# Solves the problem by trying MSTs rooted at each vertex, and improves
# using other local optimization algorithms above.
def solveMST(G):
    """
    Args:
        G: networkx.Graph

    Returns:
        T: networkx.Graph
    """
    g = nx.minimum_spanning_tree(G)
    d = average_pairwise_distance_fast(g)
    for v in G:
        t = msdt(G, v)
        if len(G) <= 50:
            t = localSearch(G, t)
        t = contract(G, t, True)
        t = expand(G, t)
        dt = average_pairwise_distance_fast(t)
        if dt == 0:
            return t
        if dt < d:
            g = t
            d = dt 
    return g if len(G) <= 50 else localSearch(G, g)

# Given a graph G and a vertex r, compute the SPT starting from r.
# Yields when the partial SPT is a dominating tree. 
# For small and medium graphs, also yields the complete SPT. 
# Change h for different heuristics (not necessarily yieding SPTs).
def spdt(G, root, h=0):
    flag = True
    g = nx.Graph()
    g.add_node(root)
    heap = []
    heuristics = []
    heuristics.append(lambda d, w, v: d + w)
    heuristics.append(lambda d, w, v: d + w / G.degree[v])
    heuristics.append(lambda d, w, v: d + w / (sum([nb in g for nb in G[v]]) + 1))
    for v, d in G.adj[root].items():
        heappush(heap, (heuristics[h](0, d['weight'], v), d['weight'], v, root))
    while len(heap) > 0:
        if flag and nx.is_dominating_set(G, g):
            if len(G) > 50: break
            flag = False
            yield g
        dist, w, u, prev = heappop(heap)
        if u not in g:
            g.add_node(u)
            g.add_edge(prev, u, weight=w)
            for v, d in G.adj[u].items():
                if v not in g:
                    heappush(heap, (heuristics[h](dist, d['weight'], v), d['weight'], v, u))
    yield g

# Solves the problem by trying SPTs rooted at each vertex, and improves
# using other local optimization algorithms above.
def solveSPDT(G):
    g = None
    d = float('inf')
    for v in G:
        for t in spdt(G, v):
            if len(G) <= 50:
                t = localSearch(G, t)
            t = contract(G, t, True)
            t = expand(G, t)
            dt = average_pairwise_distance_fast(t)
            if dt == 0:
                return t
            if dt < d:
                g = t
                d = dt
    return g if len(G) <= 50 else localSearch(G, g)

# Attempts to solve the problem by connecting vertices in the approaximate
# minimum dominating set with a approaximate steiner tree. Not very effective.
def solveDominating(G, dsalg = nx.dominating_set):
    for u in G:
        if len(G.adj[u]) == len(G) - 1:
            g = nx.Graph()
            g.add_node(u)
            return g
    ds = dsalg(G)
    if len(ds) == 1:
        g = nx.Graph()
        g.add_nodes_from(ds)
        return g
    g = steiner_tree(G, ds)
    return expand(G, g)
solveDominatingAlt = lambda G: solveDominating(G, min_weighted_dominating_set)

# Attempts to solve the problem by connecting vertices in the approaximate
# minimum dominating set with a approaximate steiner tree. The dominating set 
# is certered at the vertex with minimum conbined distances to all other vertices.
# Not very effective.
def solveDijkstra(G):
    center = None
    dist = float('inf')
    for u, d in nx.shortest_path_length(G):
        s = sum(d.values())
        if s < dist:
            center = u
            dist = s
    ds = nx.dominating_set(G, center)
    if len(ds) == 1:
        g = nx.Graph()
        g.add_nodes_from(ds)
        return g
    g = steiner_tree(G, ds)
    return expand(G, g)

# Aggregates the results from each algorithm above and selects the best one.
def solve(G):
    algs = [solveMST, solveDominating, solveDominatingAlt, solveDijkstra, solveSPDT]
    trees = [alg(G) for alg in algs]
    dists = [average_pairwise_distance_fast(t) for t in trees]
    return trees[dists.index(min(dists))], dists

# Usage: python3 solver.py test.in

if __name__ == '__main__':
    if len(sys.argv) == 2:
        path = sys.argv[1]
        G = read_input_file(path)
        T, dists = solve(G)
        assert is_valid_network(G, T)
        for d in dists: print(d)
        outputfile = path.replace("in", "out")
        current = read_output_file(outputfile, G)
        if min(dists) < average_pairwise_distance_fast(current):
            print("Improved.")
            write_output_file(T, outputfile)
        write_output_file(T, 'outputs/test.out')
    else:
        inputs = os.fsencode("inputs/")
        timestamp = datetime.datetime.now().timestamp()
        with open('log.csv', 'w') as logfile:
            logwriter = csv.writer(logfile)
            for file in os.listdir(inputs):
                filename = os.fsdecode(file)
                G = read_input_file(os.path.join("inputs/", filename))
                if len(G) > 25: continue
                print("Solving " + filename)
                T, dists = solve(G)
                for d in dists: print(d)
                outputfile = os.path.join("outputs/", filename.replace("in", "out"))
                current = read_output_file(outputfile, G)
                if min(dists) < average_pairwise_distance_fast(current):
                    print("Improved output for " + filename)
                    logwriter.writerow(dists + [filename, timestamp])
                    write_output_file(T, outputfile)