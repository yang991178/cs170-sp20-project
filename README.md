# CS 170 Project Spring 2020

Files:
- `parse.py`: functions to read/write inputs and outputs
- `solver.py`: main solver program for the problem
- `solver.ipynb`: helper file that invokes solver.py on an array of inputs
- `utils.py`: contains functions to compute cost and validate NetworkX graphs

Run `python solver.py` to execute the solver on all files after extracting them to `/inputs`. 
Change the `h=` parameters of `spdt` and `msdt` functions to try different heuristics.

Computing all outputs should take a few hours (tested on an i7-9750H CPU). 
Ranked #6 overall on the [project leaderboard](https://berkeley-cs170.github.io/project-leaderboard/).

## Reflection

The algorithm can be roughly divided into two main phases: initial tree generation and local optimization. 

For tree generation, a modified version of shortest path tree is used. A SPT is generated for every vertex in the graph and the algorithm yields the partial SPT when the vertices in the tree constitute a dominating set for the graph. For small and medium graphs, the complete SPT is also yielded. Different heuristics are also optionally included in the SPT generation algorithm to diversify the initial trees. These heuristics divide the weight of edges by the degree of the connected vertices when pushing vertices onto the heap and result in (incorrect) SPTs that may have smaller average routing costs. MSTs are also considered in the tree generation by running Prim’s from every vertex, although their costs are usually higher than SPTs. Heuristics considering degrees of vertices are also provided for MSTs.

After obtaining the initial trees, three local optimization algorithms are applied to them. The first one considers every edge from the tree and try replacing them with another edge from the graph to reduce the average routing cost, taken from [Appendix A of this paper](http://www.orstw.org.tw/ijor/vol10no4/ijor_vol10_no4_153_160.pdf). The second one is a combination of two contraction algorithms. The first contraction algorithm greedily prunes leaves if removing them reduces the average routing cost while keeping the tree dominating. The order of pruning may be randomized for smaller graphs. The second contraction algorithm prune all leaves regardless of the cost and only keep vertices such that they still constitute a dominating set. The contracted tree with the lowest cost among randomized runs of the greedy algorithm and the cost-blind algorithm is kept. The third optimizing algorithm considers all vertices directly neighboring the tree in the graph and adds them to the tree if their inclusion can lead to a lower average cost. 

After applying the optimization algorithms to the initial trees, the tree with the lowest average routing cost is outputted.
