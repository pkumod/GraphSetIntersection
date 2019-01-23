# GraphSetIntersection
This repository is about the source code of our paper "Speeding Up Set Intersections in Graph Algorithms using SIMD Instructions".

To compile our code, just type `make` in the `src` directory, and all executables (including `tc, mc, sm, reorder`) will be ok.

The software can be run as follows:

```
./tc [graph_file_name]
```

For example, count the number of triangles in youtube_cont.txt (in original order):

```
./tc ../data/youtube_cont.txt
```

Take the reordered graph file as the input:

```
./tc ../data/youtube_cont_GRO.txt
```

Executables `mc` (maximal clique algorithm) and `sm` (subgraph matching algorithm) can be run in the same way.

To reorder a graph, run:

```
./reorder <graph_file_name> [-order OPT]
```

Optional argument OPT can be: `gro, hybrid, mloggapa, metis, slashburn, bfsr, dfs`. Please refer to our paper for the different graph orderings.


The format of the input graph file is shown as follows.

```
0       1
0       2
1       3
...
```

Each line is an edge of the graph and there should be M lines in the input file if the graph has M edges. The first integer of each line denotes the start node of the edge and the second integer denotes the end node of the edge. The node IDs should be continuous and start with 0.

Edges are directed in default. For undirected graphs, double the edges in both directions.

```
0       1
0       2
1       0
1       3
2       0
3       1
...
```

Real graph datasets can be obtained from SNAP ("http://snap.stanford.edu/data/") and KONECT ("http://konect.uni-koblenz.de/").

Our code are tested on Linux system using GCC 4.8.5 and GCC 5.3.0.

If you use our code for research work, please kindly cite the following paper:

Shuo Han, Lei Zou, and Jeffrey Xu Yu. Speeding Up Set Intersections in Graph Algorithms using SIMD Instructions. Proceedings of SIGMOD, 2018.
