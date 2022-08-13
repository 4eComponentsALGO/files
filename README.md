This repository is an accompaniment to ESA22 Track B paper 60, titled "Computing the 4-Edge-Connected Components of a Graph: An Experimental Study".

It contains implementations of various linear-time algorithms for computing 3-edge-cuts in 3-edge-connected undirected graphs, and 4-edge-connected components in general undirected graphs.

The programs get3cutsGIK, get3cutsLinear, get3cutsNRSS and get3cutsNRSSsort, compute all 3-edge cuts of a 3-edge-connected graph.

 * get3cutsGIK implements the algorithm of Georgiadis, Italiano, and Kosinas. (Full version at: https://arxiv.org/abs/2105.02910.)

 * get3cutsLinear implements our improved algorithm which runs in linear time in pointer-machines. (See https://arxiv.org/abs/2108.08558.)

 * get3cutsNRSS and get3cutsNRSSsort implement two different variations of (an improvement of) the randomized algorithm of Nadara, Radecki, Smulewicz, and Sokołowski. (Full version at: https://arxiv.org/abs/2105.01699.) 

The respective get4eComponents* algorithms compute the 4-edge-connected components in general graphs. Specifically, they compute the number of k-edge-connected components, for every k<=4, the number of minimal k-edge cuts, for every k<=3, and they output in a file a label for every vertex signifying the 4-edge-connected component in which it belongs.

All programs output the time they need to perform their computations. For this purpose we use the "high_resolution_clock" function of the standard library "chrono".

All programs take as input a text file representing a graph (undirected multigraph).
The format of a graph-file is as follows.
 * The first line contains two integers, n and m, representing the number of vertices and edges, respectively, of the graph.
 * The next m lines have the form "x y", meaning that there is an edge (x,y).
The vertices are numbered from 0 to n-1 (inclusive).
Multiple edges are allowed.

For example, the following file (K4.txt) represents the complete graph on four vertices.

4 6  
0 1  
0 2  
0 3  
1 2  
1 3  
2 3  

A larger example of a graph file (used in our experiments) is "artist.txt" and "artist-SC.txt" (its sparse certificate for 4-edge connectivity).

Every 3-cut or 4-component program P is called as ./P <input_graph> <output_file>.

For the get3cuts* programs, the output_file contains the 3-edge cuts of the (3-edge-connected) input graph, and has the following format.
  * The first line contains the number k of 3-edge cuts.
  * The next k lines have the form "x1 y1 x2 y2 x3 y3", signifying the existence of 3-edge cut {(x1,y1),(x2,y2),(x3,y3)}.
  * The last line contains the time it took for the total computation.

For example, a sample run of ./get3cutsGIK K4.txt out.txt gives as output

4  
1 0 0 3 0 2  
3 2 0 3 1 3  
3 2 2 1 2 0  
2 1 1 0 3 1  
0.000008


For the get4eComponents* programs, the output_file provides labels to the vertices, signifying the 4-edge-connected component in which they belong, and has the following format.
  * The first line contains the number n of vertices.
  * The next n lines contain the labels. Line i is the label for vertex i-1.
  * The second to last line contains the time it took for the total computation.
  * The last line contains the total time it took for computing the 4-edge-connected components of the auxiliary 3-edge-connected graphs.

For example, a sample run of ./get4eComponentsGIK K4.txt out.txt gives as output:

4  
0  
2  
3  
1  
0.000103  
0.000045

(This graph has no non-singleton 4-edge-connected components.)


We also provide the following additional programs.

genrand.cpp generates a random multigraph, with a specified number of vertices and edges.
It is called as ./genrand n m seed <output_graph>, where n is the number of vertices, m is the number of edges, and seed is a number to initialize the rand function.

make3eConnected.cpp make a general graph 3-edge-connected by introducing a few more edges to the graph.  
It is called as ./make3eConnected <input_graph> <output_graph>.

sparsify.cpp creates a sparse certificate of a graph (for 4-edge connectivity). It is called as ./sparsify <input_graph> <output_graph>.

