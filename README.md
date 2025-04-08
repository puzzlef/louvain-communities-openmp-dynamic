Design of OpenMP-based Parallel Dynamic [Louvain algorithm] for [community detection].

<br>

Community detection is the problem of recognizing natural divisions in networks. A relevant challenge in this problem is to find communities on rapidly evolving graphs. In this report we present our Parallel Dynamic Frontier (DF) Louvain algorithm, which given a batch update of edge deletions and insertions, incrementally identifies and processes an approximate set of affected vertices in the graph with minimal overhead, while using a novel approach of incrementally updating weighted-degrees of vertices and total edge weights of communities. We also present our parallel implementations of Naive-dynamic (ND) and Delta-screening (DS) Louvain. On a server with a 64-core AMD EPYC-7742 processor, our experiments show that DF Louvain obtains speedups of 179x, 7.2x, and 5.3x on real-world dynamic graphs, compared to Static, ND, and DS Louvain, respectively, and is 183x, 13.8x, and 8.7x faster, respectively, on large graphs with random batch updates. Moreover, DF Louvain improves its performance by 1.6x for every doubling of threads.

<br>

Below we illustrate the mean runtime and modularity of communities obtained with our parallel implementation of Static, Naive-dynamic (ND), Delta-screening (DS), and Dynamic Frontier (DF) Louvain on real-world dynamic graphs on batch updates of size `10^-5|Eᴛ|` to `10^-3|Eᴛ|` (where `|Eᴛ|` is the number of temporal edges). In (a), the speedup of each approach with respect to Static Louvain is labeled.

[![](https://i.imgur.com/cdB6Aq0.png)][sheets-o1]

Next, we plot the average time taken by Static, ND, DS, and DF Louvain on large graphs with random batch updates of size `10^-7|E|` to `0.1|E|`. In (a), the speedup of each approach with respect to Static Louvain is labeled.

[![](https://i.imgur.com/mw4zoeE.png)][sheets-o2]

Refer to our technical report for more details: \
[DF Louvain: Fast Incrementally Expanding Approach for Community Detection on Dynamic Graphs][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.

[Louvain algorithm]: https://en.wikipedia.org/wiki/Louvain_method
[community detection]: https://en.wikipedia.org/wiki/Community_search
[modularity]: https://en.wikipedia.org/wiki/Modularity_(networks)
[Delta-screening]: https://ieeexplore.ieee.org/document/9384277
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[sheets-o1]: https://docs.google.com/spreadsheets/d/1Ad-D7n9bZCOEIGNhfHobUIJDQNR3CLZnrGGc8AcnDc0/edit?usp=sharing
[sheets-o2]: https://docs.google.com/spreadsheets/d/1zg93KjbVUEsjy-U1lg-CmDm5UZTSVEG93uPTPrwf9XA/edit?usp=sharing
[report]: https://arxiv.org/abs/2404.19634

<br>
<br>


### Code structure

The code structure of our multicore implementation of Dynamic Frontier (DF) Louvain is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/louvain.hxx: Louvain algorithm functions
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetrize.hxx: Graph Symmetrization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [Fetch-and-add using OpenMP atomic operations](https://stackoverflow.com/a/7918281/1413259)

<br>
<br>


[![](https://i.imgur.com/3abceEx.png)](https://www.youtube.com/watch?v=yqO7wVBTuLw&pp)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/louvain-communities-openmp-dynamic)
