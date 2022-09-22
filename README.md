Comparing naive dynamic approaches of the Louvain algorithm for community
detection.

`TODO`

[Louvain] is an algorithm for **detecting communities in graphs**. *Community*
*detection* helps us understand the *natural divisions in a network* in an
**unsupervised manner**. It is used in **e-commerce** for *customer*
*segmentation* and *advertising*, in **communication networks** for *multicast*
*routing* and setting up of *mobile networks*, and in **healthcare** for
*epidemic causality*, setting up *health programmes*, and *fraud detection* is
hospitals. *Community detection* is an **NP-hard** problem, but heuristics exist
to solve it (such as this). **Louvain algorithm** is an **agglomerative-hierarchical**
community detection method that **greedily optimizes** for **modularity**
(**iteratively**).

**Modularity** is a score that measures *relative density of edges inside* vs
*outside* communities. Its value lies between `−0.5` (*non-modular clustering*)
and `1.0` (*fully modular clustering*). Optimizing for modularity *theoretically*
results in the best possible grouping of nodes in a graph.

Given an *undirected weighted graph*, all vertices are first considered to be
*their own communities*. In the **first phase**, each vertex greedily decides to
move to the community of one of its neighbors which gives greatest increase in
modularity. If moving to no neighbor's community leads to an increase in
modularity, the vertex chooses to stay with its own community. This is done
sequentially for all the vertices. If the total change in modularity is more
than a certain threshold, this phase is repeated. Once this **local-moving**
**phase** is complete, all vertices have formed their first hierarchy of
communities. The **next phase** is called the **aggregation phase**, where all
the *vertices belonging to a community* are *collapsed* into a single
**super-vertex**, such that edges between communities are represented as edges
between respective super-vertices (edge weights are combined), and edges within
each community are represented as self-loops in respective super-vertices
(again, edge weights are combined). Together, the local-moving and the
aggregation phases constitute a **pass**. This super-vertex graph is then used
as input for the next pass. This process continues until the increase in
modularity is below a certain threshold. As a result from each pass, we have a
*hierarchy of community memberships* for each vertex as a **dendrogram**. We
generally consider the *top-level hierarchy* as the *final result* of community
detection process.

*Louvain* algorithm is a hierarchical algorithm, and thus has two different
tolerance parameters: `tolerance` and `passTolerance`. **tolerance** defines the
minimum amount of increase in modularity expected, until the local-moving phase
of the algorithm is considered to have converged. We compare the increase in
modularity in each iteration of the local-moving phase to see if it is below
`tolerance`. **passTolerance** defines the minimum amount of increase in
modularity expected, until the entire algorithm is considered to have converged.
We compare the increase in modularity across all iterations of the local-moving
phase in the current pass to see if it is below `passTolerance`. `passTolerance`
is normally set to `0` (we want to maximize our modularity gain), but the same
thing does not apply for `tolerance`. Adjusting values of `tolerance` between
each pass have been observed to impact the runtime of the algorithm, without
significantly affecting the modularity of obtained communities.

In this experiment, we compare the performance of *two different types* of **naive**
**dynamic Louvain** with respect to the *static* version. The **last** approach
(`louvainSeqDynamicLast`) considers the community membership of each vertex
*after* the Louvain algorithm has *converged* (community membership from the "last"
pass) and then performs the Louvain algorithm upon the new (updated) graph. This
is *similar* to naive dynamic approaches with other algorithms. On the other hand,
the **first** approach (`louvainSeqDynamicFirst`) considers the community
membership of each vertex right after the *first pass* of the Louvain algorithm
(this is the first community membership hierarchy) and then performs the Louvain
algorithm upon the updated graph. With this approach, we allow the affected
vertices to choose their community membership from the first pass itself, which
which to my intuition would lead to better communities.

First, we compute the community membership of each vertex using the static
Louvain algorithm (`louvainSeqLast`). We also run the static Louvain algorithm
for only one pass (`louvainSeqFirst`). We then generate *batches* of *insertions*
*(+)* and *deletions (-)* of edges of sizes 500, 1000, 5000, ... 100000. For each
batch size, we generate *five* different batches for the purpose of *averaging*.
Each batch of edges (insertion / deletion) is generated randomly such that the
selection of each vertex (as endpoint) is *equally probable*. We choose the
Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for local-moving
phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDeclineFactor = 10`, and a `passTolerance = 0.0` (when passes stop).
In addition we limit the maximum number of iterations in a single local-moving
phase with `maxIterations = 500`, and limit the maximum number of passes with
`maxPasses = 500`. We run the Louvain algorithm until convergence (or until the
maximum limits are exceeded), and measure the **time** **taken** for the
*computation* (performed 5 times for averaging), the **modularity score**, the
**total number of iterations** (in the *local-moving* *phase*), and the number
of **passes**. This is repeated for *seventeen* different graphs.

From the results, we make make the following observations. The performance of
dynamic approaches upon a batch of deletions appears to *increase* with *increasing*
batch size*. This makes sense since, as the graph keeps getting smaller, the
computation would complete *sooner*. Next, the `first` naive dynamic approach is
found to be *significantly slower* (~0.3x speedup) than the `last` approach.
However, the `first` approach is *still faster* than the static approach upto a
batch size of `50000`. On the other hand, the `last` approach is *faster* than the
static approach for all batch sizes. A similar behavior is observed with the
total number of iterations. The `first` approach seems to have a *slightly higher*
modularity with respect to the `last` approach. Since the modularity between the
two dynamic approaches are almost the same, the **last** approach is clearly the
**best choice**.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 3985272 [directed] {} (symmetricize)
# [-0.000497 modularity] noop
# [0e+00 batch_size; 0 batch_count; 00792.589 ms; 0025 iters.; 009 passes; 0.923382580 modularity] louvainSeqLast
# [0e+00 batch_size; 0 batch_count; 00573.409 ms; 0004 iters.; 001 passes; 0.766543329 modularity] louvainSeqFirst
# [5e+02 batch_size; 1 batch_count; 00224.703 ms; 0004 iters.; 004 passes; 0.914939582 modularity] louvainSeqDynamicLast
# [5e+02 batch_size; 1 batch_count; 00480.665 ms; 0025 iters.; 009 passes; 0.923243225 modularity] louvainSeqDynamicFirst
# [5e+02 batch_size; 2 batch_count; 00227.230 ms; 0004 iters.; 004 passes; 0.914955676 modularity] louvainSeqDynamicLast
# ...
# [-1e+05 batch_size; 4 batch_count; 00397.172 ms; 0011 iters.; 006 passes; 0.876155496 modularity] louvainSeqDynamicFirst
# [-1e+05 batch_size; 5 batch_count; 00206.570 ms; 0004 iters.; 004 passes; 0.869377553 modularity] louvainSeqDynamicLast
# [-1e+05 batch_size; 5 batch_count; 00406.253 ms; 0012 iters.; 006 passes; 0.876216054 modularity] louvainSeqDynamicFirst
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 13298940 [directed] {} (symmetricize)
# [-0.000316 modularity] noop
# [0e+00 batch_size; 0 batch_count; 01248.049 ms; 0028 iters.; 009 passes; 0.935839474 modularity] louvainSeqLast
# [0e+00 batch_size; 0 batch_count; 00926.964 ms; 0005 iters.; 001 passes; 0.798873782 modularity] louvainSeqFirst
# [5e+02 batch_size; 1 batch_count; 00527.935 ms; 0003 iters.; 003 passes; 0.932618558 modularity] louvainSeqDynamicLast
# [5e+02 batch_size; 1 batch_count; 00884.497 ms; 0025 iters.; 009 passes; 0.935987055 modularity] louvainSeqDynamicFirst
# [5e+02 batch_size; 2 batch_count; 00529.566 ms; 0003 iters.; 003 passes; 0.932615817 modularity] louvainSeqDynamicLast
# ...
```

[![](https://i.imgur.com/HJAS3Di.png)][sheetp]
[![](https://i.imgur.com/4iQ7CzY.png)][sheetp]
[![](https://i.imgur.com/E9nDrAI.png)][sheetp]
[![](https://i.imgur.com/BZGF6Yt.png)][sheetp]

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>

[![](https://i.imgur.com/9HITKSz.jpg)](https://www.youtube.com/watch?v=wCUV6N4Qtew&t=447s)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/519984922.svg)](https://zenodo.org/badge/latestdoi/519984922)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[Louvain]: https://en.wikipedia.org/wiki/Louvain_method
[gist]: https://gist.github.com/wolfram77/9c1bff3cc327acd80c9e2479ef7c4e57
[charts]: https://imgur.com/a/3vhRU3c
[sheets]: https://docs.google.com/spreadsheets/d/189GRfvpTxSMWLrqafHMvHyV7ddaXFZQ0u056MfnO2uU/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRrGpRtzagVZCmtPuIUaD03I8SGY2PEZGusNV90ojCgntRbiEg0r8wCp-YiT8A7e8ZqzqQqAJveqGOD/pubhtml
