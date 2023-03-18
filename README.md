Comparing static vs dynamic approaches of OpenMP-based [Louvain algorithm]
for [community detection].

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
significantly affecting the modularity of obtained communities. In this
experiment, we compare the performance of *three different types* of OpenMP-based
**dynamic Louvain** with respect to the *static* version. We also compare with
sequential static approach.

**Naive dynamic**:
- We start with previous community membership of each vertex (instead of each vertex its own community).

**Delta screening**:
- All edge batches are undirected, and sorted by source vertex-id.
- For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
  `i`'s neighbors and `j*`'s community is marked as affected.
- For edge deletions within the same community `i` and `j`,
  `i`'s neighbors and `j`'s community is marked as affected.

**Frontier-based**:
- All source and destination vertices are marked as affected for insertions and deletions.
- For edge additions across communities with source vertex `i` and destination vertex `j`,
  `i` is marked as affected.
- For edge deletions within the same community `i` and `j`,
  `i` is marked as affected.
- Vertices whose communities change in local-moving phase have their neighbors marked as affected.

First, we compute the community membership of each vertex using the static
Louvain algorithm. We then generate *batches* of *insertions (+)* and
*deletions (-)* of edges of sizes 500, 1000, 5000, ... 100000. For each batch
size, we generate *five* different batches for the purpose of *averaging*. Each
batch of edges (insertion / deletion) is generated randomly such that the
selection of each vertex (as endpoint) is *equally probable*. We choose the
Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for local-moving
phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDeclineFactor = 10`, and a `passTolerance = 0.0` (when passes stop).
In addition we limit the maximum number of iterations in a single local-moving
phase with `maxIterations = 20`, and limit the maximum number of passes with
`maxPasses = 20`. We run the Louvain algorithm until convergence (or until the
maximum limits are exceeded), and measure the **time taken** for the
*computation* (performed 5 times for averaging), the **modularity score**, the
**total number of iterations** (in the *local-moving phase*), and the number
of **passes**. This is repeated for *seventeen* different graphs.

From the results, we make make the following observations. The frontier-based
dynamic approach converges the fastest, which obtaining communities with only
slightly lower modularity than other approaches. We also observe that
delta-screening based dynamic Louvain algorithm has the same performance as that
of the naive dynamic approach. Therefore, **frontier-based dynamic Louvain**
would be the **best choice**. However, one of the most interesting things we
note is that sequential static approach is only slightly slower than
OpenMP-based approach with 12 threads. This is indeed suprising, and is likely
due to higher pressure on cache coherence system as well as the algorithm
becoming closes to an unordered approach, which is inherently slower than an
ordered approach. Trying to avoid community swaps with OpenMP-based approach
does not seem to improve performance by any significant amount. However, it is
possible that if unordered approach is used with OpenMP, then its performance
may be a bit better.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].

[![](https://i.imgur.com/v08cdyB.png)][sheetp]
[![](https://i.imgur.com/pAUh0Ql.png)][sheetp]
[![](https://i.imgur.com/bA7CMBF.png)][sheetp]
[![](https://i.imgur.com/VlsiC93.png)][sheetp]
[![](https://i.imgur.com/oO9fXcX.png)][sheetp]
[![](https://i.imgur.com/7NGU3wh.png)][sheetp]

[Louvain algorithm]: https://en.wikipedia.org/wiki/Louvain_method
[community detection]: https://en.wikipedia.org/wiki/Community_search

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

[![](https://i.imgur.com/UGB0g2L.jpg)](https://www.youtube.com/watch?v=pIF3wOet-zw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/540003159.svg)](https://zenodo.org/badge/latestdoi/540003159)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[Louvain]: https://en.wikipedia.org/wiki/Louvain_method
[gist]: https://gist.github.com/wolfram77/2d64f933f6524ba15ee7593f7e3b10f5
[charts]: https://imgur.com/a/QvH0m2i
[sheets]: https://docs.google.com/spreadsheets/d/1F6Z-lWNDYynm6m2PTsIN_nxMu8Y9CrkIQagCU0Nr2LU/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vS5LH03ALzgcv6QNV9I9Wl1_000Vl9BNZKnMdF04d4qeG5dqQ60fFHL4xynG_8LnVFbsyaJAucWuen6/pubhtml
