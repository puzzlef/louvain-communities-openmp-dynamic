Design of OpenMP-based Dynamic [Louvain algorithm] for [community detection].

This is an implementation of the popular [Louvain algorithm]. It is an
**agglomerative-hierarchical** community detection method that *greedily*
*optimizes* for [modularity]. Given an undirected weighted graph, all
vertices are first considered to be their own communities. In the
**local-moving phase**, each vertex greedily decides to move to the community of
one of its neighbors which gives greatest increase in modularity (multiple
iterations). In the **aggregation phase**, all vertices belonging to a
community* are *collapsed* into a single *super-vertex*. This super-vertex
graph is then used as input for the next pass. The process continues until the
increase in modularity is below a certain threshold. As a result from each pass,
we have a *hierarchy of community memberships* for each vertex as a
**dendrogram**. See [extended report] for details. For HIPC2023 submission, see
[submission-hipc23].

There are three different dynamic approaches we are trying out:
- **Naive-dynamic**: We simply use the previous community memberships and perform the algorithm.
- **Dynamic Delta-screening**: We find a set of affected vertices as per the [Delta-screening] paper.
- **Dyanmic Frontier**: We mark endpoints of each vertex as affected, and expand out iteratively.

The input data used for below experiments is available from the [SuiteSparse Matrix Collection].
The experiments were done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].

[Louvain algorithm]: https://en.wikipedia.org/wiki/Louvain_method
[extended report]: https://gist.github.com/wolfram77/91b2d2ac50b9aba6b203e88b291c7671
[submission-hipc23]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/submission-hipc23
[community detection]: https://en.wikipedia.org/wiki/Community_search
[modularity]: https://en.wikipedia.org/wiki/Modularity_(networks)
[Delta-screening]: https://ieeexplore.ieee.org/document/9384277
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>


### Comparision on large graphs

In this experiment ([input-large]), we first compute the community membership of
each vertex using the static Louvain algorithm. We then generate random batch
updates consisting of an equal mix of *deletions (-)* and  *insertions (+)* of
edges of size `10^-7 |E|` to `0.1 |E|` in multiples of `10` (where `|E|` is the
number of edges in the original graph after making it undirected). For each
batch size, we generate *five* different batches for the purpose of *averaging*.
Each batch of edges (insertion / deletion) is generated randomly such that the
selection of each vertex (as endpoint) is *equally probable*. We choose the
Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for local-moving
phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDrop = 10`, an `aggregationTolerance = 0.8` which considers
communities to have converged when a few communities merge together, and a
`passTolerance = 0.0` (when passes stop). In addition we limit the maximum
number of iterations in a single local-moving phase with `maxIterations = 20`,
and limit the maximum number of passes with `maxPasses = 10`. We run the Louvain
algorithm until convergence (or until the maximum limits are exceeded), and
measure the **time taken** for the *computation* and *pre-processing* (for
dynamic approaches), the **modularity** **score**, the **total number of**
**iterations** (in the *local-moving phase*), and the number of **passes**. This
is repeated for each input graph.

From the results, we make make the following observations. **Dynamic Frontier**
based **Louvain** converges the fastest, which obtaining communities with
equivalent modularity. We also observe that **Dynamic Delta-screening** based
**Louvain** has poorer performance than the Naive-dynamic approach (due to its
high pre-processing cost/overhead). Therefore, **Dynamic Frontier based**
**Louvain** would be the **best choice**. We also not that **Louvain** algorithm
does not scale too well with an increase in the number of threads. This is
likely due to higher pressure on cache coherence system as well as the algorithm
becoming closer to a synchronous approach, which is inherently slower than an
asynchronous approach. Trying to avoid community swaps with parallel approach
does not seem to improve performance by any significant amount. However, it is
possible that if synchronous approach is used with OpenMP, then its performance
may be a bit better.

> See
> [code](https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/input-large),
> [output](https://gist.github.com/wolfram77/adbef451db5bf46f1a7243349121a860), or
> [sheets].


[![](https://i.imgur.com/qCYVeh4.png)][sheets]
[![](https://i.imgur.com/4PBroEt.png)][sheets]
[![](https://i.imgur.com/T2LG2RJ.png)][sheets]
[![](https://i.imgur.com/WdOu1ON.png)][sheets]
[![](https://i.imgur.com/fGrM5an.png)][sheets]
[![](https://i.imgur.com/3xfJBsD.png)][sheets]

We also measure the pass-wise and phase-wise split of time taken for each
approach. **Dynamic Frontier** approach only performs more than one pass when
batch size is large enough.

[![](https://i.imgur.com/qccysjM.png)][sheets]
[![](https://i.imgur.com/gK9NS6b.png)][sheets]
[![](https://i.imgur.com/p8439Ry.png)][sheets]
[![](https://i.imgur.com/kCyhyMf.png)][sheets]

Each pass of Louvain is divided into two main phases: *local-moving*, and
*aggregation*. In addition to the two, *other work* also need to be
performed. These include:
- Initializing vertex weights, community weights, and community memberships
- Re-numbering community IDs
- Flattening dendrogram (lookups)
- Obtaining vertices beloning to each community
- Clearing various buffers

With **Dynamic Frontier**, most time is spent doing this *other work*, followed by
*local-moving* phase, and then the *aggregation* phase.

[![](https://i.imgur.com/DrHWAgO.png)][sheets]
[![](https://i.imgur.com/2EFYN1X.png)][sheets]
[![](https://i.imgur.com/5eSzJzf.png)][sheets]
[![](https://i.imgur.com/xSREmHM.png)][sheets]

[sheets]: https://docs.google.com/spreadsheets/d/1F6Z-lWNDYynm6m2PTsIN_nxMu8Y9CrkIQagCU0Nr2LU/edit?usp=sharing

<br>


### Measure communities

In this experiment ([measure-communities]), we **measure** the **properties of**
**communities obtained** with *Static Louvain* algorithm. These include the
*number of communities*, the *size distribution of communities* (*gini*
*coefficient*), and the *overall modularity score*.

[measure-communities]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/measure-communities

<br>


### Measure affected vertices

In this experiment ([measure-affected]), we **measure** the number of **affected**
**vertices** with *Dynamic* *Delta-screening* and *Dynamic Frontier* based
*Louvain* for random batch updates consisting of edge insertions, with the size
of batch update varying from `10^-6 |E|` to `0.1 |E|`.

Results show that *Dynamic Delta-screening* marks `15000x`, `2000x`, `440x`,
`44x`, `6.4x`, and `1.7x` the number of affected vertices as *Dynamic Frontier*
based approach on batch updates of size `10^-6 |E|` to `0.1 |E|`.

[measure-affected]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/measure-affected

<br>


### Multi-batch updates

In this experiment ([multi-batch]), we generate `5000` random **multi-batch updates** consisting
of *edge insertions* of size `10^-3 |E|` one after the other on graphs
`web-Stanford` and `web-BerkStan` and observe the performance and modularity of
communities obtained with *Static*, *Naive-dynamic*, *Dynamic Delta-screening*,
and *Dynamic Frontier* based *Louvain*. We do this to measure after how many
batch updates do we need to re-run the static algorithm.

Our results indicate that we need to rerun the static algorithm after `~1300`
batch updates with *Dynamic Delta-screening* based *Louvain*, and after `~2800`
batch updates with *Dynamic Frontier* based *Louvain*.

[multi-batch]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/multi-batch

<br>
<br>


## Build instructions

To run the [input-large] experiment, download this repository and run the
following. Note that input graphs must be placed in `~/Data` directory, and
output logs will be written to `~/Logs` directory.

```bash
# Perform comparision on large graphs
$ DOWNLOAD=0 ./mains.sh

# Perform comparision on large graphs with custom number of threads
$ DOWNLOAD=0 MAX_THREADS=4 ./mains.sh
```

[input-large]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/input-large


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


[![](https://i.imgur.com/UGB0g2L.jpg)](https://www.youtube.com/watch?v=pIF3wOet-zw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
