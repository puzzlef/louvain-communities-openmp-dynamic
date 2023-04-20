In this experiment, we generate `5000` random **multi-batch updates** consisting
of *edge insertions* of size `10^-3 |E|` one after the other on graphs
`web-Stanford` and `web-BerkStan` and observe the performance and modularity of
communities obtained with *Static*, *Naive-dynamic*, *Dynamic Delta-screening*,
and *Dynamic Frontier* based *Louvain*. We do this to measure after how many
batch updates do we need to re-run the static algorithm.

Our results indicate that we need to rerun the static algorithm after `~1300`
batch updates with *Dynamic Delta-screening* based *Louvain*, and after `~2800`
batch updates with *Dynamic Frontier* based *Louvain*. Some [charts] are also
included below, generated from [sheets].

<br>

[![](https://i.imgur.com/DyYKDdD.png)][sheetp]
[![](https://i.imgur.com/sbLogMq.png)][sheetp]
[![](https://i.imgur.com/vVEjSDZ.png)][sheetp]
[![](https://i.imgur.com/pwbLHJb.png)][sheetp]
[![](https://i.imgur.com/RCLYopU.png)][sheetp]
[![](https://i.imgur.com/VU0wTq1.png)][sheetp]
[![](https://i.imgur.com/Jn9uPBF.png)][sheetp]

<br>
<br>


[![](https://i.imgur.com/UGB0g2L.jpg)](https://www.youtube.com/watch?v=pIF3wOet-zw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)


[charts]: https://imgur.com/a/PNoOpuV
[sheets]: https://docs.google.com/spreadsheets/d/1hN5dJ9VWAAORWwSdC-L3Ww4Svfnhs8yokjdWppFRp4A/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRHsslN9ghcDMO-O1bv8GBDqo7mNdA5Lm9u7s-ftRmgbyguRxsUEF1yCXoovzxVabInR3j0G9MJxJcY/pubhtml
