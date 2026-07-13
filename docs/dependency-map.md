# Core dependency map

This map was generated before the core refactor by parsing the public entry
points and tracing package-local calls. Dynamic calls and S4 dispatch were
checked manually.

## `CreatScDNSobject()`

`scDNS` S4 class; input/network validation; `getDoubleDropout()` ->
`preFilterNet()` -> `detect_Duplicate_edge4net()`; `sub2ind()`;
`CalGeneInfAmount()`; inherited normalization; directed-bipartite attribute
preservation.

## `scDNS_1_CalDivs()`

`getKLD_cKLDnetwork()` -> KNN/KDE density estimation, radius and raw-density
utilities, coarse graining, technology-noise removal, joint and conditional
density construction, JSD/KLD matrix calculations, parallel batching and
matrix/index utilities.

## `scDNS_2_creatNEAModel_v2()`

Dropout/support matrices; directed bipartite edge swaps; candidate edge pools;
`creatNEAModel_test()`; shuffled labels; inherited random-network divergence;
distribution fitting; degree/support adjustment; `GeneVariability`; parallel
divergence calculation and statistical utilities.

## `scDNS_3_GeneZscore_v2()`

`getZscore_v2()` -> edge-score aggregation, degree normalization, random-model
fits, label-shuffle and random-distribution models; Stouffer combination;
P-value/FDR/rank/accuracy utilities. The sf_event adaptation appends
`node_type` without changing numeric score calculations.

## `scDNS_4_scContribution()`

The retained implementation delegates to the numerically guarded inherited
cell-contribution implementation: network subsetting; point-wise divergence;
joint/conditional density indices; node-edge contribution weights; matrix
aggregation and cell-level `scZscore`. In sf_event mode the default node set is
restricted to network sources.

## Removed non-dependencies

Automatic SF-to-gene-to-event expansion, standalone SF activity scoring,
activity differential testing and correlation-based rewiring were not called
by any of the five entry points or their recursive package-local dependencies.
