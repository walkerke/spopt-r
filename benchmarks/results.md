# Stage timings: before (branch fix/azp-maxp-search) vs after (perf/regionalization-hotpaths)

Seconds, single run, hexagonal grid fixture with 4 attributes, neighbour list supplied. Speedup = before / after.

| n | stage | before | after | speedup |
|---|---|---:|---:|---:|
| 374 | weights: prepare_weights(bridge) | 0.117 | 0.040 | 2.9x |
| 500 | azp basic from skater, nb given | 0.039 | 0.028 | 1.4x |
| 500 | azp sa from skater, nb given | 0.309 | 0.062 | 5.0x |
| 500 | max_p 20k floor, nb given | 0.068 | 0.041 | 1.7x |
| 500 | max_p construction only | 0.022 | 0.018 | 1.2x |
| 500 | skater, floor, nb given | 0.088 | 0.007 | 12.6x |
| 500 | skater, nb given | 0.094 | 0.007 | 13.4x |
| 500 | spenc, nb given | 0.075 | 0.062 | 1.2x |
| 500 | ward_spatial, nb given | 0.021 | 0.007 | 3.0x |
| 500 | weights: poly2nb queen | 0.566 | 0.449 | 1.3x |
| 500 | weights: prepare_weights() | 0.056 | 0.041 | 1.4x |
| 2235 | weights: prepare_weights(bridge) | 2.260 | 0.460 | 4.9x |
| 2500 | azp basic from skater, nb given | 0.185 | 0.130 | 1.4x |
| 2500 | azp sa from skater, nb given | 6.627 | 0.940 | 7.0x |
| 2500 | max_p 20k floor, nb given | 0.436 | 0.321 | 1.4x |
| 2500 | max_p construction only | 0.193 | 0.180 | 1.1x |
| 2500 | skater, floor, nb given | 1.932 | 0.064 | 30.2x |
| 2500 | skater, nb given | 2.203 | 0.067 | 32.9x |
| 2500 | spenc, nb given | 7.601 | 6.698 | 1.1x |
| 2500 | ward_spatial, nb given | 0.392 | 0.080 | 4.9x |
| 2500 | weights: poly2nb queen | 0.179 | 0.141 | 1.3x |
| 2500 | weights: prepare_weights() | 0.203 | 0.165 | 1.2x |
| 10000 | azp basic from skater, nb given | 1.815 | 1.001 | 1.8x |
| 10000 | azp sa from skater, nb given | 107.164 | 5.138 | 20.9x |
| 10000 | max_p 20k floor, nb given | 2.992 | 2.233 | 1.3x |
| 10000 | max_p construction only | 1.188 | 1.309 | 0.9x |
| 10000 | skater, floor, nb given | 30.994 | 0.801 | 38.7x |
| 10000 | skater, nb given | 30.967 | 0.780 | 39.7x |
| 10000 | ward_spatial, nb given | 7.294 | 0.703 | 10.4x |
| 10000 | weights: poly2nb queen | 1.211 | 1.286 | 0.9x |
| 10000 | weights: prepare_weights() | 0.842 | 0.944 | 0.9x |
