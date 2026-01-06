# Vignette Screenshots

This folder contains screenshots for the spopt vignettes. The following images need to be generated and placed here:

## Regionalization vignette (regionalization.qmd)

- `dallas-income.png` - Dallas County Census tracts colored by median household income
- `maxp-result.png` - Max-P regionalization results showing custom regions
- `skater-result.png` - SKATER clustering results with 15 regions
- `azp-result.png` - AZP results with 20 regions using tabu search
- `spenc-result.png` - SPENC spectral clustering results
- `ward-result.png` - Ward spatial clustering results

## Facility location vignette (facility-location.qmd)

- `facility-setup.png` - Problem setup showing demand points (blue) and candidate facilities (black)
- `pmedian-result.png` - P-Median solution with 5 facilities
- `pcenter-comparison.png` - Comparing P-Median (blue) and P-Center (red) facility locations
- `pdispersion-result.png` - P-Dispersion solution with 10 maximally spread facilities

## Huff model vignette (huff-model.qmd)

- `huff-market-areas.png` - Market areas by primary store assignment
- `huff-entropy.png` - Competition intensity (entropy) across the market

## Travel-time matrices vignette (travel-time-matrices.qmd)

- `ttm-comparison.png` - Comparing P-Median solutions: travel time (red) vs Euclidean (blue)

## Generating screenshots

Run each vignette code block, then use mapgl's screenshot functionality or take manual screenshots of the resulting maps. Images should be approximately 800-1000px wide for optimal display.
