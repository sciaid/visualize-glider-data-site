# Seafloor and land grids

Subsets of the GEBCO_2026 Grid (15 arc-second), downloaded from https://download.gebco.net
(2D netCDF), used by the daily GitHub build for seafloor and land under each deployment.

| File | North | South | West | East |
|---|---|---|---|---|
| gebco_2026_antarctic_peninsula.nc | -60.5 | -65.5 | -66 | -55 |
| gebco_2026_california.nc | 50 | 30 | -135 | -114 |
| gebco_2026_hawaii.nc | 21.6 | 20.0 | -158.7 | -157.3 |

A deployment uses the grid covering at least 95% of its positions. The Hawaii grid only covers
the area around Oʻahu; replace it with a larger box (e.g. N 26, S 15, W -165, E -150) before
Hawaii deployments farther afield are added. Keep "gebco" in the file names.

GEBCO Compilation Group (2026) GEBCO_2026 Grid, doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa
