# visualize-glider-data-site

AMLR glider deployments – website

A [Quarto](https://quarto.org) website for NOAA SWFSC AMLR Slocum glider deployments:

- **Overview** – map of all deployments, filterable by start year and name
- **Deployments** – one page per deployment with a 3D curtain (temperature, salinity,
  chlorophyll, oxygen), a depth–time section with seafloor, σθ contours, mixed-layer and
  chlorophyll-maximum depth, a day/night strip, and a T–S diagram
- **Methods** – how everything is computed, with references

GitHub builds and publishes the site automatically every time you push to `main`
(see `.github/workflows/publish.yml`). You do not need Quarto installed.

## What's in the repository

| Path | What it is |
|---|---|
| `index.qmd`, `deployments.qmd`, `methods.qmd`, `deployments/*.qmd` | site pages (markdown) |
| `_quarto.yml` | site settings and navigation bar |
| `map/index.html` | overview map (generated) |
| `viewer/index.html`, `viewer/plotly.min.js` | shared 3D curtain viewer (generated) |
| `viewer/data/*.json` | one data file per deployment (generated) |
| `scripts/` | Python that generates the files above from the NetCDF files |

The NetCDF files and the GEBCO grid are **not** stored here; only the generated data is.

## Publishing

The site is published at **https://sciaid.github.io/visualize-glider-data-site/** by the
workflow in `.github/workflows/publish.yml`, which runs on every push to `main`.

One-time GitHub setting: **Settings → Pages → Build and deployment → Source: GitHub Actions**.

## Adding or updating deployments

On the computer with the NetCDF files (Python 3 needed):

```
pip install -r requirements.txt
python scripts/build_site.py "C:\Claude_work\PROJECTS\Glider L-3 ncfiles"
```

Then commit and push (or open a pull request). The site rebuilds itself.

- **Only new or changed files are processed.** Each NetCDF file's fingerprint is stored in
  `viewer/data/index.json`; unchanged deployments keep byte-identical data files, so git does not
  store new copies of them. Add `--force` to rebuild everything (e.g. after changing the processing).
- **Folders:** give one or more folders; they are searched recursively, so a layout like
  `ncfiles\Antarctic\`, `ncfiles\Hawaii\`, `ncfiles\California\` works. The overview map is
  built from all files found, so keep every deployment's file available when rebuilding.
  `--prune` removes deployments whose file is no longer in the folders.

### Regions

`regions.json` lists named lat/lon boxes (Antarctic Peninsula, Hawaii, California). Each
deployment is assigned by its median position; anything outside every box is "Other". When
more than one region is present, the map gets a Region filter (opening on the region with the
most deployments), the Deployments page is grouped by region, and the viewer's picker is
grouped too. Edit the boxes or names freely, then rerun the build.

### Bathymetry

Put one GEBCO grid per region (download at [download.gebco.net](https://download.gebco.net),
2D netCDF, ice-surface elevation) in the NetCDF folder or one folder below it, keeping
"gebco" in the file name. Each deployment uses whichever grid covers at least 95% of its
positions (the smallest such grid if several do); deployments with no covering grid are
built without seafloor and the build prints a note.

## Size

About 120 MB of data for 16 deployments (~7.5 MB each). GitHub Pages allows published sites up to
1 GB, so roughly 130 deployments at the current resolution. Beyond that, reduce precision or host
the data files outside the repository.

## Data and credits

Glider data: NOAA SWFSC Ecosystem Science Division Glider Team, public domain (CC0-1.0),
processed with pyglider and esdglider. Bathymetry: GEBCO_2026 Grid,
doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa.
