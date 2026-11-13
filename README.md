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

Then commit and push in GitHub Desktop. The site rebuilds itself.

Keep the GEBCO grid in (or one folder below) the NetCDF folder, with "gebco" in its file name,
so seafloor and land are included. New deployments outside 60.5–65.5°S, 55–66°W need a larger
GEBCO download.

To preview locally (optional), install Quarto and run `quarto preview`.

## Size

About 120 MB of data for 16 deployments. GitHub Pages allows published sites up to 1 GB.

## Data and credits

Glider data: NOAA SWFSC Ecosystem Science Division Glider Team, public domain (CC0-1.0),
processed with pyglider and esdglider. Bathymetry: GEBCO_2026 Grid,
doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa.
