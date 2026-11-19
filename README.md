# visualize-glider-data-site

ESD glider deployments – website

A [Quarto](https://quarto.org) website for glider deployments by the NOAA SWFSC Ecosystem Science
Division (ESD) Glider Team, including the Antarctic Marine Living Resources (AMLR) program:

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
  `ncfiles\Antarctica\`, `ncfiles\Hawaii\`, `ncfiles\California\` works. The overview map is
  built from all files found, so keep every deployment's file available when rebuilding.
  `--prune` removes deployments whose file is no longer in the folders.

### Daily build from ERDDAP (automatic)

`.github/workflows/publish.yml` runs every day (and on every push to `main`, and on demand from
**Actions → Publish site → Run workflow**). It downloads new or updated datasets from the IOOS
Glider DAC ERDDAP, builds their pages and publishes the site. Nothing needs to run on your computer.

- **ERDDAP deployments are built on GitHub and not committed**, so the repository does not grow;
  the gridded downloads are kept in GitHub's cache between runs.
- **Deployments published from local L3 files are never replaced** by their ERDDAP versions. To
  switch one over, add its name to `replace_local_with_erddap` in `erddap.json` (or `"*"` for all).
- **Failures don't break the site**: a dataset that fails keeps its previous version; if ERDDAP is
  down, nothing changes. The run is then marked failed (GitHub emails you) and the home page says so.
- The home and **Real-time** pages show when ERDDAP was last checked.
- Seafloor grids for the GitHub build are in `bathymetry/` (see its README).
- Pull requests run the same build, including a real ERDDAP download, without publishing.

### Deployments from ERDDAP

`scripts/fetch_erddap.py` downloads deployments from the
[IOOS Glider DAC ERDDAP](https://gliders.ioos.us/erddap/) and grids them into the same 5 m
layout as the local L3 files, into `<folder>\ERDDAP\`:

```
python scripts/fetch_erddap.py "C:\Claude_work\PROJECTS\Glider L-3 ncfiles" --dry-run   # see the plan
python scripts/fetch_erddap.py "C:\Claude_work\PROJECTS\Glider L-3 ncfiles"
python scripts/build_site.py  "C:\Claude_work\PROJECTS\Glider L-3 ncfiles"
```

- Only datasets from *NOAA National Marine Fisheries Service Southwest Fisheries Science Center*
  (the older *Antarctic Ecosystem Research Division* datasets are skipped).
- Delayed-mode data are used when published; otherwise real-time data, labelled **real-time
  (provisional)** on the site. When the delayed dataset appears, the real-time one is removed.
- A deployment with a local L3 file is not downloaded: local files win.
- Only new or changed datasets are downloaded (`ERDDAP\manifest.json`); raw downloads are
  deleted after gridding (~1–2 MB per gridded deployment).
- Checked against calanus-20260824: the gridded ERDDAP data matched the local L3 file in
  73,826 of 73,828 cells; the other 2 differ by one sample on a bin edge (< 0.001 °C).

Files under the folder that are not depth × profile grids (e.g. a raw ERDDAP download) are
skipped by the build with a note.

### Regions

`regions.json` lists named lat/lon boxes (Antarctic Peninsula, Hawaii, Northern and Southern
California, split at 36°N). Each deployment is assigned by its median position; anything
outside every box is "Other". When more than one region is present, the overview page shows
region buttons above the map (`_regionbar.md`, opening on the region with the most recent data;
`?region=<name>` or `?region=All` in the address overrides that), the map's year and project
chips list only the selected region's, the Deployments page is grouped by region, and the
viewer's picker is grouped too. Edit the boxes or names freely, then rerun the build.

### Project labels

`projects.json` maps the start of each file's `project` attribute to the short label shown on
the site, e.g. "CalCOFI Survey Augmentation" → "CalCOFI Augmentation". The full name stays on
each deployment page. Unlisted projects are shortened to the text before any colon. Edit the
list and rerun the build; the labels update without reprocessing any data.

### Bathymetry

Put one GEBCO grid per region (download at [download.gebco.net](https://download.gebco.net),
2D netCDF, ice-surface elevation) in the NetCDF folder or one folder below it, keeping
"gebco" in the file name. Each deployment uses whichever grid covers at least 95% of its
positions (the smallest such grid if several do); deployments with no covering grid are
built without seafloor and the build prints a note.

## Ocean currents on the real-time page

`scripts/fetch_currents.py` runs in the daily build. For each glider sending real-time data it
writes `currents/currents.json` (not committed): Copernicus Marine model currents (yesterday to
+5 days, every model level down to the glider's dive depth), on the US West Coast NOAA WCOFS
model currents and HF radar surface currents for the last 24 h, and the glider's own
depth-averaged currents for the last 7 days (`u`, `v` from the real-time dataset, which
`fetch_erddap.py` keeps in each real-time grid as `dac_u`, `dac_v`, ...; values that cannot be
real, such as `u` equal to `v` for every dive, are left off with a note). `realtime_map/index.html` (written by `build_site.py` from
`scripts/currents_map.py`) draws them, with day and depth sliders.

The model layer needs a free Copernicus Marine account
(<https://data.marine.copernicus.eu/register>). Add its user name and password as repository
secrets named `COPERNICUSMARINE_SERVICE_USERNAME` and `COPERNICUSMARINE_SERVICE_PASSWORD`
(Settings → Secrets and variables → Actions → New repository secret). Without them the map
shows the gliders and HF radar only, with a note.

## Size

About 120 MB of data for 16 deployments (~7.5 MB each). GitHub Pages allows published sites up to
1 GB, so roughly 130 deployments at the current resolution. Beyond that, reduce precision or host
the data files outside the repository.

## Data and credits

Glider data: NOAA SWFSC Ecosystem Science Division Glider Team, public domain (CC0-1.0),
processed with pyglider and esdglider. Bathymetry: GEBCO_2026 Grid,
doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa.
