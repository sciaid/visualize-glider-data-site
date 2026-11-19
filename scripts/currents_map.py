"""
currents_map.py - The real-time map page (realtime_map/index.html): gliders that are out now,
with model current forecasts (depth slider, day slider) and observed HF radar currents.

The page reads ../currents/currents.json, written each day by fetch_currents.py. It is a static
template; build_site.py writes it with write_map().
"""
from pathlib import Path

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Real-time gliders and currents</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  :root { --surface:#fcfcfb; --surface-2:#f3f2ef; --border:#e2e0db; --ink:#0b0b0b; --ink-2:#52514e; --muted:#7a7873; --accent:#2a78d6; }
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; font:13px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif; color:var(--ink); background:var(--surface); }
  #app { display:grid; grid-template-columns:300px 1fr; height:100%; }
  aside { border-right:1px solid var(--border); overflow:auto; padding:12px 14px; }
  #map { min-width:0; }
  h2 { font-size:11px; text-transform:uppercase; letter-spacing:.05em; color:var(--ink-2); margin:14px 0 6px; font-weight:600; }
  h2:first-child { margin-top:0; }
  .muted { color:var(--muted); font-size:12px; }
  label.line { display:flex; align-items:center; gap:8px; margin:4px 0; }
  input[type=range] { width:100%; }
  .val { font-weight:600; }
  .gl { display:block; width:100%; text-align:left; padding:5px 8px; margin:3px 0; border:1px solid var(--border);
        border-radius:6px; background:var(--surface); cursor:pointer; font:inherit; color:inherit; }
  .gl:hover { border-color:var(--accent); }
  .gl small { display:block; color:var(--muted); }
  .legend { background:var(--surface); padding:6px 9px; border:1px solid var(--border); border-radius:6px; font-size:12px; line-height:1.5; }
  .ramp { display:flex; height:8px; width:150px; border-radius:2px; overflow:hidden; margin:3px 0 1px; }
  .ramp span { flex:1; }
  .ticks { display:flex; justify-content:space-between; width:150px; color:var(--ink-2); font-size:11px; }
  .note { background:#fff4e5; border:1px solid #f3d3a6; border-radius:6px; padding:6px 8px; margin:6px 0; font-size:12px; }
  .pop td { padding:0 8px 0 0; vertical-align:top; } .pop td:first-child { color:#666; white-space:nowrap; }
</style>
</head>
<body><div id="app">
<aside>
  <h2>Gliders out now</h2>
  <div id="gliders" class="muted">Loading…</div>

  <h2>Model current (forecast)</h2>
  <label class="line"><input type="checkbox" id="mOn" checked> Show</label>
  <div>Day: <span class="val" id="dayLabel"></span></div>
  <input type="range" id="day" min="0" max="0" value="0">
  <div>Mean from the surface to <span class="val" id="depthLabel"></span></div>
  <input type="range" id="depth" min="0" max="100" value="100" step="5">
  <div class="muted" id="depthNote"></div>

  <h2>Observed current</h2>
  <label class="line"><input type="checkbox" id="hOn" checked> HF radar, mean of the last 24 h</label>
  <div class="muted" id="hfrNote"></div>

  <h2>Track</h2>
  <label class="line"><input type="checkbox" id="tOn" checked> Glider tracks</label>

  <h2>Sources</h2>
  <div class="muted" id="sources"></div>
  <div id="notes"></div>
</aside>
<div id="map"></div>
</div>
<script>
(async function () {
  // Speed colour ramp (one hue, light -> dark = slow -> fast), 0..MAXS m/s
  const RAMP = ['#cfe1f6', '#a8c9ef', '#7fb0e7', '#5598e7', '#3987e5', '#2a78d6', '#1c5cab', '#184f95', '#104281', '#0d366b'];
  const MAXS = 0.6;
  const colorFor = s => RAMP[Math.min(RAMP.length - 1, Math.floor(s / MAXS * RAMP.length))];
  const HFR_COLOR = '#d95926';

  const map = L.map('map', { worldCopyJump: true }).setView([0, -120], 2);
  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}',
    { maxZoom: 13, attribution: 'Tiles &copy; Esri, GEBCO, NOAA, NGS' }).addTo(map);
  L.control.scale({ imperial: false }).addTo(map);
  const arrowsModel = L.layerGroup().addTo(map), arrowsHfr = L.layerGroup().addTo(map), tracks = L.layerGroup().addTo(map);
  const renderer = L.canvas({ padding: 0.3 });

  let C;
  try { C = await fetch('../currents/currents.json', { cache: 'no-cache' }).then(r => { if (!r.ok) throw 0; return r.json(); }); }
  catch (e) { C = { gliders: [], notes: ['Current data are not available yet (they are made by the daily site build).'] }; }
  const G = C.gliders || [];
  const fmtT = s => new Date(s * 1000).toISOString().replace('T', ' ').slice(0, 16) + ' UTC';
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const compass = (u, v) => { const d = (Math.atan2(u, v) * 180 / Math.PI + 360) % 360;
    return ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'][Math.round(d / 22.5) % 16]; };
  const fmtCur = uv => !uv || uv[0] == null ? '—' : `${Math.hypot(uv[0], uv[1]).toFixed(0)} cm/s toward ${compass(uv[0], uv[1])}`;

  (C.notes || []).forEach(n => { const d = document.createElement('div'); d.className = 'note'; d.textContent = n; document.getElementById('notes').appendChild(d); });
  if (!G.length) {
    document.getElementById('gliders').textContent = 'No glider is sending real-time data at the moment.';
    map.setView([0, -120], 2);
  }

  // ---- day and depth controls
  const models = G.filter(g => g.model);
  const dates = models.length ? models[0].model.days.map(d => d.date) : [];
  const today = models.length ? models[0].model.today : new Date().toISOString().slice(0, 10);
  const dayEl = document.getElementById('day'), depthEl = document.getElementById('depth');
  dayEl.max = Math.max(0, dates.length - 1);
  dayEl.value = Math.max(0, dates.indexOf(today));
  const maxDepth = Math.max(100, ...models.map(g => g.model.maxDepth));
  depthEl.max = maxDepth; depthEl.value = Math.min(100, maxDepth);
  const dayText = i => {
    if (!dates.length) return 'no model data';
    const d = dates[i], n = Math.round((Date.parse(d) - Date.parse(today)) / 864e5);
    const wd = new Date(d + 'T12:00:00Z').toLocaleDateString('en-US', { weekday: 'short', month: 'short', day: 'numeric', timeZone: 'UTC' });
    return `${wd} (${n < 0 ? 'analysis, yesterday' : n === 0 ? 'today' : '+' + n + ' day' + (n > 1 ? 's' : '') + ', forecast'})`;
  };
  const depthText = D => D <= 0 ? 'the surface (top model level only)' : `${D} m`;

  // ---- depth mean on the model grid: thickness-weighted over levels down to D
  function weights(m, D) {
    if (D <= 0) return m.levels.map((_, k) => (k === 0 ? 1 : 0));
    return m.levels.map((_, k) => Math.max(0, Math.min(m.edges[k + 1], D) - Math.min(m.edges[k], D)));
  }
  function meanCol(us, vs, w, D) {           // us, vs: per-level values (cm/s, may be null)
    let su = 0, sv = 0, sw = 0;
    for (let k = 0; k < w.length; k++) if (w[k] > 0 && us[k] != null && vs[k] != null) { su += w[k] * us[k]; sv += w[k] * vs[k]; sw += w[k]; }
    const need = D <= 0 ? 1 : 0.5 * w.reduce((a, b) => a + b, 0);
    return sw >= need && sw > 0 ? [su / sw, sv / sw] : null;
  }

  // ---- arrows
  function arrow(group, lat, lon, u, v, spacing, color) {   // u, v in cm/s; spacing in degrees latitude
    const s = Math.hypot(u, v) / 100; if (!(s > 0.005)) return;
    const len = spacing * 0.9 * Math.min(1.6, s / 0.3);    // 0.3 m/s ~ one grid spacing
    const k = Math.cos(lat * Math.PI / 180);
    const dy = (v / 100 / s) * len, dx = (u / 100 / s) * len / Math.max(0.2, k);
    const tip = [lat + dy, lon + dx];
    const a = Math.atan2(dy, dx * k), h = len * 0.35;
    const left = [tip[0] - h * Math.sin(a - 0.45), tip[1] - h * Math.cos(a - 0.45) / Math.max(0.2, k)];
    const right = [tip[0] - h * Math.sin(a + 0.45), tip[1] - h * Math.cos(a + 0.45) / Math.max(0.2, k)];
    const c = color || colorFor(s);
    L.polyline([[lat, lon], tip], { color: c, weight: 1.6, opacity: 0.95, renderer, interactive: false }).addTo(group);
    L.polyline([left, tip, right], { color: c, weight: 1.6, opacity: 0.95, renderer, interactive: false }).addTo(group);
  }

  function drawModel() {
    arrowsModel.clearLayers();
    if (!document.getElementById('mOn').checked) return;
    const di = +dayEl.value, D = +depthEl.value;
    models.forEach(g => {
      const m = g.model, day = m.days[Math.min(di, m.days.length - 1)], gr = m.grid;
      const Dg = Math.min(D, m.maxDepth), w = weights(m, Dg);
      for (let r = 0; r < gr.ny; r++) for (let q = 0; q < gr.nx; q++) {
        const c = r * gr.nx + q;
        const uv = meanCol(day.u.map(lv => lv[c]), day.v.map(lv => lv[c]), w, Dg);
        if (uv) arrow(arrowsModel, gr.lat0 + r * gr.dlat, gr.lon0 + q * gr.dlon, uv[0], uv[1], Math.abs(gr.dlat));
      }
    });
    document.getElementById('dayLabel').textContent = dayText(di);
    document.getElementById('depthLabel').textContent = depthText(D);
    const capped = models.filter(g => D > g.model.maxDepth).map(g => `${g.name} (${g.model.maxDepth} m)`);
    document.getElementById('depthNote').textContent = capped.length ? `Shown to each glider's dive depth: ${capped.join(', ')}.` : '';
  }

  function drawHfr() {
    arrowsHfr.clearLayers();
    const on = document.getElementById('hOn').checked, notes = [];
    G.filter(g => g.hfr).forEach(g => {
      notes.push(`${g.name}: ${g.hfr.start} to ${g.hfr.end}`);
      if (!on) return;
      const gr = g.hfr.grid;
      for (let r = 0; r < gr.ny; r++) for (let q = 0; q < gr.nx; q++) {
        const c = r * gr.nx + q;
        if (gr.u[c] != null && gr.v[c] != null) arrow(arrowsHfr, gr.lat0 + r * gr.dlat, gr.lon0 + q * gr.dlon, gr.u[c], gr.v[c], Math.abs(gr.dlat), HFR_COLOR);
      }
    });
    document.getElementById('hfrNote').textContent = notes.length ? notes.join(' · ') :
      'None near the gliders that are out (HF radar covers the US West Coast only).';
  }

  // ---- gliders
  const list = document.getElementById('gliders');
  if (G.length) list.textContent = '';
  const bounds = L.latLngBounds([]);
  G.forEach(g => {
    const tr = g.track, recent = Date.now() / 1000 - 7 * 86400;
    const all = tr.lat.map((la, i) => [la, tr.lon[i]]);
    L.polyline(all, { color: '#7a7873', weight: 1.5, opacity: 0.7 }).addTo(tracks);
    L.polyline(all.filter((_, i) => tr.t[i] >= recent), { color: '#0b0b0b', weight: 2.5 }).addTo(tracks);
    const mk = L.circleMarker([g.lat, g.lon], { radius: 7, color: '#fff', weight: 2, fillColor: '#eb6834', fillOpacity: 1 }).addTo(map);
    mk.bindTooltip(g.name, { permanent: true, direction: 'right', offset: [8, 0] });
    mk.on('click', () => mk.bindPopup(popup(g)).openPopup());
    bounds.extend([g.lat, g.lon]);
    if (g.model) { const gr = g.model.grid; bounds.extend([gr.lat0, gr.lon0]); bounds.extend([gr.lat0 + gr.dlat * (gr.ny - 1), gr.lon0 + gr.dlon * (gr.nx - 1)]); }
    const b = document.createElement('button'); b.className = 'gl';
    const ageH = (Date.now() / 1000 - g.time) / 3600;
    b.innerHTML = `<b>${esc(g.name)}</b><small>Last fix ${fmtT(g.time)} (${ageH < 48 ? ageH.toFixed(0) + ' h' : (ageH / 24).toFixed(1) + ' days'} ago)` +
                  `${g.diveDepth ? ' · dives to ~' + g.diveDepth + ' m' : ''}</small>`;
    b.onclick = () => { map.setView([g.lat, g.lon], 9, { animate: false }); mk.bindPopup(popup(g)).openPopup(); };
    list.appendChild(b);
  });
  if (bounds.isValid()) map.fitBounds(bounds.pad(0.1));
  const wantG = new URLSearchParams(location.search).get('g');   // ?g=<deployment> zooms to that glider

  function popup(g) {
    let rows = `<tr><td>Last fix</td><td>${fmtT(g.time)}</td></tr><tr><td>Position</td><td>${g.lat.toFixed(4)}, ${g.lon.toFixed(4)}</td></tr>` +
               (g.diveDepth ? `<tr><td>Dive depth</td><td>~${g.diveDepth} m (recent dives)</td></tr>` : '');
    if (g.model) {
      const m = g.model, di = Math.min(+dayEl.value, m.days.length - 1), day = m.days[di], D = Math.min(+depthEl.value, m.maxDepth);
      const col = D2 => meanCol(day.colU, day.colV, weights(m, D2), D2);
      rows += `<tr><td colspan="2" style="padding-top:6px;color:inherit"><b>Model, ${esc(dayText(di))}</b></td></tr>` +
              `<tr><td>Surface</td><td>${fmtCur(col(0))}</td></tr>` +
              `<tr><td>0–${D} m</td><td>${fmtCur(col(D))}</td></tr>` +
              (g.diveDepth ? `<tr><td>0–${m.maxDepth} m (dive)</td><td>${fmtCur(col(m.maxDepth))}</td></tr>` : '');
    }
    if (g.hfr) {
      const gr = g.hfr.grid, r = Math.round((g.lat - gr.lat0) / gr.dlat), q = Math.round((g.lon - gr.lon0) / gr.dlon);
      const c = r >= 0 && r < gr.ny && q >= 0 && q < gr.nx ? r * gr.nx + q : -1;
      rows += `<tr><td colspan="2" style="padding-top:6px;color:inherit"><b>HF radar, last 24 h</b></td></tr>` +
              `<tr><td>Surface</td><td>${c >= 0 ? fmtCur([gr.u[c], gr.v[c]]) : '—'}</td></tr>`;
    }
    return `<div class="pop"><b>${esc(g.name)}</b><table>${rows}</table></div>`;
  }

  // ---- legend
  const legend = L.control({ position: 'bottomright' });
  legend.onAdd = () => {
    const d = L.DomUtil.create('div', 'legend');
    d.innerHTML = `<b>Model current</b><div class="ramp">${RAMP.map(c => `<span style="background:${c}"></span>`).join('')}</div>` +
      `<div class="ticks"><span>0</span><span>${(MAXS / 2).toFixed(1)}</span><span>${MAXS}+ m/s</span></div>` +
      `<div><span style="color:${HFR_COLOR};font-weight:700">→</span> HF radar (observed)</div>` +
      `<div><span style="display:inline-block;width:16px;height:3px;background:#0b0b0b;vertical-align:middle"></span> last 7 days of track</div>`;
    return d;
  };
  legend.addTo(map);

  // ---- sources
  const src = [];
  if (models.length) src.push(`Model: <a href="${models[0].model.doi}" target="_blank">${esc(models[0].model.source)}</a>, run of ${esc(models[0].model.today)}.`);
  if (G.some(g => g.hfr)) src.push(`Observed: ${esc(G.find(g => g.hfr).hfr.source)}.`);
  if (C.made) src.push(`Updated ${esc(C.made)} by the daily site build.`);
  document.getElementById('sources').innerHTML = src.join('<br>');

  dayEl.oninput = depthEl.oninput = drawModel;
  document.getElementById('mOn').onchange = drawModel;
  document.getElementById('hOn').onchange = drawHfr;
  document.getElementById('tOn').onchange = e => e.target.checked ? tracks.addTo(map) : map.removeLayer(tracks);
  if (!models.length) { document.getElementById('dayLabel').textContent = 'no model data'; document.getElementById('depthLabel').textContent = '—'; }
  drawModel(); drawHfr();
  if (wantG) { const k = G.findIndex(g => g.name === wantG); if (k >= 0) list.children[k].click(); }
})();
</script>
</body></html>
"""


def write_map(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    (out / "index.html").write_text(HTML, encoding="utf-8")
