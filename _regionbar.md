```{=html}
<div id="regionbar" class="regionbar" data-default="Antarctic Peninsula" hidden>
  <span class="rb-label">Region</span>
  <button type="button" data-region="Antarctic Peninsula" aria-pressed="false"><span class="rb-name">Antarctic Peninsula</span><span class="rb-meta">16 deployments · 2018–2026</span></button>
</div>
<iframe id="mapframe" title="Map of glider deployments"
        style="width:100%; height:72vh; border:1px solid #e2e0db; border-radius:6px;"></iframe>
<style>
  .regionbar { display:flex; flex-wrap:wrap; gap:8px; align-items:stretch; margin:0 0 10px; }
  .regionbar .rb-label { align-self:center; font-weight:600; margin-right:4px; }
  .regionbar button { text-align:left; border:1px solid #d6d3cc; background:#fff; border-radius:8px;
                      padding:6px 12px; cursor:pointer; line-height:1.25; font:inherit; color:inherit; }
  .regionbar button:hover { border-color:#2a78d6; }
  .regionbar button[aria-pressed="true"] { border-color:#2a78d6; background:#eaf2fc; box-shadow:inset 0 0 0 1px #2a78d6; }
  .regionbar .rb-name { display:block; font-weight:600; }
  .regionbar .rb-meta { display:block; font-size:.82em; color:#5d5b56; }
  .regionbar .rb-rt { color:#b4501a; font-weight:600; }
</style>
<script>
(function () {
  const bar = document.getElementById('regionbar'), frame = document.getElementById('mapframe');
  const buttons = [...bar.querySelectorAll('button')];
  const asked = new URLSearchParams(location.search).get('region');
  const known = buttons.map(b => b.dataset.region);
  let region = known.includes(asked) ? asked : bar.dataset.default;
  function mark(r) {
    region = r;
    buttons.forEach(b => b.setAttribute('aria-pressed', b.dataset.region === r));
    const u = new URL(location.href);
    if (r === bar.dataset.default) u.searchParams.delete('region'); else u.searchParams.set('region', r);
    history.replaceState(null, '', u);
  }
  frame.src = 'map/index.html?regionbar=' + (bar.hidden ? '0' : '1') + '&region=' + encodeURIComponent(region);
  mark(region);
  buttons.forEach(b => b.onclick = () => {
    mark(b.dataset.region);
    try { frame.contentWindow.setRegion(b.dataset.region); }
    catch (e) { frame.src = 'map/index.html?regionbar=1&region=' + encodeURIComponent(b.dataset.region); }
  });
  window.addEventListener('message', e => {           // keep the buttons in step with the map
    if (e.source === frame.contentWindow && e.data && e.data.glidermapRegion) mark(e.data.glidermapRegion);
  });
})();
</script>
```
