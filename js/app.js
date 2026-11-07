/* Interfaz: estado, entradas, renderizado reactivo. Requiere linalg.js, pcg.js, model.js, qa.js */
(function () {
  'use strict';
  const F = window.FinCG;
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const KEYS = ['ideal', 'stress', 'ill'];
  const TITLES = { ideal: 'Ideal', stress: 'Bajo estrés financiero', ill: 'Mal condicionado' };

  /* ---- formato numérico ---- */
  const fx = v => (v === null || v === undefined || Number.isNaN(v)) ? 'n/d' : (v === 0 ? '0' : (Math.abs(v) >= 1e5 || Math.abs(v) < 1e-4 ? v.toExponential(4) : v.toFixed(6)));
  const fe = v => (v === null || v === undefined || Number.isNaN(v)) ? 'no definido' : (v === 0 ? '0' : v.toExponential(3));
  const fk = v => !isFinite(v) ? '∞' : (v >= 1e4 ? v.toExponential(3) : v.toFixed(3));
  const fv = v => '[' + v.map(a => (Math.abs(a) >= 1e-3 && Math.abs(a) < 1e5 ? a.toFixed(4) : a.toExponential(2))).join(', ') + ']';

  /* ---- estado ---- */
  const state = {
    names: F.ASSETS.map(a => a.name),
    scenarios: F.defaultScenarios(),
    active: 'ideal',
    cfg: { tol: 1e-6, maxIter: 1000, pre: 'jacobi', omega: 1, direct: 'cholesky', x0: [] },
    results: {}
  };
  const n = () => state.scenarios.ideal.b.length;
  const act = () => state.scenarios[state.active];
  const res = () => state.results[state.active];

  function readCfg() {
    const c = state.cfg, tol = parseFloat($('#cfgTol').value), mx = parseInt($('#cfgMax').value, 10), om = parseFloat($('#cfgOmega').value);
    c.tol = tol > 0 ? tol : 1e-6; c.maxIter = mx >= 1 ? mx : 1000; c.pre = $('#cfgPre').value; c.omega = om > 0 && om < 2 ? om : 1; c.direct = $('#cfgDirect').value;
    $('#omegaWrap').hidden = c.pre !== 'ssor';
    c.x0 = $$('#x0Row input').map(i => { const v = parseFloat(i.value); return isFinite(v) ? v : 0; });
    if (c.x0.length !== n()) c.x0 = F.zeros(n());
  }

  /* ---- cálculo ---- */
  function recompute() {
    readCfg();
    for (const k of KEYS) {
      try { state.results[k] = F.analyze(state.scenarios[k], state.cfg); }
      catch (e) { state.results[k] = { key: k, validation: { ok: false, errors: ['Error inesperado: ' + e.message], eig: [], minors: [] } }; }
    }
    renderOutputs();
  }

  /* ================= ENTRADA ================= */
  function renderInputs() {
    const s = act(), N = n();
    let h = '<thead><tr><th></th>' + s.b.map((_, j) => '<th>x<sub>' + (j + 1) + '</sub></th>').join('') + '<th>b</th></tr></thead><tbody>';
    for (let i = 0; i < N; i++) {
      h += '<tr><th>f<sub>' + (i + 1) + '</sub></th>';
      for (let j = 0; j < N; j++) h += '<td><input data-i="' + i + '" data-j="' + j + '" type="number" step="any" aria-label="a' + (i + 1) + (j + 1) + '" value="' + s.A[i][j] + '"></td>';
      h += '<td class="bcol"><input data-b="' + i + '" type="number" step="any" aria-label="b' + (i + 1) + '" value="' + s.b[i] + '"></td></tr>';
    }
    $('#matrixInput').innerHTML = h + '</tbody>';
    const old = state.cfg.x0.length === N ? state.cfg.x0 : F.zeros(N);
    $('#x0Row').innerHTML = old.map((v, i) => '<label>x₀' + (i + 1) + '<input type="number" step="any" value="' + v + '"></label>').join('');
    $$('.chip').forEach(b => b.setAttribute('aria-pressed', b.dataset.scn === state.active));
  }

  function onCellInput(e) {
    const t = e.target; if (t.tagName !== 'INPUT') return;
    const s = act(), raw = t.value.replace(',', '.'), v = raw.trim() === '' ? NaN : parseFloat(raw);
    t.classList.toggle('bad', Number.isNaN(v));
    if (t.dataset.b !== undefined) s.b[+t.dataset.b] = v;
    else if (t.dataset.i !== undefined) {
      const i = +t.dataset.i, j = +t.dataset.j; s.A[i][j] = v;
      if ($('#mirror').checked && i !== j) {
        s.A[j][i] = v;
        const m = $('#matrixInput input[data-i="' + j + '"][data-j="' + i + '"]'); if (m) m.value = t.value;
      }
    } else return;
    recompute();
  }

  function structural(fn) { fn(); renderInputs(); recompute(); }
  function resetActive() {
    const d = F.defaultScenarios();
    if (act().b.length === d.ideal.b.length && state.names.length === d.ideal.b.length) state.scenarios[state.active] = d[state.active];
    else { state.scenarios = d; state.names = F.ASSETS.map(a => a.name); state.cfg.x0 = []; }
  }
  function addAsset() {
    if (n() >= 12) return;
    for (const k of KEYS) {
      const s = state.scenarios[k], N = s.b.length, avg = s.A.reduce((a, r, i) => a + r[i], 0) / N;
      s.A.forEach(r => r.push(0)); const row = F.zeros(N + 1); row[N] = avg; s.A.push(row);
      s.b.push(s.b.reduce((a, v) => a + v, 0) / N);
    }
    state.names.push('Unidad ' + (n())); state.cfg.x0 = [];
  }
  function delAsset() {
    if (n() <= 4) { showTransient('Se requieren al menos n = 4 unidades.'); return; }
    for (const k of KEYS) { const s = state.scenarios[k]; s.A.pop(); s.A.forEach(r => r.pop()); s.b.pop(); }
    state.names.pop(); state.cfg.x0 = [];
  }
  function showTransient(msg) { const a = $('#alerts'); a.insertAdjacentHTML('afterbegin', '<div class="alert warn transient">' + esc(msg) + '</div>'); setTimeout(() => $$('.transient', a).forEach(x => x.remove()), 2500); }
  function setActive(k) { state.active = k; renderInputs(); renderOutputs(); }

  /* ================= SALIDAS ================= */
  function renderOutputs() {
    $$('.chip').forEach(b => b.setAttribute('aria-pressed', b.dataset.scn === state.active));
    renderAlerts(); renderForm(); renderCards(); renderRef(); renderGCP(); renderErr(); renderCompare(); renderConclusions();
  }

  function renderAlerts() {
    const r = res(), box = $('#alerts'), pill = $('#statusPill');
    const keep = $$('.transient', box).map(x => x.outerHTML).join('');
    let h = '';
    if (!r.validation.ok) {
      h += '<div class="alert bad"><b>Sistema no válido — no se ejecutan los algoritmos.</b><ul>' + r.validation.errors.map(m => '<li>' + esc(m) + '</li>').join('') + '</ul>Corrige la matriz o pulsa «Restaurar escenario».</div>';
      pill.className = 'pill bad'; pill.textContent = '✗ A no válida';
    } else {
      pill.className = 'pill ok'; pill.textContent = '✓ A simétrica y definida positiva';
      if (r.preWarning) h += '<div class="alert warn">' + esc(r.preWarning) + '</div>';
      if (!r.gcp.converged) h += '<div class="alert warn">GCP no convergió: ' + esc(r.gcp.note || '') + '</div>';
    }
    box.innerHTML = keep + h;
  }

  function invalidMsg() { return '<div class="alert bad">El escenario activo no es válido (ver sección 4). Corrige los datos para ver resultados.</div>'; }

  function renderForm() {
    const N = n(), s = act(), r = res();
    let h = '<h3>Variables</h3><div class="scroll"><table class="tbl"><tr><th>Variable</th><th>Unidad (editable)</th><th>Significado financiero</th></tr>';
    for (let i = 0; i < N; i++) h += '<tr><td>x<sub>' + (i + 1) + '</sub></td><td><input class="nameIn" data-i="' + i + '" value="' + esc(state.names[i]) + '" aria-label="Nombre unidad ' + (i + 1) + '"></td><td>Reasignación presupuestaria a esta unidad, como fracción del presupuesto base (positiva = refuerzo, negativa = recorte)</td></tr>';
    $('#formVars').innerHTML = h + '</table></div>' +
      '<h3>Coeficientes</h3><ul><li><b>a<sub>ii</sub></b> = s<sub>i</sub>²: costo de escala de la unidad <i>i</i> (rendimientos decrecientes: cuanto más presupuesto recibe, menos valor marginal aporta cada unidad adicional).</li><li><b>a<sub>ij</sub></b> = s<sub>i</sub>s<sub>j</sub>ρ<sub>ij</sub>: interdependencia entre las unidades <i>i</i> y <i>j</i> (recursos, personal y sistemas compartidos); ρ<sub>ij</sub> alto significa que reforzar una unidad encarece reforzar la otra.</li><li><b>b<sub>i</sub></b> = (r<sub>i</sub> − c)/γ: valor neto marginal de la unidad <i>i</i>, es decir, su retorno esperado por unidad presupuestada menos el costo de capital c, ponderado por el factor de curvatura γ = ' + F.GAMMA + '.</li></ul>';
    // sistema
    let t = '<h3>Sistema A x = b — escenario «' + TITLES[state.active] + '»</h3><p class="muted">' + esc(s.why) + '</p><div class="scroll"><table class="tbl"><tr><th></th>' + s.b.map((_, j) => '<th class="num">x' + (j + 1) + '</th>').join('') + '<th>=</th><th class="num">b</th></tr>';
    for (let i = 0; i < N; i++) t += '<tr><th>f' + (i + 1) + '</th>' + s.A[i].map(v => '<td class="num">' + (Number.isFinite(v) ? v.toPrecision(5) : 'NaN') + '</td>').join('') + '<td>=</td><td class="num">' + (Number.isFinite(s.b[i]) ? s.b[i].toPrecision(5) : 'NaN') + '</td></tr>';
    $('#formSystem').innerHTML = t + '</table></div>';
    // SPD
    let p = '<h3>Verificación formal: A = Aᵀ y A definida positiva</h3><p>Se comprueban tres criterios equivalentes: (i) simetría <code>max|aᵢⱼ − aⱼᵢ|</code>; (ii) autovalores λ &gt; 0; (iii) Criterio de Sylvester: todos los menores principales líderes Δ<sub>k</sub> &gt; 0; y (iv) existe la factorización de Cholesky.</p><div class="scroll"><table class="tbl"><tr><th>Escenario</th><th>Simetría (defecto)</th><th class="num">λ<sub>mín</sub></th><th class="num">λ<sub>máx</sub></th><th>Δ<sub>k</sub> &gt; 0 (k = 1..n)</th><th>Cholesky</th><th>¿SPD?</th></tr>';
    for (const k of KEYS) {
      const v = state.results[k].validation;
      const dm = v.minors && v.minors.length ? Math.min(...v.minors) : NaN;
      p += '<tr><td>' + TITLES[k] + '</td><td>' + (v.symmetric ? '✓ ' : '✗ ') + (v.symDefect !== undefined ? fe(v.symDefect) : '—') + '</td><td class="num">' + (v.eig.length ? fe(v.eig[0]) : '—') + '</td><td class="num">' + (v.eig.length ? fe(v.eig[v.eig.length - 1]) : '—') + '</td><td>' + (v.minors && v.minors.length ? (v.minors.every(m => m > 0) ? '✓ todos (mín Δ = ' + fe(dm) + ')' : '✗ alguno ≤ 0') : '—') + '</td><td>' + (v.symmetric && v.eig.length && F.cholesky(state.scenarios[k].A) ? '✓' : '✗') + '</td><td class="' + (v.ok ? 'ok' : 'bad') + '">' + (v.ok ? 'Sí' : 'No') + '</td></tr>';
    }
    p += '</table></div>';
    if (r.validation.eig.length) p += '<details class="card"><summary>Detalle del escenario activo: autovalores y menores principales</summary><p><b>Autovalores</b> (ascendentes): ' + r.validation.eig.map(fe).join(', ') + '</p><p><b>Menores Δ<sub>1</sub>…Δ<sub>n</sub></b>: ' + r.validation.minors.map(fe).join(', ') + '</p></details>';
    $('#formSPD').innerHTML = p;
  }

  function kpi(label, val) { return '<div class="kpi"><b>' + val + '</b><span>' + label + '</span></div>'; }

  function renderCards() {
    $('#scnCards').innerHTML = KEYS.map(k => {
      const r = state.results[k], ok = r.validation.ok;
      return '<div class="card scncard' + (k === state.active ? ' active' : '') + '"><h3>' + TITLES[k] + '<span class="' + (ok ? 'ok' : 'bad') + '">' + (ok ? '✓ válido' : '✗ inválido') + '</span></h3><p>' + esc(state.scenarios[k].why) + '</p>' +
        (ok ? '<div class="kpis">' + kpi('κ(A)', fk(r.condA)) + kpi('κ(M⁻¹A)', fk(r.condMA)) + kpi('Iter. GC', r.gc.iters) + kpi('Iter. GCP', r.gcp.iters) + '</div>' : '') +
        '<button class="btn primary" data-load="' + k + '">Cargar y ejecutar</button></div>';
    }).join('');
  }

  function renderRef() {
    const r = res(), N = n();
    if (!r.validation.ok) { $('#refOut').innerHTML = invalidMsg(); return; }
    let h = '<p>Método directo: <b>' + (state.cfg.direct === 'lu' ? 'Eliminación LU con pivoteo parcial (PA = LU)' : 'Factorización de Cholesky (A = L Lᵀ)') + '</b>. Se considera la <b>solución numérica de referencia</b> (puede clasificarse como solución exacta salvo errores de redondeo).</p>';
    h += '<div class="kpis">' + kpi('Residuo relativo ‖b−Ax_ref‖/‖b‖', fe(r.refResidual)) + kpi('Reasignación neta total Σxᵢ',fx(r.xref.reduce((a, v) => a + v, 0))) + '</div>';
    h += '<h3>Escenario «' + TITLES[state.active] + '»</h3><div class="scroll"><table class="tbl"><tr><th>Variable</th><th>Unidad</th><th class="num">x_ref</th></tr>' + r.xref.map((v, i) => '<tr><td>x<sub>' + (i + 1) + '</sub></td><td>' + esc(state.names[i]) + '</td><td class="num">' + fx(v) + '</td></tr>').join('') + '</table></div>';
    h += '<h3>Soluciones de referencia de los tres escenarios</h3><div class="scroll"><table class="tbl"><tr><th>Variable</th>' + KEYS.map(k => '<th class="num">' + TITLES[k] + '</th>').join('') + '</tr>';
    for (let i = 0; i < N; i++) h += '<tr><td>x<sub>' + (i + 1) + '</sub> · ' + esc(state.names[i]) + '</td>' + KEYS.map(k => '<td class="num">' + (state.results[k].xref ? fx(state.results[k].xref[i]) : '—') + '</td>').join('') + '</tr>';
    $('#refOut').innerHTML = h + '</table></div><p class="muted">En el escenario Mal condicionado las componentes están en unidades reescaladas: x′ᵢ = xᵢ/sᵢ con S = diag(1, 2, 5, 10, 20, 50, 100, 200).</p>';
  }

  function chart(series, tol) {
    const W = 680, H = 300, L = 54, R = 14, T = 12, B = 34, floor = 1e-18;
    const maxK = Math.max(...series.map(s => s.v.length - 1), 1);
    const ys = series.flatMap(s => s.v.map(y => Math.log10(Math.max(y, floor))));
    let lo = Math.floor(Math.min(...ys, Math.log10(tol)) - 0.5), hi = Math.ceil(Math.max(...ys, 0) + 0.2);
    if (hi <= lo) hi = lo + 1;
    const X = k => L + (W - L - R) * k / maxK, Y = y => T + (H - T - B) * (hi - Math.log10(Math.max(y, floor))) / (hi - lo);
    let g = '';
    for (let p = Math.ceil(lo); p <= hi; p += Math.max(1, Math.round((hi - lo) / 8))) g += '<line class="grid" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(10 ** p) + '" y2="' + Y(10 ** p) + '"/><text x="' + (L - 6) + '" y="' + (Y(10 ** p) + 4) + '" text-anchor="end">1e' + p + '</text>';
    const step = Math.max(1, Math.ceil(maxK / 10));
    for (let k = 0; k <= maxK; k += step) g += '<text x="' + X(k) + '" y="' + (H - 14) + '" text-anchor="middle">' + k + '</text>';
    g += '<text x="' + (W / 2) + '" y="' + (H - 1) + '" text-anchor="middle">iteración k</text>';
    g += '<line x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(tol) + '" y2="' + Y(tol) + '" stroke="var(--bad)" stroke-dasharray="5 4"/><text x="' + (W - R) + '" y="' + (Y(tol) - 4) + '" text-anchor="end" style="fill:var(--bad)">tol</text>';
    const lines = series.map(s => '<polyline fill="none" stroke="' + s.color + '" stroke-width="2.2" points="' + s.v.map((y, k) => X(k) + ',' + Y(y)).join(' ') + '"/>' + s.v.map((y, k) => '<circle cx="' + X(k) + '" cy="' + Y(y) + '" r="3.2" fill="' + s.color + '"><title>' + s.name + ' k=' + k + ': ' + y.toExponential(2) + '</title></circle>').join('')).join('');
    return '<svg class="chart" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Residuo relativo por iteración (escala logarítmica)">' + g + lines + '</svg>';
  }

  function renderGCP() {
    const r = res(), N = n();
    if (!r.validation.ok) { $('#gcpOut').innerHTML = invalidMsg(); return; }
    const conv = x => x ? '<span class="ok">Sí</span>' : '<span class="bad">No</span>';
    let h = '<p>Precondicionador: <b>' + esc(r.pre.label) + '</b> · x₀ = ' + fv(state.cfg.x0) + ' · criterio de parada ‖r<sub>k</sub>‖₂/‖b‖₂ &lt; ' + state.cfg.tol + ' · N<sub>máx</sub> = ' + state.cfg.maxIter + '.</p>';
    h += '<div class="kpis">' + kpi('Iteraciones GC (sin precond.)', r.gc.iters) + kpi('Iteraciones GCP', r.gcp.iters) + kpi('Residuo final GC', fe(r.gc.residual)) + kpi('Residuo final GCP', fe(r.gcp.residual)) + kpi('κ(A)', fk(r.condA)) + kpi('κ(M⁻¹A)', fk(r.condMA)) + '</div>';
    h += '<div class="scroll"><table class="tbl"><tr><th>Diagnóstico</th><th>GC simple</th><th>GCP</th></tr><tr><td>¿Convergió?</td><td>' + conv(r.gc.converged) + '</td><td>' + conv(r.gcp.converged) + '</td></tr><tr><td>Observaciones</td><td>' + esc(r.gc.note || 'Tolerancia alcanzada') + '</td><td>' + esc(r.gcp.note || 'Tolerancia alcanzada') + '</td></tr></table></div>';
    h += '<h3>Vectores aproximados</h3><div class="scroll"><table class="tbl"><tr><th>Variable</th><th class="num">x_ref</th><th class="num">x GC</th><th class="num">x GCP</th></tr>';
    for (let i = 0; i < N; i++) h += '<tr><td>x<sub>' + (i + 1) + '</sub> · ' + esc(state.names[i]) + '</td><td class="num">' + fx(r.xref[i]) + '</td><td class="num">' + fx(r.gc.x[i]) + '</td><td class="num">' + fx(r.gcp.x[i]) + '</td></tr>';
    h += '</table></div><h3>Convergencia</h3><div class="legend"><span><i style="background:var(--c-gc)"></i>GC simple</span><span><i style="background:var(--c-gcp)"></i>GCP</span><span><i style="background:var(--bad)"></i>tolerancia</span></div>';
    h += chart([{ name: 'GC', color: 'var(--c-gc)', v: r.gc.history }, { name: 'GCP', color: 'var(--c-gcp)', v: r.gcp.history }], state.cfg.tol);
    const m = Math.max(r.gc.history.length, r.gcp.history.length);
    h += '<details class="card"><summary>Historial de residuos por iteración</summary><div class="scroll"><table class="tbl"><tr><th>k</th><th class="num">‖r‖/‖b‖ GC</th><th class="num">‖r‖/‖b‖ GCP</th></tr>';
    for (let k = 0; k < m; k++) h += '<tr><td>' + k + '</td><td class="num">' + (r.gc.history[k] !== undefined ? fe(r.gc.history[k]) : '—') + '</td><td class="num">' + (r.gcp.history[k] !== undefined ? fe(r.gcp.history[k]) : '—') + '</td></tr>';
    $('#gcpOut').innerHTML = h + '</table></div></details>' +
      '<p class="muted">Algoritmo (Suñagua, 2020): r₀ = b − Ax₀; resolver Mz₀ = r₀; p₀ = z₀. Para k = 0,1,…: α<sub>k</sub> = r<sub>k</sub>ᵀz<sub>k</sub> / p<sub>k</sub>ᵀAp<sub>k</sub>; x<sub>k+1</sub> = x<sub>k</sub> + α<sub>k</sub>p<sub>k</sub>; r<sub>k+1</sub> = r<sub>k</sub> − α<sub>k</sub>Ap<sub>k</sub>; resolver Mz<sub>k+1</sub> = r<sub>k+1</sub>; β<sub>k</sub> = r<sub>k+1</sub>ᵀz<sub>k+1</sub> / r<sub>k</sub>ᵀz<sub>k</sub>; p<sub>k+1</sub> = z<sub>k+1</sub> + β<sub>k</sub>p<sub>k</sub>. El GC simple es el mismo algoritmo con M = I. El residuo del gráfico es el recursivo; el «residuo final» es el verdadero ‖b − Ax‖/‖b‖.</p>';
  }

  function renderErr() {
    const r = res(), N = n();
    if (!r.validation.ok) { $('#errOut').innerHTML = invalidMsg(); return; }
    const which = ($('#errWhich') || {}).value || 'gcp';
    const E = which === 'gc' ? r.errorsGC : r.errors, X = which === 'gc' ? r.gc : r.gcp;
    let h = '<label class="inline">Aproximación analizada: <select id="errWhich"><option value="gcp"' + (which === 'gcp' ? ' selected' : '') + '>GCP</option><option value="gc"' + (which === 'gc' ? ' selected' : '') + '>GC simple</option></select></label>';
    h += '<div class="kpis">' + kpi('Residuo relativo global ‖b−Ax‖/‖b‖', fe(X.residual)) + kpi('Máx. error absoluto', fe(Math.max(...E.map(e => e.ea)))) + kpi('κ(A)', fk(r.condA)) + kpi('κ(M⁻¹A)', fk(r.condMA)) + '</div>';
    h += '<div class="scroll"><table class="tbl"><tr><th>Variable</th><th class="num">x_ref</th><th class="num">x aprox.</th><th class="num">E<sub>a</sub> = |x_ref − x|</th><th class="num">E<sub>r</sub> = E<sub>a</sub>/|x_ref|</th><th class="num">E<sub>p</sub> = E<sub>r</sub>·100 %</th></tr>';
    for (let i = 0; i < N; i++) h += '<tr><td>x<sub>' + (i + 1) + '</sub> · ' + esc(state.names[i]) + '</td><td class="num">' + fx(r.xref[i]) + '</td><td class="num">' + fx(X.x[i]) + '</td><td class="num">' + fe(E[i].ea) + '</td><td class="num">' + fe(E[i].er) + '</td><td class="num">' + (E[i].ep === null ? 'no definido' : fe(E[i].ep) + ' %') + '</td></tr>';
    $('#errOut').innerHTML = h + '</table></div><p class="muted">Si x_ref,i = 0 el error relativo y porcentual no están definidos para esa componente.</p>';
    $('#errWhich').addEventListener('change', renderErr);
  }

  function renderCompare() {
    const rows = [
      ['Solución numérica de referencia', r => r.xref ? fv(r.xref) : '—'],
      ['Solución GCP', r => r.gcp ? fv(r.gcp.x) : '—'],
      ['Error absoluto (E<sub>a</sub>) por variable', r => r.errors ? fv(r.errors.map(e => e.ea)) : '—'],
      ['Error porcentual (E<sub>p</sub>) por variable (%)', r => r.errors ? '[' + r.errors.map(e => e.ep === null ? 'n/d' : e.ep.toExponential(1)).join(', ') + ']' : '—'],
      ['Iteraciones GC (sin precondicionar)', r => r.gc ? r.gc.iters : '—'],
      ['Iteraciones GCP (precondicionado)', r => r.gcp ? r.gcp.iters : '—'],
      ['Número de condición κ(A)', r => r.gc ? fk(r.condA) : '—'],
      ['Número de condición κ(M⁻¹A)', r => r.gc ? fk(r.condMA) : '—'],
      ['Residuo relativo final GCP', r => r.gcp ? fe(r.gcp.residual) : '—'],
      ['¿Convergió el algoritmo?', r => r.gcp ? (r.gcp.converged ? 'Sí' : 'No') : '—'],
      ['Interpretación financiera', (r, k) => r.xref ? esc(F.interpretation(state.scenarios[k], r, state.names.map(nm => nm))) : 'Sistema no válido.']
    ];
    let h = '<table class="tbl"><tr><th>Indicador</th>' + KEYS.map(k => '<th>' + TITLES[k] + '</th>').join('') + '</tr>';
    for (const [lab, fn] of rows) h += '<tr><td>' + lab + '</td>' + KEYS.map(k => '<td' + (lab.startsWith('Interp') ? '' : ' class="num"') + '>' + fn(state.results[k], k) + '</td>').join('') + '</tr>';
    $('#compareOut').innerHTML = h + '</table>';
  }

  function renderConclusions() {
    const c = {}; // resultados canónicos: datos originales, Jacobi, tol 1e-6
    const D = F.defaultScenarios();
    for (const k of KEYS) c[k] = F.analyze(D[k], { tol: 1e-6, maxIter: 1000, pre: 'jacobi', omega: 1, direct: 'cholesky' });
    const tot = k => c[k].xref.reduce((a, v) => a + v, 0);
    $('#concl').innerHTML =
      '<p><b>Interpretación matemática.</b> En los tres escenarios A es simétrica y definida positiva, por lo que GC/GCP están bien definidos. Con los datos originales (n = 8, tol = 10⁻⁶) el GCP-Jacobi converge en <b>' + c.ideal.gcp.iters + ', ' + c.stress.gcp.iters + ' y ' + c.ill.gcp.iters + '</b> iteraciones (Ideal, Estrés, Mal cond.), con residuos relativos ≤ ' + Math.max(...KEYS.map(k => c[k].gcp.residual)).toExponential(1) + ' y errores absolutos ≤ ' + Math.max(...KEYS.map(k => Math.max(...c[k].errors.map(e => e.ea)))).toExponential(1) + ' respecto de la solución directa.</p>' +
      '<p><b>Condicionamiento.</b> κ(A) pasa de ' + fk(c.ideal.condA) + ' (Ideal) a ' + fk(c.stress.condA) + ' (Estrés) y a ' + fk(c.ill.condA) + ' (Mal cond.). El estrés encarece el sistema moderadamente, pero lo esencial es que <i>desplaza el punto de equilibrio</i> (reasignación neta total de ' + (tot('ideal') * 100).toFixed(0) + ' % a ' + (tot('stress') * 100).toFixed(0) + ' % del presupuesto base) sin degradar el algoritmo. En el escenario mal escalado, el GC simple necesita <b>' + c.ill.gc.iters + '</b> iteraciones frente a <b>' + c.ill.gcp.iters + '</b> del GCP: Jacobi reduce κ de ' + fk(c.ill.condA) + ' a ' + fk(c.ill.condMA) + ' (el mismo valor que en el escenario Ideal, porque el escalamiento diagonal se elimina por completo).</p>' +
      '<p><b>Interpretación financiera y organizacional.</b> La distribución óptima refuerza las unidades con mayor valor neto marginal y recorta las que comparten recursos con unidades ya reforzadas (interdependencia alta). Bajo presión financiera (mayor costo de capital, menores retornos y más acoplamiento) la organización reduce la reasignación neta y se vuelve más conservadora. Un mal escalamiento de las unidades contables no cambia la economía del problema (x′ = S⁻¹x), pero sí puede arruinar la convergencia o la estabilidad numérica de un solver sin precondicionar: de ahí la importancia de normalizar los datos al integrar sistemas.</p>' +
      '<p><b>Limitaciones.</b> (1) Los parámetros son ilustrativos, no tomados de los registros contables de una organización real. (2) Con n = 8 el GC termina teóricamente en ≤ 8 pasos, por lo que las diferencias en iteraciones son modestas en valor absoluto; la ventaja crece con n y κ. (3) El modelo no impone restricciones de signo ni de presupuesto total (permite recortes y financiamiento adicional). (4) Con A densa, IC(0) coincide con Cholesky completo y converge en una iteración: solo ilustra el límite ideal del precondicionamiento.</p>';
  }

  /* ================= QA ================= */
  function runQA() {
    const env = {
      simulateEdit() {
        const prev = state.active, prevMirror = $('#mirror').checked;
        try {
          setActive('ideal'); recompute();
          const before = state.results.ideal.xref.slice(), cellRef = $('#refOut td.num').textContent;
          const a11 = $('#matrixInput input[data-i="0"][data-j="0"]'), b2 = $('#matrixInput input[data-b="1"]');
          a11.value = String(parseFloat(a11.value) * 1.5); a11.dispatchEvent(new Event('input', { bubbles: true }));
          b2.value = String(parseFloat(b2.value) * 1.2); b2.dispatchEvent(new Event('input', { bubbles: true }));
          const r = state.results.ideal, cellNew = $('#refOut td.num').textContent;
          const dx = F.norm2(F.sub(r.xref, before)), s = state.scenarios.ideal;
          const res0 = F.norm2(F.sub(s.b, F.matvec(s.A, r.xref))) / F.norm2(s.b);
          const maxEa = Math.max(...r.errors.map(e => e.ea));
          const pass = r.validation.ok && dx > 1e-6 && cellNew !== cellRef && res0 < 1e-12 && maxEa < 1e-5;
          return { pass, got: 'Δx_ref = ' + fe(dx) + '; tabla de referencia actualizada: ' + (cellNew !== cellRef ? 'sí' : 'no') + '; ‖b−Ax_ref‖/‖b‖ = ' + fe(res0) + '; máx E_a(GCP) = ' + fe(maxEa) };
        } finally { state.scenarios.ideal = F.defaultScenarios().ideal; $('#mirror').checked = prevMirror; setActive(prev); recompute(); }
      }
    };
    const rows = F.QA.run(env);
    $('#qaTable').innerHTML = '<tr><th>ID · Prueba</th><th>Entrada</th><th>Resultado esperado</th><th>Resultado obtenido</th><th>Conclusión</th></tr>' +
      rows.map(r => '<tr><td>' + r.id + '<br>' + esc(r.name) + '</td><td>' + esc(r.input) + '</td><td>' + esc(r.expected) + '</td><td>' + esc(r.got) + '</td><td class="' + (r.pass ? 'res-pass' : 'res-fail') + '">' + (r.pass ? '✓ ' : '✗ ') + esc(r.conclusion) + '</td></tr>').join('');
    const ok = rows.filter(r => r.pass).length;
    $('#qaSummary').textContent = ok + ' de ' + rows.length + ' pruebas aprobadas';
    $('#qaSummary').className = ok === rows.length ? 'ok' : 'bad';
  }

  /* ================= INIT ================= */
  function init() {
    $('#year').textContent = new Date().getFullYear();
    try { const st = JSON.parse(localStorage.getItem('gcp-student') || '{}'); ['stName', 'stCareer', 'stMention'].forEach(id => { if (st[id]) $('#' + id).value = st[id]; }); } catch (e) { /* sin almacenamiento */ }
    $$('.student input').forEach(i => i.addEventListener('input', () => { try { const o = {}; $$('.student input').forEach(x => o[x.id] = x.value); localStorage.setItem('gcp-student', JSON.stringify(o)); } catch (e) { /* ignorar */ } }));
    $('#matrixInput').addEventListener('input', onCellInput);
    $('#x0Row').addEventListener('input', recompute);
    ['cfgTol', 'cfgMax', 'cfgPre', 'cfgOmega', 'cfgDirect'].forEach(id => $('#' + id).addEventListener('input', recompute));
    $('#btnReset').onclick = () => structural(resetActive);
    $('#btnAdd').onclick = () => structural(addAsset);
    $('#btnDel').onclick = () => structural(delAsset);
    $('#btnX0').onclick = () => { $$('#x0Row input').forEach(i => i.value = 0); recompute(); };
    $('#btnBadSym').onclick = () => structural(() => { act().A[0][1] += 0.05 * Math.max(1e-3, Math.abs(act().A[0][1]) + 1e-3); });
    $('#btnBadPD').onclick = () => structural(() => { act().A[0][0] = -Math.abs(act().A[0][0]) - 0.01; });
    $$('.chip').forEach(b => b.addEventListener('click', () => setActive(b.dataset.scn)));
    $('#scnCards').addEventListener('click', e => { const k = e.target.dataset.load; if (k) setActive(k); });
    $('#formVars').addEventListener('input', e => { if (e.target.classList.contains('nameIn')) { state.names[+e.target.dataset.i] = e.target.value; renderRef(); renderGCP(); renderErr(); renderCompare(); } });
    $('#btnQA').onclick = runQA;
    renderInputs(); recompute();
    window.__gcpState = state; // útil para depurar en consola
  }
  document.addEventListener('DOMContentLoaded', init);
})();
