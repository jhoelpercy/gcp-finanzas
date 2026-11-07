/* Modelo de distribución de recursos entre unidades de una organización, escenarios y análisis.
   Valor organizacional neto:  V(x) = (r − c·1)ᵀx − (γ/2) xᵀQx      (Q = matriz de costos de escala e interdependencia)
   ∇V = 0  ⇒  Q x = (r − c·1)/γ      (A = Q, b = (r − c)/γ, x = reasignación como fracción del presupuesto base) */
(function (root) {
  'use strict';
  const F = root.FinCG;

  const ASSETS = [
    { id: 'x1', name: 'Desarrollo de Software', short: 'Software' },
    { id: 'x2', name: 'Infraestructura TI', short: 'Infra TI' },
    { id: 'x3', name: 'Ciberseguridad', short: 'Ciberseg.' },
    { id: 'x4', name: 'Datos y Analítica', short: 'Datos' },
    { id: 'x5', name: 'Soporte al Cliente', short: 'Soporte' },
    { id: 'x6', name: 'Marketing Digital', short: 'Marketing' },
    { id: 'x7', name: 'Investigación e Innovación (I+D)', short: 'I+D' },
    { id: 'x8', name: 'Capacitación y Talento', short: 'Talento' }
  ];
  const GROUP = [0, 0, 1, 1, 2, 2, 3, 3]; // bloques de procesos (mayor acoplamiento intra-bloque)
  const corr = (base, bump) => F.A_from(8, (i, j) => (i === j ? 1 : base + (GROUP[i] === GROUP[j] ? bump : 0)));

  const r6 = v => Math.round(v * 1e6) / 1e6;

  /* Q_ij = s_i s_j ρ_ij: s_i = escala del costo propio de la unidad i, ρ_ij = acoplamiento entre unidades i y j
     (simétrica por construcción; sigma = s, rho = ρ en los parámetros) */
  function covariance(sigma, rho) {
    const n = sigma.length, A = F.A_from(n, () => 0);
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) { A[i][j] = A[j][i] = r6(sigma[i] * sigma[j] * (i === j ? 1 : rho[i][j])); }
    return A;
  }

  const GAMMA = 3;

  function build(params) {
    const A = covariance(params.sigma, params.rho);
    const b = params.mu.map(m => r6((m - params.rf) / GAMMA));
    return { A, b };
  }

  const P_IDEAL = { sigma: [.14, .20, .16, .17, .18, .19, .22, .24], rho: corr(.45, .25), mu: [.080, .095, .085, .075, .070, .090, .105, .110], rf: .030 };
  const P_STRESS = { sigma: [.24, .34, .27, .29, .31, .32, .37, .41], rho: corr(.78, .10), mu: [.065, .060, .060, .055, .050, .060, .070, .070], rf: .045 };
  const SCALE = [1, 2, 5, 10, 20, 50, 100, 200]; // unidades heterogéneas de las series (escenario 3)

  function scenarioIdeal() {
    const { A, b } = build(P_IDEAL);
    return { key: 'ideal', title: 'Ideal', A, b, params: P_IDEAL,
      why: 'Organización estable: costos de escala de magnitud similar entre unidades (0,14–0,24), acoplamiento interdepartamental moderado (0,45 entre procesos distintos, 0,70 dentro de un mismo bloque de procesos), retornos esperados por unidad presupuestada de 7–11 % y costo de capital del 3 %. La matriz Q está equilibrada, de modo que κ(A) es bajo.' };
  }
  function scenarioStress() {
    const { A, b } = build(P_STRESS);
    return { key: 'stress', title: 'Bajo estrés financiero', A, b, params: P_STRESS,
      why: 'Presión financiera: la inflación de costos eleva ≈ 70 % los costos de escala (0,24–0,41), la escasez de recursos compartidos aumenta el acoplamiento a 0,78 (0,88 intra-bloque), los retornos esperados caen y el costo de capital sube de 3 % a 4,5 %. Cambian A y b, pero sin cambiar la escala relativa: se desplaza el punto de equilibrio (menor reasignación neta) y κ(A) crece de forma moderada.' };
  }
  function scenarioIll() {
    const base = build(P_IDEAL);
    const A = F.A_from(8, (i, j) => r6(SCALE[i] * SCALE[j] * base.A[i][j])); // A' = S Σ S
    const Asym = F.A_from(8, (i, j) => (i <= j ? A[i][j] : A[j][i]));
    const b = base.b.map((v, i) => r6(SCALE[i] * v));
    return { key: 'ill', title: 'Mal condicionado', A: Asym, b, params: { ...P_IDEAL, scale: SCALE },
      why: 'Integración de sistemas con unidades heterogéneas: cada unidad reporta en su propia unidad contable (Bs, miles de Bs, USD, horas-hombre valorizadas, etc.). Con x′ = S⁻¹x se obtiene A′ = S Q S y b′ = S b, con S = diag(1, 2, 5, 10, 20, 50, 100, 200). Es el mismo problema del escenario Ideal, pero con un fuerte mal escalamiento diagonal: κ(A′) explota y el precondicionador de Jacobi lo corrige.' };
  }
  const defaultScenarios = () => ({ ideal: scenarioIdeal(), stress: scenarioStress(), ill: scenarioIll() });

  /* ---------- Validación de A (cuadrada, simétrica, definida positiva) ---------- */
  function validate(A, b, minN = 4) {
    const res = { ok: false, errors: [], n: A ? A.length : 0, symmetric: false, spd: false, eig: [], minors: [], cond: NaN };
    if (!Array.isArray(A) || A.length === 0) { res.errors.push('La matriz A está vacía.'); return res; }
    const n = A.length;
    if (A.some(r => !Array.isArray(r) || r.length !== n)) { res.errors.push('A debe ser cuadrada (n×n).'); return res; }
    if (A.some(r => r.some(v => typeof v !== 'number' || !isFinite(v))) || (b || []).some(v => typeof v !== 'number' || !isFinite(v))) {
      res.errors.push('A y b solo pueden contener números reales finitos.'); return res;
    }
    if (!b || b.length !== n) { res.errors.push('El vector b debe tener n componentes.'); return res; }
    if (n < minN) res.errors.push('Se requieren al menos ' + minN + ' ecuaciones e incógnitas (n ≥ ' + minN + ').');
    res.symmetric = F.isSymmetric(A);
    res.symDefect = F.symmetryDefect(A);
    if (!res.symmetric) res.errors.push('A no es simétrica: máx |aᵢⱼ − aⱼᵢ| = ' + res.symDefect.toExponential(3) + '.');
    res.minors = F.leadingMinors(A);
    if (res.symmetric) {
      res.eig = F.eigenvaluesSym(A);
      res.cond = F.cond2FromEigs(res.eig);
      res.spd = res.eig[0] > 0 && res.minors.every(m => m > 0) && F.cholesky(A) !== null;
      if (!res.spd) res.errors.push('A no es definida positiva: λ_mín = ' + res.eig[0].toExponential(3) + '.');
    }
    res.ok = res.errors.length === 0;
    return res;
  }

  /* ---------- Errores por componente ---------- */
  function errorsVs(xref, x) {
    return xref.map((r, i) => {
      const ea = Math.abs(r - x[i]);
      const defined = Math.abs(r) > 0;
      const er = defined ? ea / Math.abs(r) : null;
      return { ea, er, ep: er === null ? null : er * 100 };
    });
  }

  /* ---------- Análisis completo de un escenario ---------- */
  function analyze(scn, cfg) {
    const A = scn.A, b = scn.b, n = b.length;
    const out = { key: scn.key, validation: validate(A, b), cfg };
    if (!out.validation.ok) return out;
    // 1) Solución de referencia (método directo)
    let xref;
    if (cfg.direct === 'lu') xref = F.luSolve(F.luDecompose(A), b);
    else xref = F.choleskySolve(F.cholesky(A), b);
    out.xref = xref;
    out.refResidual = F.norm2(F.sub(b, F.matvec(A, xref))) / (F.norm2(b) || 1);
    // 2) GC simple y GCP
    const x0 = cfg.x0 && cfg.x0.length === n ? cfg.x0 : F.zeros(n);
    out.gc = F.pcg(A, b, x0, cfg.tol, cfg.maxIter, F.makePrecond('none', A));
    try {
      out.pre = F.makePrecond(cfg.pre, A, { omega: cfg.omega });
    } catch (e) {
      out.preWarning = e.message + ' Se usa Jacobi.';
      out.pre = F.makePrecond('jacobi', A);
    }
    if (out.pre.warning) out.preWarning = (out.preWarning ? out.preWarning + ' ' : '') + out.pre.warning;
    out.gcp = F.pcg(A, b, x0, cfg.tol, cfg.maxIter, out.pre);
    // 3) Errores, condicionamiento
    out.errors = errorsVs(xref, out.gcp.x);
    out.errorsGC = errorsVs(xref, out.gc.x);
    out.condA = out.validation.cond;
    out.condMA = F.condPrecond(A, out.pre.M);
    return out;
  }

  /* ---------- Texto de interpretación financiera (dinámico) ---------- */
  const pct = v => (v * 100).toFixed(1).replace('.', ',') + ' %';
  function interpretation(scn, res, names) {
    if (!res.xref) return 'Sistema no válido: no hay interpretación.';
    const x = res.xref, tot = x.reduce((a, v) => a + v, 0);
    const parts = x.map((v, i) => (names[i] || 'x' + (i + 1)) + ' ' + pct(v));
    const neg = x.some(v => v < 0);
    let t = 'Reasignación óptima (fracción del presupuesto base): ' + parts.join('; ') + '. Reasignación neta total ' + pct(tot) + ' ('
      + (tot >= 0 ? 'requiere financiar ' + pct(tot) + ' adicional del presupuesto base al costo de capital' : 'libera ' + pct(-tot) + ' del presupuesto base')
      + ')' + (neg ? '; los valores negativos son recortes a esa unidad (modelo sin restricciones de signo).' : '.') + ' ';
    if (scn.key === 'ideal') t += 'Con κ(A) ≈ ' + res.condA.toFixed(1) + ' el problema está bien planteado: pequeñas variaciones en Q o en los retornos producen cambios proporcionales en la distribución y ambos métodos iterativos convergen casi de inmediato.';
    else if (scn.key === 'stress') t += 'La presión financiera reduce la reasignación neta (retornos menores, costos de escala, acoplamiento y costo de capital mayores): la organización pasa a una posición más conservadora. El condicionamiento (κ(A) ≈ ' + res.condA.toFixed(1) + ') sigue siendo moderado: el algoritmo no se degrada, solo cambia el punto de equilibrio.';
    else t += 'Las componentes de x′ están en unidades contables distintas (x′ᵢ = xᵢ/sᵢ); al deshacer la escala se recupera la distribución del escenario Ideal. Con κ(A) ≈ ' + res.condA.toExponential(2) + ' el GC simple sufre; Jacobi reduce κ(M⁻¹A) a ≈ ' + res.condMA.toFixed(1) + ', recuperando la convergencia.';
    return t;
  }

  Object.assign(F, { ASSETS, GAMMA, defaultScenarios, covariance, validate, errorsVs, analyze, interpretation });
})(typeof window !== 'undefined' ? window : globalThis);
