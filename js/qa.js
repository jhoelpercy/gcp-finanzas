/* Módulo QA: pruebas ejecutables (Entrada, Resultado esperado, Resultado obtenido, Conclusión).
   Funciona igual en el navegador y en Node (tests/run-qa.js). */
(function (root) {
  'use strict';
  const F = root.FinCG;
  const CFG = { tol: 1e-6, maxIter: 1000, pre: 'jacobi', omega: 1, direct: 'cholesky' };
  const e = v => (typeof v === 'number' ? v.toExponential(2) : String(v));
  const T = (id, name, input, expected, fn) => ({ id, name, input, expected, fn });

  function tests(env = {}) {
    return [
      T('QA-1', 'Carga de escenarios',
        'Cargar Ideal, Bajo estrés y Mal condicionado (n = 8, tol = 1e-6, Jacobi).',
        'Los 3 escenarios son válidos (A simétrica y definida positiva) y el GCP converge con residuo relativo < 1e-6.',
        () => {
          const S = F.defaultScenarios(), rows = [];
          let ok = true;
          for (const k of ['ideal', 'stress', 'ill']) {
            const r = F.analyze(S[k], CFG);
            const good = r.validation.ok && r.gcp.converged && r.gcp.residual < CFG.tol;
            ok = ok && good;
            rows.push(S[k].title + ': ' + (good ? 'OK' : 'FALLÓ') + ' (GCP ' + r.gcp.iters + ' it., res ' + e(r.gcp.residual) + ')');
          }
          return { pass: ok, got: rows.join(' · ') };
        }),
      T('QA-2a', 'Validación: matriz no simétrica',
        'A del escenario Ideal con a₁₂ modificado (a₁₂ ≠ a₂₁).',
        'El sistema captura el error, informa que A no es simétrica y no ejecuta los algoritmos (sin excepciones).',
        () => {
          const s = F.defaultScenarios().ideal; s.A[0][1] += 0.01;
          let r, threw = false;
          try { r = F.analyze(s, CFG); } catch (err) { threw = true; }
          const pass = !threw && r.validation.ok === false && !r.validation.symmetric && !r.gcp;
          return { pass, got: threw ? 'Lanzó excepción' : 'ok=' + r.validation.ok + '; ' + r.validation.errors.join(' ') };
        }),
      T('QA-2b', 'Validación: matriz simétrica no definida positiva',
        'A del escenario Ideal con a₁₁ = −0,05 (autovalor negativo).',
        'El sistema detecta que A no es definida positiva, notifica al usuario y no colapsa.',
        () => {
          const s = F.defaultScenarios().ideal; s.A[0][0] = -0.05;
          let r, threw = false;
          try { r = F.analyze(s, CFG); } catch (err) { threw = true; }
          const pass = !threw && r.validation.ok === false && r.validation.symmetric && !r.validation.spd && !r.gcp;
          return { pass, got: threw ? 'Lanzó excepción' : 'ok=' + r.validation.ok + '; ' + r.validation.errors.join(' ') };
        }),
      T('QA-2c', 'Validación: entrada no numérica / dimensiones',
        'b con un valor NaN y, aparte, una matriz de 3×3 (n < 4).',
        'Ambas entradas se rechazan con un mensaje claro.',
        () => {
          const s = F.defaultScenarios().ideal; s.b[2] = NaN;
          const r1 = F.validate(s.A, s.b);
          const A3 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]], r2 = F.validate(A3, [1, 1, 1]);
          const pass = !r1.ok && !r2.ok && r2.errors.some(m => /n ≥ 4/.test(m));
          return { pass, got: 'NaN → ' + r1.errors[0] + ' · 3×3 → ' + r2.errors[0] };
        }),
      T('QA-3', 'Edición dinámica de coeficientes',
        env.simulateEdit ? 'Se escribe un nuevo valor en la celda a₁₁ y en b₂ de la interfaz (evento real "input").' : 'a₁₁ ← 1,5·a₁₁ y b₂ ← 1,2·b₂ (simétrico) sobre el escenario Ideal.',
        'La solución de referencia, el GCP y los errores se recalculan al instante, cambian respecto de la anterior y cumplen A·x = b.',
        () => {
          if (env.simulateEdit) return env.simulateEdit();
          const s = F.defaultScenarios().ideal, before = F.analyze(s, CFG);
          s.A[0][0] *= 1.5; s.b[1] *= 1.2;
          const after = F.analyze(s, CFG);
          const moved = F.norm2(F.sub(after.xref, before.xref)) > 1e-6;
          const res = F.norm2(F.sub(s.b, F.matvec(s.A, after.xref))) / F.norm2(s.b);
          const gcpOk = after.errors.every(x => x.ea < 1e-5);
          return { pass: after.validation.ok && moved && res < 1e-12 && gcpOk, got: 'Δx_ref = ' + e(F.norm2(F.sub(after.xref, before.xref))) + '; ‖b−Ax_ref‖/‖b‖ = ' + e(res) + '; máx Ea(GCP) = ' + e(Math.max(...after.errors.map(x => x.ea))) };
        }),
      T('QA-4', 'Comparativa de precondicionamiento',
        'Escenario Mal condicionado, tol = 1e-6, GC simple vs GCP-Jacobi.',
        'Iteraciones(GCP) ≤ iteraciones(GC) y κ(M⁻¹A) < κ(A).',
        () => {
          const r = F.analyze(F.defaultScenarios().ill, CFG);
          const pass = r.gc.converged && r.gcp.converged && r.gcp.iters <= r.gc.iters && r.condMA < r.condA;
          return { pass, got: 'GC = ' + r.gc.iters + ' it.; GCP = ' + r.gcp.iters + ' it.; κ(A) = ' + e(r.condA) + '; κ(M⁻¹A) = ' + r.condMA.toFixed(2) };
        }),
      T('QA-5', 'Extra: exactitud del método directo (solución conocida)',
        'x* = (1, −2, 3, …) y b = A·x* en el escenario Ideal; se resuelve con Cholesky y con LU.',
        'Ambos métodos directos recuperan x* con error absoluto < 1e-10.',
        () => {
          const s = F.defaultScenarios().ideal, xs = s.b.map((_, i) => (i % 2 ? -1 : 1) * (i + 1));
          const b = F.matvec(s.A, xs);
          const a = F.choleskySolve(F.cholesky(s.A), b), l = F.luSolve(F.luDecompose(s.A), b);
          const ea = Math.max(...F.errorsVs(xs, a).map(x => x.ea)), el = Math.max(...F.errorsVs(xs, l).map(x => x.ea));
          return { pass: ea < 1e-10 && el < 1e-10, got: 'máx Ea Cholesky = ' + e(ea) + '; máx Ea LU = ' + e(el) };
        }),
      T('QA-6', 'Extra: criterio de parada y N_máx',
        'Escenario Mal condicionado con N_máx = 2 (insuficiente) y con N_máx = 1000.',
        'Con N_máx = 2 reporta «no convergió»; con 1000 converge. Nunca queda en bucle infinito.',
        () => {
          const s = F.defaultScenarios().ill;
          const a = F.analyze(s, { ...CFG, maxIter: 2, pre: 'none' }), b = F.analyze(s, CFG);
          return { pass: !a.gcp.converged && a.gcp.iters === 2 && b.gcp.converged, got: 'N_máx=2 → converged=' + a.gcp.converged + ' (' + a.gcp.iters + ' it.); N_máx=1000 → converged=' + b.gcp.converged };
        }),
      T('QA-7', 'Extra: error relativo con componente nula',
        'x_ref = (0, 2, 4, 5) y x_GCP = (0.1, 2, 4, 5).',
        'Para la componente con x_ref = 0 el error relativo es «no definido» (null) y no aparece Infinity/NaN.',
        () => {
          const er = F.errorsVs([0, 2, 4, 5], [0.1, 2, 4, 5]);
          return { pass: er[0].er === null && er[0].ep === null && Math.abs(er[0].ea - 0.1) < 1e-15 && er[1].er === 0, got: 'Er₁ = ' + er[0].er + ', Ea₁ = ' + er[0].ea + ', Er₂ = ' + er[1].er };
        })
    ];
  }

  function run(env = {}) {
    return tests(env).map(t => {
      let r;
      try { r = t.fn(); } catch (err) { r = { pass: false, got: 'Excepción: ' + err.message }; }
      return { id: t.id, name: t.name, input: t.input, expected: t.expected, got: r.got, pass: !!r.pass,
        conclusion: r.pass ? 'APROBADA: el comportamiento coincide con lo esperado.' : 'FALLIDA: revisar el resultado obtenido.' };
    });
  }
  F.QA = { run };
})(typeof window !== 'undefined' ? window : globalThis);
