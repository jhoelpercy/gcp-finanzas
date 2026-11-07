/* Gradientes Conjugados (Precondicionado) y precondicionadores. Requiere linalg.js */
(function (root) {
  'use strict';
  const F = root.FinCG;
  const { dot, norm2, matvec, sub, zeros, A_from } = F;

  /* ---------- Precondicionadores ----------
     Devuelven { type, label, solve(r) -> z (resuelve M z = r), M (matriz densa), warning? } */
  function makePrecond(type, A, opts = {}) {
    const n = A.length;
    const omega = opts.omega ?? 1;
    if (type === 'none') {
      return { type, label: 'Ninguno (M = I)', solve: r => r.slice(), M: F.identity(n) };
    }
    if (type === 'jacobi') {
      const d = A.map((r, i) => r[i]);
      if (d.some(v => !(v > 0))) throw new Error('Jacobi requiere diagonal positiva.');
      return { type, label: 'Jacobi (M = diag(A))', solve: r => r.map((v, i) => v / d[i]), M: A_from(n, (i, j) => (i === j ? d[i] : 0)) };
    }
    if (type === 'ssor') {
      if (!(omega > 0 && omega < 2)) throw new Error('SSOR requiere 0 < ω < 2.');
      const d = A.map((r, i) => r[i]);
      if (d.some(v => !(v > 0))) throw new Error('SSOR requiere diagonal positiva.');
      // M = ω/(2-ω) (D/ω + L) (D/ω)^-1 (D/ω + L^T)
      const T = A_from(n, (i, j) => (i === j ? d[i] / omega : (i > j ? A[i][j] : 0)));
      const solve = r => {
        const y = F.forwardSolve(T, r);                 // (D/ω+L) y = r
        const w = y.map((v, i) => v * d[i] / omega);    // (D/ω) y
        const z = F.backSolveT(T, w);                   // (D/ω+L)^T z = w
        return z.map(v => v * (2 - omega) / omega);
      };
      const Dw = A_from(n, (i, j) => (i === j ? omega / d[i] : 0));
      const M = F.matmul(F.matmul(T, Dw), F.transpose(T)).map(r => r.map(v => v * omega / (2 - omega)));
      return { type, label: 'SSOR (ω = ' + omega + ')', solve, M };
    }
    if (type === 'ic0') {
      // Cholesky incompleto IC(0): L conserva el patrón de ceros de A
      const L = A_from(n, () => 0);
      for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
        if (A[i][j] === 0) continue;
        let s = A[i][j];
        for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
        if (i === j) { if (!(s > 0)) throw new Error('IC(0): pivote no positivo (la factorización incompleta falló).'); L[i][i] = Math.sqrt(s); }
        else L[i][j] = s / L[j][j];
      }
      const M = F.matmul(L, F.transpose(L));
      const dense = A.every(r => r.every(v => v !== 0));
      return { type, label: 'Cholesky incompleto IC(0)', solve: r => F.choleskySolve(L, r), M,
        warning: dense ? 'Con A densa el patrón de IC(0) es completo, así que M = A y GCP converge en ~1 iteración.' : undefined };
    }
    throw new Error('Precondicionador desconocido: ' + type);
  }

  /* κ2(M^-1 A) = cond2 de L^-1 A L^-T, con M = L L^T */
  function condPrecond(A, M) {
    const L = F.cholesky(M);
    if (!L) return NaN;
    const n = A.length;
    const cols = [];
    for (let j = 0; j < n; j++) cols.push(F.forwardSolve(L, A.map(r => r[j]))); // columnas de L^-1 A
    const Y = A_from(n, (i, j) => cols[j][i]);
    const B = A_from(n, (i, j) => F.forwardSolve(L, Y[j])[i]);                   // L^-1 (L^-1 A)^T
    const S = A_from(n, (i, j) => 0.5 * (B[i][j] + B[j][i]));
    return F.cond2FromEigs(F.eigenvaluesSym(S));
  }

  /* ---------- Gradientes Conjugados Precondicionado ----------
     r0 = b - A x0;  M z0 = r0;  p0 = z0
     α_k = (r_k·z_k)/(p_k·A p_k);  x_{k+1} = x_k + α_k p_k;  r_{k+1} = r_k - α_k A p_k
     M z_{k+1} = r_{k+1};  β_k = (r_{k+1}·z_{k+1})/(r_k·z_k);  p_{k+1} = z_{k+1} + β_k p_k
     Con M = I se obtiene el GC simple. */
  function pcg(A, b, x0, tol, maxIter, pre) {
    const n = b.length, nb = norm2(b);
    const x = x0 ? x0.slice() : zeros(n);
    const r = sub(b, matvec(A, x));
    const history = [];
    const rel = v => (nb > 0 ? norm2(v) / nb : norm2(v));
    const finish = (converged, iters, note) => {
      const trueRes = rel(sub(b, matvec(A, x)));
      return { x, iters, converged, history, residual: trueRes, note };
    };
    history.push(rel(r));
    if (history[0] < tol) return finish(true, 0, null);
    let z = pre.solve(r);
    const p = z.slice();
    let rz = dot(r, z);
    for (let k = 1; k <= maxIter; k++) {
      const Ap = matvec(A, p), pAp = dot(p, Ap);
      if (!(pAp > 0)) return finish(false, k - 1, 'p·Ap ≤ 0: A no es definida positiva (ruptura)');
      const alpha = rz / pAp;
      for (let i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
      const rr = rel(r);
      history.push(rr);
      if (!isFinite(rr)) return finish(false, k, 'Divergencia numérica');
      if (rr < tol) return finish(true, k, null);
      z = pre.solve(r);
      const rzNew = dot(r, z), beta = rzNew / rz; rz = rzNew;
      for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    }
    return finish(false, maxIter, 'Se alcanzó N_máx sin cumplir la tolerancia');
  }

  Object.assign(F, { makePrecond, condPrecond, pcg });
})(typeof window !== 'undefined' ? window : globalThis);
