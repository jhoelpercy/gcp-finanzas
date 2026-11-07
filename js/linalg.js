/* Álgebra lineal densa mínima (sin dependencias). Namespace global: FinCG */
(function (root) {
  'use strict';
  const FinCG = root.FinCG || (root.FinCG = {});

  const zeros = n => new Array(n).fill(0);
  const copyM = A => A.map(r => r.slice());
  function A_from(n, f) { const M = []; for (let i = 0; i < n; i++) { const r = []; for (let j = 0; j < n; j++) r.push(f(i, j)); M.push(r); } return M; }
  const identity = n => A_from(n, (i, j) => (i === j ? 1 : 0));
  const dot = (u, v) => { let s = 0; for (let i = 0; i < u.length; i++) s += u[i] * v[i]; return s; };
  const norm2 = v => Math.sqrt(dot(v, v));
  const matvec = (A, x) => A.map(r => dot(r, x));
  const sub = (u, v) => u.map((a, i) => a - v[i]);
  function matmul(A, B) {
    const n = A.length, k = B.length;
    return A_from(n, (i, j) => { let s = 0; for (let t = 0; t < k; t++) s += A[i][t] * B[t][j]; return s; });
  }
  const transpose = A => A_from(A.length, (i, j) => A[j][i]);
  const maxAbs = A => A.reduce((m, r) => Math.max(m, ...r.map(Math.abs)), 0);

  /* Mayor |a_ij - a_ji| */
  function symmetryDefect(A) {
    let d = 0; const n = A.length;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) d = Math.max(d, Math.abs(A[i][j] - A[j][i]));
    return d;
  }
  function isSymmetric(A, tol = 1e-10) { return symmetryDefect(A) <= tol * Math.max(1, maxAbs(A)); }

  /* Cholesky A = L L^T. Devuelve null si A no es definida positiva. */
  function cholesky(A) {
    const n = A.length, L = A_from(n, () => 0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let s = A[i][j];
        for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
        if (i === j) { if (!(s > 0) || !isFinite(s)) return null; L[i][i] = Math.sqrt(s); }
        else L[i][j] = s / L[j][j];
      }
    }
    return L;
  }
  function forwardSolve(L, b) { // L y = b
    const n = b.length, y = zeros(n);
    for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; }
    return y;
  }
  function backSolveT(L, y) { // L^T x = y
    const n = y.length, x = zeros(n);
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; }
    return x;
  }
  const choleskySolve = (L, b) => backSolveT(L, forwardSolve(L, b));

  /* LU con pivoteo parcial: PA = LU */
  function luDecompose(A) {
    const n = A.length, M = copyM(A), piv = [...Array(n).keys()]; let sign = 1, singular = false;
    for (let k = 0; k < n; k++) {
      let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
      if (Math.abs(M[p][k]) < 1e-300) { singular = true; continue; }
      if (p !== k) { [M[p], M[k]] = [M[k], M[p]]; [piv[p], piv[k]] = [piv[k], piv[p]]; sign = -sign; }
      for (let i = k + 1; i < n; i++) {
        M[i][k] /= M[k][k];
        for (let j = k + 1; j < n; j++) M[i][j] -= M[i][k] * M[k][j];
      }
    }
    return { LU: M, piv, sign, singular };
  }
  function luSolve(Fc, b) {
    const n = b.length, y = zeros(n), x = zeros(n);
    for (let i = 0; i < n; i++) { let s = b[Fc.piv[i]]; for (let k = 0; k < i; k++) s -= Fc.LU[i][k] * y[k]; y[i] = s; }
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= Fc.LU[i][k] * x[k]; x[i] = s / Fc.LU[i][i]; }
    return x;
  }
  function det(A) { const Fc = luDecompose(A); if (Fc.singular) return 0; let d = Fc.sign; for (let i = 0; i < A.length; i++) d *= Fc.LU[i][i]; return d; }
  /* Menores principales líderes (Criterio de Sylvester) */
  function leadingMinors(A) {
    const out = [];
    for (let k = 1; k <= A.length; k++) out.push(det(A.slice(0, k).map(r => r.slice(0, k))));
    return out;
  }

  /* Autovalores de una matriz simétrica (Jacobi cíclico), orden ascendente */
  function eigenvaluesSym(A0) {
    const n = A0.length, A = copyM(A0);
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0, dg = 0;
      for (let i = 0; i < n; i++) { dg += A[i][i] * A[i][i]; for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j]; }
      if (off <= 1e-30 * Math.max(dg, 1e-300)) break;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
      }
    }
    return A.map((r, i) => r[i]).sort((a, b) => a - b);
  }
  function cond2FromEigs(ev) {
    const a = ev.map(Math.abs), mn = Math.min(...a), mx = Math.max(...a);
    return mn === 0 ? Infinity : mx / mn;
  }

  Object.assign(FinCG, { zeros, copyM, identity, A_from, dot, norm2, matvec, sub, matmul, transpose, maxAbs,
    symmetryDefect, isSymmetric, cholesky, forwardSolve, backSolveT, choleskySolve, luDecompose, luSolve, det,
    leadingMinors, eigenvaluesSym, cond2FromEigs });
})(typeof window !== 'undefined' ? window : globalThis);
