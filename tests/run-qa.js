// Ejecuta el módulo QA en Node:  node tests/run-qa.js
require('../js/linalg.js'); require('../js/pcg.js'); require('../js/model.js'); require('../js/qa.js');
const rows = globalThis.FinCG.QA.run();
let fail = 0;
for (const r of rows) { console.log((r.pass ? 'PASS ' : 'FAIL ') + r.id + ' ' + r.name + '\n       ' + r.got); if (!r.pass) fail++; }
console.log(fail ? fail + ' prueba(s) fallida(s)' : 'Todas las pruebas pasaron (' + rows.length + ')');
process.exit(fail ? 1 : 0);
