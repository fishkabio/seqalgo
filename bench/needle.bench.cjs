/*
 * Zero-dependency performance suite for needleAlign.
 *
 * Benchmarks the built package (dist/) so numbers reflect what consumers run.
 * Run with: npm run bench   (rebuilds first)
 *
 * Two scenarios:
 *   1. Realistic mtDNA HV zones (HV3 ~150, HV2 ~270, HV1 ~340).
 *   2. Scaling sweep (100 → 2000 bp) to expose the O(n*m) growth and memory cost.
 */

const { needleAlign } = require('../dist/index.js');

const BASES = 'ACGT';

/** Deterministic LCG so runs are reproducible. */
function makeRng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function randomSeq(len, rng) {
  let out = '';
  for (let i = 0; i < len; i++) out += BASES[(rng() * 4) | 0];
  return out;
}

/**
 * Builds a read from a reference: ~1% SNPs, one homopolymer C-tract with a +1C
 * insertion, and one single-base deletion — exercises the match, gap and
 * tie-break paths the way a real sequencing read does.
 */
function makeReadPair(len, rng) {
  let ref = randomSeq(len, rng);
  // Plant a homopolymer C-tract near the middle to stress the tie-break path.
  const mid = Math.floor(len / 2);
  ref = ref.slice(0, mid) + 'CCCCCCC' + ref.slice(mid + 7);

  const chars = ref.split('');
  for (let i = 0; i < chars.length; i++) {
    if (rng() < 0.01) chars[i] = BASES[(rng() * 4) | 0]; // SNP
  }
  let read = chars.join('');
  // +1C insertion inside the C-tract.
  read = read.slice(0, mid + 3) + 'C' + read.slice(mid + 3);
  // single-base deletion near 3/4.
  const del = Math.floor(len * 0.75);
  read = read.slice(0, del) + read.slice(del + 1);
  return { ref, read };
}

/** Times fn: warms up, then runs for ~minTime collecting per-call samples. */
function bench(fn, { warmup = 30, minTimeMs = 1500 } = {}) {
  for (let i = 0; i < warmup; i++) fn();
  const samples = [];
  const start = process.hrtime.bigint();
  do {
    const t0 = process.hrtime.bigint();
    fn();
    samples.push(Number(process.hrtime.bigint() - t0));
  } while (Number(process.hrtime.bigint() - start) / 1e6 < minTimeMs);

  samples.sort((a, b) => a - b);
  const mean = samples.reduce((a, b) => a + b, 0) / samples.length / 1e6;
  const median = samples[samples.length >> 1] / 1e6;
  const p95 = samples[Math.floor(samples.length * 0.95)] / 1e6;
  return { iters: samples.length, meanMs: mean, medianMs: median, p95Ms: p95, opsSec: 1000 / mean };
}

function row(label, n, m, r) {
  return {
    case: label,
    'ref×read': `${n}×${m}`,
    'mean ms': r.meanMs.toFixed(3),
    'median ms': r.medianMs.toFixed(3),
    'p95 ms': r.p95Ms.toFixed(3),
    'ops/sec': Math.round(r.opsSec).toLocaleString('en-US'),
    iters: r.iters,
  };
}

function peakMemoryFor(ref, read) {
  if (!global.gc) return null;
  global.gc();
  const before = process.memoryUsage().heapUsed;
  let peak = before;
  for (let i = 0; i < 10; i++) {
    needleAlign(ref, read);
    peak = Math.max(peak, process.memoryUsage().heapUsed);
  }
  return ((peak - before) / 1024 / 1024).toFixed(1);
}

console.log(`\nnode ${process.version} — needleAlign performance\n`);

// Scenario 1: realistic mtDNA HV zones.
const rng1 = makeRng(42);
const hv = [
  ['HV3', 150],
  ['HV2', 270],
  ['HV1', 340],
];
const table1 = [];
for (const [name, len] of hv) {
  const { ref, read } = makeReadPair(len, rng1);
  table1.push(
    row(
      `${name} (5')`,
      ref.length,
      read.length,
      bench(() => needleAlign(ref, read)),
    ),
  );
  table1.push(
    row(
      `${name} (3')`,
      ref.length,
      read.length,
      bench(() => needleAlign(ref, read, { gapAnchor: '3-prime' })),
    ),
  );
}
console.log('Scenario 1 — realistic mtDNA HV zones');
console.table(table1);

// Scenario 2: scaling sweep.
const rng2 = makeRng(7);
const sizes = [100, 250, 500, 1000, 2000];
const table2 = [];
for (const len of sizes) {
  const { ref, read } = makeReadPair(len, rng2);
  const r = bench(() => needleAlign(ref, read), { minTimeMs: len >= 1000 ? 2500 : 1500 });
  const mem = peakMemoryFor(ref, read);
  table2.push({ ...row(`${len} bp`, ref.length, read.length, r), ...(mem ? { 'heap MB': mem } : {}) });
}
console.log('Scenario 2 — scaling sweep (quadratic O(n*m))');
console.table(table2);

// Quadratic-scaling sanity: ms should grow ~4x per length doubling.
const base = Number(table2[0]['mean ms']);
console.log('\nScaling factor vs 100bp (ideal quadratic in parentheses):');
for (let i = 0; i < sizes.length; i++) {
  const factor = Number(table2[i]['mean ms']) / base;
  const ideal = (sizes[i] / sizes[0]) ** 2;
  console.log(`  ${String(sizes[i]).padStart(4)}bp: ${factor.toFixed(1)}x  (${ideal.toFixed(0)}x)`);
}
console.log('');
