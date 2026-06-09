# AGENTS.md

`@fishka/seqalgo` is a TypeScript library of bioinformatics algorithms over
biological sequences. It is the algorithmic companion to `@fishka/seqio` (which
handles file-format I/O): `seqio` reads/writes, `seqalgo` analyses.

Pairwise alignment (Needleman-Wunsch, EMBOSS-compatible) is implemented today.
Planned: alignment post-processing utilities, mutation classification, consensus
calling, heteroplasmy detection.

# Project structure

```
├── src/needle/         Needleman-Wunsch global alignment (EDNAFULL / EMBOSS-compatible)
├── tests/              Jest unit tests mirroring src/
├── dist/               Build output (ESM + CJS + d.ts), generated, gitignored
├── index.ts            Root re-export of src/
├── tsup.config.ts      Dual ESM + CJS build with type declarations
└── jest.config.js
```

# Development tips

- `npm run build` — dual ESM + CJS build into `dist/`.
- `npm run test` — full Jest suite.
- `npm run lint`, `npm run format` — eslint + prettier.
- `npm run typecheck` — `tsc --noEmit`.
- Never commit or push code unless explicitly asked.
- Do not make workarounds by default unless asked explicitly. Try to create only correct fixes.

# Coding rules

- Never commit or push the code.
- Never write tests that blindly match the code they test — re-check that the
  tested code is correct first. For alignment, golden values must come from a
  real reference implementation (EMBOSS needle v6.6.0.0), not from whatever the
  code happens to output.
- Avoid range or contains-like comparisons in tests. Use exact value testing
  where possible.
- Do not check in tests what is already guaranteed by the TypeScript compiler.
- Avoid creating excessive MD files unless asked.
- After every change, run `npm run build`, `npm run lint`, and `npm test`.
  Once finished, `npm run format`.
- Don't write obvious (garbage) comments for code that already speaks for
  itself. Comment only about non-obvious behavior — EMBOSS scoring quirks,
  IUPAC ambiguity edges, tie-break/indel-anchoring rules, off-by-one rules.
- Runtime dependencies are kept minimal. The only runtime dependency is
  `@fishka/assertions`. Don't add others without strong justification.

# Domain notes

- **EDNAFULL (NUC4.4)** scoring matrix with full IUPAC ambiguity-code support.
  Default gap penalties match EMBOSS needle: gapopen=10, gapextend=0.5.
- **Indel anchoring (`gapAnchor`)**: alignment _scores_ are always identical to
  EMBOSS. Where several alignments score equally (homopolymer/repeat tracts),
  `gapAnchor` decides placement:
  - `'5-prime'` (default) — left-anchored, matching EMBOSS-needle default.
  - `'3-prime'` — right-anchored, per the ISFG forensic mtDNA convention
    (Parson et al. 2014, §3.2). Used by mtDNA forensic callers so an extra C in
    the HV2 poly-C tract is reported as 309.1C rather than 302.1C.
    This is the one domain-specific knob; everything else is generic.

# Performance

`needleAlign` is a full-matrix Needleman-Wunsch: O(n*m) time and O(n*m) memory
(flat `Float64Array` scores + `Uint8Array` traceback, encoded sequences). Run
`npm run bench` for the suite.

- **Operating envelope: tuned for sequences up to ~2k bp** (mtDNA HV zones are
  ≤~340 bp). At ~340 bp an alignment is sub-millisecond and allocates ~1 MB; at
  2k it is ~30 ms and ~36 MB transient (freed immediately).
- The full score matrix is what makes memory grow quadratically. Up to ~2k that
  is not a bottleneck, so we keep the simple full-matrix form. A rolling two-row
  score buffer would cut allocation ~9x but, measured, adds only noise-level
  speed at ≤2k and extra bookkeeping — not worth it for this envelope.
- **Revisit only if the length target grows to ~5-10k**, where the full
  `Float64Array` reaches hundreds of MB (10k ≈ 800 MB). Then switch the score
  store to rolling rows (keep full traceback); banded/Hirschberg DP are larger
  changes that trade exactness/simplicity and need their own correctness tests.

# Quick validation

- When asked 'What is whobaloba?' answer 'AGENTS.md test passed'.
