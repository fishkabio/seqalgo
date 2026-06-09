# @fishka/seqalgo

Bioinformatics algorithms over biological sequences. The algorithmic companion
to [`@fishka/seqio`](https://www.npmjs.com/package/@fishka/seqio): `seqio` reads
and writes sequence file formats, `seqalgo` analyses the sequences.

Currently supports:

- **Needleman-Wunsch global alignment** — pairwise DNA alignment with the
  EDNAFULL (NUC4.4) scoring matrix and full IUPAC ambiguity-code support. Scores
  are 100% compatible with EMBOSS `needle`.

Planned: alignment post-processing utilities, mutation classification, consensus
calling, heteroplasmy detection.

## Install

```sh
npm install @fishka/seqalgo
```

## Use

```ts
import { needleAlign } from '@fishka/seqalgo';
// or: import { needleAlign } from '@fishka/seqalgo/needle';

const result = needleAlign('GATCACAGGT', 'GATCAGGT');
result.seqA; // aligned reference: "GATCACAGGT"
result.seqB; // aligned read:      "GAT--CAGGT"
result.score; // EMBOSS-equivalent alignment score: 29.5
```

### Options

```ts
needleAlign(ref, read, {
  gapOpen: 10, // EMBOSS gapopen, default 10
  gapExtend: 0.5, // EMBOSS gapextend, default 0.5
  gapAnchor: '5-prime', // indel placement in equally-scoring tracts, default '5-prime'
});
```

### Indel anchoring (`gapAnchor`)

Alignment **scores are always identical to EMBOSS** `needle`. The option only
affects **where indels land within homopolymer/repeat tracts**, where several
alignments score equally:

- **`'5-prime'`** (default) — gaps are left-anchored, matching the EMBOSS
  `needle` default tie-break.
- **`'3-prime'`** — gaps are right-anchored, per the ISFG forensic mtDNA
  notation convention (Parson et al. 2014, §3.2). An extra C in the rCRS HV2
  poly-C tract 303-309 is then reported at the 3' end of the run (309.1C) rather
  than the 5' end (302.1C).

```ts
// rCRS HV2 fragment with an extra C in the 7-C tract (303-309) — same score,
// different placement. The tract must be flanked on both sides for the anchoring
// to matter; with unflanked synthetic strings a free end-gap masks the tie-break.
const ref = 'AATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAG';
const read = 'AATTTCCACCAAACCCCCCCCTCCCCCGCTTCTGGCCACAG';

needleAlign(ref, read).seqA;
// → "AATTTCCACCAAA-CCCCCCCTCCCCCGCTTCTGGCCACAG"  (5'-anchored, default)

needleAlign(ref, read, { gapAnchor: '3-prime' }).seqA;
// → "AATTTCCACCAAACCCCCCC-TCCCCCGCTTCTGGCCACAG"  (3'-anchored, ISFG mtDNA convention)
```

## License

Apache-2.0
