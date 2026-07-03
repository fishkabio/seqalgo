# @fishka/seqalgo

Part of [fishka.bio](https://fishka.bio) — free browser-based bioinformatics tools.

Sequence algorithms for browser and Node bioinformatics apps. It is the
algorithmic companion to [`@fishka/seqio`](https://www.npmjs.com/package/@fishka/seqio):
`seqio` reads and writes sequence files, `seqalgo` analyses the sequences.

Current API:

- **Needleman-Wunsch global alignment** with EDNAFULL scoring and IUPAC
  ambiguity-code support. Scores match EMBOSS `needle`.
- **Sequence utilities**: gap removal/counting and IUPAC-aware reverse
  complement.
- **Quality trimming**: modified Mott and sliding-window trim ranges.
- **Alignment utilities**: coordinate conversion, gap-column cleanup, pairwise
  alignment combination, and gap-filled transforms.
- **Chromatogram utilities**: reverse-complement traces and inject alignment
  gaps into chromatograms.

Planned: mutation classification, consensus calling, heteroplasmy detection.

## Install

```sh
npm install @fishka/seqalgo
```

## Needleman-Wunsch

```ts
import { needleAlign } from '@fishka/seqalgo';

const result = needleAlign('GATCACAGGT', 'GATCAGGT');

result.seqA; // aligned reference: "GATCACAGGT"
result.seqB; // aligned read:      "GAT--CAGGT"
result.score; // EMBOSS-equivalent score
```

```ts
needleAlign(ref, read, {
  gapOpen: 10,
  gapExtend: 0.5,
  gapAnchor: '5-prime',
});
```

`gapAnchor` only affects where equally-scoring indels land inside
homopolymer/repeat tracts:

- `'5-prime'` (default) matches EMBOSS `needle` tie-breaking.
- `'3-prime'` right-anchors gaps for ISFG forensic mtDNA notation.

## Sequence and quality helpers

```ts
import { mottTrim, reverseComplement, slidingWindowTrim } from '@fishka/seqalgo';

reverseComplement('ACGTRY'); // "RYACGT"

mottTrim([8, 12, 30, 31, 29, 10], { cutoff: 20 }); // { start, end }
slidingWindowTrim([8, 12, 30, 31, 29, 10], { threshold: 20, windowSize: 3 });
```

## Alignment and chromatogram helpers

```ts
import {
  applyAlignmentGapsToChromatogram,
  combineAlignments,
  removeColumnsOfGaps,
  reverseComplementChromatogram,
} from '@fishka/seqalgo';

const cleaned = removeColumnsOfGaps(['A-C-', '--C-']);
const combined = combineAlignments(pairwiseAlignments);
const rc = reverseComplementChromatogram(chromatogram);
const gappedTrace = applyAlignmentGapsToChromatogram(rc, 'AC-GT');
```

Subpath imports are available for `@fishka/seqalgo/needle`,
`@fishka/seqalgo/sequence`, and `@fishka/seqalgo/alignment`.

## License

Apache-2.0
