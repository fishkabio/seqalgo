import { assertTruthy } from '@fishka/assertions';

/**
 * Needleman-Wunsch global alignment with EMBOSS-equivalent scoring (EDNAFULL
 * matrix, gapopen=10, gapextend=0.5; gap penalty formula: gapopen + (n-1) * gapextend).
 *
 * The EDNAFULL (NUC4.4) matrix covers all IUPAC ambiguity codes:
 * A, T, G, C, R (A/G), Y (C/T), K (G/T), M (A/C), S (G/C), W (A/T),
 * B (C/G/T), V (A/C/G), H (A/C/T), D (A/G/T), N (any).
 *
 * Scores are always identical to EMBOSS needle. The placement of indels within
 * homopolymer/repeat tracts (where several alignments score equally) is governed
 * by the `gapAnchor` option:
 *
 * - '5-prime' (default): left-anchored gaps, matching EMBOSS-needle default
 *   tie-breaking. An extra C in the rCRS HV2 poly-C tract 303-309 lands at the
 *   5' end of the run.
 * - '3-prime': right-anchored gaps, per the ISFG forensic mtDNA convention
 *   (Parson et al. 2014, §3.2 and rule #3 on p. 137). Within homopolymer/repeat
 *   tracts the gap is anchored to the 3' end of the run, so the same extra C is
 *   reported as 309.1C, not 302.1C.
 *
 * Implementation notes: sequences are encoded to integer codes once (O(n+m)),
 * the scoring matrix is a flat Int8Array, and the DP tables are flat typed
 * arrays (Float64Array scores + Uint8Array traceback). This keeps memory at one
 * contiguous (n+1)*(m+1) block per table and avoids per-cell string/object work.
 */

/** Where to anchor indels within homopolymer/repeat tracts at score ties. */
export type GapAnchor = '5-prime' | '3-prime';

/** Options for {@link needleAlign}. */
export interface NeedleOptions {
  /** Gap opening penalty (EMBOSS gapopen). Default: 10. */
  gapOpen?: number;
  /** Gap extension penalty (EMBOSS gapextend). Default: 0.5. */
  gapExtend?: number;
  /**
   * Indel placement within equally-scoring alignments. Default: '5-prime'
   * (EMBOSS-needle default). Use '3-prime' for the ISFG forensic mtDNA
   * notation convention.
   */
  gapAnchor?: GapAnchor;
}

/** Alignment result */
export interface NeedleAlignment {
  /** Aligned reference sequence */
  seqA: string;
  /** Aligned read sequence */
  seqB: string;
  /** Alignment score */
  score: number;
}

const DEFAULT_GAP_OPEN = 10;
const DEFAULT_GAP_EXTEND = 0.5;
const DEFAULT_GAP_ANCHOR: GapAnchor = '5-prime';

/** IUPAC nucleotide codes supported by EDNAFULL matrix */
const IUPAC_CODES = 'ATGCSWRYKMBVHDN' as const;

/** Characters that create gaps (non-IUPAC, non-X, non-Z, non-whitespace) */
const GAP_CHARS = new Set('!@#$%^&*()_+-=[]{};\':",./<>?0123456789');

/** Characters that cause errors in EMBOSS (Z and lowercase z) */
const ERROR_CHARS = new Set('Zz');

type IupacCode = (typeof IUPAC_CODES)[number];

/**
 * EDNAFULL (NUC4.4) scoring matrix.
 * Rows and columns ordered: A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N
 *
 * Values derived from EMBOSS EDNAFULL matrix:
 * - Exact match (A-A): +5
 * - Partial match (A-W, where W=A/T): +1
 * - N vs any: -2
 * - Mismatch (A-G): -4
 */
const EDNAFULL_MATRIX: Record<IupacCode, Record<IupacCode, number>> = {
  A: { A: 5, T: -4, G: -4, C: -4, S: -4, W: 1, R: 1, Y: -4, K: -4, M: 1, B: -4, V: -1, H: -1, D: -1, N: -2 },
  T: { A: -4, T: 5, G: -4, C: -4, S: -4, W: 1, R: -4, Y: 1, K: 1, M: -4, B: -1, V: -4, H: -1, D: -1, N: -2 },
  G: { A: -4, T: -4, G: 5, C: -4, S: 1, W: -4, R: 1, Y: -4, K: 1, M: -4, B: -1, V: -1, H: -4, D: -1, N: -2 },
  C: { A: -4, T: -4, G: -4, C: 5, S: 1, W: -4, R: -4, Y: 1, K: -4, M: 1, B: -1, V: -1, H: -1, D: -4, N: -2 },
  S: { A: -4, T: -4, G: 1, C: 1, S: -1, W: -4, R: -2, Y: -2, K: -2, M: -2, B: -1, V: -1, H: -3, D: -3, N: -1 },
  W: { A: 1, T: 1, G: -4, C: -4, S: -4, W: -1, R: -2, Y: -2, K: -2, M: -2, B: -3, V: -3, H: -1, D: -1, N: -1 },
  R: { A: 1, T: -4, G: 1, C: -4, S: -2, W: -2, R: -1, Y: -4, K: -2, M: -2, B: -3, V: -1, H: -3, D: -1, N: -1 },
  Y: { A: -4, T: 1, G: -4, C: 1, S: -2, W: -2, R: -4, Y: -1, K: -2, M: -2, B: -1, V: -3, H: -1, D: -3, N: -1 },
  K: { A: -4, T: 1, G: 1, C: -4, S: -2, W: -2, R: -2, Y: -2, K: -1, M: -4, B: -1, V: -3, H: -3, D: -1, N: -1 },
  M: { A: 1, T: -4, G: -4, C: 1, S: -2, W: -2, R: -2, Y: -2, K: -4, M: -1, B: -3, V: -1, H: -1, D: -3, N: -1 },
  B: { A: -4, T: -1, G: -1, C: -1, S: -1, W: -3, R: -3, Y: -1, K: -1, M: -3, B: -1, V: -2, H: -2, D: -2, N: -1 },
  V: { A: -1, T: -4, G: -1, C: -1, S: -1, W: -3, R: -1, Y: -3, K: -3, M: -1, B: -2, V: -1, H: -2, D: -2, N: -1 },
  H: { A: -1, T: -1, G: -4, C: -1, S: -3, W: -1, R: -3, Y: -1, K: -3, M: -1, B: -2, V: -2, H: -1, D: -2, N: -1 },
  D: { A: -1, T: -1, G: -1, C: -4, S: -3, W: -1, R: -1, Y: -3, K: -1, M: -3, B: -2, V: -2, H: -2, D: -1, N: -1 },
  N: { A: -2, T: -2, G: -2, C: -2, S: -1, W: -1, R: -1, Y: -1, K: -1, M: -1, B: -1, V: -1, H: -1, D: -1, N: -1 },
};

/** Number of IUPAC codes (15). 'N' is the last index and the fallback code. */
const N_CODES = IUPAC_CODES.length;

/**
 * Flat EDNAFULL scoring matrix: score(a, b) = FLAT_SCORE[a * N_CODES + b], where
 * a and b are the integer codes assigned by CHAR_TO_CODE. Built once from
 * EDNAFULL_MATRIX so the inner DP loop is a single typed-array read.
 */
const FLAT_SCORE: Int8Array = ((): Int8Array => {
  const flat = new Int8Array(N_CODES * N_CODES);
  for (let a = 0; a < N_CODES; a++) {
    const rowName = IUPAC_CODES[a] as IupacCode;
    for (let b = 0; b < N_CODES; b++) {
      flat[a * N_CODES + b] = EDNAFULL_MATRIX[rowName][IUPAC_CODES[b] as IupacCode];
    }
  }
  return flat;
})();

/**
 * Maps an (uppercase) ASCII char code to its IUPAC integer code. Anything not in
 * IUPAC_CODES maps to 'N' — matching the original "unknown is treated as N"
 * behavior (dead in practice, since preprocessSequence emits only IUPAC chars).
 */
const CHAR_TO_CODE: Int8Array = ((): Int8Array => {
  const map = new Int8Array(256).fill(N_CODES - 1); // default to 'N'
  for (let i = 0; i < N_CODES; i++) map[IUPAC_CODES.charCodeAt(i)] = i;
  return map;
})();

/** Encodes a preprocessed sequence to IUPAC integer codes, O(length). */
function encodeSequence(seq: string): Int8Array {
  const len = seq.length;
  const out = new Int8Array(len);
  for (let i = 0; i < len; i++) out[i] = CHAR_TO_CODE[seq.charCodeAt(i) & 0xff];
  return out;
}

/**
 * Preprocesses sequence for EMBOSS compatibility:
 * - X/x → N
 * - Z/z → throws error
 * - Whitespace → removed
 * - Other special chars → removed (gaps will be inserted by alignment algorithm)
 * - Lowercase → uppercase
 */
function preprocessSequence(seq: string): string {
  let result = '';
  for (const char of seq) {
    // Skip whitespace (spaces, tabs, newlines)
    if (/\s/.test(char)) {
      continue;
    }

    const upper = char.toUpperCase();

    // Z causes error
    if (ERROR_CHARS.has(char)) {
      throw new Error(`Sequence contains invalid character: ${char}`);
    }

    // X becomes N
    if (upper === 'X') {
      result += 'N';
      continue;
    }

    // Gap characters and other special chars are removed
    // (gaps will be inserted by the alignment algorithm if needed)
    if (GAP_CHARS.has(char) || GAP_CHARS.has(upper)) {
      continue;
    }

    // Valid IUPAC codes
    if (IUPAC_CODES.includes(upper as IupacCode)) {
      result += upper;
      continue;
    }

    // Any other character is removed
    continue;
  }
  return result;
}

/** Traceback directions (stored in a Uint8Array). 0 is the unset/origin value. */
const DIAG = 1; // 'M' — match/mismatch
const LEFT = 2; // 'X' — gap in seqA (insertion in seqB)
const UP = 3; // 'Y' — gap in seqB (deletion in seqA)

const GAP_CHAR = 45; // '-'

/**
 * Performs Needleman-Wunsch global alignment.
 * Scores are 100% compatible with EMBOSS needle using the EDNAFULL matrix; the
 * `gapAnchor` option controls indel placement within equally-scoring tracts.
 *
 * @param seqA Reference sequence
 * @param seqB Read sequence
 * @param options Gap penalties and indel anchoring (see {@link NeedleOptions})
 * @returns Alignment result with gapped sequences
 */
export function needleAlign(seqA: string, seqB: string, options: NeedleOptions = {}): NeedleAlignment {
  const { gapOpen = DEFAULT_GAP_OPEN, gapExtend = DEFAULT_GAP_EXTEND, gapAnchor = DEFAULT_GAP_ANCHOR } = options;

  assertTruthy(seqA && seqA.length > 0, 'Reference sequence is empty');
  assertTruthy(seqB && seqB.length > 0, 'Read sequence is empty');

  // Preprocess for EMBOSS compatibility, then encode once for scoring.
  const processedA = preprocessSequence(seqA);
  const processedB = preprocessSequence(seqB);
  const codeA = encodeSequence(processedA);
  const codeB = encodeSequence(processedB);

  const lenA = codeA.length;
  const lenB = codeB.length;
  const width = lenB + 1;

  // Flat DP tables. Float64Array starts zero-filled, matching the no-end-gap
  // initialization (first row/column scores are 0).
  const score = new Float64Array((lenA + 1) * width);
  const traceback = new Uint8Array((lenA + 1) * width);

  // First column: gaps at the beginning of seqB (UP). First row: gaps at the
  // beginning of seqA (LEFT). Origin stays NONE.
  for (let i = 1; i <= lenA; i++) traceback[i * width] = UP;
  for (let j = 1; j <= lenB; j++) traceback[j] = LEFT;

  const prefer3Prime = gapAnchor === '3-prime';

  // Fill the matrices. Tie-break order decides indel placement:
  //   3-prime → LEFT > UP > DIAG (gaps pushed to the 3' end of repeats).
  //   5-prime → DIAG > UP > LEFT (gaps anchored to the 5' end, EMBOSS default).
  for (let i = 1; i <= lenA; i++) {
    const rowBase = i * width;
    const prevBase = rowBase - width;
    const aCode = codeA[i - 1] * N_CODES;
    for (let j = 1; j <= lenB; j++) {
      const diagScore = score[prevBase + j - 1] + FLAT_SCORE[aCode + codeB[j - 1]];
      const upScore = score[prevBase + j] - (traceback[prevBase + j] === UP ? gapExtend : gapOpen);
      const leftScore = score[rowBase + j - 1] - (traceback[rowBase + j - 1] === LEFT ? gapExtend : gapOpen);

      let cellScore: number;
      let dir: number;
      if (prefer3Prime) {
        if (leftScore >= diagScore && leftScore >= upScore) {
          cellScore = leftScore;
          dir = LEFT;
        } else if (upScore >= diagScore) {
          cellScore = upScore;
          dir = UP;
        } else {
          cellScore = diagScore;
          dir = DIAG;
        }
      } else if (diagScore >= leftScore && diagScore >= upScore) {
        cellScore = diagScore;
        dir = DIAG;
      } else if (upScore >= leftScore) {
        cellScore = upScore;
        dir = UP;
      } else {
        cellScore = leftScore;
        dir = LEFT;
      }

      score[rowBase + j] = cellScore;
      traceback[rowBase + j] = dir;
    }
  }

  // Optimal end position, allowing end gaps without penalty: best of the corner,
  // the last row (end gaps in seqA) and the last column (end gaps in seqB).
  let maxScore = score[lenA * width + lenB];
  let endI = lenA;
  let endJ = lenB;
  const lastRow = lenA * width;
  for (let j = 0; j <= lenB; j++) {
    if (score[lastRow + j] > maxScore) {
      maxScore = score[lastRow + j];
      endI = lenA;
      endJ = j;
    }
  }
  for (let i = 0; i <= lenA; i++) {
    if (score[i * width + lenB] > maxScore) {
      maxScore = score[i * width + lenB];
      endI = i;
      endJ = lenB;
    }
  }

  const { seqA: alignedA, seqB: alignedB } = buildAlignedSequences(
    traceback,
    processedA,
    processedB,
    width,
    endI,
    endJ,
  );
  return { seqA: alignedA, seqB: alignedB, score: maxScore };
}

/**
 * Reconstructs the aligned strings from the traceback table without quadratic
 * string concatenation: columns are pushed into arrays of char codes in 3'→5'
 * order, reversed once, then the trailing end-gap columns are appended.
 */
function buildAlignedSequences(
  traceback: Uint8Array,
  seqA: string,
  seqB: string,
  width: number,
  endI: number,
  endJ: number,
): { seqA: string; seqB: string } {
  const revA: number[] = [];
  const revB: number[] = [];
  let i = endI;
  let j = endJ;

  // Walk from (endI, endJ) back to the origin, recording columns in reverse.
  while (i > 0 || j > 0) {
    const dir = traceback[i * width + j];
    if (dir === DIAG && i > 0 && j > 0) {
      revA.push(seqA.charCodeAt(i - 1));
      revB.push(seqB.charCodeAt(j - 1));
      i--;
      j--;
    } else if (dir === UP && i > 0) {
      revA.push(seqA.charCodeAt(i - 1));
      revB.push(GAP_CHAR);
      i--;
    } else if (dir === LEFT && j > 0) {
      revA.push(GAP_CHAR);
      revB.push(seqB.charCodeAt(j - 1));
      j--;
    } else if (i > 0) {
      // Fallback: gap in seqB.
      revA.push(seqA.charCodeAt(i - 1));
      revB.push(GAP_CHAR);
      i--;
    } else if (j > 0) {
      // Fallback: gap in seqA.
      revA.push(GAP_CHAR);
      revB.push(seqB.charCodeAt(j - 1));
      j--;
    } else {
      break;
    }
  }

  revA.reverse();
  revB.reverse();

  // Trailing end gaps (if we ended before the full length), appended at the 3' end.
  for (let k = endI; k < seqA.length; k++) {
    revA.push(seqA.charCodeAt(k));
    revB.push(GAP_CHAR);
  }
  for (let k = endJ; k < seqB.length; k++) {
    revA.push(GAP_CHAR);
    revB.push(seqB.charCodeAt(k));
  }

  return { seqA: codesToString(revA), seqB: codesToString(revB) };
}

/** Joins char codes into a string in chunks (avoids apply() arg-count limits). */
function codesToString(codes: number[]): string {
  const CHUNK = 8192;
  if (codes.length <= CHUNK) return String.fromCharCode(...codes);
  let out = '';
  for (let i = 0; i < codes.length; i += CHUNK) {
    out += String.fromCharCode(...codes.slice(i, i + CHUNK));
  }
  return out;
}
