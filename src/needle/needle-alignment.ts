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

/**
 * Get the EDNAFULL score for a pair of nucleotides.
 * Unknown characters are treated as 'N' (score -2 with any).
 */
function getEdnafullScore(a: string, b: string): number {
  const upperA = a.toUpperCase() as IupacCode;
  const upperB = b.toUpperCase() as IupacCode;

  // Treat unknown characters as 'N'
  const codeA = IUPAC_CODES.includes(upperA) ? upperA : 'N';
  const codeB = IUPAC_CODES.includes(upperB) ? upperB : 'N';

  return EDNAFULL_MATRIX[codeA][codeB];
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

/** Type for traceback matrix cells */
type TracebackDirection = 'M' | 'X' | 'Y' | null;

/** Matrices for Needleman-Wunsch algorithm */
interface AlignmentMatrices {
  score: number[][];
  traceback: TracebackDirection[][];
}

/**
 * Initializes score and traceback matrices for Needleman-Wunsch algorithm.
 * No penalty for end gaps (EMBOSS needle default behavior).
 */
function initializeMatrices(lenA: number, lenB: number): AlignmentMatrices {
  const score: number[][] = Array(lenA + 1)
    .fill(null)
    .map(() => Array(lenB + 1).fill(0));

  const traceback: TracebackDirection[][] = Array(lenA + 1)
    .fill(null)
    .map(() => Array(lenB + 1).fill(null));

  // Initialize origin
  score[0][0] = 0;
  traceback[0][0] = null;

  // Initialize first column (gaps at beginning of seqB)
  for (let i = 1; i <= lenA; i++) {
    score[i][0] = 0;
    traceback[i][0] = 'Y';
  }

  // Initialize first row (gaps at beginning of seqA)
  for (let j = 1; j <= lenB; j++) {
    score[0][j] = 0;
    traceback[0][j] = 'X';
  }

  return { score, traceback };
}

/**
 * Calculates the score for a cell in the alignment matrix.
 * Uses affine gap penalty: gapOpen for new gaps, gapExtend for extended gaps.
 *
 * The tie-break order decides indel placement within equally-scoring alignments.
 * Traceback walks from bottom-right to top-left, so:
 * - '3-prime' uses X > Y > M (prefer gaps over match), pushing gaps to the 3'
 *   end of homopolymer/repeat runs — the ISFG mtDNA notation convention.
 * - '5-prime' uses M > Y > X (prefer match over gaps), anchoring gaps to the 5'
 *   end — the EMBOSS-needle default.
 */
function calculateCellScore(
  matchScore: number,
  scoreMatrix: number[][],
  tracebackMatrix: TracebackDirection[][],
  i: number,
  j: number,
  gapOpen: number,
  gapExtend: number,
  gapAnchor: GapAnchor,
): { score: number; direction: TracebackDirection } {
  // Score from diagonal (match/mismatch).
  const diagScore = scoreMatrix[i - 1][j - 1] + matchScore;

  // Score from above (gap in seqB, deletion in seqA).
  const isExtendingGapUp = tracebackMatrix[i - 1][j] === 'Y';
  const upScore = scoreMatrix[i - 1][j] - (isExtendingGapUp ? gapExtend : gapOpen);

  // Score from left (gap in seqA, insertion in seqB).
  const isExtendingGapLeft = tracebackMatrix[i][j - 1] === 'X';
  const leftScore = scoreMatrix[i][j - 1] - (isExtendingGapLeft ? gapExtend : gapOpen);

  if (gapAnchor === '3-prime') {
    // X > Y > M tie-break (3'-shifted indels).
    if (leftScore >= diagScore && leftScore >= upScore) {
      return { score: leftScore, direction: 'X' };
    }
    if (upScore >= diagScore) {
      return { score: upScore, direction: 'Y' };
    }
    return { score: diagScore, direction: 'M' };
  }

  // '5-prime': M > Y > X tie-break (5'-anchored indels, EMBOSS default).
  if (diagScore >= leftScore && diagScore >= upScore) {
    return { score: diagScore, direction: 'M' };
  }
  if (upScore >= leftScore) {
    return { score: upScore, direction: 'Y' };
  }
  return { score: leftScore, direction: 'X' };
}

/**
 * Finds the optimal ending position allowing for end gaps without penalty.
 */
function findOptimalEndPosition(
  scoreMatrix: number[][],
  lenA: number,
  lenB: number,
): { maxScore: number; endI: number; endJ: number } {
  let maxScore = scoreMatrix[lenA][lenB];
  let endI = lenA;
  let endJ = lenB;

  // Check last row (end gaps in seqA)
  for (let j = 0; j <= lenB; j++) {
    if (scoreMatrix[lenA][j] > maxScore) {
      maxScore = scoreMatrix[lenA][j];
      endI = lenA;
      endJ = j;
    }
  }

  // Check last column (end gaps in seqB)
  for (let i = 0; i <= lenA; i++) {
    if (scoreMatrix[i][lenB] > maxScore) {
      maxScore = scoreMatrix[i][lenB];
      endI = i;
      endJ = lenB;
    }
  }

  return { maxScore, endI, endJ };
}

/**
 * Performs traceback to construct aligned sequences.
 */
function performTraceback(
  tracebackMatrix: TracebackDirection[][],
  seqA: string,
  seqB: string,
  endI: number,
  endJ: number,
): { alignedA: string; alignedB: string } {
  let alignedA = '';
  let alignedB = '';
  let i = endI;
  let j = endJ;

  // Traceback from end to beginning
  while (i > 0 || j > 0) {
    const direction = tracebackMatrix[i][j];

    if (direction === 'M' && i > 0 && j > 0) {
      alignedA = seqA[i - 1] + alignedA;
      alignedB = seqB[j - 1] + alignedB;
      i--;
      j--;
    } else if (direction === 'Y' && i > 0) {
      alignedA = seqA[i - 1] + alignedA;
      alignedB = '-' + alignedB;
      i--;
    } else if (direction === 'X' && j > 0) {
      alignedA = '-' + alignedA;
      alignedB = seqB[j - 1] + alignedB;
      j--;
    } else if (i > 0) {
      // Fallback: gap in seqB
      alignedA = seqA[i - 1] + alignedA;
      alignedB = '-' + alignedB;
      i--;
    } else if (j > 0) {
      // Fallback: gap in seqA
      alignedA = '-' + alignedA;
      alignedB = seqB[j - 1] + alignedB;
      j--;
    } else {
      break;
    }
  }

  // Add leading gaps (if we started after position 0,0)
  while (i > 0) {
    alignedA = seqA[i - 1] + alignedA;
    alignedB = '-' + alignedB;
    i--;
  }
  while (j > 0) {
    alignedA = '-' + alignedA;
    alignedB = seqB[j - 1] + alignedB;
    j--;
  }

  // Add trailing gaps (if we ended before full length)
  for (let k = endI; k < seqA.length; k++) {
    alignedA = alignedA + seqA[k];
    alignedB = alignedB + '-';
  }
  for (let k = endJ; k < seqB.length; k++) {
    alignedA = alignedA + '-';
    alignedB = alignedB + seqB[k];
  }

  return { alignedA, alignedB };
}

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

  // Preprocess sequences for EMBOSS compatibility
  const processedA = preprocessSequence(seqA);
  const processedB = preprocessSequence(seqB);

  const lenA = processedA.length;
  const lenB = processedB.length;

  // Initialize alignment matrices
  const { score: scoreMatrix, traceback: tracebackMatrix } = initializeMatrices(lenA, lenB);

  // Fill matrices using dynamic programming
  for (let i = 1; i <= lenA; i++) {
    for (let j = 1; j <= lenB; j++) {
      const matchScore = getEdnafullScore(processedA[i - 1], processedB[j - 1]);
      const { score, direction } = calculateCellScore(
        matchScore,
        scoreMatrix,
        tracebackMatrix,
        i,
        j,
        gapOpen,
        gapExtend,
        gapAnchor,
      );
      scoreMatrix[i][j] = score;
      tracebackMatrix[i][j] = direction;
    }
  }

  // Find optimal end position (allowing end gaps)
  const { maxScore, endI, endJ } = findOptimalEndPosition(scoreMatrix, lenA, lenB);

  // Traceback to construct aligned sequences
  const { alignedA, alignedB } = performTraceback(tracebackMatrix, processedA, processedB, endI, endJ);

  return {
    seqA: alignedA,
    seqB: alignedB,
    score: maxScore,
  };
}
