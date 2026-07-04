/** IUPAC-aware motif search over a nucleotide sequence (primer/motif locator). */

import { IUPAC_NUCLEOTIDES_LIST } from './iupac.utils';
import { reverseComplement } from './sequence.utils';

/** A merged hit of a motif on a sequence, in forward-sequence coordinates (half-open `[start, end)`). */
export interface MotifMatch {
  /** 0-based index of the first matched base. */
  start: number;
  /** 0-based index one past the last matched base. */
  end: number;
  /** Whether the hit is on the reverse-complement strand (query annealing to the reverse strand). */
  isReverseComplement: boolean;
  /** Number of raw overlapping hits merged into this span (≥1; >1 only for periodic motifs). */
  count: number;
}

export interface FindMotifOptions {
  /** Also search the reverse-complement strand (default: forward strand only). */
  searchReverseComplement?: boolean;
}

/** Tells whether `query` is a non-empty string of IUPAC nucleotide codes (case-insensitive). */
export function isIupacMotif(query: string): boolean {
  if (query.length === 0) return false;
  for (const ch of query.toUpperCase()) {
    if (!(ch in IUPAC_NUCLEOTIDES_LIST)) return false;
  }
  return true;
}

/** Compiles a query into a per-position set of the bases each IUPAC code admits. Throws on an unknown code. */
function compileMotif(query: string): Set<string>[] {
  return [...query.toUpperCase()].map(ch => {
    const bases = IUPAC_NUCLEOTIDES_LIST[ch];
    if (!bases) throw new Error(`Invalid IUPAC code: ${ch}`);
    return new Set(bases);
  });
}

/**
 * Tells whether a read base matches an allowed (query) set. Subset semantics: the query may be degenerate,
 * but an ambiguous read base (e.g. `N` = A/C/G/T) does NOT match a more specific query position — a locator
 * must not report unresolved bases as hits. A concrete read base matches iff the query admits it.
 */
function baseMatches(base: string, allowed: Set<string>): boolean {
  const bases = IUPAC_NUCLEOTIDES_LIST[base.toUpperCase()];
  if (!bases) return false; // Gaps and unknown symbols match nothing.
  for (const b of bases) {
    if (!allowed.has(b)) return false;
  }
  return true;
}

/** Collects every start index where the compiled motif matches, allowing overlaps. */
function scanStarts(sequence: string, motif: Set<string>[]): number[] {
  const starts: number[] = [];
  const last = sequence.length - motif.length;
  for (let i = 0; i <= last; i++) {
    let hit = true;
    for (let j = 0; j < motif.length; j++) {
      if (!baseMatches(sequence[i + j], motif[j])) {
        hit = false;
        break;
      }
    }
    if (hit) starts.push(i);
  }
  return starts;
}

/**
 * Merges overlapping hits (of one strand) into single spans; adjacent-but-not-overlapping hits stay
 * separate. Overlap only happens for periodic motifs, so a poly-A run collapses to one span while two
 * back-to-back distinct sites remain two matches. `starts` must be ascending.
 */
function mergeOverlapping(starts: number[], width: number, isReverseComplement: boolean): MotifMatch[] {
  const matches: MotifMatch[] = [];
  for (const start of starts) {
    const prev = matches[matches.length - 1];
    if (prev && start < prev.end) {
      prev.end = start + width;
      prev.count += 1;
    } else {
      matches.push({ start, end: start + width, isReverseComplement, count: 1 });
    }
  }
  return matches;
}

/**
 * Finds every occurrence of an IUPAC motif in `sequence` (forward strand, plus the reverse-complement
 * strand when {@link FindMotifOptions.searchReverseComplement} is set). Matching is exact per position
 * (no mismatch budget) and fixed-length. Overlapping hits of the same strand are merged into one span
 * (see {@link mergeOverlapping}). Results are in forward-sequence coordinates, sorted by `start` then
 * forward-before-reverse. Throws when `query` contains a non-IUPAC character.
 */
export function findMotif(sequence: string, query: string, options: FindMotifOptions = {}): MotifMatch[] {
  // Compile (and thereby validate the IUPAC codes) before any length short-circuit, so an invalid
  // query throws even when it is longer than the sequence — matching the documented contract.
  const forwardMotif = compileMotif(query);
  if (query.length === 0 || query.length > sequence.length) return [];

  const forward = mergeOverlapping(scanStarts(sequence, forwardMotif), query.length, false);
  if (!options.searchReverseComplement) return forward;

  const reverse = mergeOverlapping(scanStarts(sequence, compileMotif(reverseComplement(query))), query.length, true);
  return [...forward, ...reverse].sort(
    (a, b) => a.start - b.start || Number(a.isReverseComplement) - Number(b.isReverseComplement),
  );
}
