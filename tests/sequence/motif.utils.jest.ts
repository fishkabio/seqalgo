import { findMotif, isIupacMotif, MotifMatch } from '../../src/sequence/motif.utils';

/** A forward-strand match, for terse expectations. */
function fwd(start: number, end: number, count = 1): MotifMatch {
  return { start, end, isReverseComplement: false, count };
}

/** A reverse-complement-strand match. */
function rev(start: number, end: number, count = 1): MotifMatch {
  return { start, end, isReverseComplement: true, count };
}

describe('isIupacMotif', () => {
  it('accepts non-empty IUPAC strings, case-insensitively', () => {
    expect(isIupacMotif('ACGTN')).toBe(true);
    expect(isIupacMotif('acgt')).toBe(true);
    expect(isIupacMotif('RYSWKMBDHVN')).toBe(true);
  });

  it('rejects an empty string and non-IUPAC characters', () => {
    expect(isIupacMotif('')).toBe(false);
    expect(isIupacMotif('AU')).toBe(false); // U is not a DNA IUPAC code.
    expect(isIupacMotif('AC GT')).toBe(false);
  });
});

describe('findMotif — forward strand', () => {
  it('finds a single exact occurrence in forward coordinates', () => {
    expect(findMotif('AACGTAA', 'CGT')).toEqual([fwd(2, 5)]);
  });

  it('keeps back-to-back non-overlapping sites separate', () => {
    expect(findMotif('GAATTCGAATTC', 'GAATTC')).toEqual([fwd(0, 6), fwd(6, 12)]);
  });

  it('matches IUPAC degeneracy in the query (N matches any base)', () => {
    expect(findMotif('GGAACC', 'GGNNCC')).toEqual([fwd(0, 6)]);
  });

  it('matches a two-base ambiguity code (R = A/G)', () => {
    expect(findMotif('ATAGA', 'AR')).toEqual([fwd(2, 4)]);
  });

  it('normalizes a lowercase (soft-masked) read to uppercase', () => {
    expect(findMotif('acgt', 'ACGT')).toEqual([fwd(0, 4)]);
  });

  it('never matches across a gap character', () => {
    expect(findMotif('A-C', 'ANC')).toEqual([]);
  });

  it('lets a degenerate query match a concrete read, but not vice-versa', () => {
    // Query N matches any concrete read base.
    expect(findMotif('ACGT', 'ANGT')).toEqual([fwd(0, 4)]);
    // An ambiguous read base (N) does NOT match a concrete query position — no locator hit over unknowns.
    expect(findMotif('ANGT', 'ACGT')).toEqual([]);
    // N in both is fine (query N admits everything the read N could be).
    expect(findMotif('ANGT', 'ANGT')).toEqual([fwd(0, 4)]);
  });
});

describe('findMotif — homopolymer / periodic merging', () => {
  it('collapses overlapping hits of a homopolymer into one span with a count', () => {
    // AAA over AAAAAA hits at starts 0,1,2,3 — all overlap → a single 1..6 span.
    expect(findMotif('AAAAAA', 'AAA')).toEqual([fwd(0, 6, 4)]);
  });

  it('collapses an overlapping dinucleotide repeat', () => {
    // AT over ATATAT hits at 0,2,4 (non-overlapping) — stays three separate matches.
    expect(findMotif('ATATAT', 'AT')).toEqual([fwd(0, 2), fwd(2, 4), fwd(4, 6)]);
    // ATA over ATATAT hits at 0,2 which overlap (0..3 and 2..5) → merges 0..5, then 4? no room.
    expect(findMotif('ATATAT', 'ATA')).toEqual([fwd(0, 5, 2)]);
  });
});

describe('findMotif — reverse complement', () => {
  it('ignores the reverse strand by default', () => {
    // AAC has no forward hit here; its RC (GTT) is present but not searched.
    expect(findMotif('AACGTT', 'AAC')).toEqual([fwd(0, 3)]);
  });

  it('finds reverse-strand hits (as revcomp of the query) in forward coordinates', () => {
    // Query AAC: forward hit at 0..3; RC pattern GTT at 3..6, flagged reverse.
    expect(findMotif('AACGTT', 'AAC', { searchReverseComplement: true })).toEqual([fwd(0, 3), rev(3, 6)]);
  });

  it('merges overlapping reverse-strand hits within the reverse strand only', () => {
    // Query AA: forward none on TTTT; RC pattern TT hits 0,1,2 → one reverse span 0..4.
    expect(findMotif('TTTT', 'AA', { searchReverseComplement: true })).toEqual([rev(0, 4, 3)]);
  });
});

describe('findMotif — edge cases', () => {
  it('returns nothing for an empty query', () => {
    expect(findMotif('ACGT', '')).toEqual([]);
  });

  it('returns nothing when the query is longer than the sequence', () => {
    expect(findMotif('AC', 'ACGT')).toEqual([]);
  });

  it('throws on a non-IUPAC character in the query', () => {
    expect(() => findMotif('ACGT', 'AXG')).toThrow('Invalid IUPAC code: X');
  });

  it('validates the query before the length short-circuit (invalid + too long still throws)', () => {
    expect(() => findMotif('AC', 'AXG')).toThrow('Invalid IUPAC code: X');
  });
});
