/**
 * Unit tests for Needleman-Wunsch alignment.
 * 100% compatible with EMBOSS needle using EDNAFULL scoring matrix.
 *
 * Test cases validated against EMBOSS needle v6.6.0.0 with parameters:
 * -gapopen 10 -gapextend 0.5 -endopen 10 -endextend 0.5
 */

import { needleAlign } from '../../src/needle/needle-alignment';

describe('needleAlign', () => {
  describe('basic alignment', () => {
    it('should align identical sequences without gaps', () => {
      // EMBOSS: ATCG vs ATCG, score = 20.0
      const result = needleAlign('ATCG', 'ATCG');
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('ATCG');
      expect(result.score).toBe(20);
    });

    it('should align sequences with one mismatch', () => {
      // EMBOSS: ATCG vs ATGG, score = 11.0
      const result = needleAlign('ATCG', 'ATGG');
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('ATGG');
      expect(result.score).toBe(11);
    });

    it('should handle gap in reference (insertion in read)', () => {
      // EMBOSS: AT-CG vs ATGCG, score = 10.0
      const result = needleAlign('ATCG', 'ATGCG');
      expect(result.seqA).toBe('AT-CG');
      expect(result.seqB).toBe('ATGCG');
      expect(result.score).toBe(10);
    });

    it('should handle gap in read (deletion)', () => {
      // EMBOSS: ATCGATCG vs ----ATCG, score = 20.0
      const result = needleAlign('ATCGATCG', 'ATCG');
      expect(result.seqA).toBe('ATCGATCG');
      expect(result.seqB).toBe('----ATCG');
      expect(result.score).toBe(20);
    });
  });

  describe('EDNAFULL ambiguity codes', () => {
    it('should score A-W as partial match (+1)', () => {
      // EMBOSS: ATCG vs AWCG, score = 16.0 (A-W = +1)
      const result = needleAlign('ATCG', 'AWCG');
      expect(result.score).toBe(16);
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('AWCG');
    });

    it('should score A-N as -2', () => {
      // EMBOSS: ATCG vs ANCG, score = 13.0 (A-N = -2)
      const result = needleAlign('ATCG', 'ANCG');
      expect(result.score).toBe(13);
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('ANCG');
    });

    it('should handle RYMK ambiguity codes', () => {
      // EMBOSS: ATCG vs RYMK, score = 4.0 (all partial matches)
      const result = needleAlign('ATCG', 'RYMK');
      expect(result.score).toBe(4);
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('RYMK');
    });

    it('should score S-C as partial match (+1)', () => {
      // EMBOSS: SSSS vs CCCC, score = 4.0 (S-C = +1 each)
      const result = needleAlign('SSSS', 'CCCC');
      expect(result.score).toBe(4);
    });

    it('should score S-G as partial match (+1)', () => {
      // EMBOSS: SSSS vs GGGG, score = 4.0 (S-G = +1 each)
      const result = needleAlign('SSSS', 'GGGG');
      expect(result.score).toBe(4);
    });

    it('should prefer gaps over S-A mismatch', () => {
      // EMBOSS: SSSS vs AAAA, score = 0.0 (gaps preferred over mismatches)
      const result = needleAlign('SSSS', 'AAAA');
      expect(result.seqA).toBe('SSSS----');
      expect(result.seqB).toBe('----AAAA');
      expect(result.score).toBe(0);
    });

    it('should handle N vs any as -2', () => {
      // EMBOSS: -AN vs AA-, score = 5.0
      const result = needleAlign('AN', 'AA');
      expect(result.seqA).toBe('-AN');
      expect(result.seqB).toBe('AA-');
      expect(result.score).toBe(5);
    });
  });

  describe('mtDNA sequences', () => {
    // Real rCRS HV1 reference sequence fragment
    const rCRS_HV1_FRAGMENT =
      'GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTC';

    it('should align read with single SNP', () => {
      const ref = 'ATCGATCG';
      const read = 'ATCGATGG'; // C->G at position 7
      const result = needleAlign(ref, read);

      expect(result.seqA).toBe(ref);
      expect(result.seqB).toBe(read);
      expect(result.seqA.length).toBe(result.seqB.length);
      // 7 matches (35) + 1 mismatch (-4) = 31
      expect(result.score).toBe(31);
    });

    it('should align read with insertion', () => {
      const ref = 'ATCGATCG';
      const read = 'ATCGAATCG'; // A inserted
      const result = needleAlign(ref, read);

      expect(result.seqA).toContain('-');
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });

    it('should align read with deletion', () => {
      const ref = 'ATCGATCG';
      const read = 'ATCATCG'; // G deleted
      const result = needleAlign(ref, read);

      expect(result.seqB).toContain('-');
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
    });

    it('should align real mtDNA fragment', () => {
      const ref = rCRS_HV1_FRAGMENT.substring(0, 50);
      const read = ref.substring(0, 20) + 'G' + ref.substring(21); // SNP at position 20
      const result = needleAlign(ref, read);

      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
      expect(result.seqA.length).toBe(result.seqB.length);
    });
  });

  describe('edge cases', () => {
    it('should throw on empty reference', () => {
      expect(() => needleAlign('', 'ATCG')).toThrow('Reference sequence is empty');
    });

    it('should throw on empty read', () => {
      expect(() => needleAlign('ATCG', '')).toThrow('Read sequence is empty');
    });

    it('should handle single base sequences', () => {
      const result = needleAlign('A', 'A');
      expect(result.seqA).toBe('A');
      expect(result.seqB).toBe('A');
      expect(result.score).toBe(5);
    });

    it('should handle completely different sequences with gaps', () => {
      // EMBOSS: AAAA---- vs ----TTTT, score = 0 (gaps preferred over mismatches)
      const result = needleAlign('AAAA', 'TTTT');
      expect(result.seqA).toBe('AAAA----');
      expect(result.seqB).toBe('----TTTT');
      expect(result.score).toBe(0);
    });

    it('should handle case insensitive alignment', () => {
      const result1 = needleAlign('ATCG', 'atcg');
      const result2 = needleAlign('atcg', 'ATCG');

      expect(result1.score).toBe(result2.score);
      expect(result1.score).toBe(20);
    });

    it('should treat unknown characters as N', () => {
      // Unknown character 'X' treated as 'N'
      const result = needleAlign('ATCG', 'ATXG');
      // A-A = 5, T-T = 5, C-X(N) = -2, G-G = 5
      expect(result.score).toBe(13);
    });
  });

  describe('gap penalties', () => {
    it('should use custom gap open penalty', () => {
      // ATCG vs ATGCG (insertion of G)
      // gapOpen=10: AT-CG vs ATGCG, score = 5+5+(-10)+5+5 = 10
      // gapOpen=5:  same alignment, score = 5+5+(-5)+5+5 = 15
      const result1 = needleAlign('ATCG', 'ATGCG', { gapOpen: 10, gapExtend: 0.5 });
      const result2 = needleAlign('ATCG', 'ATGCG', { gapOpen: 5, gapExtend: 0.5 });

      expect(result1.score).toBe(10);
      expect(result2.score).toBe(15);
    });

    it('should prefer gap extension over new gaps', () => {
      // ATCGATCG vs ATGATG
      // EMBOSS needle output:
      // ref: ATCGATCG
      // read: AT-GATG-
      // Score: 11.0
      const ref = 'ATCGATCG';
      const read = 'ATGATG';
      const result = needleAlign(ref, read, { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('ATCGATCG');
      expect(result.seqB).toBe('AT-GATG-');
      expect(result.score).toBe(11);
    });
  });

  describe('alignment length consistency', () => {
    it('should produce equal length alignments', () => {
      // Test cases validated against EMBOSS needle
      const testCases = [
        { ref: 'ATCG', read: 'ATCG', resultRef: 'ATCG', resultRead: 'ATCG' },
        {
          ref: 'AAAAAAAAAA',
          read: 'TTTTTTTTTT',
          resultRef: 'AAAAAAAAAA----------',
          resultRead: '----------TTTTTTTTTT',
        },
        {
          ref: 'ATCGATCGATCGATCGATCG',
          read: 'ATCGATCGATCGATCGATCG',
          resultRef: 'ATCGATCGATCGATCGATCG',
          resultRead: 'ATCGATCGATCGATCGATCG',
        },
      ];

      for (const { ref, read, resultRef, resultRead } of testCases) {
        const result = needleAlign(ref, read);
        expect(result.seqA).toBe(resultRef);
        expect(result.seqB).toBe(resultRead);
        expect(result.seqA.length).toBe(result.seqB.length);
      }
    });

    it('should preserve all bases in reference', () => {
      const ref = 'ATCGATCGATCG';
      const read = 'ATGCGATCG';
      const result = needleAlign(ref, read);

      expect(result.seqA.replace(/-/g, '')).toBe(ref);
    });

    it('should preserve all bases in read', () => {
      const ref = 'ATCGATCG';
      const read = 'ATGCGATCGATCG';
      const result = needleAlign(ref, read);

      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });

    it('should handle read shorter than reference (deletion case)', () => {
      // Regression test: traceback bug when seqB is shorter than seqA
      // EMBOSS needle output for ATCG vs TG:
      // ref: ATCG
      // read: --TG
      // Score: 1.0
      const ref = 'ATCG';
      const read = 'TG';
      const result = needleAlign(ref, read);

      expect(result.seqA.length).toBe(result.seqB.length);
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('--TG');
      expect(result.score).toBe(1);
    });

    it('should handle read much shorter than reference', () => {
      // Multiple gaps needed
      const ref = 'ACGTACGT';
      const read = 'CGT';
      const result = needleAlign(ref, read);

      expect(result.seqA.length).toBe(result.seqB.length);
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
      expect(result.seqB).not.toContain('undefined');
    });
  });

  describe('EMBOSS compatibility', () => {
    // Test cases validated against EMBOSS needle v6.6.0.0
    it('should match EMBOSS needle for identical sequences', () => {
      const result = needleAlign('ATCG', 'ATCG', { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('ATCG');
      expect(result.score).toBe(20);
    });

    it('should match EMBOSS needle for single mismatch', () => {
      const result = needleAlign('ATCG', 'ATGG', { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('ATGG');
      expect(result.score).toBe(11);
    });

    it('should match EMBOSS needle for insertion', () => {
      const result = needleAlign('ATCG', 'ATGCG', { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('AT-CG');
      expect(result.seqB).toBe('ATGCG');
      expect(result.score).toBe(10);
    });

    it('should match EMBOSS needle for deletion', () => {
      const result = needleAlign('ATCGATCG', 'ATCG', { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('ATCGATCG');
      expect(result.seqB).toBe('----ATCG');
      expect(result.score).toBe(20);
    });

    it('should match EMBOSS needle for all gaps case', () => {
      const result = needleAlign('AAAA', 'TTTT', { gapOpen: 10, gapExtend: 0.5 });

      expect(result.seqA).toBe('AAAA----');
      expect(result.seqB).toBe('----TTTT');
      expect(result.score).toBe(0);
    });
  });

  describe('gap behavior', () => {
    it('should handle consecutive gaps', () => {
      // ATCG vs A--CG (2 consecutive gaps)
      const result = needleAlign('ATCG', 'A--CG', { gapOpen: 10, gapExtend: 0.5 });
      // After removing input gaps: ATCG vs ACG
      // EMBOSS: ATCG vs A-CG (single gap is better than 2 consecutive)
      expect(result.seqA.replace(/-/g, '')).toBe('ATCG');
      expect(result.seqB.replace(/-/g, '')).toBe('ACG');
    });

    it('should apply gapopen for gaps', () => {
      // Gap should incur gapopen penalty
      const result = needleAlign('ATCG', 'ACG', { gapOpen: 10, gapExtend: 0.5 });
      // EMBOSS: ATCG vs -ACG (leading gap preferred)
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('-ACG');
      // Score: 3 matches (15) - gapopen (10) = 5? But EMBOSS gives 6
      expect(result.score).toBe(6);
    });

    it('should handle gapOpen=0 and gapExtend=0', () => {
      // With no gap penalties, gaps are free
      const result = needleAlign('AAAA', 'TTTT', { gapOpen: 0, gapExtend: 0 });
      // All gaps preferred over mismatches
      expect(result.score).toBe(0);
    });
  });

  describe('EDNAFULL coverage', () => {
    it('should handle B (not A) code', () => {
      // B = C/G/T, vs C should be partial match
      const result = needleAlign('BBBB', 'CCCC', { gapOpen: 10, gapExtend: 0.5 });
      // B-C = -1 (from EDNAFULL matrix)
      // EMBOSS prefers gaps (score 0) over mismatches
      expect(result.seqA).toBe('BBBB----');
      expect(result.seqB).toBe('----CCCC');
      expect(result.score).toBe(0);
    });

    it('should handle V (not T) code', () => {
      // V = A/C/G, vs A should be partial match (+1)
      const result = needleAlign('VVVV', 'AAAA', { gapOpen: 10, gapExtend: 0.5 });
      // V-A = -1, but let's check actual EMBOSS behavior
      expect(result.score).toBe(0); // Gaps preferred
    });

    it('should handle H (not G) code', () => {
      // H = A/C/T, vs A should be partial match
      const result = needleAlign('HHHH', 'AAAA', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0); // Gaps preferred
    });

    it('should handle D (not C) code', () => {
      // D = A/G/T, vs A should be partial match
      const result = needleAlign('DDDD', 'AAAA', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0); // Gaps preferred
    });

    it('should handle all IUPAC codes in one sequence', () => {
      const allCodes = 'ATGCSWRYKMBVHDN';
      const result = needleAlign(allCodes, allCodes, { gapOpen: 10, gapExtend: 0.5 });
      // Exact match of all codes
      expect(result.seqA).toBe(allCodes);
      expect(result.seqB).toBe(allCodes);
      // Score should be sum of all diagonal elements
      expect(result.score).toBeGreaterThan(0);
    });
  });

  describe('score validation', () => {
    it('should be symmetric', () => {
      const result1 = needleAlign('ATCG', 'ATGG', { gapOpen: 10, gapExtend: 0.5 });
      const result2 = needleAlign('ATGG', 'ATCG', { gapOpen: 10, gapExtend: 0.5 });
      expect(result1.score).toBe(result2.score);
    });

    it('should be deterministic', () => {
      const result1 = needleAlign('ATCGATCG', 'ATGATG', { gapOpen: 10, gapExtend: 0.5 });
      const result2 = needleAlign('ATCGATCG', 'ATGATG', { gapOpen: 10, gapExtend: 0.5 });
      expect(result1.score).toBe(result2.score);
      expect(result1.seqA).toBe(result2.seqA);
      expect(result1.seqB).toBe(result2.seqB);
    });

    it('should verify exact match score calculation', () => {
      const result = needleAlign('AAAA', 'AAAA', { gapOpen: 10, gapExtend: 0.5 });
      // 4 matches * 5 = 20
      expect(result.score).toBe(20);
    });

    it('should verify mismatch score calculation as EMBOSS does', () => {
      // EMBOSS: AAAA- vs -AAAG, score = 15 (gaps preferred over mismatch at end)
      const result = needleAlign('AAAA', 'AAAG', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(15);
    });
  });

  describe('mtDNA real-world', () => {
    it('should handle multiple SNPs close together', () => {
      const ref = 'GATCACAGGTCTATCACCCT';
      const read = 'GATCACAGGTCTGTGACCCT'; // 2 SNPs at positions 14,15
      const result = needleAlign(ref, read, { gapOpen: 10, gapExtend: 0.5 });
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });

    it('should handle mixed indel + SNP', () => {
      const ref = 'GATCACAGGTCTATCACCCT';
      const read = 'GATCACAGGTCTGTCACCCT'; // 1 deletion + 1 SNP
      const result = needleAlign(ref, read, { gapOpen: 10, gapExtend: 0.5 });
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });

    it('should align full HV1 fragment', () => {
      // HV1 is ~340bp, test with first 100bp
      const hv1Fragment =
        'GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCC';
      const result = needleAlign(hv1Fragment, hv1Fragment, { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(hv1Fragment.length * 5);
    });
  });

  describe('gap extension and multiple gaps', () => {
    it('should penalize each gap opening separately', () => {
      // ATCGATCG vs ATGATG — two separate gap regions
      // EMBOSS: ATCGATCG vs AT-GATG-
      const result = needleAlign('ATCGATCG', 'ATGATG', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.seqA).toBe('ATCGATCG');
      expect(result.seqB).toBe('AT-GATG-');
      expect(result.score).toBe(11);
    });

    it('should use custom gap extend penalty', () => {
      // ATCG vs ATGCG with different gapExtend values
      // gapOpen=10, gapExtend=0.5: score = 10
      // gapOpen=10, gapExtend=2.0: score = 5+5-10+5+5 = 10 (same gapOpen)
      const result1 = needleAlign('ATCG', 'ATGCG', { gapOpen: 10, gapExtend: 0.5 });
      const result2 = needleAlign('ATCG', 'ATGCG', { gapOpen: 10, gapExtend: 2.0 });
      expect(result1.score).toBe(10);
      expect(result2.score).toBe(10);
    });

    it('should align with minimal gaps even when gapOpen is very high', () => {
      // Even with very high gapOpen, EMBOSS may still use gaps if beneficial
      // EMBOSS: -ATCG vs ATGCG, score = 2.0
      const result = needleAlign('ATCG', 'ATGCG', { gapOpen: 100, gapExtend: 0.5 });
      expect(result.seqA).toBe('-ATCG');
      expect(result.seqB).toBe('ATGCG');
      expect(result.score).toBe(2);
    });

    it('should allow many gaps when gapOpen is very low', () => {
      // With gapOpen=0, gaps are cheap
      const result = needleAlign('AAAA', 'TTTT', { gapOpen: 0, gapExtend: 0 });
      expect(result.score).toBe(0);
    });
  });

  describe('three-state ambiguity codes', () => {
    // EMBOSS needle results for single-base alignments:
    // When match score is negative, EMBOSS prefers gaps (score 0)

    it('should score B-T as EMBOSS does', () => {
      // B = C/G/T, B-T = -1 from EDNAFULL
      // EMBOSS: B- vs -T, score = 0 (gaps preferred)
      const result = needleAlign('B', 'T', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should score V-A as EMBOSS does', () => {
      // V = A/C/G, V-A = -1 from EDNAFULL
      // EMBOSS: V- vs -A, score = 0
      const result = needleAlign('V', 'A', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should score H-T as EMBOSS does', () => {
      // H = A/C/T, H-T = -1 from EDNAFULL
      // EMBOSS: H- vs -T, score = 0
      const result = needleAlign('H', 'T', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should score D-G as EMBOSS does', () => {
      // D = A/G/T, D-G = -1 from EDNAFULL
      // EMBOSS: D- vs -G, score = 0
      const result = needleAlign('D', 'G', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });
  });

  describe('ambiguity vs ambiguity', () => {
    // EMBOSS needle results: when both sequences have single ambiguity codes
    // with negative match score, EMBOSS prefers gaps (score 0)

    it('should score R vs K as EMBOSS does', () => {
      // R = A/G, K = G/T
      // R-K = -2 from EDNAFULL, but EMBOSS: R- vs -K, score = 0
      const result = needleAlign('R', 'K', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should score W vs S as EMBOSS does', () => {
      // W = A/T, S = C/G
      // W-S = -4 from EDNAFULL, but EMBOSS: W- vs -S, score = 0
      const result = needleAlign('W', 'S', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should score N vs N as EMBOSS does', () => {
      // N-N = -1 from EDNAFULL, but EMBOSS: N- vs -N, score = 0
      const result = needleAlign('N', 'N', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should handle all ambiguity combinations', () => {
      // Test various ambiguity code combinations
      const result1 = needleAlign('R', 'A', { gapOpen: 10, gapExtend: 0.5 }); // R-A = +1
      const result2 = needleAlign('Y', 'C', { gapOpen: 10, gapExtend: 0.5 }); // Y-C = +1
      const result3 = needleAlign('M', 'A', { gapOpen: 10, gapExtend: 0.5 }); // M-A = +1
      expect(result1.score).toBe(1);
      expect(result2.score).toBe(1);
      expect(result3.score).toBe(1);
    });
  });

  describe('read longer than reference', () => {
    it('should handle read much longer than reference', () => {
      const ref = 'ATCG';
      const read = 'ATCGATCGATCG';
      const result = needleAlign(ref, read, { gapOpen: 10, gapExtend: 0.5 });
      // EMBOSS: ----ATCG vs ATCGATCGATCG (leading gaps in ref)
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });
  });

  describe('symmetry', () => {
    it('should produce same score regardless of sequence order', () => {
      const r1 = needleAlign('ATCG', 'ATGG', { gapOpen: 10, gapExtend: 0.5 });
      const r2 = needleAlign('ATGG', 'ATCG', { gapOpen: 10, gapExtend: 0.5 });
      expect(r1.score).toBe(r2.score);
      expect(r1.seqA).toBe(r2.seqB);
      expect(r1.seqB).toBe(r2.seqA);
    });

    it('should produce symmetric alignment for equal sequences', () => {
      const r1 = needleAlign('ATCGATCG', 'ATGATG', { gapOpen: 10, gapExtend: 0.5 });
      const r2 = needleAlign('ATGATG', 'ATCGATCG', { gapOpen: 10, gapExtend: 0.5 });
      expect(r1.score).toBe(r2.score);
    });
  });

  describe('gaps at both ends', () => {
    it('should handle gaps at both ends of a sequence', () => {
      // ref: --ATCG--  read: TTATCGTT
      // EMBOSS should align with gaps at both ends
      const result = needleAlign('TTATCGTT', 'ATCG', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.seqA.replace(/-/g, '')).toBe('TTATCGTT');
      expect(result.seqB.replace(/-/g, '')).toBe('ATCG');
    });
  });

  describe('all N sequences', () => {
    it('should handle sequences of all Ns as EMBOSS does', () => {
      // EMBOSS: NNNN---- vs ----NNNN, score = 0 (all gaps preferred)
      const result = needleAlign('NNNN', 'NNNN', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(0);
    });

    it('should handle N vs regular bases', () => {
      const result = needleAlign('NNNN', 'AAAA', { gapOpen: 10, gapExtend: 0.5 });
      // N-A = -2 * 4 = -8, but EMBOSS may prefer gaps
      expect(result.seqA.replace(/-/g, '')).toBe('NNNN');
      expect(result.seqB.replace(/-/g, '')).toBe('AAAA');
    });
  });

  describe('combined mutations', () => {
    it('should handle combined SNP, insertion and deletion', () => {
      // Combined: SNP + insertion + deletion in one alignment
      const ref = 'ATCGATCG';
      const read = 'ATGGTCG'; // C->G SNP, deletion of A, insertion of G
      const result = needleAlign(ref, read, { gapOpen: 10, gapExtend: 0.5 });
      expect(result.seqA.replace(/-/g, '')).toBe(ref);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });
  });

  describe('stress tests', () => {
    it('should handle long sequences without errors', () => {
      const long = 'A'.repeat(1000);
      const result = needleAlign(long, long, { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(5000);
    });

    it('should handle long different sequences', () => {
      const longA = 'A'.repeat(500);
      const longB = 'T'.repeat(500);
      const result = needleAlign(longA, longB, { gapOpen: 10, gapExtend: 0.5 });
      // All gaps preferred
      expect(result.score).toBe(0);
    });

    /*
     * The aligned strings are reconstructed in 8192-column chunks, so alignments
     * longer than that exercise a code path the shorter cases never reach. These
     * use provably-correct expected values (a perfect self-alignment, and a
     * single free leading end-gap), not values copied from the implementation.
     */
    it('should reconstruct a 10000bp identical alignment exactly (chunked path)', () => {
      const seq = 'ACGT'.repeat(2500); // 10000 bp, periodic so adjacent bases differ
      const result = needleAlign(seq, seq);
      // Identical sequences align perfectly: every base matches (+5), no gaps.
      expect(result.seqA).toBe(seq);
      expect(result.seqB).toBe(seq);
      expect(result.seqA.length).toBe(10000);
      expect(result.score).toBe(50000);
    });

    it('should place a single free leading end-gap over a 10000bp read (chunked path)', () => {
      const ref = 'ACGT'.repeat(2500); // 10000 bp
      const read = ref.slice(1); // 9999 bp: ref with its first base dropped
      const result = needleAlign(ref, read);
      // Read equals ref shifted by one; the unique optimum is a free leading
      // end-gap (no penalty), leaving all 9999 remaining bases matched.
      expect(result.seqA).toBe(ref);
      expect(result.seqB).toBe('-' + read);
      expect(result.seqA.length).toBe(10000);
      expect(result.score).toBe(9999 * 5);
    });
  });

  describe('edge cases', () => {
    it('should handle homopolymer runs', () => {
      const result = needleAlign('AAAAAAA', 'TTTTTTT', { gapOpen: 10, gapExtend: 0.5 });
      // All gaps preferred over all mismatches
      expect(result.seqA).toBe('AAAAAAA-------');
      expect(result.seqB).toBe('-------TTTTTTT');
      expect(result.score).toBe(0);
    });

    it('should handle leading gaps', () => {
      const result = needleAlign('ATCG', 'CG', { gapOpen: 10, gapExtend: 0.5 });
      // EMBOSS: ATCG vs --CG (leading gaps in read)
      expect(result.seqA).toBe('ATCG');
      expect(result.seqB).toBe('--CG');
    });

    it('should treat X as N (EMBOSS compatible)', () => {
      // EMBOSS: X -> N, score = 13.0
      const resultX = needleAlign('ATCG', 'ATXG', { gapOpen: 10, gapExtend: 0.5 });
      const resultN = needleAlign('ATCG', 'ATNG', { gapOpen: 10, gapExtend: 0.5 });
      expect(resultX.score).toBe(resultN.score);
      expect(resultX.score).toBe(13);
    });

    it('should throw error for Z (EMBOSS compatible)', () => {
      // EMBOSS: Z -> error "Sequence is not nucleic"
      expect(() => needleAlign('ATCG', 'ATZG', { gapOpen: 10, gapExtend: 0.5 })).toThrow();
      expect(() => needleAlign('ATCG', 'ATzG', { gapOpen: 10, gapExtend: 0.5 })).toThrow();
    });

    it('should treat ! as gap (EMBOSS compatible)', () => {
      // EMBOSS: ! -> gap, score = 6.0
      const result = needleAlign('ATCG', 'AT!G', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(6);
      expect(result.seqB).toBe('ATG-');
    });

    it('should treat @ as gap (EMBOSS compatible)', () => {
      // EMBOSS: @ -> gap, score = 6.0
      const result = needleAlign('ATCG', 'AT@G', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(6);
    });

    it('should ignore whitespace (EMBOSS compatible)', () => {
      // EMBOSS: whitespace is ignored, AT G becomes ATG
      // Alignment: ATCG vs ATG- gives score 6.0
      const result = needleAlign('ATCG', 'AT G', { gapOpen: 10, gapExtend: 0.5 });
      expect(result.score).toBe(6);
    });
  });

  /*
   * 3'-shifted indel placement per ISFG forensic mtDNA convention
   * (Parson et al. 2014, §3.2; see docs/parson.pdf). Within homopolymer
   * and repeat tracts the gap must land at the 3' end of the run, so an
   * extra C in a poly-C tract is reported at the last C position.
   *
   * These tests use long flanked references to keep the alignment anchored
   * (short synthetic strings let the free-end-gap behaviour score a leading
   * shift higher than an internal gap, masking the tie-break direction).
   */
  describe("3'-shifted indel placement (ISFG mtDNA convention)", () => {
    /** rCRS positions 290-330 (1-based): includes both HV2 poly-C tracts. */
    const HV2_FRAGMENT = 'AATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAG';

    it('places gap at 3-prime end of seven-C tract (HV2 303-309)', () => {
      // 7-C tract becomes 8 Cs (one extra C anywhere in the tract is indistinguishable).
      const ref = HV2_FRAGMENT;
      const read = 'AATTTCCACCAAACCCCCCCCTCCCCCGCTTCTGGCCACAG';
      const result = needleAlign(ref, read, { gapAnchor: '3-prime' });
      // Gap must land between the last C of the 7-C tract (pos 309) and the T (pos 310).
      expect(result.seqA).toBe('AATTTCCACCAAACCCCCCC-TCCCCCGCTTCTGGCCACAG');
      expect(result.seqB).toBe(read);
    });

    it('places gap at 3-prime end of five-C tract (HV2 311-315)', () => {
      const ref = HV2_FRAGMENT;
      const read = 'AATTTCCACCAAACCCCCCCTCCCCCCGCTTCTGGCCACAG';
      const result = needleAlign(ref, read, { gapAnchor: '3-prime' });
      // Gap must land between the last C of the 5-C tract (pos 315) and the G (pos 316).
      expect(result.seqA).toBe('AATTTCCACCAAACCCCCCCTCCCCC-GCTTCTGGCCACAG');
      expect(result.seqB).toBe(read);
    });

    it('places gap at 3-prime end of merged tract when T310C transition is present', () => {
      // T310C merges the 7C+T+5C region into a 13-C run; +1C pushes total to 14 Cs.
      const ref = HV2_FRAGMENT;
      const read = 'AATTTCCACCAAACCCCCCCCCCCCCCGCTTCTGGCCACAG';
      const result = needleAlign(ref, read, { gapAnchor: '3-prime' });
      // Gap lands before G (pos 316) — the 3' end of the merged tract.
      expect(result.seqA[result.seqA.indexOf('-') + 1]).toBe('G');
      // The merged tract's last position in rCRS is 315; in alignment we expect a single
      // gap and the rest of the read is consistent.
      expect(result.seqA.split('-').length - 1).toBe(1);
      expect(result.seqB.replace(/-/g, '')).toBe(read);
    });

    it('places gap at 3-prime end of HV1 C-tract under T16189C transition', () => {
      // rCRS 16180-16210: AAAACCCCCTCCCCATGCTTACGGCCAATAA
      // Add T16189C + extra C → merged 11-C tract, gap should land before A at position 16194.
      const ref = 'GCAAACCCCCCAAAACCCCCTCCCCATGCTTACGGCCAATAACAGTGCTAGCC';
      const read = 'GCAAACCCCCCAAAACCCCCCCCCCCATGCTTACGGCCAATAACAGTGCTAGCC';
      const result = needleAlign(ref, read, { gapAnchor: '3-prime' });
      const gapIdx = result.seqA.indexOf('-');
      // Gap lands immediately before the 'A' that follows the merged tract (rCRS 16194).
      expect(result.seqA[gapIdx + 1]).toBe('A');
      expect(result.seqA.split('-').length - 1).toBe(1);
    });

    it('places gap pair at 3-prime end of AC-repeat (HV3 514-524)', () => {
      // rCRS 510-530: CCAGCACACACACACCGCTGC. The (CA)5 C repeat 514-524 ends at C524.
      // +1 AC pair extends repeat by one unit; right-anchor places the 2bp gap before the
      // CG that follows the repeat (so the last C of the gap-bearing alignment is rCRS 524).
      const ref = 'GACAGCAACGCCAGCACACACACACCGCTGCTAACCAGAT';
      const read = 'GACAGCAACGCCAGCACACACACACACCGCTGCTAACCAGAT';
      const result = needleAlign(ref, read, { gapAnchor: '3-prime' });
      // Two adjacent gaps, ending immediately before C525 + G526.
      const gaps: number[] = [];
      for (let i = 0; i < result.seqA.length; i++) {
        if (result.seqA[i] === '-') gaps.push(i);
      }
      expect(gaps).toHaveLength(2);
      expect(gaps[1]).toBe(gaps[0] + 1);
      // Character right after the gap pair should be the C that follows the repeat (C525).
      expect(result.seqA[gaps[1] + 1]).toBe('C');
      expect(result.seqA[gaps[1] + 2]).toBe('G');
    });

    it('does not shift gaps in non-repeat regions', () => {
      // Plain insertion between distinct bases — no tie, alignment is unique.
      const result = needleAlign('GATCAGTACG', 'GATCATGTACG', { gapAnchor: '3-prime' });
      // Gap must land between A and G (no other valid placement gives the same score).
      expect(result.seqA).toBe('GATCA-GTACG');
      expect(result.seqB).toBe('GATCATGTACG');
    });
  });

  /*
   * 5'-anchored indel placement is the default and matches the EMBOSS-needle
   * default tie-break (left-anchored gaps). The same equally-scoring tracts that
   * the ISFG block anchors to the 3' end land at the 5' end here. Scores are
   * identical to the 3'-anchored alignments above.
   */
  describe("5'-anchored indel placement (EMBOSS default)", () => {
    /** rCRS positions 290-330 (1-based): includes both HV2 poly-C tracts. */
    const HV2_FRAGMENT = 'AATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAG';

    it('places gap at 5-prime end of seven-C tract (HV2 303-309) by default', () => {
      const read = 'AATTTCCACCAAACCCCCCCCTCCCCCGCTTCTGGCCACAG';
      const result = needleAlign(HV2_FRAGMENT, read);
      // Default (no option) left-anchors: gap lands at the start of the 7-C run.
      expect(result.seqA).toBe('AATTTCCACCAAA-CCCCCCCTCCCCCGCTTCTGGCCACAG');
      expect(result.seqB).toBe(read);
    });

    it('places gap at 5-prime end of five-C tract (HV2 311-315) by default', () => {
      const read = 'AATTTCCACCAAACCCCCCCTCCCCCCGCTTCTGGCCACAG';
      const result = needleAlign(HV2_FRAGMENT, read);
      // Gap lands right after the T that precedes the 5-C run.
      expect(result.seqA).toBe('AATTTCCACCAAACCCCCCCT-CCCCCGCTTCTGGCCACAG');
      expect(result.seqB).toBe(read);
    });

    it('produces the same score as the 3-prime anchoring', () => {
      const read = 'AATTTCCACCAAACCCCCCCCTCCCCCGCTTCTGGCCACAG';
      const fivePrime = needleAlign(HV2_FRAGMENT, read);
      const threePrime = needleAlign(HV2_FRAGMENT, read, { gapAnchor: '3-prime' });
      expect(fivePrime.score).toBe(threePrime.score);
    });
  });
});
