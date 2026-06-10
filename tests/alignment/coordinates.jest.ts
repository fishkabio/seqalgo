import {
  calculateAlignmentPosition,
  calculateBiologicalPosition,
  calculateGapRegions,
  calculateReadsCoverageRange,
  isPositionInGap,
} from '../../src/alignment/coordinates';

describe('calculateGapRegions', () => {
  it('returns half-open [start, end) regions for each gap run, shifted by offset', () => {
    // 'AT--CG-': gaps at columns 2,3 (run -> [2,4)) and column 6 (-> [6,7)).
    expect(calculateGapRegions('AT--CG-', 0)).toEqual([
      { start: 2, end: 4 },
      { start: 6, end: 7 },
    ]);
  });

  it('returns no regions when there are no gaps', () => {
    expect(calculateGapRegions('ATCG', 0)).toEqual([]);
  });

  it('treats an all-gap sequence as a single region', () => {
    expect(calculateGapRegions('---', 0)).toEqual([{ start: 0, end: 3 }]);
  });

  it('applies the offset to region coordinates', () => {
    expect(calculateGapRegions('-A-', 5)).toEqual([
      { start: 5, end: 6 },
      { start: 7, end: 8 },
    ]);
  });
});

describe('isPositionInGap', () => {
  const regions = [{ start: 2, end: 4 }]; // gap columns 2 and 3 (end exclusive)

  it('is true inside the gap run and false outside', () => {
    expect(isPositionInGap(1, regions)).toBe(false);
    expect(isPositionInGap(2, regions)).toBe(true);
    expect(isPositionInGap(3, regions)).toBe(true);
    expect(isPositionInGap(4, regions)).toBe(false);
  });
});

describe('calculateBiologicalPosition', () => {
  // 'AT--CG' -> gap columns 2,3 -> regions [2,4).
  const regions = [{ start: 2, end: 4 }];

  it('is the count of non-gap columns before the alignment position', () => {
    expect(calculateBiologicalPosition(0, regions)).toBe(0); // before A
    expect(calculateBiologicalPosition(1, regions)).toBe(1); // before T
    expect(calculateBiologicalPosition(4, regions)).toBe(2); // C is the 2nd base (0-based)
    expect(calculateBiologicalPosition(5, regions)).toBe(3); // G is the 3rd base
  });

  it('inside a gap, counts the non-gap columns strictly before it', () => {
    expect(calculateBiologicalPosition(3, regions)).toBe(2);
  });

  it('subtracts gaps from multiple regions', () => {
    // gaps at columns 2,3 and 6 -> 3 gap columns before column 8.
    const two = [
      { start: 2, end: 4 },
      { start: 6, end: 7 },
    ];
    expect(calculateBiologicalPosition(8, two)).toBe(5);
  });
});

describe('calculateAlignmentPosition', () => {
  // gaps at columns 2,3 -> the non-gap columns are 0,1,4,5,...
  const regions = [{ start: 2, end: 4 }];

  it('returns the alignment column of the Nth non-gap base (N is 1-based)', () => {
    expect(calculateAlignmentPosition(1, regions, 0)).toBe(0); // 1st non-gap base
    expect(calculateAlignmentPosition(2, regions, 0)).toBe(1); // 2nd non-gap base
    expect(calculateAlignmentPosition(3, regions, 0)).toBe(4); // 3rd non-gap base (cols 2,3 are gaps)
  });

  it('with no gaps the column is N-1', () => {
    expect(calculateAlignmentPosition(3, [], 0)).toBe(2);
  });

  it('passes positions before the reference start through unchanged', () => {
    expect(calculateAlignmentPosition(2, [], 5)).toBe(2);
  });
});

describe('calculateGapRegions + coordinate mapping (end-to-end on a real reference)', () => {
  it('maps non-gap columns to their 0-based base index', () => {
    const ref = 'AT--CG'; // bases A,T,C,G at columns 0,1,4,5
    const regions = calculateGapRegions(ref, 0);
    expect(regions).toEqual([{ start: 2, end: 4 }]);
    expect(isPositionInGap(2, regions)).toBe(true);
    expect(isPositionInGap(4, regions)).toBe(false);
    expect(calculateBiologicalPosition(4, regions)).toBe(2); // C is base index 2
    expect(calculateBiologicalPosition(5, regions)).toBe(3); // G is base index 3
  });
});

describe('calculateReadsCoverageRange', () => {
  function withReads(reads: Array<{ offset: number; sequence: string }>): {
    reads: Array<{ offset: number; sequence: string }>;
  } {
    return { reads };
  }

  it('should return undefined when there are no reads', () => {
    expect(calculateReadsCoverageRange(withReads([]))).toBeUndefined();
  });

  it('should return single read range', () => {
    expect(calculateReadsCoverageRange(withReads([{ offset: 20, sequence: 'T'.repeat(60) }]))).toEqual({
      start: 20,
      end: 80,
    });
  });

  it('should return union when reads cover different ranges', () => {
    expect(
      calculateReadsCoverageRange(
        withReads([
          { offset: 10, sequence: 'T'.repeat(50) },
          { offset: 20, sequence: 'G'.repeat(60) },
        ]),
      ),
    ).toEqual({ start: 10, end: 80 });
  });

  it('should handle a trimmed read on the left - other reads extend coverage', () => {
    expect(
      calculateReadsCoverageRange(
        withReads([
          { offset: 30, sequence: 'T'.repeat(50) },
          { offset: 10, sequence: 'G'.repeat(80) },
        ]),
      ),
    ).toEqual({ start: 10, end: 90 });
  });

  it('should handle multiple non-overlapping reads', () => {
    expect(
      calculateReadsCoverageRange(
        withReads([
          { offset: 0, sequence: 'T'.repeat(30) },
          { offset: 70, sequence: 'G'.repeat(30) },
        ]),
      ),
    ).toEqual({ start: 0, end: 100 });
  });

  it('should apply biologicalOffset correctly', () => {
    expect(calculateReadsCoverageRange(withReads([{ offset: 20, sequence: 'T'.repeat(60) }]), 100)).toEqual({
      start: 120,
      end: 180,
    });
  });

  it('should count non-gapped length only (internal gaps do not extend coverage)', () => {
    // sequence has 4 non-gap bases -> covers offset..offset+4.
    expect(calculateReadsCoverageRange(withReads([{ offset: 10, sequence: 'A-T-C-G' }]))).toEqual({
      start: 10,
      end: 14,
    });
  });
});
