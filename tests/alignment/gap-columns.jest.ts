import { removeColumnsOfGaps, removeLeadingAndTrailingGapColumns } from '../../src/alignment/gap-columns';

describe('removeLeadingAndTrailingGapColumns', () => {
  it('should handle empty array', () => {
    expect(removeLeadingAndTrailingGapColumns([])).toEqual([]);
  });

  it('should remove leading gap columns where all sequences have gaps', () => {
    expect(removeLeadingAndTrailingGapColumns(['--ATCG', '--AACG'])).toEqual(['ATCG', 'AACG']);
  });

  it('should remove trailing gap columns where all sequences have gaps', () => {
    expect(removeLeadingAndTrailingGapColumns(['ATCG--', 'AACG--'])).toEqual(['ATCG', 'AACG']);
  });

  it('should remove both leading and trailing gap columns', () => {
    expect(removeLeadingAndTrailingGapColumns(['--ATCG--', '--AACG--'])).toEqual(['ATCG', 'AACG']);
  });

  it('should preserve leading gaps if not all sequences have gaps', () => {
    expect(removeLeadingAndTrailingGapColumns(['ATCG', '-ACG'])).toEqual(['ATCG', '-ACG']);
  });

  it('should preserve trailing gaps if not all sequences have gaps', () => {
    expect(removeLeadingAndTrailingGapColumns(['ATCG', 'ACG-'])).toEqual(['ATCG', 'ACG-']);
  });

  it('should preserve internal gaps', () => {
    expect(removeLeadingAndTrailingGapColumns(['--AT-CG--', '--AA-CG--'])).toEqual(['AT-CG', 'AA-CG']);
  });

  it('should handle sequences of different lengths', () => {
    expect(removeLeadingAndTrailingGapColumns(['--ATCG', '--AA'])).toEqual(['ATCG', 'AA']);
  });

  it('should handle all gaps sequence', () => {
    expect(removeLeadingAndTrailingGapColumns(['----', '----'])).toEqual(['', '']);
  });

  it('should handle mixed gap patterns', () => {
    expect(removeLeadingAndTrailingGapColumns(['--ATCG---', '---ACG--', '--TTCG---'])).toEqual([
      'ATCG',
      '-ACG',
      'TTCG',
    ]);
  });
});

describe('removeColumnsOfGaps', () => {
  it('should handle empty array', () => {
    expect(removeColumnsOfGaps([])).toEqual([]);
  });

  it('should remove all gap columns', () => {
    expect(removeColumnsOfGaps(['A-TC-G', 'A-AC-G'])).toEqual(['ATCG', 'AACG']);
  });

  it('should preserve columns with at least one non-gap', () => {
    // Column 1 has '-' in first seq but 'A' in second, so it is preserved.
    expect(removeColumnsOfGaps(['A-TC', 'AAAG'])).toEqual(['A-TC', 'AAAG']);
  });

  it('should handle all gaps sequence', () => {
    expect(removeColumnsOfGaps(['----', '----'])).toEqual(['', '']);
  });

  it('should handle mixed patterns', () => {
    // Only columns where ALL seqs have gaps are removed.
    expect(removeColumnsOfGaps(['-A-T-C-', '--AT-C-', '-AAT-C-'])).toEqual(['A-TC', '-ATC', 'AATC']);
  });
});
