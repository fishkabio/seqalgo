import { combineAlignments, decomposeToPairwiseAlignments, shiftFlanksInRead } from '../../src/alignment/combine';
import { removeColumnsOfGaps } from '../../src/alignment/gap-columns';
import { Alignment, CombinedAlignment, PairwiseAlignment } from '../../src/alignment/types';
import { removeGaps } from '../../src/sequence';

function createPairwise(pairs: { rf: string; rd: string }[]): PairwiseAlignment[] {
  return pairs.map(({ rf, rd }, i) => ({
    reference: { name: 'ref1', sequence: rf },
    read: { name: `read${i + 1}`, sequence: rd },
    isReverseComplement: false,
  }));
}

function checkResultAlignment(result: CombinedAlignment, reference: string, ...reads: string[]): void {
  expect(result.reference.sequence).toBe(reference);
  for (let i = 0; i < reads.length; i++) {
    expect(result.reads[i].sequence).toBe(reads[i]);
  }
}

function checkResultCorrespondsToInput(input: PairwiseAlignment[], result: CombinedAlignment): void {
  expect(removeGaps(result.reference.sequence)).toBe(removeGaps(input[0].reference.sequence));
  expect(result.reads.length).toBe(input.length);

  for (let i = 0; i < input.length; i++) {
    expect(result.reference.name).toBe(input[i].reference.name);
    expect(removeGaps(result.reads[i].sequence)).toBe(removeGaps(input[i].read.sequence));
  }

  // Removing all-gap columns from the result reproduces each original pairwise alignment.
  for (let i = 0; i < input.length; i++) {
    const ungappedOriginalPair = removeColumnsOfGaps([input[i].reference.sequence, input[i].read.sequence]);
    const ungappedResultPair = removeColumnsOfGaps([result.reference.sequence, result.reads[i].sequence]);
    expect(ungappedResultPair).toEqual(ungappedOriginalPair);
  }
}

describe('combineAlignments', () => {
  it('should handle single alignment without changes', () => {
    const input = createPairwise([{ rf: 'ATCG', rd: 'ATCG' }]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    expect(result.reference.sequence).toBe('ATCG');
    expect(result.reads[0].sequence).toBe('ATCG');
  });

  it('should throw error for empty alignment array', () => {
    expect(() => combineAlignments([])).toThrow('List of alignments is empty.');
  });

  it('should throw error for different reference names', () => {
    const input = createPairwise([
      { rf: 'ATCG', rd: 'ATCG' },
      { rf: 'ATCG', rd: 'AACG' },
    ]);
    input[1].reference.name = 'ref2';
    expect(() => combineAlignments(input)).toThrow('Reference sequences have different names: ref1 vs ref2');
  });

  it('should throw error for different reference sequences when gaps removed', () => {
    const input = createPairwise([
      { rf: 'ATCG', rd: 'ATCG' },
      { rf: 'TTCG', rd: 'AACG' },
    ]);
    expect(() => combineAlignments(input)).toThrow('Reference sequences are different when gaps are removed');
  });

  it('should handle different internal gaps in reads', () => {
    const input = createPairwise([
      { rf: 'ATACG', rd: 'AT-CG' },
      { rf: 'ATACG', rd: 'A-ACG' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'ATACG', 'AT-CG', 'A-ACG');
  });

  it('should handle same internal gaps in reads', () => {
    const input = createPairwise([
      { rf: 'ATACG', rd: 'AT-CG' },
      { rf: 'ATACG', rd: 'AT-CG' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'ATACG', 'AT-CG', 'AT-CG');
  });

  it('should handle different internal gaps in reference', () => {
    const input = createPairwise([
      { rf: 'AT-CG', rd: 'ATGCG' },
      { rf: 'ATC-G', rd: 'ATCAG' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'AT-C-G', 'ATGC-G', 'AT-CAG');
  });

  it('should handle same internal gaps in reference', () => {
    const input = createPairwise([
      { rf: 'AT-CG', rd: 'ATGCG' },
      { rf: 'AT-CG', rd: 'ATACG' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'AT-CG', 'ATGCG', 'ATACG');
  });

  it('should correctly combine alignments with different leading/trailing gap patterns', () => {
    const input = createPairwise([
      { rf: 'ATCG--', rd: 'ATCGAA' },
      { rf: '--ATCG', rd: 'AAATCG' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, '--ATCG--', '--ATCGAA', 'AAATCG--');
  });

  it('should handle a read that is much shorter than the reference', () => {
    const input = createPairwise([{ rf: 'ATGCGT', rd: 'A--C--' }]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'ATGCGT', 'A--C--');
  });

  it('should handle a read that is longer than the reference (insertion)', () => {
    const input = createPairwise([{ rf: 'AT--GT', rd: 'ATCCGT' }]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'AT--GT', 'ATCCGT');
  });

  it('should handle a complex mix of matching and short reads', () => {
    const input = createPairwise([
      { rf: 'A-TCG', rd: 'ACTCG' },
      { rf: 'AT-CG', rd: 'A--C-' },
      { rf: 'AT-CG', rd: 'A---G' },
    ]);
    const result = combineAlignments(input);
    checkResultCorrespondsToInput(input, result);
    checkResultAlignment(result, 'A-T-CG', 'ACT-CG', 'A---C-', 'A----G');
  });
});

describe('shiftFlanksInRead', () => {
  it('should handle empty or null sequences', () => {
    const bounds = { start: 0, end: 3 };
    expect(shiftFlanksInRead('', bounds)).toBe('');
    expect(shiftFlanksInRead(null as unknown as string, bounds)).toBe(null);
    expect(shiftFlanksInRead(undefined as unknown as string, bounds)).toBe(undefined);
  });

  it('should not modify sequences without flank gaps', () => {
    expect(shiftFlanksInRead('ATCG', { start: 0, end: 3 })).toBe('ATCG');
    expect(shiftFlanksInRead('A-T-', { start: 0, end: 2 })).toBe('A-T-');
  });

  it('should shift leading flank gaps to the beginning', () => {
    expect(shiftFlanksInRead('GGG-ATCG', { start: 4, end: 7 })).toBe('-GGGATCG');
    expect(shiftFlanksInRead('AT--CGGT', { start: 4, end: 7 })).toBe('--ATCGGT');
  });

  it('should shift trailing flank gaps to the end', () => {
    expect(shiftFlanksInRead('ATCG-GGG', { start: 0, end: 3 })).toBe('ATCGGGG-');
    expect(shiftFlanksInRead('ATCG--TT', { start: 0, end: 3 })).toBe('ATCGTT--');
  });

  it('should shift both leading and trailing flank gaps', () => {
    expect(shiftFlanksInRead('GGG-ATCG-TTT', { start: 4, end: 7 })).toBe('-GGGATCGTTT-');
    expect(shiftFlanksInRead('AT--CG--GT', { start: 4, end: 5 })).toBe('--ATCGGT--');
  });

  it('should handle gaps in reference region without moving them', () => {
    expect(shiftFlanksInRead('GG-A-T-C-G-TT', { start: 2, end: 7 })).toBe('GG-A-T-CGTT--');
  });

  it('should preserve sequences that are all within reference bounds', () => {
    expect(shiftFlanksInRead('ATCGGGGTTT', { start: 0, end: 9 })).toBe('ATCGGGGTTT');
    expect(shiftFlanksInRead('----------', { start: 0, end: 9 })).toBe('----------');
  });

  it('should handle single character flanks with gaps', () => {
    expect(shiftFlanksInRead('A-TCGG', { start: 2, end: 5 })).toBe('-ATCGG');
    expect(shiftFlanksInRead('TCGG-A', { start: 0, end: 3 })).toBe('TCGGA-');
  });

  it('should handle sequences where flanks are already correctly positioned', () => {
    expect(shiftFlanksInRead('-ATCG', { start: 1, end: 4 })).toBe('-ATCG');
    expect(shiftFlanksInRead('ATCG-', { start: 0, end: 3 })).toBe('ATCG-');
    expect(shiftFlanksInRead('-ATCG-', { start: 1, end: 4 })).toBe('-ATCG-');
  });

  it('should handle edge cases with reference bounds', () => {
    expect(shiftFlanksInRead('A-T', { start: 1, end: 1 })).toBe('A-T');
    expect(shiftFlanksInRead('AT-GG', { start: 0, end: 1 })).toBe('ATGG-');
    expect(shiftFlanksInRead('GG-AT', { start: 2, end: 3 })).toBe('GG-AT');
  });
});

describe('decomposeToPairwiseAlignments', () => {
  it('pairs each read with the gap-filled reference and carries the read id', () => {
    const alignment: Alignment = {
      name: 'a',
      reference: { name: 'Ref', sequence: 'ATCG', offset: 1 },
      reads: [
        { id: 'r1', name: 'read1', sequence: 'TC', offset: 2, isReverseComplement: false },
        { id: 'r2', name: 'read2', sequence: 'CG', offset: 3, isReverseComplement: true },
      ],
    };

    // gap-filled (total length 5): ref '-ATCG', read1 '--TC-', read2 '---CG'.
    expect(decomposeToPairwiseAlignments(alignment)).toEqual([
      {
        reference: { name: 'Ref', sequence: '-ATCG' },
        read: { name: 'read1', sequence: '--TC-' },
        isReverseComplement: false,
        id: 'r1',
      },
      {
        reference: { name: 'Ref', sequence: '-ATCG' },
        read: { name: 'read2', sequence: '---CG' },
        isReverseComplement: true,
        id: 'r2',
      },
    ]);
  });

  it('returns an empty array for an alignment with no reads', () => {
    const alignment: Alignment = {
      name: 'a',
      reference: { name: 'Ref', sequence: 'ATCG', offset: 0 },
      reads: [],
    };
    expect(decomposeToPairwiseAlignments(alignment)).toEqual([]);
  });
});
