import { fromGapFilledAlignment, generateUniqueName, toGapFilledAlignment } from '../../src/alignment/transforms';

describe('toGapFilledAlignment', () => {
  it('should convert alignment data to gap-filled format', () => {
    const alignment = {
      reference: { name: 'Reference', sequence: 'ATCG', offset: 2 },
      reads: [
        { sequence: 'AT', offset: 0 },
        { sequence: 'CGA', offset: 3 },
      ],
    };

    expect(toGapFilledAlignment(alignment)).toEqual([
      '--ATCG', // reference
      'AT----', // read 1
      '---CGA', // read 2
    ]);
  });

  it('should handle alignment with no reads', () => {
    const alignment = { reference: { name: 'Reference', sequence: 'ATCG', offset: 0 }, reads: [] };
    expect(toGapFilledAlignment(alignment)).toEqual(['ATCG']);
  });
});

describe('fromGapFilledAlignment', () => {
  it('should convert gap-filled alignment back to offset-based format', () => {
    const result = fromGapFilledAlignment(['--ATCG-', 'AT-----', '---CGA-']);
    expect(result.reference.sequence).toBe('ATCG');
    expect(result.reference.offset).toBe(2);
    expect(result.reads).toEqual([
      { sequence: 'AT', offset: 0 },
      { sequence: 'CGA', offset: 3 },
    ]);
  });

  it('should handle empty array', () => {
    const result = fromGapFilledAlignment([]);
    expect(result.reference.sequence).toBe('');
    expect(result.reference.offset).toBe(0);
    expect(result.reads).toEqual([]);
  });

  it('should handle all-gap sequences', () => {
    const result = fromGapFilledAlignment(['---', '---', '---']);
    expect(result.reference.sequence).toBe('');
    expect(result.reference.offset).toBe(0);
    expect(result.reads).toEqual([
      { sequence: '', offset: 0 },
      { sequence: '', offset: 0 },
    ]);
  });

  it('should preserve internal gaps in reads', () => {
    const result = fromGapFilledAlignment(['AT-GCG', 'A--GCG']);
    expect(result.reference.sequence).toBe('AT-GCG');
    expect(result.reference.offset).toBe(0);
    expect(result.reads[0].sequence).toBe('A--GCG');
    expect(result.reads[0].offset).toBe(0);
  });

  it('should preserve internal gaps with leading gaps', () => {
    const result = fromGapFilledAlignment(['--AT-GCG--', '---A--CG--']);
    expect(result.reference.sequence).toBe('AT-GCG');
    expect(result.reference.offset).toBe(2);
    expect(result.reads[0].sequence).toBe('A--CG');
    expect(result.reads[0].offset).toBe(3);
  });

  it('round-trips with toGapFilledAlignment', () => {
    const gapFilled = ['--ATCG', 'AT----', '---CGA'];
    const offsetBased = fromGapFilledAlignment(gapFilled);
    expect(toGapFilledAlignment(offsetBased)).toEqual(gapFilled);
  });
});

describe('generateUniqueName', () => {
  it('should return original name if no conflict exists', () => {
    expect(generateUniqueName('read3', new Set(['read1', 'read2']))).toBe('read3');
  });

  it('should add _1 suffix when name conflicts', () => {
    expect(generateUniqueName('read1', new Set(['read1', 'read2']))).toBe('read1_1');
  });

  it('should increment suffix for multiple conflicts', () => {
    expect(generateUniqueName('read1', new Set(['read1', 'read1_1', 'read1_2']))).toBe('read1_3');
  });

  it('should find the first available suffix when there are gaps in the sequence', () => {
    expect(generateUniqueName('read1', new Set(['read1', 'read1_1', 'read1_3']))).toBe('read1_2');
  });

  it('should work with an empty set', () => {
    expect(generateUniqueName('read1', new Set<string>())).toBe('read1');
  });

  it('should handle names with underscores', () => {
    expect(generateUniqueName('my_read', new Set(['my_read', 'my_read_1']))).toBe('my_read_2');
  });

  it('should handle names with numeric suffixes already', () => {
    expect(generateUniqueName('read1_5', new Set(['read1_5', 'read1_5_1']))).toBe('read1_5_2');
  });

  it('should generate sequential unique names when called multiple times', () => {
    const existingNames = new Set(['read1']);
    const r1 = generateUniqueName('read1', existingNames);
    expect(r1).toBe('read1_1');
    existingNames.add(r1);
    const r2 = generateUniqueName('read1', existingNames);
    expect(r2).toBe('read1_2');
    existingNames.add(r2);
    expect(generateUniqueName('read1', existingNames)).toBe('read1_3');
  });
});
