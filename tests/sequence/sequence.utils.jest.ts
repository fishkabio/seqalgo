import {
  complement,
  countLeadingGaps,
  countNonGappedLength,
  removeGaps,
  removeLeadingAndTrailingGaps,
  removeLeadingGaps,
  removeTrailingGaps,
  reverseComplement,
} from '../../src/sequence/sequence.utils';

describe('removeGaps', () => {
  it('should remove all gaps from sequence', () => {
    expect(removeGaps('AT-C-G')).toBe('ATCG');
    expect(removeGaps('---ATCG---')).toBe('ATCG');
    expect(removeGaps('ATCG')).toBe('ATCG');
    expect(removeGaps('----')).toBe('');
    expect(removeGaps('')).toBe('');
  });
});

describe('removeTrailingGaps', () => {
  it('should remove trailing gaps from sequence', () => {
    expect(removeTrailingGaps('ATCG---')).toBe('ATCG');
    expect(removeTrailingGaps('ATCG')).toBe('ATCG');
    expect(removeTrailingGaps('---')).toBe('');
    expect(removeTrailingGaps('AT-CG--')).toBe('AT-CG');
  });

  it('should handle empty string', () => {
    expect(removeTrailingGaps('')).toBe('');
  });

  it('should preserve internal gaps', () => {
    expect(removeTrailingGaps('AT--CG---')).toBe('AT--CG');
  });
});

describe('removeLeadingGaps', () => {
  it('should remove leading gaps from sequence', () => {
    expect(removeLeadingGaps('---ATCG')).toBe('ATCG');
    expect(removeLeadingGaps('ATCG')).toBe('ATCG');
    expect(removeLeadingGaps('----')).toBe('');
    expect(removeLeadingGaps('--AT-CG')).toBe('AT-CG');
  });

  it('should handle empty string', () => {
    expect(removeLeadingGaps('')).toBe('');
  });

  it('should preserve internal and trailing gaps', () => {
    expect(removeLeadingGaps('--AT--CG---')).toBe('AT--CG---');
  });
});

describe('removeLeadingAndTrailingGaps', () => {
  it('should remove both leading and trailing gaps', () => {
    expect(removeLeadingAndTrailingGaps('---ATCG---')).toBe('ATCG');
    expect(removeLeadingAndTrailingGaps('--AT-CG--')).toBe('AT-CG');
    expect(removeLeadingAndTrailingGaps('ATCG')).toBe('ATCG');
    expect(removeLeadingAndTrailingGaps('----')).toBe('');
    expect(removeLeadingAndTrailingGaps('')).toBe('');
  });

  it('should preserve internal gaps', () => {
    expect(removeLeadingAndTrailingGaps('--AT--CG--')).toBe('AT--CG');
  });
});

describe('countLeadingGaps', () => {
  it('should calculate correct offset for leading gaps', () => {
    expect(countLeadingGaps('---ATCG')).toBe(3);
    expect(countLeadingGaps('ATCG')).toBe(0);
    expect(countLeadingGaps('----')).toBe(4);
    expect(countLeadingGaps('-A-T-CG')).toBe(1);
  });

  it('should handle empty string', () => {
    expect(countLeadingGaps('')).toBe(0);
  });
});

describe('countNonGappedLength', () => {
  it('should count non-gap characters correctly', () => {
    expect(countNonGappedLength('ATCG')).toBe(4);
    expect(countNonGappedLength('AT-CG')).toBe(4);
    expect(countNonGappedLength('---ATCG---')).toBe(4);
    expect(countNonGappedLength('----')).toBe(0);
    expect(countNonGappedLength('')).toBe(0);
  });
});

describe('reverseComplement', () => {
  it('reverse-complements A/T/G/C', () => {
    expect(reverseComplement('ATCG')).toBe('CGAT');
    expect(reverseComplement('AAAA')).toBe('TTTT');
    expect(reverseComplement('GGCC')).toBe('GGCC');
  });

  it('complements all IUPAC ambiguity codes (R↔Y, K↔M, B↔V, D↔H; S/W/N self)', () => {
    expect(reverseComplement('R')).toBe('Y');
    expect(reverseComplement('Y')).toBe('R');
    // RYKM -> reverse MKYR -> complement K,M,R,Y
    expect(reverseComplement('RYKM')).toBe('KMRY');
    // BVDH -> reverse HDVB -> complement D,H,B,V
    expect(reverseComplement('BVDH')).toBe('DHBV');
    // S, W, N are self-complementary
    expect(reverseComplement('SWN')).toBe('NWS');
  });

  it('preserves case (soft-masking survives)', () => {
    expect(reverseComplement('atcg')).toBe('cgat');
    // aTcG -> reverse GcTa -> complement C,g,A,t
    expect(reverseComplement('aTcG')).toBe('CgAt');
  });

  it('preserves N and gaps', () => {
    expect(reverseComplement('N')).toBe('N');
    expect(reverseComplement('ATN')).toBe('NAT');
    expect(reverseComplement('AT-CG')).toBe('CG-AT');
  });

  it('passes characters not in the table through unchanged', () => {
    // '*' is not a nucleotide code, so it is preserved as-is.
    expect(reverseComplement('A*')).toBe('*T');
  });

  it('handles empty string', () => {
    expect(reverseComplement('')).toBe('');
  });
});

describe('complement', () => {
  it('complements A/T/G/C without reversing order', () => {
    expect(complement('ATCG')).toBe('TAGC');
    expect(complement('AAAA')).toBe('TTTT');
    expect(complement('GGCC')).toBe('CCGG');
  });

  it('complements all IUPAC ambiguity codes (R↔Y, K↔M, B↔V, D↔H; S/W/N self)', () => {
    expect(complement('R')).toBe('Y');
    expect(complement('Y')).toBe('R');
    expect(complement('RYKM')).toBe('YRMK');
    expect(complement('BVDH')).toBe('VBHD');
    expect(complement('SWN')).toBe('SWN');
  });

  it('preserves case (soft-masking survives)', () => {
    expect(complement('atcg')).toBe('tagc');
    expect(complement('aTcG')).toBe('tAgC');
  });

  it('preserves N and gaps', () => {
    expect(complement('N')).toBe('N');
    expect(complement('ATN')).toBe('TAN');
    expect(complement('AT-CG')).toBe('TA-GC');
  });

  it('passes characters not in the table through unchanged', () => {
    expect(complement('A*')).toBe('T*');
  });

  it('handles empty string', () => {
    expect(complement('')).toBe('');
  });
});
