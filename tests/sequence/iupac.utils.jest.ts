import { IUPAC_NUCLEOTIDES_LIST, NUCLEOTIDES_TO_IUPAC } from '../../src/sequence/iupac.utils';

describe('IUPAC_NUCLEOTIDES_LIST', () => {
  it('maps each code to the bases it represents', () => {
    expect(IUPAC_NUCLEOTIDES_LIST.A).toEqual(['A']);
    expect(IUPAC_NUCLEOTIDES_LIST.R).toEqual(['A', 'G']);
    expect(IUPAC_NUCLEOTIDES_LIST.Y).toEqual(['C', 'T']);
    expect(IUPAC_NUCLEOTIDES_LIST.B).toEqual(['C', 'G', 'T']);
    expect(IUPAC_NUCLEOTIDES_LIST.N).toEqual(['A', 'C', 'G', 'T']);
  });
});

describe('NUCLEOTIDES_TO_IUPAC', () => {
  it('maps a sorted base set back to its IUPAC code', () => {
    expect(NUCLEOTIDES_TO_IUPAC.A).toBe('A');
    expect(NUCLEOTIDES_TO_IUPAC.AG).toBe('R'); // A,G -> R
    expect(NUCLEOTIDES_TO_IUPAC.CT).toBe('Y'); // C,T -> Y
    expect(NUCLEOTIDES_TO_IUPAC.CGT).toBe('B'); // C,G,T -> B
    expect(NUCLEOTIDES_TO_IUPAC.ACGT).toBe('N'); // A,C,G,T -> N
  });

  it('is the inverse of IUPAC_NUCLEOTIDES_LIST under sorted-join keys', () => {
    for (const [code, bases] of Object.entries(IUPAC_NUCLEOTIDES_LIST)) {
      const key = [...bases].sort().join('');
      expect(NUCLEOTIDES_TO_IUPAC[key]).toBe(code);
    }
  });
});
