import { STANDARD_GENETIC_CODE, translate } from '../../src/translation';

describe('translate', () => {
  it('translates a simple ORF including the stop codon', () => {
    // ATG AAA TAA -> M K *
    expect(translate('ATGAAATAA')).toBe('MK*');
  });

  it('reads lowercase and RNA (U) input', () => {
    expect(translate('augaaauaa')).toBe('MK*');
  });

  it('ignores a trailing partial codon', () => {
    expect(translate('ATGAA')).toBe('M');
  });

  it('stops at the first stop codon when toStop is set', () => {
    expect(translate('ATGTAAATG', { toStop: true })).toBe('M');
    expect(translate('ATGTAAATG')).toBe('M*M');
  });

  it('marks a codon with an unknown/ambiguous base as X', () => {
    expect(translate('ATGNNNAAA')).toBe('MXK');
    expect(translate('ATGNNNAAA', { unknown: '-' })).toBe('M-K');
  });

  it('covers every codon in the standard code', () => {
    const bases = ['T', 'C', 'A', 'G'];
    for (const a of bases) {
      for (const b of bases) {
        for (const c of bases) {
          const codon = a + b + c;
          expect(translate(codon)).toBe(STANDARD_GENETIC_CODE[codon]);
        }
      }
    }
  });
});
