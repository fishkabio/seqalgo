import { COMMON_RESTRICTION_ENZYMES } from '../../src/restriction/enzymes';
import { findRestrictionSites } from '../../src/restriction/find-sites';
import { RestrictionEnzyme } from '../../src/restriction/types';

const EcoRI: RestrictionEnzyme = { name: 'EcoRI', site: 'GAATTC', cutTop: 1, cutBottom: 5 };
const EcoRV: RestrictionEnzyme = { name: 'EcoRV', site: 'GATATC', cutTop: 3, cutBottom: 3 };
const KpnI: RestrictionEnzyme = { name: 'KpnI', site: 'GGTACC', cutTop: 5, cutBottom: 1 };
// A synthetic non-palindromic enzyme so forward/reverse strand math can be checked in isolation.
const NP: RestrictionEnzyme = { name: 'NP', site: 'GGGAAA', cutTop: 1, cutBottom: 4 };

describe('findRestrictionSites — forward strand, cut coordinates', () => {
  it('locates a palindromic site and its top/bottom cuts (5′ overhang)', () => {
    expect(findRestrictionSites('AAGAATTCAA', [EcoRI])).toEqual([
      {
        enzyme: 'EcoRI',
        site: 'GAATTC',
        start: 2,
        length: 6,
        strand: 1,
        position: 3, // 2 + cutTop(1)
        bottomCut: 7, // 2 + cutBottom(5)
        overhang: 4,
        crossesOrigin: false,
      },
    ]);
  });

  it('reports back-to-back sites separately, sorted by cut position', () => {
    const sites = findRestrictionSites('GAATTCGAATTC', [EcoRI]);
    expect(sites.map(s => s.position)).toEqual([1, 7]);
    expect(sites.map(s => s.start)).toEqual([0, 6]);
  });

  it('reports a blunt cutter with a zero overhang and coincident cuts', () => {
    const [site] = findRestrictionSites('GATATC', [EcoRV]);
    expect(site).toMatchObject({ position: 3, bottomCut: 3, overhang: 0 });
  });

  it('reports a 3′ overhang cutter (cutTop > cutBottom) with a negative overhang', () => {
    const [site] = findRestrictionSites('GGTACC', [KpnI]);
    expect(site).toMatchObject({ position: 5, bottomCut: 1, overhang: -4 });
  });

  it('reports a palindromic site only once (does not double-count the reverse strand)', () => {
    expect(findRestrictionSites('GAATTC', [EcoRI])).toHaveLength(1);
  });
});

describe('findRestrictionSites — reverse strand (non-palindromic)', () => {
  it('finds the recognition site on the forward strand', () => {
    const sites = findRestrictionSites('CCGGGAAACC', [NP]);
    expect(sites).toEqual([
      {
        enzyme: 'NP',
        site: 'GGGAAA',
        start: 2,
        length: 6,
        strand: 1,
        position: 3, // 2 + cutTop(1)
        bottomCut: 6, // 2 + cutBottom(4)
        overhang: 3,
        crossesOrigin: false,
      },
    ]);
  });

  it('finds the reverse occurrence and mirrors the cuts onto the correct strands', () => {
    // 'TTTCCC' = revcomp('GGGAAA') at index 1; enzyme binds in reverse orientation.
    const sites = findRestrictionSites('ATTTCCCA', [NP]);
    expect(sites).toEqual([
      {
        enzyme: 'NP',
        site: 'GGGAAA',
        start: 1,
        length: 6,
        strand: -1,
        position: 3, // start + L - cutBottom = 1 + 6 - 4
        bottomCut: 6, // start + L - cutTop = 1 + 6 - 1
        overhang: 3, // orientation-independent
        crossesOrigin: false,
      },
    ]);
  });
});

describe('findRestrictionSites — circular molecules', () => {
  it('finds a site spanning the origin and wraps its cut coordinates', () => {
    // last base 'G' + first five 'AATTC' = GAATTC, starting at index 9 of a length-10 circle.
    const seq = 'AATTCGGGGG';
    const [site] = findRestrictionSites(seq, [EcoRI], { circular: true });
    expect(site).toMatchObject({
      start: 9,
      strand: 1,
      position: 0, // (9 + 1) mod 10
      bottomCut: 4, // (9 + 5) mod 10
      crossesOrigin: true,
    });
  });

  it('does not find an origin-spanning site when the molecule is linear', () => {
    expect(findRestrictionSites('AATTCGGGGG', [EcoRI])).toEqual([]);
  });
});

describe('findRestrictionSites — matching semantics and validation', () => {
  it('does not report a site where the sequence has an ambiguous base', () => {
    expect(findRestrictionSites('GANTTC', [EcoRI])).toEqual([]);
  });

  it('matches IUPAC degeneracy in the recognition sequence', () => {
    const HinfI: RestrictionEnzyme = { name: 'HinfI', site: 'GANTC', cutTop: 1, cutBottom: 4 };
    expect(findRestrictionSites('GACTC', [HinfI]).map(s => s.position)).toEqual([1]);
    expect(findRestrictionSites('GAGTC', [HinfI]).map(s => s.position)).toEqual([1]);
  });

  it('throws when a recognition sequence carries a non-IUPAC character', () => {
    expect(() => findRestrictionSites('GAATTC', [{ name: 'Bad', site: 'GAZTTC', cutTop: 1, cutBottom: 5 }])).toThrow();
  });

  it('sorts a multi-enzyme digest by cut position', () => {
    const sites = findRestrictionSites('GAATTCGGGGATATC', COMMON_RESTRICTION_ENZYMES);
    const names = sites.map(s => s.enzyme);
    expect(names).toContain('EcoRI');
    expect(names).toContain('EcoRV');
    const positions = sites.map(s => s.position);
    expect([...positions]).toEqual([...positions].sort((a, b) => a - b));
  });
});
