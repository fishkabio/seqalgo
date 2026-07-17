import {
  GENETIC_CODES,
  STANDARD_GENETIC_CODE,
  getGeneticCode,
  getGeneticCodeTable,
  translate,
} from '../../src/translation';

/**
 * Golden codon->amino-acid assignments taken from known biology / NCBI's published tables — NOT from
 * this module's own output. Each is a codon whose meaning is a well-documented difference of that
 * table from the Standard Code (or a Standard-Code anchor), so the set independently pins both the
 * data and the codon ordering used to build the maps.
 */
const GOLDEN: ReadonlyArray<readonly [number, string, string]> = [
  // 1 — Standard: anchors for order + the three stops.
  [1, 'ATG', 'M'],
  [1, 'TGG', 'W'],
  [1, 'TAA', '*'],
  [1, 'TAG', '*'],
  [1, 'TGA', '*'],
  [1, 'TTT', 'F'],
  [1, 'AAA', 'K'],
  [1, 'AGA', 'R'],
  [1, 'ATA', 'I'],
  [1, 'CTG', 'L'],
  // 2 — Vertebrate Mitochondrial: TGA=W, ATA=M, AGA/AGG=stop.
  [2, 'TGA', 'W'],
  [2, 'ATA', 'M'],
  [2, 'AGA', '*'],
  [2, 'AGG', '*'],
  // 3 — Yeast Mitochondrial: all four CTN read as Thr, TGA=W.
  [3, 'CTT', 'T'],
  [3, 'CTC', 'T'],
  [3, 'CTA', 'T'],
  [3, 'CTG', 'T'],
  [3, 'TGA', 'W'],
  // 5 — Invertebrate Mitochondrial: AGA/AGG=Ser, ATA=Met, TGA=W.
  [5, 'ATA', 'M'],
  [5, 'AGA', 'S'],
  [5, 'AGG', 'S'],
  [5, 'TGA', 'W'],
  // 6 — Ciliate Nuclear: TAA/TAG read as Gln, TGA still stop.
  [6, 'TAA', 'Q'],
  [6, 'TAG', 'Q'],
  [6, 'TGA', '*'],
  // 9 — Echinoderm/Flatworm Mito: AAA=Asn, AGA/AGG=Ser, TGA=W.
  [9, 'AAA', 'N'],
  [9, 'AGA', 'S'],
  [9, 'AGG', 'S'],
  [9, 'TGA', 'W'],
  // 10 — Euplotid Nuclear: TGA=Cys.
  [10, 'TGA', 'C'],
  // 13 — Ascidian Mitochondrial: AGA/AGG=Gly.
  [13, 'AGA', 'G'],
  [13, 'AGG', 'G'],
  // 16 — Chlorophycean Mito: TAG=Leu.
  [16, 'TAG', 'L'],
  // 22 — Scenedesmus obliquus Mito: TCA=stop, TAG=Leu.
  [22, 'TCA', '*'],
  [22, 'TAG', 'L'],
  // 24 — Rhabdopleuridae Mito: AGA=Ser, AGG=Lys.
  [24, 'AGA', 'S'],
  [24, 'AGG', 'K'],
  // 25 — Candidate Division SR1: TGA=Gly.
  [25, 'TGA', 'G'],
];

describe('genetic codes — NCBI translation tables', () => {
  it('exposes exactly the 27 NCBI tables, in id order, with no retired ids', () => {
    const ids = GENETIC_CODES.map(t => t.id);
    expect(ids).toEqual([
      1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
    ]);
    // Retired / never-assigned ids resolve to nothing.
    for (const gap of [0, 7, 8, 17, 18, 19, 20, 34, 99]) expect(getGeneticCode(gap)).toBeUndefined();
    expect(getGeneticCodeTable(7)).toBeUndefined();
  });

  it('gives every table all 64 codons, each an amino-acid letter or a stop', () => {
    for (const table of GENETIC_CODES) {
      const codons = Object.keys(table.code);
      expect(codons).toHaveLength(64);
      for (const aa of Object.values(table.code)) expect(aa).toMatch(/^[A-Z*]$/);
    }
  });

  it('matches golden codon assignments from known biology (data + codon order)', () => {
    for (const [id, codon, aa] of GOLDEN) {
      expect(getGeneticCode(id)?.[codon]).toBe(aa);
    }
  });

  it('carries the official NCBI names', () => {
    expect(getGeneticCodeTable(1)?.name).toBe('Standard');
    expect(getGeneticCodeTable(2)?.name).toBe('Vertebrate Mitochondrial');
    expect(getGeneticCodeTable(11)?.name).toBe('Bacterial, Archaeal and Plant Plastid');
  });

  it('keeps STANDARD_GENETIC_CODE identical to table 1', () => {
    expect(STANDARD_GENETIC_CODE).toEqual(getGeneticCode(1));
    // Table 11 (Bacterial) shares the standard amino-acid map (differs only in start codons).
    expect(getGeneticCode(11)).toEqual(STANDARD_GENETIC_CODE);
  });

  it('drives translate() with a chosen table — vertebrate mito reads TGA as Trp, AGA as stop', () => {
    const vertMito = getGeneticCode(2);
    // ATG TGA AGA -> M W * under table 2 (vs M * R under the standard code).
    expect(translate('ATGTGAAGA', { geneticCode: vertMito })).toBe('MW*');
    expect(translate('ATGTGAAGA')).toBe('M*R');
  });
});
