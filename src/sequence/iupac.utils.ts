/**
 * IUPAC nucleotide code mapping.
 * Maps each IUPAC code to the list of nucleotides it represents.
 */
export const IUPAC_NUCLEOTIDES_LIST: Record<string, string[]> = {
  A: ['A'],
  C: ['C'],
  G: ['G'],
  T: ['T'],
  R: ['A', 'G'], // puRine
  Y: ['C', 'T'], // pYrimidine
  S: ['G', 'C'], // Strong
  W: ['A', 'T'], // Weak
  K: ['G', 'T'], // Keto
  M: ['A', 'C'], // aMino
  B: ['C', 'G', 'T'], // not A
  D: ['A', 'G', 'T'], // not C
  H: ['A', 'C', 'T'], // not G
  V: ['A', 'C', 'G'], // not T
  N: ['A', 'C', 'G', 'T'], // aNy
};

/** Maps a sorted-and-joined set of nucleotides to its IUPAC code (e.g. 'AG' -> 'R'). */
export const NUCLEOTIDES_TO_IUPAC: Record<string, string> = ((): Record<string, string> => {
  const map: Record<string, string> = {};
  for (const [code, bases] of Object.entries(IUPAC_NUCLEOTIDES_LIST)) {
    const key = [...bases].sort().join('');
    map[key] = code;
  }
  return map;
})();
