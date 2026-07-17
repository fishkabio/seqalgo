import { GeneticCode, STANDARD_GENETIC_CODE } from './genetic-code';

/** Options for {@link translate}. */
export interface TranslateOptions {
  /** Genetic code to use; defaults to the Standard Code. */
  geneticCode?: GeneticCode;
  /** Stop at the first stop codon instead of translating through it. */
  toStop?: boolean;
  /** Amino-acid symbol for a codon absent from the code (e.g. one with an IUPAC ambiguity base). */
  unknown?: string;
}

/**
 * Translate a nucleotide sequence to a one-letter amino-acid string, reading frame 0.
 * A trailing partial codon is ignored; unknown codons become {@link TranslateOptions.unknown} (default `X`);
 * `U`/`u` are read as `T`. Slice the input yourself to translate a different frame.
 */
export function translate(sequence: string, options: TranslateOptions = {}): string {
  const code = options.geneticCode ?? STANDARD_GENETIC_CODE;
  const unknown = options.unknown ?? 'X';
  // Drop whitespace so formatted input (newlines, spaces) never silently shifts the reading frame.
  const dna = sequence.toUpperCase().replace(/\s/g, '').replace(/U/g, 'T');
  const codons = Math.floor(dna.length / 3);
  let protein = '';
  for (let i = 0; i < codons; i++) {
    const aminoAcid = code[dna.substr(i * 3, 3)] ?? unknown;
    if (options.toStop && aminoAcid === '*') break;
    protein += aminoAcid;
  }
  return protein;
}
