/** Primitive operations on gapped nucleotide sequences (gap char is '-'). */

/** Removes all gaps from the sequence. */
export function removeGaps(sequence: string): string {
  return sequence.replace(/-/g, '');
}

/** Removes leading gaps from sequence. */
export function removeLeadingGaps(sequence: string): string {
  return sequence.replace(/^-+/, '');
}

/** Removes trailing gaps from sequence. */
export function removeTrailingGaps(sequence: string): string {
  return sequence.replace(/-+$/, '');
}

/** Removes both leading and trailing gaps from sequence. */
export function removeLeadingAndTrailingGaps(sequence: string): string {
  return removeLeadingGaps(removeTrailingGaps(sequence));
}

/**
 * Counts the number of non-gap characters in a sequence.
 * Optimized version that doesn't create new strings.
 */
export function countNonGappedLength(sequence: string): number {
  let count = 0;
  for (let i = 0; i < sequence.length; i++) {
    if (sequence[i] !== '-') {
      count++;
    }
  }
  return count;
}

/** Counts the number of leading gaps in a sequence. */
export function countLeadingGaps(sequence: string): number {
  let count = 0;
  for (let i = 0; i < sequence.length; i++) {
    if (sequence[i] === '-') {
      count++;
    } else {
      break;
    }
  }
  return count;
}

/**
 * Complement table covering all IUPAC nucleotide codes in both cases, plus the
 * gap character. Built so an ambiguity code maps to the code of the complements
 * of its bases (e.g. R = A/G -> Y = T/C). Case is preserved; '-' maps to itself.
 */
const COMPLEMENT: ReadonlyMap<string, string> = ((): Map<string, string> => {
  const map = new Map<string, string>();
  // Each pair is mutually complementary; S/W/N are self-complementary.
  for (const [a, b] of ['AT', 'GC', 'RY', 'KM', 'BV', 'DH', 'SS', 'WW', 'NN']) {
    map.set(a, b);
    map.set(b, a);
    map.set(a.toLowerCase(), b.toLowerCase());
    map.set(b.toLowerCase(), a.toLowerCase());
  }
  map.set('-', '-');
  return map;
})();

/**
 * Creates the reverse complement of a DNA sequence.
 *
 * All IUPAC codes are complemented (A↔T, G↔C, R↔Y, K↔M, B↔V, D↔H; S/W/N map to
 * themselves) and case is preserved (so soft-masking survives). The gap '-' is
 * kept; any character not in the table passes through unchanged.
 */
export function reverseComplement(sequence: string): string {
  const out = new Array<string>(sequence.length);
  for (let i = sequence.length - 1, j = 0; i >= 0; i--, j++) {
    const ch = sequence[i];
    out[j] = COMPLEMENT.get(ch) ?? ch;
  }
  return out.join('');
}

/**
 * Complements a DNA sequence WITHOUT reversing it (order preserved). Uses the same table as
 * {@link reverseComplement}: all IUPAC codes complemented, case preserved, gap kept, unknown
 * characters pass through unchanged.
 */
export function complement(sequence: string): string {
  const out = new Array<string>(sequence.length);
  for (let i = 0; i < sequence.length; i++) {
    const ch = sequence[i];
    out[i] = COMPLEMENT.get(ch) ?? ch;
  }
  return out.join('');
}
