import { reverseComplement } from '../sequence/sequence.utils';
import { GeneticCode, STANDARD_GENETIC_CODE } from '../translation/genetic-code';

/** An open reading frame in 0-based half-open genomic coordinates ([start, end), strand ±1). */
export interface Orf {
  /** 0-based start on the plus strand (inclusive). For a minus-strand ORF this is still the lower coordinate. */
  readonly start: number;
  /** 0-based end on the plus strand (exclusive). Includes the stop codon when `hasStop`. */
  readonly end: number;
  /** Reading direction: +1 plus strand, -1 minus strand. */
  readonly strand: 1 | -1;
  /** Reading frame 0..2 within the ORF's own strand (offset of the start codon along that strand). */
  readonly frame: number;
  /** Number of amino acids in the ORF, excluding the stop codon. */
  readonly aaLength: number;
  /** The start codon (5'->3' on the ORF's strand). */
  readonly startCodon: string;
  /** The stop codon (5'->3'), or '' when the ORF runs to the sequence end without a stop. */
  readonly stopCodon: string;
  /** Whether the ORF is terminated by a stop codon (false only when it runs off the end). */
  readonly hasStop: boolean;
}

export interface FindOrfsOptions {
  /** Minimum ORF length in amino acids (excluding the stop). Default 30. */
  readonly minAminoAcids?: number;
  /** Genetic code used to identify start/stop codons. Default {@link STANDARD_GENETIC_CODE}. */
  readonly geneticCode?: GeneticCode;
  /** Also scan the reverse strand. Default true. */
  readonly includeReverse?: boolean;
  /** Codons that open an ORF (5'->3'). Default `['ATG']`. */
  readonly startCodons?: readonly string[];
  /**
   * Require the ORF to end in a stop codon. Default true. When false, a start codon with no
   * downstream in-frame stop yields an ORF running to the end of the (in-frame) sequence.
   */
  readonly requireStop?: boolean;
}

/** One ORF found in a single strand's own 5'->3' coordinate system, before mapping to genomic. */
interface StrandOrf {
  readonly seqStart: number;
  readonly seqEnd: number;
  readonly frame: number;
  readonly aaLength: number;
  readonly startCodon: string;
  readonly stopCodon: string;
  readonly hasStop: boolean;
}

/** Scan one strand (already oriented 5'->3') for ORFs in all three frames. */
function scanStrand(
  seq: string,
  code: GeneticCode,
  startCodons: ReadonlySet<string>,
  minAa: number,
  requireStop: boolean,
): StrandOrf[] {
  const orfs: StrandOrf[] = [];
  for (let frame = 0; frame < 3; frame++) {
    let startIdx = -1;
    let startCodon = '';
    for (let pos = frame; pos + 3 <= seq.length; pos += 3) {
      const codon = seq.slice(pos, pos + 3); // `seq` is already upper-case DNA (see findOrfs)
      if (startIdx === -1 && startCodons.has(codon)) {
        startIdx = pos;
        startCodon = codon;
      }
      if (code[codon] === '*') {
        if (startIdx !== -1) {
          const aaLength = (pos - startIdx) / 3; // codons before the stop = amino acids
          if (aaLength >= minAa) {
            orfs.push({ seqStart: startIdx, seqEnd: pos + 3, frame, aaLength, startCodon, stopCodon: codon, hasStop: true });
          }
        }
        startIdx = -1;
        startCodon = '';
      }
    }
    // An open ORF with no in-frame stop before the end of this frame's codons.
    if (!requireStop && startIdx !== -1) {
      const lastCodon = frame + Math.floor((seq.length - frame) / 3) * 3; // one past the last full codon
      const aaLength = (lastCodon - startIdx) / 3;
      if (aaLength >= minAa) {
        orfs.push({ seqStart: startIdx, seqEnd: lastCodon, frame, aaLength, startCodon, stopCodon: '', hasStop: false });
      }
    }
  }
  return orfs;
}

/**
 * Find open reading frames on both strands of a linear DNA sequence.
 *
 * An ORF is the span from the first start codon in a stop-delimited region to (and including) that
 * region's stop codon, per reading frame. Coordinates are 0-based half-open on the plus strand; a
 * minus-strand ORF's `[start, end)` still denotes plus-strand coordinates with `strand: -1`, and its
 * codons read 5'->3' along the minus strand (i.e. right-to-left in plus-strand coordinates).
 *
 * Circular origin-spanning ORFs are NOT modelled here — pass a linearised sequence.
 *
 * Input is normalised before scanning: whitespace is removed and RNA is read as DNA (U→T), so the
 * returned coordinates are relative to the WHITESPACE-STRIPPED sequence (for a `GenbankDocument`'s
 * `sequence`, which carries no layout whitespace, they equal the record's own coordinates). The
 * genetic code is expected to be keyed by upper-case DNA codons.
 */
export function findOrfs(sequence: string, options: FindOrfsOptions = {}): Orf[] {
  const clean = sequence.replace(/\s/g, '').toUpperCase().replace(/U/g, 'T');
  const length = clean.length;
  const code = options.geneticCode ?? STANDARD_GENETIC_CODE;
  const minAa = Math.max(1, options.minAminoAcids ?? 30);
  const includeReverse = options.includeReverse ?? true;
  const requireStop = options.requireStop ?? true;
  const startCodons = new Set((options.startCodons ?? ['ATG']).map(c => c.toUpperCase().replace(/U/g, 'T')));

  const orfs: Orf[] = [];

  for (const plus of scanStrand(clean, code, startCodons, minAa, requireStop)) {
    orfs.push({
      start: plus.seqStart,
      end: plus.seqEnd,
      strand: 1,
      frame: plus.frame,
      aaLength: plus.aaLength,
      startCodon: plus.startCodon,
      stopCodon: plus.stopCodon,
      hasStop: plus.hasStop,
    });
  }

  if (includeReverse) {
    const rc = reverseComplement(clean);
    for (const minus of scanStrand(rc, code, startCodons, minAa, requireStop)) {
      // rc coordinate p maps to plus-strand coordinate (length - p); an rc span [s, e) becomes the
      // plus-strand span [length - e, length - s).
      orfs.push({
        start: length - minus.seqEnd,
        end: length - minus.seqStart,
        strand: -1,
        frame: minus.frame,
        aaLength: minus.aaLength,
        startCodon: minus.startCodon,
        stopCodon: minus.stopCodon,
        hasStop: minus.hasStop,
      });
    }
  }

  // Deterministic order: by start, then strand (plus first), then longer first.
  orfs.sort((a, b) => a.start - b.start || b.strand - a.strand || b.aaLength - a.aaLength);
  return orfs;
}
