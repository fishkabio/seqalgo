import { AlignedSequence } from './types';

/** The minimal `{ offset, sequence }` shape the gap-fill transforms operate on. */
type OffsetSequence = Pick<AlignedSequence, 'sequence' | 'offset'>;

/**
 * Converts offset-based alignment data to fully gap-filled rows (all the same
 * length): `[reference, ...reads]`, each padded with leading gaps for its offset
 * and trailing gaps up to the alignment length.
 * @param alignment Reference plus reads, by `{ sequence, offset }`
 * @returns Array of gap-filled sequences `[reference, ...reads]`
 */
export function toGapFilledAlignment(alignment: { reference: OffsetSequence; reads: OffsetSequence[] }): string[] {
  const sequences = [
    { sequence: alignment.reference.sequence, offset: alignment.reference.offset },
    ...alignment.reads.map(r => ({ sequence: r.sequence, offset: r.offset })),
  ];

  if (sequences.length === 0) return [];

  const totalLength = Math.max(...sequences.map(s => s.offset + s.sequence.length));

  return sequences.map(
    seq =>
      '-'.repeat(seq.offset) + seq.sequence + '-'.repeat(Math.max(0, totalLength - seq.offset - seq.sequence.length)),
  );
}

/** Extracts offset and sequence from a gap-filled string, preserving internal gaps. */
function extractOffsetAndSequence(gapFilledSequence: string): OffsetSequence {
  const offset = Math.max(0, gapFilledSequence.search(/[^-]/));
  const sequence = gapFilledSequence.substring(offset).replace(/-+$/, '');
  return { offset, sequence };
}

/**
 * Converts gap-filled sequences back to offset-based format (the inverse of
 * {@link toGapFilledAlignment}). Reads are returned as the minimal `{ sequence,
 * offset }` shape — callers re-attach any per-read metadata (id, chromatogram, …).
 * @param gapFilledSequences Array where the first element is the reference, the rest are reads
 * @returns Reference (named "Reference") and reads in offset-based form
 */
export function fromGapFilledAlignment(gapFilledSequences: string[]): {
  reference: AlignedSequence;
  reads: OffsetSequence[];
} {
  if (gapFilledSequences.length === 0) {
    return { reference: { name: 'Reference', sequence: '', offset: 0 }, reads: [] };
  }

  const alignedReference = gapFilledSequences[0];
  const alignedReads = gapFilledSequences.slice(1);

  const { offset: referenceOffset, sequence: referenceSequence } = extractOffsetAndSequence(alignedReference);
  const reads = alignedReads.map(alignedRead => extractOffsetAndSequence(alignedRead));

  return {
    reference: { name: 'Reference', sequence: referenceSequence, offset: referenceOffset },
    reads,
  };
}

/**
 * Generates a unique name by adding a numeric suffix if the name already exists.
 * If `baseName` exists, tries `baseName_1`, `baseName_2`, … until an unused one is found.
 * @param baseName The original name to make unique
 * @param existingNames Set of existing names to avoid conflicts
 * @returns A unique name
 */
export function generateUniqueName(baseName: string, existingNames: Set<string>): string {
  if (!existingNames.has(baseName)) {
    return baseName;
  }

  let counter = 1;
  let candidateName = `${baseName}_${counter}`;
  while (existingNames.has(candidateName)) {
    counter++;
    candidateName = `${baseName}_${counter}`;
  }
  return candidateName;
}
