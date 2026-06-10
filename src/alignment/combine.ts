import { assertTruthy } from '@fishka/assertions';
import { removeGaps } from '../sequence';
import { removeLeadingAndTrailingGapColumns } from './gap-columns';
import { toGapFilledAlignment } from './transforms';
import { Alignment, CombinedAlignment, PairwiseAlignment, Region } from './types';

/**
 * Combines multiple pairwise anchor alignments into a single, square-like alignment block.
 *
 * 1. Builds a 'consensus reference' = the common ungapped reference with a gap inserted
 *    at every position where a gap exists in *any* input reference.
 * 2. Projects each read onto that consensus reference, preserving its original alignment
 *    within the new gapping pattern. Reads are NOT aligned against each other.
 *
 * @param pairwiseAlignments Pairwise alignments that share the same underlying (ungapped)
 *   reference sequence and name. Must not be empty.
 * @returns A `CombinedAlignment` with the consensus reference and the re-projected reads.
 */
export function combineAlignments(pairwiseAlignments: PairwiseAlignment[]): CombinedAlignment {
  assertTruthy(pairwiseAlignments.length > 0, 'List of alignments is empty.');

  const ref0Name = pairwiseAlignments[0].reference.name;
  const ungappedRef0 = removeGaps(pairwiseAlignments[0].reference.sequence);

  // Validate that all alignments share the same underlying reference sequence and name.
  for (let i = 1; i < pairwiseAlignments.length; i++) {
    const refI = pairwiseAlignments[i].reference;
    const ungappedRefI = removeGaps(refI.sequence);
    assertTruthy(refI.name === ref0Name, `Reference sequences have different names: ${ref0Name} vs ${refI.name}`);
    assertTruthy(ungappedRefI === ungappedRef0, 'Reference sequences are different when gaps are removed');
  }

  // If there's only one alignment, no combination is needed.
  if (pairwiseAlignments.length === 1) {
    return { reference: pairwiseAlignments[0].reference, reads: [pairwiseAlignments[0].read] };
  }

  // Build a consensus reference sequence that includes all gaps from all input references.
  const consensusReference = buildConsensusReference(pairwiseAlignments.map(a => a.reference.sequence));

  // Project each read onto the new consensus reference.
  const alignedReadSequences = pairwiseAlignments.map(({ read, reference }) =>
    alignReadToConsensus(read.sequence, reference.sequence, consensusReference),
  );

  // Pad all sequences to the same length to create a square alignment.
  const maxLength = Math.max(consensusReference.length, ...alignedReadSequences.map(s => s.length));
  const finalConsensus = consensusReference.padEnd(maxLength, '-');
  const finalReads = alignedReadSequences.map(seq => seq.padEnd(maxLength, '-'));

  // Remove any columns that contain only gaps from the start and end of the alignment.
  const allSequences = [finalConsensus, ...finalReads];
  const trimmedSequences = removeLeadingAndTrailingGapColumns(allSequences);

  // Find reference boundaries in the final trimmed alignment and shift flank gaps.
  const trimmedReference = trimmedSequences[0];
  const referenceBounds = findNonGapBounds(trimmedReference);
  const shiftedReads = trimmedSequences.slice(1).map(seq => shiftFlanksInRead(seq, referenceBounds));

  return {
    reference: {
      name: ref0Name,
      sequence: trimmedReference,
    },
    reads: pairwiseAlignments.map((pairwise, i) => ({
      ...pairwise.read,
      sequence: shiftedReads[i],
    })),
  };
}

/** Finds the bounds of non-gap characters in a sequence (inclusive `end`). */
function findNonGapBounds(sequence: string): Region {
  let start = 0;
  while (start < sequence.length && sequence[start] === '-') {
    start++;
  }
  let end = sequence.length - 1;
  while (end >= start && sequence[end] === '-') {
    end--;
  }
  return { start, end };
}

/** Counts gaps in a sequence region (inclusive indices). */
function countGapsInRegion(sequence: string, startIndex: number, endIndex: number): number {
  if (startIndex < 0 || startIndex >= sequence.length) return 0;
  const endPos = Math.min(endIndex, sequence.length - 1);
  const region = sequence.slice(startIndex, endPos + 1);
  return (region.match(/-/g) || []).length;
}

/**
 * Post-processes a read sequence to fix "torn flanks": a gap inserted between an
 * unaligned flank and the main alignment body. Using reference bounds, gaps are moved
 * from the flanks to the edges, e.g. 'GGG-A...' becomes '-GGGA...' and '...A-GGG'
 * becomes '...AGGG-'.
 */
export function shiftFlanksInRead(sequence: string, referenceBounds: Region): string {
  if (!sequence || sequence.length === 0) {
    return sequence;
  }

  // Count gaps in flanks.
  const leftFlankGaps = countGapsInRegion(sequence, 0, referenceBounds.start - 1);
  const rightFlankGaps = countGapsInRegion(sequence, referenceBounds.end + 1, sequence.length - 1);

  // Extract components using substrings.
  const leftFlankRegion = sequence.slice(0, referenceBounds.start);
  const alignmentRegion = sequence.slice(referenceBounds.start, referenceBounds.end + 1);
  const rightFlankRegion = sequence.slice(referenceBounds.end + 1);

  // Reconstruct: leading gaps + left flank chars + alignment + right flank chars + trailing gaps.
  return (
    '-'.repeat(leftFlankGaps) +
    removeGaps(leftFlankRegion) +
    alignmentRegion +
    removeGaps(rightFlankRegion) +
    '-'.repeat(rightFlankGaps)
  );
}

/**
 * Builds a consensus reference by merging all gap patterns from a set of reference
 * sequences: the result has a gap wherever any input had a gap.
 * @param referenceSequences Gapped reference sequence strings (same underlying sequence).
 * @returns The consensus reference sequence string.
 */
function buildConsensusReference(referenceSequences: string[]): string {
  const ungappedRef = removeGaps(referenceSequences[0]);
  // `insertions` stores the longest gap chunk found at each position. Size is
  // `ungappedRef.length + 1` to cover insertions before the first and after the last base.
  const insertions = Array(ungappedRef.length + 1)
    .fill(0)
    .map(() => '');

  for (const ref of referenceSequences) {
    let refCharCount = 0;
    let i = 0;
    while (i < ref.length) {
      if (ref[i] !== '-') {
        refCharCount++;
        i++;
      } else {
        // Found a gap in the reference, representing an insertion point.
        let insertionChunk = '';
        while (i < ref.length && ref[i] === '-') {
          insertionChunk += '-';
          i++;
        }
        // If this insertion is the longest we've seen at this position, record it.
        if (insertionChunk.length > insertions[refCharCount].length) {
          insertions[refCharCount] = insertionChunk;
        }
      }
    }
  }

  // Reconstruct the consensus string from the ungapped reference and the collected insertions.
  let consensus = insertions[0];
  for (let i = 0; i < ungappedRef.length; i++) {
    consensus += ungappedRef[i] + insertions[i + 1];
  }
  return consensus;
}

/**
 * Projects a read's alignment from its original reference onto a new consensus reference,
 * walking the consensus and, for each position, placing the right character from the read.
 * @param readSequence The original gapped read sequence.
 * @param referenceSequence The original gapped reference the read was aligned to.
 * @param consensusSequence The new consensus reference to align to.
 * @returns The new gapped read sequence, aligned to the consensus.
 */
function alignReadToConsensus(readSequence: string, referenceSequence: string, consensusSequence: string): string {
  let newRead = '';
  let originalReferencePointer = 0;
  let originalReadPointer = 0;

  for (let consensusPointer = 0; consensusPointer < consensusSequence.length; consensusPointer++) {
    // If we have exhausted the original alignment, pad with gaps.
    if (originalReferencePointer >= referenceSequence.length) {
      newRead += '-';
      continue;
    }

    const consensusChar = consensusSequence[consensusPointer];
    const refChar = referenceSequence[originalReferencePointer];

    if (consensusChar === refChar) {
      // Same character: copy the corresponding character from the original read.
      newRead += readSequence[originalReadPointer];
      originalReferencePointer++;
      originalReadPointer++;
    } else if (refChar === '-') {
      // Original reference has a gap (read insertion): preserve it, re-evaluate this consensus column.
      newRead += readSequence[originalReadPointer];
      originalReferencePointer++;
      originalReadPointer++;
      consensusPointer--;
    } else {
      // Consensus has a gap but the original reference has a base: introduce a gap in the read.
      newRead += '-';
    }
  }
  return newRead;
}

/**
 * Decomposes an alignment into pairwise alignments, pairing each read with the reference.
 * Each read's stable `id` is carried onto the resulting `PairwiseAlignment`.
 * Note: the resulting pairs may contain pairs of gaps.
 * @param alignment The alignment to decompose
 * @returns Array of `PairwiseAlignment` objects
 */
export function decomposeToPairwiseAlignments(alignment: Alignment): PairwiseAlignment[] {
  // Convert alignment to gap-filled format to work with full sequence positions.
  const gapFilledSequences = toGapFilledAlignment(alignment);
  if (gapFilledSequences.length === 0) {
    return [];
  }
  const gapFilledReference = gapFilledSequences[0];
  const reads = gapFilledSequences.slice(1);

  return alignment.reads.map((read, index) => ({
    reference: {
      name: alignment.reference.name,
      sequence: gapFilledReference,
    },
    read: {
      name: read.name,
      sequence: reads[index],
    },
    isReverseComplement: read.isReverseComplement,
    id: read.id,
  }));
}
