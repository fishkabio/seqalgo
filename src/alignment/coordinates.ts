import { countNonGappedLength } from '../sequence';
import { Region } from './types';

/**
 * Calculates gap regions from a sequence.
 * Finds all stretches of '-' characters and returns them as a Region array, where
 * each region is half-open `[start, end)` in alignment coordinates (plus `offset`).
 * @param sequence The sequence to scan for gaps
 * @param offset The offset to add to all region positions (alignment coordinates)
 * @returns Array of gap regions in alignment coordinates
 */
export function calculateGapRegions(sequence: string, offset: number): Region[] {
  const gapRegions: Region[] = [];
  for (let i = 0; i < sequence.length; i++) {
    if (sequence[i] === '-') {
      const start = i;
      while (i < sequence.length && sequence[i] === '-') {
        i++;
      }
      gapRegions.push({
        start: offset + start,
        end: offset + i,
      });
      i--; // Compensate for loop increment
    }
  }
  return gapRegions;
}

/**
 * Converts alignment position to biological position (subtracts gaps).
 * @param alignmentPos Position in alignment coordinates
 * @param gapRegions Gap regions
 * @returns Biological position (without gaps)
 */
export function calculateBiologicalPosition(alignmentPos: number, gapRegions: Region[]): number {
  let gapCount = 0;
  for (const region of gapRegions) {
    if (region.end <= alignmentPos) {
      gapCount += region.end - region.start;
    } else if (region.start < alignmentPos) {
      gapCount += alignmentPos - region.start;
      break;
    }
  }
  return alignmentPos - gapCount;
}

/**
 * Converts biological position to alignment position (accounting for gaps in reference).
 * Finds the alignment position that corresponds to the Nth non-gap character.
 * @param biologicalPosInRef 0-based position relative to reference start (without gaps).
 * @param gapRegions Gap regions in alignment coordinates.
 * @param referenceStart Reference start in alignment coordinates.
 * @returns Alignment position (with gaps). 0-based coordinates.
 */
export function calculateAlignmentPosition(
  biologicalPosInRef: number,
  gapRegions: Region[],
  referenceStart: number,
): number {
  if (biologicalPosInRef < referenceStart) return biologicalPosInRef;
  let alignmentPos = referenceStart;
  let biologicalCount = 0;

  // Walk through alignment, counting non-gap positions.
  while (biologicalCount < biologicalPosInRef) {
    if (isPositionInGap(alignmentPos, gapRegions)) {
      alignmentPos++;
    } else {
      biologicalCount++;
      alignmentPos++;
    }
  }

  if (biologicalCount > 0) {
    alignmentPos--;
  }

  return alignmentPos;
}

/**
 * Checks if a position is within a gap region (half-open `[start, end)`).
 * @param alignmentPos Position to check
 * @param gapRegions Array of gap regions
 * @returns true if position is a gap
 */
export function isPositionInGap(alignmentPos: number, gapRegions: Region[]): boolean {
  return gapRegions.some(gap => alignmentPos >= gap.start && alignmentPos < gap.end);
}

/**
 * Calculates the coverage range based on reads only (union of all reads).
 * The reference is NOT considered — only reads determine the coverage.
 *
 * @param alignment Object with the reads to analyze (only offset/sequence are read)
 * @param biologicalOffset Optional offset to add to the returned range
 * @returns Region covered by at least one read, or undefined if no reads
 *
 * @example
 * // Read1: offset=10, length=50 (covers 10-60)
 * // Read2: offset=20, length=60 (covers 20-80)
 * // Result: start=10, end=80 (union of both reads)
 */
export function calculateReadsCoverageRange(
  alignment: { reads: ReadonlyArray<{ offset: number; sequence: string }> },
  biologicalOffset = 0,
): Region | undefined {
  if (alignment.reads.length === 0) return undefined;

  // Union of reads coverage (min start, max end)
  let minStart = Infinity;
  let maxEnd = -Infinity;

  for (const read of alignment.reads) {
    const readStart = read.offset;
    const readEnd = read.offset + countNonGappedLength(read.sequence);
    minStart = Math.min(minStart, readStart);
    maxEnd = Math.max(maxEnd, readEnd);
  }

  if (minStart >= maxEnd) return undefined;

  return { start: minStart + biologicalOffset, end: maxEnd + biologicalOffset };
}
