import { assertTruthy } from '@fishka/assertions';
import { countNonGappedLength } from '../sequence/sequence.utils';
import { ChannelSignals, Chromatogram } from './types';

/** Returns a reversed copy, leaving the input untouched. */
function reversedCopy(values: number[]): number[] {
  return values.slice().reverse();
}

/**
 * Reverse-complement a chromatogram so it reads as the opposite strand.
 *
 * Three coordinated transforms, all reversible (applying this twice is the
 * identity), and a faithful re-view of the existing trace — never a re-basecall:
 *  - signals: the A and T channels swap, and C and G swap (a peak called A on
 *    this strand is a T on its complement), and every channel is reversed so the
 *    trace mirrors left-to-right.
 *  - positions: each peak's sample-point index is mirrored against the signal
 *    length (`L-1-p`); the array is then reversed so positions stay ascending.
 *  - confidences: reversed to track the now-reversed base order.
 */
export function reverseComplementChromatogram(chromatogram: Chromatogram): Chromatogram {
  const { signals } = chromatogram;
  const reversedSignals: ChannelSignals = {
    A: reversedCopy(signals.T),
    T: reversedCopy(signals.A),
    C: reversedCopy(signals.G),
    G: reversedCopy(signals.C),
  };
  // The four channels share one length (ABIF DATA blocks are equal-length); mirror against it.
  const signalLength = reversedSignals.A.length;
  const positions = chromatogram.positions.map(p => signalLength - 1 - p).reverse();
  const confidences = chromatogram.confidences === undefined ? undefined : reversedCopy(chromatogram.confidences);
  return {
    positions,
    signals: reversedSignals,
    confidences,
    samplingRate: chromatogram.samplingRate,
  };
}

/** One run of `count` gap columns to insert at clean-sequence base index `position`. */
interface GapInsertion {
  position: number;
  count: number;
}

/** Signal samples added per inserted gap column — one flat (zero) sample, as in the source model. */
const GAP_SAMPLE_COUNT = 1;

/**
 * Inject alignment gaps into a chromatogram so it lines up, column-for-column, with a gapped
 * sequence from an alignment. `gappedSequence` is the read's bases with '-' gaps; the chromatogram
 * must correspond to the same sequence *without* gaps. Each gap run adds `count` flat columns — new
 * peak positions plus one zero sample per gap in every channel — and shifts the following peaks, so
 * the result indexes directly against the alignment columns. Returns the input unchanged when there
 * are no gaps.
 */
export function applyAlignmentGapsToChromatogram(chromatogram: Chromatogram, gappedSequence: string): Chromatogram {
  const nonGappedLength = countNonGappedLength(gappedSequence);
  assertTruthy(
    chromatogram.positions.length === nonGappedLength,
    `Chromatogram length (${chromatogram.positions.length}) doesn't match clean sequence length (${nonGappedLength}).`,
  );

  if (!gappedSequence.includes('-')) {
    return chromatogram;
  }

  // Collect gap runs as (clean-sequence base index, length).
  const gapsToInsert: GapInsertion[] = [];
  let cleanPos = 0;
  for (let i = 0; i < gappedSequence.length; i++) {
    if (gappedSequence[i] === '-') {
      let count = 0;
      while (i < gappedSequence.length && gappedSequence[i] === '-') {
        count++;
        i++;
      }
      i--; // adjust for the loop increment
      gapsToInsert.push({ position: cleanPos, count });
    } else {
      cleanPos++;
    }
  }

  if (gapsToInsert.length === 0) {
    return chromatogram;
  }

  const newPositions = [...chromatogram.positions];
  const newConfidences = chromatogram.confidences ? [...chromatogram.confidences] : undefined;
  const newSignals: ChannelSignals = {
    A: [...chromatogram.signals.A],
    C: [...chromatogram.signals.C],
    G: [...chromatogram.signals.G],
    T: [...chromatogram.signals.T],
  };

  // Insert from the highest base index down so earlier insertions don't shift the indices of later ones.
  const sortedGaps = [...gapsToInsert].sort((a, b) => b.position - a.position);
  for (const gap of sortedGaps) {
    const { position, count } = gap;
    assertTruthy(
      position >= 0 && position <= newPositions.length,
      `Invalid insertion position ${position} for a sequence of length ${newPositions.length}.`,
    );

    let insertSignalIndex: number;
    if (position === 0) {
      insertSignalIndex = 0;
    } else if (position >= newPositions.length) {
      insertSignalIndex = newSignals.A.length;
    } else {
      insertSignalIndex = newPositions[position];
    }

    const samplesToInsert = count * GAP_SAMPLE_COUNT;
    const zeroSignals = new Array<number>(samplesToInsert).fill(0);
    newSignals.A.splice(insertSignalIndex, 0, ...zeroSignals);
    newSignals.C.splice(insertSignalIndex, 0, ...zeroSignals);
    newSignals.G.splice(insertSignalIndex, 0, ...zeroSignals);
    newSignals.T.splice(insertSignalIndex, 0, ...zeroSignals);

    for (let i = 0; i < newPositions.length; i++) {
      if (newPositions[i] >= insertSignalIndex) {
        newPositions[i] += samplesToInsert;
      }
    }

    // One new peak position per gap column: insertSignalIndex, +1, +2, ...
    const gapPositions = Array.from({ length: count }, (_, i) => insertSignalIndex + i);
    newPositions.splice(position, 0, ...gapPositions);

    if (newConfidences) {
      newConfidences.splice(position, 0, ...new Array<number>(count).fill(0));
    }
  }

  return {
    ...chromatogram,
    positions: newPositions,
    signals: newSignals,
    confidences: newConfidences,
  };
}
