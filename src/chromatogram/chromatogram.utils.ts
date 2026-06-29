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
